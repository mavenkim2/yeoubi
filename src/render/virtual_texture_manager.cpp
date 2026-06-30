#include "render/virtual_texture_manager.h"
#include "core/assert.h"
#include "core/bit.h"
#include "core/ring_buffer.h"
#include "core/thread_context.h"
#include "device/device.h"
#include "display_device/display_device.h"
#include <unordered_set>

namespace ybi
{

VirtualTextureManager::VirtualTextureManager() : uploadPool(8) {}

// DONE:
// make lru node not 1to1 with page descriptors, instead based on page pool entries
// update the lru logic when i unlink and relink lrus for resident pages
// use my new allocation to allocate resident pages, and to evict when unavailable
// create a ringbuffer

// NOT DONE:
// actually write the file in the offline tiler, with requisite header information to load at
// runtime
// write feedback to ring buffer, and read feedback on tbb task group thread
// have thread write back to another ring buffer, which is read on host
// host dispatches a small kernel that updates the page descriptors on the gpu
// dispatch threads to load tiles, upload them synchronously on independent streams to gpu
// select correct mip on gpu, flat diffuse shading for now
// write feedback to buffer, async copy to disk

bool VirtualTextureManager::AllocatePage(uint32_t textureIndex,
                                         uint32_t &poolIndex,
                                         uint32_t &layer)
{
    const TextureDescriptor &descriptor = textureDescriptors[textureIndex];
    for (size_t index = 0; index < pagePools.size(); index++)
    {
        PagePool &pool = pagePools[index];
        if (pool.texture.info.format != descriptor.info.format)
        {
            continue;
        }
        if (pool.freeList.size())
        {
            layer = pool.freeList.Pop();
            poolIndex = index;
            return true;
        }
    }

    // TODO: not hardcoded
    const int physicalPoolNumLayers = 32;
    const int physicalPoolByteSize =
        pageSize * pageSize * physicalPoolNumLayers * GetTextureFormatSize(descriptor.info.format);
    if (physicalPoolByteSize + totalPhysicalAllocation > physicalAllocationLimit)
    {
        return false;
    }

    DeviceTextureInfo info;
    info.format = descriptor.info.format;
    info.filterType = FilterType::Point;
    info.width = pageSize;
    info.height = pageSize;
    info.layers = physicalPoolNumLayers;

    DeviceTexture texture = device->CreateTexture(info);

    int pagePoolIndex = static_cast<int>(pagePools.size());

    pagePools.EmplaceBack();
    PagePool &outPool = pagePools.Last();
    outPool.texture = texture;

    for (int i = physicalPoolNumLayers - 1; i >= 0; i--)
    {
        outPool.freeList.EmplaceBack(i);
    }
    outPool.lruIndexStart = lruNodes.size();

    for (int i = 0; i < physicalPoolNumLayers; i++)
    {
        lruNodes.EmplaceBack();
        LRUNode &node = lruNodes.Last();
        node.lruNext = -1;
        node.lruPrev = -1;
        node.pagePoolIndex = pagePoolIndex;
        node.layer = i;
    }
    layer = outPool.freeList.Pop();
    return true;
}

int VirtualTextureManager::CalculatePageDescriptorIndex(uint32_t textureIndex,
                                                        uint32_t tileX,
                                                        uint32_t tileY,
                                                        uint32_t mip) const
{
    // TODO
    const int pageSize = 128;
    const int pageShift = 7;
    const TextureDescriptor &textureDescriptor = textureDescriptors[textureIndex];

    if (IsPowerOfTwo(textureDescriptor.info.width) && IsPowerOfTwo(textureDescriptor.info.height))
    {
        uint32_t offset = 0;
        uint32_t width = textureDescriptor.info.width;
        uint32_t height = textureDescriptor.info.height;

        for (uint32_t i = 0; i < mip; i++)
        {
            offset += (width >> pageShift) * (height >> pageShift);
            width >>= 1;
            height >>= 1;
        }
        offset += tileY * (width >> pageShift) + tileX;
        return textureDescriptor.pageDescriptorOffset + offset;
    }
    else
    {
        YBI_LOGFATAL("Non power of two textures not supported\n");
        return 0;
    }
}

void VirtualTextureManager::Update()
{
    const uint32_t currentFrame = displayDevice->GetCurrentFrameIndex();
    // TODO: only if async copy is ready

    ScratchArena scratch;

    FeedbackRequest *requests = hostFeedbackRequests[currentFrame];
    uint32_t numRequests = requests[0].textureIndex;

    View<FeedbackRequest> compactedFeedbackRequests =
        scratch.arena->PushArray<FeedbackRequest>(numRequests);

    feedbackJob.Wait();

    feedbackJob.Run([this, compactedFeedbackRequests, numRequests, requests]() {
        std::unordered_set<uint64_t> seenRequests;
        uint32_t numCompacted = 0;
        // Compact duplicate feedback requests
        for (uint32_t i = 1; i < numRequests + 1; i++)
        {
            const FeedbackRequest &request = requests[i];
            uint32_t textureIndex = request.textureIndex;
            uint32_t tileX = request.tileX_tileY_mip & 0xff;
            uint32_t tileY = (request.tileX_tileY_mip >> 8u) & 0xff;
            uint32_t mip = (request.tileX_tileY_mip >> 16u) & 0xf;

            uint64_t key = request.textureIndex | (uint64_t(request.tileX_tileY_mip) << 32u);

            if (seenRequests.insert(key).second)
            {
                compactedFeedbackRequests[numCompacted++] = request;
            }
        }

        // Update LRU for resident pages, compact non resident requests
        uint32_t numNonResidentFeedback = 0;
        int lruUsed = -1;
        for (uint32_t i = 0; i < numCompacted; i++)
        {
            const FeedbackRequest &request = compactedFeedbackRequests[i];
            uint32_t textureIndex = request.textureIndex;
            uint32_t tileX = request.tileX_tileY_mip & 0xff;
            uint32_t tileY = (request.tileX_tileY_mip >> 8u) & 0xff;
            uint32_t mip = (request.tileX_tileY_mip >> 16u) & 0xf;

            int pageDescriptorIndex =
                CalculatePageDescriptorIndex(textureIndex, tileX, tileY, mip);

            PageDescriptor &descriptor = pageDescriptors[pageDescriptorIndex];
            if (descriptor.bits != ~0u)
            {
                uint32_t poolIndex = descriptor.bits & 0x07ffffffu;
                uint32_t layerIndex = descriptor.bits >> 27u;

                uint32_t lruIndex = pagePools[poolIndex].lruIndexStart + layerIndex;
                LRUNode &node = lruNodes[lruIndex];
                YBI_ASSERT(node.lruPrev != -1 && node.lruNext != -1);
                lruUsed = lruUsed == -1 ? static_cast<int>(lruIndex) : lruUsed;

                UnlinkLRU(lruIndex);
                LinkLRU(lruIndex);
            }
            else
            {
                compactedFeedbackRequests[numNonResidentFeedback++] = request;
            }
        }

        Array<PageAllocationInfo> allocationInfos;
        Array<PageDescriptorUpdateRequest> evictRequests;

        // Find physical memory location for non resident pages
        for (uint32_t i = 0; i < numNonResidentFeedback; i++)
        {
            const FeedbackRequest &request = compactedFeedbackRequests[i];
            uint32_t textureIndex = request.textureIndex;
            uint32_t tileX = request.tileX_tileY_mip & 0xff;
            uint32_t tileY = (request.tileX_tileY_mip >> 8u) & 0xff;
            uint32_t mip = (request.tileX_tileY_mip >> 16u) & 0xf;

            uint32_t layer, poolIndex;
            bool allocated = AllocatePage(textureIndex, poolIndex, layer);

            while (!allocated)
            {
                LRUNode &node = lruNodes[lruTail];
                int entryToEvict = node.lruPrev;

                if (entryToEvict == lruUsed)
                {
                    break;
                }

                UnlinkLRU(entryToEvict);

                pagePools[node.pagePoolIndex].freeList.EmplaceBack(node.layer);
                evictRequests.EmplaceBack();
                PageDescriptorUpdateRequest &evictRequest = evictRequests.Last();
                evictRequest.bits = ~0u;
                evictRequest.pageDescriptorIndex = node.pageDescriptorIndex;

                pageDescriptors[node.pageDescriptorIndex].bits = ~0u;

                allocated = AllocatePage(textureIndex, poolIndex, layer);
            }

            if (!allocated)
                continue;

            PagePool &pool = pagePools[poolIndex];
            LRUNode &node = lruNodes[pool.lruIndexStart + layer];
            node.pageDescriptorIndex =
                CalculatePageDescriptorIndex(textureIndex, tileX, tileY, mip);
            LinkLRU(pool.lruIndexStart + layer);

            allocationInfos.EmplaceBack();
            PageAllocationInfo &info = allocationInfos.Last();
            info.textureIndex = textureIndex;
            info.tileX = tileX;
            info.tileY = tileY;
            info.mip = mip;
            info.pool = &pool - pagePools.data();
            info.layer = layer;

            // TODO IMPORTANT: in my previous project i write this info to a ring buffer,
            // and then on the main thread flush the ring buffer and write to another one.
            // why don't I just use one ring buffer?
        }
    });

    // Update GPU page descriptors
    uploadPool.Run([]() {
        int threadIndex = tbb::this_task_arena::current_thread_index();
        YBI_ASSERT(threadIndex >= 0);

        // Get device stream
    });
}

void VirtualTextureManager::UnlinkLRU(int entryIndex)
{
    LRUNode &node = lruNodes[entryIndex];
    YBI_ASSERT(node.lruPrev != -1 && node.lruNext != -1);

    int prevEntry = node.lruPrev;
    int nextEntry = node.lruNext;

    lruNodes[prevEntry].lruNext = nextEntry;
    lruNodes[nextEntry].lruPrev = prevEntry;
}

void VirtualTextureManager::LinkLRU(int index)
{
    int nextEntry = lruNodes[lruHead].lruNext;
    lruNodes[index].lruNext = nextEntry;
    lruNodes[index].lruPrev = lruHead;

    lruNodes[nextEntry].lruPrev = index;
    lruNodes[lruHead].lruNext = index;
}

} // namespace ybi
