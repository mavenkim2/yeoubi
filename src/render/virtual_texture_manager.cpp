#include "render/virtual_texture_manager.h"
#include "core/assert.h"
#include "core/bit.h"
#include "core/ring_buffer.h"
#include "core/thread_context.h"
#include "device/device.h"
#include "device/device_texture.h"
#include "display_device/display_device.h"
#include <cstring>
#include <fstream>
#include <unordered_set>

namespace ybi
{

static View<uint8_t>
LoadVirtualTextureTile(Arena *arena, const TextureDescriptor &descriptor, uint32_t pageEntryIndex)
{
    YBI_ASSERT(pageEntryIndex < descriptor.pageEntries.size());

    const VirtualTexturePageEntry &entry = descriptor.pageEntries[pageEntryIndex];
    YBI_ASSERT(entry.packedSize > 0);

    View<uint8_t> tileBytes = arena->PushArray<uint8_t>(entry.packedSize);

    std::ifstream input(descriptor.tileFilename, std::ios::binary);
    YBI_ERROR(input.is_open(),
              "Failed to open virtual texture file: %s\n",
              descriptor.tileFilename.c_str());

    input.seekg(static_cast<std::streamoff>(entry.offset), std::ios::beg);
    input.read(reinterpret_cast<char *>(tileBytes.data()), entry.packedSize);
    YBI_ASSERT(input.good());

    return tileBytes;
}

VirtualTextureManager::VirtualTextureManager(Arena &arena,
                                             Device *device,
                                             DisplayDevice *displayDevice)
    : device(device), displayDevice(displayDevice), totalPhysicalAllocation(0),
      physicalAllocationLimit(2ull * 1024ull * 1024ull * 1024ull), pageSize(128), uploadPool(8),
      ringBufferFeedbackRequests(&arena, 1u << 20u),
      ringBufferPageAllocationInfos(&arena, 1u << 20u),
      ringBufferPageAllocationInfos2(&arena, 1u << 20u),
      ringBufferPageUploadCompletions(&arena, 1u << 20u),
      ringBufferUpdateRequests(&arena, 1u << 20u)

{
    DeviceKernelDesc desc;
    desc.filename = YBI_VIRTUAL_TEXTURE_KERNEL_PATH;
    desc.entryName = "UpdatePageDescriptors";
    updatePageDescriptors = device->CreateDeviceKernel(desc);

    uploadStreams.reserve(uploadPool.GetMaxThreads());
    for (int streamIndex = 0; streamIndex < uploadPool.GetMaxThreads(); streamIndex++)
    {
        uploadStreams.push_back(device->CreateDeviceStream());
    }

    const uint32_t numFramesInFlight = displayDevice->GetNumFramesInFlight();
    hostPinnedUpdateRequests.resize(numFramesInFlight);
    deviceUpdateRequests.resize(numFramesInFlight);

    lruHead = 0;
    lruTail = 1;

    LRUNode node = {};
    node.lruPrev = -1;
    node.lruNext = lruTail;
    lruNodes.EmplaceBack(node);
    node.lruPrev = lruHead;
    node.lruNext = -1;
    lruNodes.EmplaceBack(node);
}

// DONE:
// make lru node not 1to1 with page descriptors, instead based on page pool entries
// update the lru logic when i unlink and relink lrus for resident pages
// use my new allocation to allocate resident pages, and to evict when unavailable
// create a ringbuffer
// have thread write back to another ring buffer, which is read on host
// host dispatches a small kernel that updates the page descriptors on the gpu
// write feedback to ring buffer, and read feedback on tbb task group thread
// actually write the file in the offline tiler, with requisite header information to load at
// runtime
// dispatch threads to load tiles, upload them synchronously on independent streams to gpu

// NOT DONE:
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
        if (pool.texture.info.format != descriptor.tileFileHeader.format)
        {
            continue;
        }
        if (pool.freeList.size())
        {
            poolIndex = index;
            layer = pool.freeList.Pop();
            return true;
        }
    }

    // TODO: not hardcoded
    const int physicalPoolNumLayers = 32;
    const int physicalPoolByteSize = pageSize * pageSize * physicalPoolNumLayers *
                                     GetTextureFormatSize(descriptor.tileFileHeader.format);
    if (physicalPoolByteSize + totalPhysicalAllocation > physicalAllocationLimit)
    {
        return false;
    }

    DeviceTextureInfo info;
    info.format = descriptor.tileFileHeader.format;
    info.filterType = FilterType::Point;
    info.width = pageSize;
    info.height = pageSize;
    info.layers = physicalPoolNumLayers;

    DeviceTexture texture = device->CreateTexture(info);
    totalPhysicalAllocation += physicalPoolByteSize;

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
    poolIndex = pagePoolIndex;
    layer = outPool.freeList.Pop();
    return true;
}

int VirtualTextureManager::CalculatePageDescriptorIndex(uint32_t textureIndex,
                                                        uint32_t tileX,
                                                        uint32_t tileY,
                                                        uint32_t mip) const
{
    // TODO
    const int pageShift = 7;
    const TextureDescriptor &textureDescriptor = textureDescriptors[textureIndex];

    uint32_t width = textureDescriptor.tileFileHeader.width;
    uint32_t height = textureDescriptor.tileFileHeader.height;
    if (IsPowerOfTwo(width) && IsPowerOfTwo(height))
    {
        uint32_t offset = 0;

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

size_t VirtualTextureManager::AddTexture(const std::string &filename)
{
    YBI_ASSERT(!filename.empty());

    textureDescriptors.EmplaceBack();
    TextureDescriptor &descriptor = textureDescriptors.Last();
    descriptor.tileFilename = filename;

    std::ifstream input(descriptor.tileFilename, std::ios::binary);
    YBI_ERROR(input.is_open(),
              "Failed to open virtual texture file: %s\n",
              descriptor.tileFilename.c_str());

    VirtualTextureFileHeader header = {};
    input.read(reinterpret_cast<char *>(&header), sizeof(header));
    YBI_ASSERT(input.good());
    YBI_ASSERT(header.magic == YBVT_MAGIC);
    YBI_ASSERT(header.version == 1);
    YBI_ASSERT(header.numPages > 0);

    descriptor.tileFileHeader = header;
    descriptor.pageEntries.Resize(header.numPages);

    input.seekg(static_cast<std::streamoff>(header.pageTableOffset), std::ios::beg);
    input.read(reinterpret_cast<char *>(descriptor.pageEntries.data()),
               sizeof(VirtualTexturePageEntry) * header.numPages);
    YBI_ASSERT(input.good());

    return textureDescriptors.size() - 1;
}

void VirtualTextureManager::Update()
{
    const uint32_t currentFrame = displayDevice->GetCurrentFrameIndex();
    // TODO: only if async copy is ready

    ScratchArena scratch;

    FeedbackRequest *requests = hostFeedbackRequests[currentFrame];
    uint32_t numRequests = requests[0].textureIndex;

    ringBufferFeedbackRequests.SynchronizedWrite(requests + 1, numRequests);

    feedbackJob.Wait();

    feedbackJob.Run([this]() {
        std::unordered_set<uint64_t> seenRequests;

        ScratchArena scratch;

        int lruUsed = -1;
        uint32_t numCompletions = 0;
        PageUploadCompletion *completions =
            ringBufferPageUploadCompletions.SynchronizedRead(scratch.arena, numCompletions);
        Array<PageDescriptorUpdateRequest> completedUploadRequests;
        for (uint32_t completionIndex = 0; completionIndex < numCompletions; completionIndex++)
        {
            const PageUploadCompletion &completion = completions[completionIndex];
            YBI_ASSERT(completion.lruIndex < lruNodes.size());

            LRUNode &node = lruNodes[completion.lruIndex];
            if (pendingPageUploads.erase(node.pageDescriptorIndex) == 0)
            {
                continue;
            }

            YBI_ASSERT(node.lruPrev == -1 && node.lruNext == -1);

            const uint32_t bits = node.pagePoolIndex | (node.layer << 27u);
            pageDescriptors[node.pageDescriptorIndex].bits = bits;

            lruUsed = lruUsed == -1 ? static_cast<int>(completion.lruIndex) : lruUsed;
            LinkLRU(completion.lruIndex);

            completedUploadRequests.EmplaceBack();
            PageDescriptorUpdateRequest &updateRequest = completedUploadRequests.Last();
            updateRequest.pageDescriptorIndex = node.pageDescriptorIndex;
            updateRequest.bits = bits;
        }

        if (completedUploadRequests.size())
        {
            ringBufferUpdateRequests.SynchronizedWrite(completedUploadRequests.data(),
                                                       completedUploadRequests.size());
        }

        uint32_t numRequests = 0;
        FeedbackRequest *requests =
            ringBufferFeedbackRequests.SynchronizedRead(scratch.arena, numRequests);
        uint32_t numCompacted = 0;
        // Compact duplicate feedback requests
        for (uint32_t i = 0; i < numRequests; i++)
        {
            const FeedbackRequest &request = requests[i];
            uint32_t textureIndex = request.textureIndex;
            uint32_t tileX = request.tileX_tileY_mip & 0xff;
            uint32_t tileY = (request.tileX_tileY_mip >> 8u) & 0xff;
            uint32_t mip = (request.tileX_tileY_mip >> 16u) & 0xf;

            uint64_t key = request.textureIndex | (uint64_t(request.tileX_tileY_mip) << 32u);

            if (seenRequests.insert(key).second)
            {
                requests[numCompacted++] = request;
            }
        }

        // Update LRU for resident pages, compact non resident requests
        uint32_t numNonResidentFeedback = 0;
        for (uint32_t i = 0; i < numCompacted; i++)
        {
            const FeedbackRequest &request = requests[i];
            uint32_t textureIndex = request.textureIndex;
            uint32_t tileX = request.tileX_tileY_mip & 0xff;
            uint32_t tileY = (request.tileX_tileY_mip >> 8u) & 0xff;
            uint32_t mip = (request.tileX_tileY_mip >> 16u) & 0xf;

            int pageDescriptorIndex =
                CalculatePageDescriptorIndex(textureIndex, tileX, tileY, mip);

            PageDescriptor &descriptor = pageDescriptors[pageDescriptorIndex];
            if (pendingPageUploads.find(pageDescriptorIndex) != pendingPageUploads.end())
            {
                continue;
            }
            else if (descriptor.bits != ~0u)
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
                requests[numNonResidentFeedback++] = request;
            }
        }

        Array<PageAllocationInfo> allocationInfos;
        Array<PageDescriptorUpdateRequest> evictRequests;

        // Find physical memory location for non resident pages
        for (uint32_t i = 0; i < numNonResidentFeedback; i++)
        {
            const FeedbackRequest &request = requests[i];
            uint32_t textureIndex = request.textureIndex;
            uint32_t tileX = request.tileX_tileY_mip & 0xff;
            uint32_t tileY = (request.tileX_tileY_mip >> 8u) & 0xff;
            uint32_t mip = (request.tileX_tileY_mip >> 16u) & 0xf;
            uint32_t pageDescriptorIndex =
                CalculatePageDescriptorIndex(textureIndex, tileX, tileY, mip);

            uint32_t layer, poolIndex;
            bool allocated = AllocatePage(textureIndex, poolIndex, layer);

            while (!allocated)
            {
                int entryToEvict = lruNodes[lruTail].lruPrev;

                if (entryToEvict == lruHead || entryToEvict == lruUsed)
                {
                    break;
                }

                UnlinkLRU(entryToEvict);

                LRUNode &node = lruNodes[entryToEvict];
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

            // TODO: very unlikely, but what happens when entry is evicted before it's
            // even uploaded? then you could potentially have two writes to the same
            // GPU page descriptor, leading to invalid memory

            PagePool &pool = pagePools[poolIndex];
            uint32_t lruIndex = pool.lruIndexStart + layer;
            LRUNode &node = lruNodes[lruIndex];
            YBI_ASSERT(node.lruPrev == -1 && node.lruNext == -1);
            node.pageDescriptorIndex = pageDescriptorIndex;
            pendingPageUploads.insert(pageDescriptorIndex);

            allocationInfos.EmplaceBack();
            PageAllocationInfo &info = allocationInfos.Last();
            info.textureIndex = textureIndex;
            info.tileX = tileX;
            info.tileY = tileY;
            info.mip = mip;
            info.texture = pool.texture;
            info.pool = &pool - pagePools.data();
            info.layer = layer;
            info.lruIndex = lruIndex;
        }
        if (allocationInfos.size())
        {
            ringBufferPageAllocationInfos2.SynchronizedWrite(allocationInfos.data(),
                                                             allocationInfos.size());
        }
        if (evictRequests.size())
        {
            ringBufferUpdateRequests.SynchronizedWrite(evictRequests.data(), evictRequests.size());
        }
    });

    uint32_t numInfos = 0;
    PageAllocationInfo *infos =
        ringBufferPageAllocationInfos2.SynchronizedRead(scratch.arena, numInfos);

    uint32_t numUpdateRequests = 0;
    PageDescriptorUpdateRequest *updateRequests =
        ringBufferUpdateRequests.SynchronizedRead(scratch.arena, numUpdateRequests);
    std::shared_ptr<DeviceEvent> pageTableUpdatesDone;

    // NOTE: this ensures evictions always occur before copying new pages to device
    if (numInfos)
    {
        ringBufferPageAllocationInfos.SynchronizedWrite(infos, numInfos);
    }

    // Update device page descriptors
    if (numUpdateRequests)
    {
        DeviceBuffer &hostPinnedUpdateRequestBuffer = hostPinnedUpdateRequests[currentFrame];
        DeviceBuffer &deviceUpdateRequestBuffer = deviceUpdateRequests[currentFrame];

        size_t updateRequestSize = sizeof(PageDescriptorUpdateRequest) * numUpdateRequests;
        if (deviceUpdateRequestBuffer.sizeInBytes < updateRequestSize)
        {
            if (deviceUpdateRequestBuffer.sizeInBytes)
            {
                device->Free(deviceUpdateRequestBuffer);
            }
            deviceUpdateRequestBuffer = device->Alloc(updateRequestSize);
        }

        if (hostPinnedUpdateRequestBuffer.sizeInBytes < updateRequestSize)
        {
            if (hostPinnedUpdateRequestBuffer.sizeInBytes)
            {
                device->FreeHost(hostPinnedUpdateRequestBuffer);
            }
            hostPinnedUpdateRequestBuffer = device->AllocHost(updateRequestSize);
        }

        memcpy(hostPinnedUpdateRequestBuffer.Data(), updateRequests, updateRequestSize);
        device->CopyToDeviceAsync(deviceUpdateRequestBuffer.Data(),
                                  hostPinnedUpdateRequestBuffer.Data(),
                                  updateRequestSize);

        void *arguments[] = {
            devicePageDescriptors.KernelArg(),
            deviceUpdateRequestBuffer.KernelArg(),
            &numUpdateRequests,
        };

        DeviceKernelLaunchDesc desc;
        desc.gridX = (numUpdateRequests + 31) / 32;
        desc.blockX = 32;
        desc.arguments = arguments;
        device->ExecuteKernel(*updatePageDescriptors, desc);

        pageTableUpdatesDone = device->CreateDeviceEvent();
        device->RecordDeviceEvent(nullptr, *pageTableUpdatesDone);
    }

    // Load pages
    uploadPool.Run([this, pageTableUpdatesDone]() {
        ScratchArena scratch;
        int threadIndex = tbb::this_task_arena::current_thread_index();
        YBI_ASSERT(threadIndex >= 0);
        YBI_ASSERT(static_cast<size_t>(threadIndex) < uploadStreams.size());
        DeviceStream &stream = *uploadStreams[threadIndex];

        uint32_t numInfos = 1024;
        PageAllocationInfo *infos =
            ringBufferPageAllocationInfos.SynchronizedRead(scratch.arena, numInfos);

        if (numInfos == 0)
            return;

        if (pageTableUpdatesDone)
        {
            device->WaitDeviceEvent(stream, *pageTableUpdatesDone);
        }

        DeviceTextureBatchCopyDesc *copyDescs =
            scratch.arena->PushArray<DeviceTextureBatchCopyDesc>(numInfos).data();

        PageUploadCompletion *completions =
            scratch.arena->PushArray<PageUploadCompletion>(numInfos).data();

        for (uint32_t infoIndex = 0; infoIndex < numInfos; infoIndex++)
        {
            PageAllocationInfo &info = infos[infoIndex];
            const TextureDescriptor &descriptor = textureDescriptors[info.textureIndex];
            // LoadVirtualTextureTileFile(descriptor);

            YBI_ASSERT(descriptor.tileFileHeader.compression == 0);

            const uint32_t pageDescriptorIndex =
                CalculatePageDescriptorIndex(info.textureIndex, info.tileX, info.tileY, info.mip);
            const uint32_t pageEntryIndex = pageDescriptorIndex - descriptor.pageDescriptorOffset;
            View<uint8_t> tileBytes =
                LoadVirtualTextureTile(scratch.arena, descriptor, pageEntryIndex);

            const uint32_t bytesPerTexel = GetTextureFormatSize(descriptor.tileFileHeader.format);
            const uint32_t tileSize = descriptor.tileFileHeader.pageSize;
            YBI_ASSERT(tileBytes.size() ==
                       tileSize * tileSize * static_cast<size_t>(bytesPerTexel));

            DeviceTextureBatchCopyDesc &copyDesc = copyDescs[infoIndex];
            copyDesc = {};
            copyDesc.texture = &info.texture;
            copyDesc.copy.src = tileBytes.data();
            copyDesc.copy.srcRowPitch = tileSize * bytesPerTexel;
            copyDesc.copy.widthInBytes = tileSize * bytesPerTexel;
            copyDesc.copy.height = tileSize;
            copyDesc.copy.dstLayer = info.layer;

            PageUploadCompletion &completion = completions[infoIndex];
            completion.lruIndex = info.lruIndex;

            // TODO: decompress if descriptor.tileFileHeader.compression != 0.
        }

        device->CopyToDevice3DBatchAsync(stream, copyDescs, numInfos);
        device->SyncDeviceStream(stream);
        ringBufferPageUploadCompletions.SynchronizedWrite(completions, numInfos);
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
    node.lruPrev = -1;
    node.lruNext = -1;
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
