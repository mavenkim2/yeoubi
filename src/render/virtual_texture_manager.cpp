#include "render/virtual_texture_manager.h"
#include "core/assert.h"
#include "core/bit.h"
#include "core/ring_buffer.h"
#include "core/thread_context.h"
#include "device/device.h"
#include "device/device_texture.h"
#include "device/kernel_params.h"
#include "display_device/display_device.h"
#include <algorithm>
#include <cctype>
#include <cstring>
#include <filesystem>
#include <fstream>
#include <limits>
#include <unordered_set>

namespace ybi
{

static constexpr uint32_t kInvalidTextureIndex = ~0u;

static std::string GetVirtualTextureFilename(const std::string &filename)
{
    std::filesystem::path path(filename);
    path.replace_extension(".ybvt");
    return path.string();
}

static View<uint8_t>
LoadVirtualTextureTile(Arena *arena, const TextureDescriptor &descriptor, uint32_t pageEntryIndex)
{
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
    hostPinnedFeedbackRequests.resize(numFramesInFlight);
    deviceUpdateRequests.resize(numFramesInFlight);
    feedbackJob.done = true;

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

VirtualTextureManager::~VirtualTextureManager()
{
    // No task may retain a manager buffer, texture, or upload stream when it is released.
    feedbackJob.Wait();
    uploadPool.Wait();

    // The feedback and page-table transfers use the default stream.  Wait for it before
    // freeing their host-pinned and device buffers.
    if (!uploadStreams.empty())
    {
        std::shared_ptr<DeviceEvent> deviceIdle = device->CreateDeviceEvent();
        device->RecordDeviceEvent(nullptr, *deviceIdle);
        device->WaitDeviceEvent(*uploadStreams.front(), *deviceIdle);
        device->SyncDeviceStream(*uploadStreams.front());
    }

    if (feedback)
        device->Free(feedback);
    if (deviceTextureDescriptors)
        device->Free(deviceTextureDescriptors);
    if (deviceTexturePools)
        device->Free(deviceTexturePools);
    if (devicePageDescriptors)
        device->Free(devicePageDescriptors);
    for (DeviceBuffer &buffer : pageDescriptorBits)
    {
        if (buffer)
            device->Free(buffer);
    }
    for (DeviceBuffer &buffer : deviceUpdateRequests)
    {
        if (buffer)
            device->Free(buffer);
    }
    for (DeviceBuffer &buffer : hostPinnedFeedbackRequests)
    {
        if (buffer)
            device->FreeHost(buffer);
    }
    for (DeviceBuffer &buffer : hostPinnedUpdateRequests)
    {
        if (buffer)
            device->FreeHost(buffer);
    }

    for (PagePool &pool : pagePools)
        device->DestroyTexture(pool.texture);

    updatePageDescriptors.reset();
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
// select correct mip on gpu, flat diffuse shading for now

// NOT DONE:
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
        node.textureObject = static_cast<uint64_t>(texture.textureObject);
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
    const int pageShift = 7;
    const TextureDescriptor &textureDescriptor = textureDescriptors[textureIndex];
    const VirtualTextureFileHeader &header = textureDescriptor.tileFileHeader;
    YBI_ASSERT(mip < header.numMips);

    const uint32_t offset =
        CalculatePageDescriptorOffset(header.width, header.height, mip, pageShift, tileX, tileY);
    YBI_ASSERT(offset < header.numPages);
    return static_cast<int>(textureDescriptor.pageDescriptorOffset + offset);
}

MaterialTextureHandle
VirtualTextureManager::AddTexture(const std::vector<std::pair<std::string, uint32_t>> &filenames)
{
    MaterialTextureHandle handle = {};
    if (filenames.empty())
    {
        return handle;
    }

    struct PendingTexture
    {
        uint32_t udim;
        uint32_t textureIndex;
    };

    ScratchArena scratch;
    View<PendingTexture> pendingTextures =
        scratch.arena->PushArray<PendingTexture>(filenames.size());
    uint32_t count = 0;

    uint32_t minUdim = std::numeric_limits<uint32_t>::max();
    uint32_t maxUdim = 0;
    for (const auto &[filename, udim] : filenames)
    {
        YBI_ASSERT(!filename.empty());

        auto itr = filenameToDescriptorMap.find(filename);
        if (itr != filenameToDescriptorMap.end())
        {
            const TextureDescriptor &descriptor = textureDescriptors[itr->second];
            YBI_ASSERT(descriptor.tileFileHeader.udim == udim);
            pendingTextures[count++] = {descriptor.tileFileHeader.udim, itr->second};
            minUdim = Min(minUdim, udim);
            maxUdim = Max(maxUdim, udim);
            continue;
        }

        textureDescriptors.EmplaceBack();
        TextureDescriptor &descriptor = textureDescriptors.Last();
        descriptor.tileFilename = GetVirtualTextureFilename(filename);

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
        YBI_ASSERT(header.width > 0);
        YBI_ASSERT(header.height > 0);
        YBI_ASSERT(header.udim == udim);
        YBI_ASSERT(header.pageSize == static_cast<uint32_t>(pageSize));
        YBI_ASSERT(header.numMips > 0 && header.numMips <= 16);
        YBI_ASSERT(CalculateVirtualTexturePageCount(header.width, 7) <= 256);
        YBI_ASSERT(CalculateVirtualTexturePageCount(header.height, 7) <= 256);

        descriptor.tileFileHeader = header;
        descriptor.pageEntries.Resize(header.numPages);

        input.seekg(static_cast<std::streamoff>(header.pageTableOffset), std::ios::beg);
        input.read(reinterpret_cast<char *>(descriptor.pageEntries.data()),
                   sizeof(VirtualTexturePageEntry) * header.numPages);
        YBI_ASSERT(input.good());

        const uint32_t textureIndex = static_cast<uint32_t>(textureDescriptors.size() - 1);
        filenameToDescriptorMap.emplace(filename, textureIndex);

        pendingTextures[count++] = {udim, textureIndex};
        minUdim = Min(minUdim, udim);
        maxUdim = Max(maxUdim, udim);
    }

    handle.udimBase = minUdim;
    handle.udimTableOffset = static_cast<uint32_t>(udimTextureIndices.size());

    const uint32_t numUdimEntries = maxUdim - minUdim + 1;
    const size_t oldSize = udimTextureIndices.size();
    udimTextureIndices.Resize(oldSize + numUdimEntries);
    for (uint32_t i = 0; i < numUdimEntries; i++)
    {
        udimTextureIndices[oldSize + i] = kInvalidTextureIndex;
    }
    for (const PendingTexture &pending : pendingTextures)
    {
        udimTextureIndices[handle.udimTableOffset + pending.udim - minUdim] = pending.textureIndex;
    }

    return handle;
}

void VirtualTextureManager::Finalize()
{
    uint32_t pageDescriptorOffset = 0;
    for (TextureDescriptor &descriptor : textureDescriptors)
    {
        descriptor.pageDescriptorOffset = pageDescriptorOffset;
        pageDescriptorOffset += static_cast<uint32_t>(descriptor.tileFileHeader.numPages);
    }

    pageDescriptors.Resize(pageDescriptorOffset);
    for (size_t i = 0; i < pageDescriptors.size(); i++)
    {
        pageDescriptors[i].textureObject = 0ull;
        pageDescriptors[i].bits = ~0u;
    }

    hostTextureDescriptors.Resize(udimTextureIndices.size());
    for (size_t i = 0; i < udimTextureIndices.size(); i++)
    {
        DeviceTextureDescriptor &deviceDescriptor = hostTextureDescriptors[i];
        deviceDescriptor = {};
        deviceDescriptor.textureIndex = kInvalidTextureIndex;
        deviceDescriptor.fallback = Vec4f(1.f, 0.f, 1.f, 1.f);

        const uint32_t textureIndex = udimTextureIndices[i];
        if (textureIndex == kInvalidTextureIndex)
        {
            continue;
        }

        const TextureDescriptor &descriptor = textureDescriptors[textureIndex];
        deviceDescriptor.textureIndex = textureIndex;
        deviceDescriptor.format = descriptor.tileFileHeader.format;
        deviceDescriptor.pageDescriptorOffset = descriptor.pageDescriptorOffset;
        deviceDescriptor.pageDescriptorCount = descriptor.tileFileHeader.numPages;
        deviceDescriptor.width = descriptor.tileFileHeader.width;
        deviceDescriptor.height = descriptor.tileFileHeader.height;
        deviceDescriptor.fallback = descriptor.tileFileHeader.fallbackColor;
    }

    if (hostTextureDescriptors.size())
    {
        deviceTextureDescriptors =
            device->Alloc(sizeof(DeviceTextureDescriptor) * hostTextureDescriptors.size());
        device->CopyToDevice(deviceTextureDescriptors.Data(),
                             hostTextureDescriptors.data(),
                             deviceTextureDescriptors.SizeInBytes());
    }

    if (pageDescriptors.size())
    {
        devicePageDescriptors = device->Alloc(sizeof(PageDescriptor) * pageDescriptors.size());
        device->CopyToDevice(devicePageDescriptors.Data(),
                             pageDescriptors.data(),
                             devicePageDescriptors.SizeInBytes());

        const size_t numBitWords = (pageDescriptors.size() + 31) / 32;
        for (DeviceBuffer &buffer : pageDescriptorBits)
        {
            buffer = device->Alloc(sizeof(uint32_t) * numBitWords);
            device->Memzero(buffer.Data(), buffer.SizeInBytes());
        }
    }

    feedbackCapacity = (1u << 20u) - 1u;
    const size_t feedbackSize =
        sizeof(FeedbackRequest) * (static_cast<size_t>(feedbackCapacity) + 1);
    feedback = device->Alloc(feedbackSize);
    device->Memzero(feedback.Data(), sizeof(FeedbackRequest));

    for (DeviceBuffer &hostFeedback : hostPinnedFeedbackRequests)
    {
        hostFeedback = device->AllocHost(feedbackSize);
        memset(hostFeedback.Data(), 0, feedbackSize);
    }
}

void VirtualTextureManager::Update(KernelParams *params)
{
    const uint32_t currentFrame = displayDevice->GetCurrentFrameIndex();
    // TODO: only if async copy is ready

    ScratchArena scratch;

    if (feedbackJob.IsDone() && pagePools.size())
    {
        const size_t texturePoolSize = sizeof(DeviceTexture) * pagePools.size();
        if (deviceTexturePools.SizeInBytes() < texturePoolSize)
        {
            if (deviceTexturePools)
            {
                device->Free(deviceTexturePools);
            }
            deviceTexturePools = device->Alloc(texturePoolSize);
        }

        View<DeviceTexture> hostTexturePools =
            scratch.arena->PushArray<DeviceTexture>(pagePools.size());
        for (size_t i = 0; i < pagePools.size(); i++)
        {
            hostTexturePools[i] = pagePools[i].texture;
        }
        device->CopyToDevice(
            deviceTexturePools.Data(), hostTexturePools.data(), hostTexturePools.numBytes());
    }

    YBI_ASSERT(params);
    params->textureDescriptors = deviceTextureDescriptors.As<DeviceTextureDescriptor>().Get();
    params->pageDescriptors = devicePageDescriptors.As<PageDescriptor>().Get();
    // params->texturePools =
    //     deviceTexturePools ? deviceTexturePools.As<DeviceTexture>().Get() : nullptr;
    params->pageDescriptorBits = pageDescriptorBits[currentFrame].As<uint32_t>().Get();
    params->feedbackRequests = feedback.As<FeedbackRequest>().Get();

    device->Memzero(pageDescriptorBits[currentFrame].Data(),
                    pageDescriptorBits[currentFrame].SizeInBytes());

    device->Memzero(feedback.Data(), sizeof(FeedbackRequest));

    FeedbackRequest *requests =
        reinterpret_cast<FeedbackRequest *>(hostPinnedFeedbackRequests[currentFrame].Data());
    uint32_t numRequests = requests[0].textureIndex;
    numRequests = std::min(numRequests, feedbackCapacity);

    ringBufferFeedbackRequests.SynchronizedWrite(requests + 1, numRequests);

    if (feedbackJob.IsDone())
    {
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

                const uint32_t bits = (node.layer << 27u) | node.pagePoolIndex;
                pageDescriptors[node.pageDescriptorIndex].textureObject = node.textureObject;
                pageDescriptors[node.pageDescriptorIndex].bits = bits;

                lruUsed = lruUsed == -1 ? static_cast<int>(completion.lruIndex) : lruUsed;
                LinkLRU(completion.lruIndex);

                completedUploadRequests.EmplaceBack();
                PageDescriptorUpdateRequest &updateRequest = completedUploadRequests.Last();
                updateRequest.pageDescriptorIndex = node.pageDescriptorIndex;
                updateRequest.textureObject = node.textureObject;
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

                YBI_ERROR(mip != 15, "%u %u %u %u\n", textureIndex, tileX, tileY, mip);

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
                    evictRequest.textureObject = 0ull;
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

            if (evictRequests.size())
            {
                printf("num evivctions: %zu\n", evictRequests.size());
            }
            if (allocationInfos.size())
            {
                ringBufferPageAllocationInfos2.SynchronizedWrite(allocationInfos.data(),
                                                                 allocationInfos.size());
            }
            if (evictRequests.size())
            {
                ringBufferUpdateRequests.SynchronizedWrite(evictRequests.data(),
                                                           evictRequests.size());
            }
        });
    }

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
        ++residencyRevision;

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

            YBI_ASSERT(descriptor.tileFileHeader.compression == 0);

            const uint32_t pageDescriptorIndex =
                CalculatePageDescriptorIndex(info.textureIndex, info.tileX, info.tileY, info.mip);
            const uint32_t pageEntryIndex = pageDescriptorIndex - descriptor.pageDescriptorOffset;

            YBI_ERROR(pageEntryIndex < descriptor.pageEntries.size(),
                      "%u %zu , %u %u %u %u\n",
                      pageEntryIndex,
                      descriptor.pageEntries.size(),
                      info.textureIndex,
                      info.tileX,
                      info.tileY,
                      info.mip);
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
            copyDesc.copy.dstZ = info.layer;

            PageUploadCompletion &completion = completions[infoIndex];
            completion.lruIndex = info.lruIndex;

            // TODO: decompress if descriptor.tileFileHeader.compression != 0.
        }

        device->CopyToDevice3DBatchAsync(stream, copyDescs, numInfos);
        device->SyncDeviceStream(stream);
        ringBufferPageUploadCompletions.SynchronizedWrite(completions, numInfos);
    });
}

void VirtualTextureManager::UpdateParams(KernelParams *params)
{
    params->textureDescriptors = deviceTextureDescriptors.As<DeviceTextureDescriptor>().Get();
    params->pageDescriptors = devicePageDescriptors.As<PageDescriptor>().Get();
    // params->texturePools =
    //     deviceTexturePools ? deviceTexturePools.As<DeviceTexture>().Get() : nullptr;
    params->pageDescriptorBits = pageDescriptorBits[0].As<uint32_t>().Get();
    params->feedbackRequests = feedback.As<FeedbackRequest>().Get();
}

void VirtualTextureManager::CopyFeedback()
{
    const uint32_t currentFrame = displayDevice->GetCurrentFrameIndex();
    YBI_ASSERT(feedback);
    YBI_ASSERT(currentFrame < hostPinnedFeedbackRequests.size());
    DeviceBuffer &hostFeedback = hostPinnedFeedbackRequests[currentFrame];
    YBI_ASSERT(hostFeedback);
    device->CopyToHostAsync(hostFeedback.Data(), feedback.Data(), feedback.SizeInBytes());
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
