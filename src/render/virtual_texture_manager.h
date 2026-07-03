#pragma once

#include "core/array.h"
#include "core/ring_buffer.h"
#include "core/tbb_helpers.h"
#include "device/device_memory.h"
#include "device/device_texture.h"
#include "render/virtual_texture_file.h"
#include "render/virtual_texture_interop.h"
#include <cstdint>
#include <memory>
#include <mutex>
#include <string>
#include <unordered_set>
#include <vector>

namespace ybi
{

struct Device;
struct DeviceKernel;
struct DeviceStream;
struct DeviceTexture;
struct DisplayDevice;

struct FeedbackRequest
{
    uint32_t textureIndex;
    // 8, 8, 4 bits
    uint32_t tileX_tileY_mip;
};

struct PageAllocationInfo
{
    uint32_t textureIndex;
    uint32_t tileX;
    uint32_t tileY;
    uint32_t mip;

    DeviceTexture texture;
    uint32_t pool;
    uint32_t layer;
    uint32_t lruIndex;
};

struct PageUploadCompletion
{
    uint32_t lruIndex;
};

// NOTE: instantiated once at startup
struct TextureDescriptor
{
    uint32_t pageDescriptorOffset;
    std::string tileFilename;
    VirtualTextureFileHeader tileFileHeader = {};
    Array<VirtualTexturePageEntry> pageEntries;
};

// NOTE: 1 to 1 with CPU page descriptors
struct LRUNode
{
    uint32_t pagePoolIndex;
    uint32_t layer;

    uint32_t pageDescriptorIndex;

    uint32_t lruPrev;
    uint32_t lruNext;
};

struct PagePool
{
    DeviceTexture texture;
    Array<uint32_t> freeList;
    int lruIndexStart;
};

class VirtualTextureManager
{
    Device *device;
    DisplayDevice *displayDevice;

    // NOTE: CUDA pinned memory
    FeedbackRequest **hostFeedbackRequests;
    DeviceBuffer feedback;

    // LRU
    Array<LRUNode> lruNodes;
    int lruHead;
    int lruTail;

    Array<PagePool> pagePools;
    Array<PageDescriptor> pageDescriptors;
    Array<TextureDescriptor> textureDescriptors;
    size_t totalPhysicalAllocation;
    size_t physicalAllocationLimit;
    int pageSize;

    // Device
    std::unique_ptr<DeviceKernel> updatePageDescriptors;
    DeviceBuffer devicePageDescriptors;
    std::vector<DeviceBuffer> hostPinnedUpdateRequests;
    std::vector<DeviceBuffer> deviceUpdateRequests;

    // Synchronization
    AsyncJob feedbackJob;
    TaskPool uploadPool;
    std::vector<std::unique_ptr<DeviceStream>> uploadStreams;
    std::unordered_set<uint32_t> pendingPageUploads;
    RingBuffer<FeedbackRequest> ringBufferFeedbackRequests;
    RingBuffer<PageAllocationInfo> ringBufferPageAllocationInfos2;
    RingBuffer<PageAllocationInfo> ringBufferPageAllocationInfos;
    RingBuffer<PageUploadCompletion> ringBufferPageUploadCompletions;
    RingBuffer<PageDescriptorUpdateRequest> ringBufferUpdateRequests;

    int CalculatePageDescriptorIndex(uint32_t textureIndex,
                                     uint32_t tileX,
                                     uint32_t tileY,
                                     uint32_t mip) const;
    void UnlinkLRU(int entryIndex);
    void LinkLRU(int index);
    bool AllocatePage(uint32_t textureIndex, uint32_t &poolIndex, uint32_t &layer);

public:
    VirtualTextureManager(Arena &arena, Device *device, DisplayDevice *displayDevice);
    size_t AddTexture(const std::string &filename);
    void Update();
};

} // namespace ybi
