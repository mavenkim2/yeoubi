#pragma once

#include "core/array.h"
#include "core/tbb_helpers.h"
#include "device/device_memory.h"
#include "device/device_texture.h"
#include <cstdint>

namespace ybi
{

struct Device;
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

    uint32_t pool;
    uint32_t layer;
};

// NOTE: instantiated once at startup
struct TextureDescriptor
{
    uint32_t pageDescriptorOffset;
    DeviceTextureInfo info;
};

struct PageDescriptor
{
    // texture, page
    uint32_t bits;
};

struct PageDescriptorUpdateRequest
{
    uint32_t pageDescriptorIndex;
    uint32_t bits;
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
    AsyncJob feedbackJob;

    Array<PagePool> pagePools;
    Array<PageDescriptor> pageDescriptors;
    Array<TextureDescriptor> textureDescriptors;
    size_t totalPhysicalAllocation;
    size_t physicalAllocationLimit;
    int pageSize;
    TaskPool uploadPool;

    int CalculatePageDescriptorIndex(uint32_t textureIndex,
                                     uint32_t tileX,
                                     uint32_t tileY,
                                     uint32_t mip) const;
    void UnlinkLRU(int entryIndex);
    void LinkLRU(int index);
    bool AllocatePage(uint32_t textureIndex, uint32_t &poolIndex, uint32_t &layer);

public:
    VirtualTextureManager();
    void Update();
};

} // namespace ybi
