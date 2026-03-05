#pragma once

#include "device/device.h"
#include "render/launch_params.h"
#include "texture/virtual_texture/tile_file.h"

#include <array>
#include <cstdint>
#include <list>
#include <string>
#include <unordered_map>
#include <vector>

namespace ybi
{
namespace texture
{

struct VirtualTextureManagerConfig
{
    uint32_t pageSize = 128u;
    uint64_t cacheBytes = 1ull << 30u;
    uint32_t tailMaxDim = 32u;
    uint32_t maxUploadsPerPass = 1024u;
};

struct VirtualTextureTailSource
{
    uint32_t udim = 1001u;
    uint32_t width = 0u;
    uint32_t height = 0u;
    const unsigned char *rgba8 = nullptr;
    size_t rgba8Bytes = 0u;
};

struct VirtualTextureUdimExtent
{
    uint32_t udim = 1001u;
    uint32_t width = 0u;
    uint32_t height = 0u;
};

struct VirtualTextureRegisterInput
{
    uint32_t textureId = 0u;
    uint32_t width = 0u;
    uint32_t height = 0u;
    std::string tileFilePath;
    std::vector<uint32_t> activeUdims;
    std::vector<VirtualTextureUdimExtent> udimExtents;
    std::vector<VirtualTextureTailSource> tailSources;
};

struct VirtualTextureUpdateStats
{
    uint32_t feedbackCount = 0u;
    uint32_t uniqueCount = 0u;
    uint32_t hits = 0u;
    uint32_t misses = 0u;
    uint32_t uploads = 0u;
    uint32_t evictions = 0u;
    uint32_t failed = 0u;
};

class VirtualTextureManager
{
public:
    VirtualTextureManager() = default;
    ~VirtualTextureManager();

    bool Initialize(Device *device, const VirtualTextureManagerConfig &config, std::string *outError);
    bool RegisterTexture(const VirtualTextureRegisterInput &input, std::string *outError);
    bool Finalize(std::string *outError);
    void BindLaunchParams(LaunchParams *params) const;
    bool ProcessFeedback(const unsigned long long *keys,
                         uint32_t keyCount,
                         VirtualTextureUpdateStats *outStats,
                         std::string *outError);
    void Shutdown();

private:
    struct TextureState
    {
        uint32_t textureId = 0u;
        uint32_t width = 0u;
        uint32_t height = 0u;
        uint32_t mipCount = 0u;
        uint32_t tailFirstMip = 0u;
        std::array<int16_t, 128> udimToLocal = {};
        std::vector<uint32_t> activeUdims;
        std::string tileFilePath;
        VirtualTextureTileFile tileFile = {};
        bool tileFileOpen = false;

        std::vector<LaunchParams::VirtualTextureMipInfo> mipInfos;
        uint32_t mipInfoOffset = 0u;
        std::vector<LaunchParams::VirtualTextureUdimInfo> udimInfos;
        uint32_t udimInfoOffset = 0u;
        std::vector<uint32_t> udimWidths;
        std::vector<uint32_t> udimHeights;

        std::vector<uint8_t> tailPixelsHost;
        DeviceMemoryView<uint8_t> tailPixelsDevice = {};
        uint32_t tailPageCountX = 0u;
        uint32_t tailPageCountY = 0u;
        std::vector<uint32_t> tailPageIndexByLocalUdim;
    };

    struct StreamSlotState
    {
        bool occupied = false;
        unsigned long long key = 0ull;
        std::list<uint32_t>::iterator lruIt;
    };

    struct KeyVirtualInfo
    {
        bool valid = false;
        bool isTail = false;
        uint32_t textureIndex = 0u;
        uint32_t textureId = 0u;
        uint32_t mip = 0u;
        uint32_t udim = 1001u;
        uint32_t udimLocal = 0u;
        uint32_t tileX = 0u;
        uint32_t tileY = 0u;
        uint32_t vaX = 0u;
        uint32_t vaY = 0u;
    };

    struct LevelRectRequest
    {
        uint32_t textureIndex = 0u;
        uint32_t mip = 0u;
        uint32_t localUdim = 0u;
        uint32_t width = 0u;
        uint32_t height = 0u;
    };

    bool BuildMipReservations(std::string *outError);
    bool BuildTailPagesForTexture(TextureState *texture,
                                  const VirtualTextureRegisterInput &input,
                                  std::string *outError);
    bool AllocateDeviceState(std::string *outError);
    bool OpenTileFileIfNeeded(TextureState *texture, std::string *outError);
    bool ResolveKey(unsigned long long key, KeyVirtualInfo *outInfo) const;
    uint32_t PackPageTableEntry(uint32_t pageX, uint32_t pageY, uint32_t pageType, uint32_t flags) const;
    void UnpackPageTableEntry(uint32_t packed,
                              uint32_t *outPageX,
                              uint32_t *outPageY,
                              uint32_t *outPageType,
                              uint32_t *outFlags) const;
    bool UpdatePageTableTexel(uint32_t mip, uint32_t x, uint32_t y, uint32_t entry, std::string *outError);
    bool UploadStreamPage(uint32_t slotIndex,
                          const std::vector<unsigned char> &rgba8,
                          uint32_t width,
                          uint32_t height,
                          std::string *outError);
    bool LoadStreamPageForKey(const KeyVirtualInfo &info,
                              std::vector<unsigned char> *outRgba8,
                              uint32_t *outWidth,
                              uint32_t *outHeight,
                              std::string *outError);
    bool TouchResidentKey(unsigned long long key);
    bool EvictOne(uint32_t *outSlotIndex, std::string *outError);
    bool AllocateStreamSlot(uint32_t *outSlotIndex, bool *outEvicted, std::string *outError);

    Device *device_ = nullptr;
    VirtualTextureManagerConfig config_ = {};
    bool finalized_ = false;

    std::vector<TextureState> textures_;
    std::unordered_map<uint32_t, uint32_t> textureIdToIndex_;

    std::vector<uint32_t> pageTableMipOffsets_;
    std::vector<uint32_t> pageTableMipWidths_;
    std::vector<uint32_t> pageTableMipHeights_;
    std::vector<uint32_t> pageTableEntriesHost_;
    DeviceMemoryView<uint8_t> pageTableEntriesDevice_ = {};
    DeviceMemoryView<uint8_t> pageTableMipOffsetsDevice_ = {};
    DeviceMemoryView<uint8_t> pageTableMipWidthsDevice_ = {};
    DeviceMemoryView<uint8_t> pageTableMipHeightsDevice_ = {};

    std::vector<LaunchParams::VirtualTextureMipInfo> mipInfosHost_;
    DeviceMemoryView<uint8_t> mipInfosDevice_ = {};
    std::vector<LaunchParams::VirtualTextureUdimInfo> udimInfosHost_;
    DeviceMemoryView<uint8_t> udimInfosDevice_ = {};
    std::vector<LaunchParams::VirtualTextureTextureMeta> textureMetaHost_;
    DeviceMemoryView<uint8_t> textureMetaDevice_ = {};

    uint32_t streamPageCountX_ = 0u;
    uint32_t streamPageCountY_ = 0u;
    uint32_t streamPageCount_ = 0u;
    std::vector<uint8_t> streamPixelsHost_;
    DeviceMemoryView<uint8_t> streamPixelsDevice_ = {};

    std::vector<StreamSlotState> streamSlots_;
    std::unordered_map<unsigned long long, uint32_t> keyToStreamSlot_;
    std::list<uint32_t> lruSlots_;
    std::vector<uint32_t> freeSlots_;
};

} // namespace texture
} // namespace ybi
