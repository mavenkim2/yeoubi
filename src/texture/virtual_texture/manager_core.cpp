#include "texture/virtual_texture/manager.h"

#include "texture/virtual_texture/feedback.h"
#include "texture/virtual_texture/key.h"

#include <algorithm>
#include <cmath>
#include <cstring>

namespace ybi
{
namespace texture
{

namespace
{
static constexpr uint32_t kPageTypeStream = 1u;
static constexpr uint32_t kPageTypeTail = 2u;
static constexpr uint32_t kUdimMin = 1001u;
static constexpr uint32_t kUdimMax = 1128u;

static uint32_t NextPow2(uint32_t v)
{
    if (v <= 1u)
    {
        return 1u;
    }
    v--;
    v |= v >> 1u;
    v |= v >> 2u;
    v |= v >> 4u;
    v |= v >> 8u;
    v |= v >> 16u;
    return v + 1u;
}

static uint32_t ComputeMipCount(uint32_t width, uint32_t height)
{
    uint32_t w = std::max(1u, width);
    uint32_t h = std::max(1u, height);
    uint32_t count = 1u;
    while (w > 1u || h > 1u)
    {
        w = std::max(1u, w >> 1u);
        h = std::max(1u, h >> 1u);
        count++;
    }
    return count;
}

static void ResizeToPage(const std::vector<uint8_t> &src,
                         uint32_t srcW,
                         uint32_t srcH,
                         uint32_t pageSize,
                         uint8_t *dstPage)
{
    std::memset(dstPage, 0, static_cast<size_t>(pageSize) * static_cast<size_t>(pageSize) * 4u);
    for (uint32_t y = 0; y < pageSize; ++y)
    {
        const uint32_t sy = std::min(srcH - 1u, (y * srcH) / pageSize);
        for (uint32_t x = 0; x < pageSize; ++x)
        {
            const uint32_t sx = std::min(srcW - 1u, (x * srcW) / pageSize);
            const size_t srcBase =
                (static_cast<size_t>(sy) * static_cast<size_t>(srcW) + static_cast<size_t>(sx)) * 4u;
            const size_t dstBase =
                (static_cast<size_t>(y) * static_cast<size_t>(pageSize) + static_cast<size_t>(x)) * 4u;
            dstPage[dstBase + 0u] = src[srcBase + 0u];
            dstPage[dstBase + 1u] = src[srcBase + 1u];
            dstPage[dstBase + 2u] = src[srcBase + 2u];
            dstPage[dstBase + 3u] = src[srcBase + 3u];
        }
    }
}

} // namespace

VirtualTextureManager::~VirtualTextureManager()
{
    Shutdown();
}

bool VirtualTextureManager::Initialize(Device *device,
                                       const VirtualTextureManagerConfig &config,
                                       std::string *outError)
{
    Shutdown();
    if (!device)
    {
        if (outError)
        {
            *outError = "VirtualTextureManager Initialize: null device";
        }
        return false;
    }
    device_ = device;
    config_ = config;
    config_.pageSize = 128u;
    if (config_.cacheBytes == 0u)
    {
        config_.cacheBytes = 1ull << 30u;
    }
    if (config_.tailMaxDim == 0u)
    {
        config_.tailMaxDim = 32u;
    }
    finalized_ = false;
    return true;
}

bool VirtualTextureManager::BuildTailPagesForTexture(TextureState *texture, std::string *outError)
{
    if (!texture)
    {
        return false;
    }

    texture->tailPixelsHost.clear();
    texture->tailPageIndexByLocalUdim.assign(texture->activeUdims.size(), UINT32_MAX);
    texture->tailPageCountX = 0u;
    texture->tailPageCountY = 0u;

    if (texture->tailFirstMip >= texture->mipCount)
    {
        return true;
    }

    if (!OpenTileFileIfNeeded(texture, outError))
    {
        return false;
    }

    const size_t pageBytes =
        static_cast<size_t>(config_.pageSize) * static_cast<size_t>(config_.pageSize) * 4u;
    uint32_t tailPageCount = 0u;
    for (size_t local = 0; local < texture->activeUdims.size(); ++local)
    {
        const uint32_t udim = texture->activeUdims[local];
        std::vector<unsigned char> tailRgba8;
        uint32_t tailW = 0u;
        uint32_t tailH = 0u;
        uint32_t tailMip = 0u;
        uint64_t tailSourceBytes = 0u;
        std::string tailError;
        if (!ReadVirtualTextureTailMip(&texture->tileFile,
                                       udim,
                                       config_.tailMaxDim,
                                       &tailRgba8,
                                       &tailW,
                                       &tailH,
                                       &tailMip,
                                       &tailSourceBytes,
                                       &tailError))
        {
            if (outError)
            {
                *outError = "VirtualTextureManager: required tail mip missing (textureId=" +
                            std::to_string(texture->textureId) +
                            " udim=" + std::to_string(udim) + "): " + tailError;
            }
            return false;
        }
        if (tailRgba8.empty() || tailW == 0u || tailH == 0u)
        {
            if (outError)
            {
                *outError = "VirtualTextureManager: invalid tail mip data (textureId=" +
                            std::to_string(texture->textureId) +
                            " udim=" + std::to_string(udim) + ")";
            }
            return false;
        }

        texture->tailPageIndexByLocalUdim[local] = tailPageCount++;
        const size_t oldSize = texture->tailPixelsHost.size();
        texture->tailPixelsHost.resize(oldSize + pageBytes, 0u);
        ResizeToPage(tailRgba8, tailW, tailH, config_.pageSize, texture->tailPixelsHost.data() + oldSize);
    }

    if (tailPageCount > 0u)
    {
        texture->tailPageCountX = tailPageCount;
        texture->tailPageCountY = 1u;
    }
    return true;
}

bool VirtualTextureManager::RegisterTexture(const VirtualTextureRegisterInput &input, std::string *outError)
{
    if (!device_)
    {
        if (outError)
        {
            *outError = "VirtualTextureManager RegisterTexture: manager not initialized";
        }
        return false;
    }
    if (finalized_)
    {
        if (outError)
        {
            *outError = "VirtualTextureManager RegisterTexture: called after Finalize";
        }
        return false;
    }
    if (input.activeUdims.empty() || input.width == 0u || input.height == 0u)
    {
        return true;
    }

    TextureState texture = {};
    texture.textureId = input.textureId;
    texture.width = input.width;
    texture.height = input.height;
    texture.tileFilePath = input.tileFilePath;
    texture.activeUdims = input.activeUdims;
    texture.mipCount = ComputeMipCount(texture.width, texture.height);
    texture.udimToLocal.fill(-1);
    texture.udimWidths.assign(texture.activeUdims.size(), texture.width);
    texture.udimHeights.assign(texture.activeUdims.size(), texture.height);
    std::unordered_map<uint32_t, VirtualTextureUdimExtent> extentByUdim;
    extentByUdim.reserve(input.udimExtents.size());
    for (const VirtualTextureUdimExtent &extent : input.udimExtents)
    {
        if (extent.width == 0u || extent.height == 0u)
        {
            continue;
        }
        extentByUdim[extent.udim] = extent;
    }
    for (size_t i = 0; i < texture.activeUdims.size(); ++i)
    {
        const uint32_t udim = texture.activeUdims[i];
        if (udim < kUdimMin || udim > kUdimMax)
        {
            continue;
        }
        texture.udimToLocal[static_cast<size_t>(udim - kUdimMin)] = static_cast<int16_t>(i);
        auto extentIt = extentByUdim.find(udim);
        if (extentIt != extentByUdim.end())
        {
            texture.udimWidths[i] = extentIt->second.width;
            texture.udimHeights[i] = extentIt->second.height;
        }
    }

    uint32_t tailFirstMip = texture.mipCount;
    for (uint32_t mip = 0; mip < texture.mipCount; ++mip)
    {
        const uint32_t mipW = std::max(1u, texture.width >> mip);
        const uint32_t mipH = std::max(1u, texture.height >> mip);
        if (std::max(mipW, mipH) <= config_.tailMaxDim)
        {
            tailFirstMip = mip;
            break;
        }
    }
    texture.tailFirstMip = tailFirstMip;

    if (!BuildTailPagesForTexture(&texture, outError))
    {
        if (outError && outError->empty())
        {
            *outError = "VirtualTextureManager RegisterTexture: failed building tail pages";
        }
        return false;
    }

    textureIdToIndex_[texture.textureId] = static_cast<uint32_t>(textures_.size());
    textures_.push_back(std::move(texture));
    return true;
}

bool VirtualTextureManager::BuildMipReservations(std::string *outError)
{
    uint32_t maxMipCount = 0u;
    for (const TextureState &texture : textures_)
    {
        maxMipCount = std::max(maxMipCount, texture.mipCount);
    }
    if (maxMipCount == 0u)
    {
        maxMipCount = 1u;
    }

    std::vector<std::vector<LevelRectRequest>> levelRequests(maxMipCount);
    for (uint32_t textureIndex = 0u; textureIndex < textures_.size(); ++textureIndex)
    {
        TextureState &texture = textures_[textureIndex];
        texture.mipInfos.resize(texture.mipCount);
        const uint32_t udimCount = static_cast<uint32_t>(texture.activeUdims.size());
        texture.udimInfos.assign(static_cast<size_t>(texture.mipCount) * static_cast<size_t>(udimCount),
                                 LaunchParams::VirtualTextureUdimInfo{});
        for (uint32_t mip = 0u; mip < texture.mipCount; ++mip)
        {
            LaunchParams::VirtualTextureMipInfo &mipInfo = texture.mipInfos[mip];
            mipInfo.level = mip;
            mipInfo.udimInfoOffset = mip * udimCount;
            mipInfo.udimInfoCount = udimCount;
            mipInfo._padding0 = 0u;
            for (uint32_t localUdim = 0u; localUdim < udimCount; ++localUdim)
            {
                const uint32_t udimW =
                    (localUdim < texture.udimWidths.size()) ? texture.udimWidths[localUdim] : texture.width;
                const uint32_t udimH =
                    (localUdim < texture.udimHeights.size()) ? texture.udimHeights[localUdim] : texture.height;
                const uint32_t mipW = std::max(1u, udimW >> mip);
                const uint32_t mipH = std::max(1u, udimH >> mip);
                const uint32_t rectW = std::max(1u, (mipW + config_.pageSize - 1u) / config_.pageSize);
                const uint32_t rectH = std::max(1u, (mipH + config_.pageSize - 1u) / config_.pageSize);
                levelRequests[mip].push_back({textureIndex, mip, localUdim, rectW, rectH});
            }
        }
    }

    pageTableMipWidths_.assign(maxMipCount, 1u);
    pageTableMipHeights_.assign(maxMipCount, 1u);
    pageTableMipOffsets_.assign(maxMipCount, 0u);

    for (uint32_t level = 0u; level < maxMipCount; ++level)
    {
        std::vector<LevelRectRequest> &requests = levelRequests[level];
        if (requests.empty())
        {
            pageTableMipWidths_[level] = 1u;
            pageTableMipHeights_[level] = 1u;
            continue;
        }

        uint64_t area = 0u;
        uint32_t maxRectW = 1u;
        for (const LevelRectRequest &request : requests)
        {
            area += static_cast<uint64_t>(request.width) * static_cast<uint64_t>(request.height);
            maxRectW = std::max(maxRectW, request.width);
        }
        uint32_t width = NextPow2(std::max(maxRectW, static_cast<uint32_t>(std::sqrt(double(area)) + 1.0)));
        width = std::max(width, 8u);

        bool packed = false;
        while (!packed)
        {
            uint32_t cursorX = 0u;
            uint32_t cursorY = 0u;
            uint32_t rowHeight = 0u;
            packed = true;
            for (LevelRectRequest &request : requests)
            {
                if (request.width > width)
                {
                    packed = false;
                    break;
                }
                if (cursorX + request.width > width)
                {
                    cursorX = 0u;
                    cursorY += rowHeight;
                    rowHeight = 0u;
                }
                TextureState &texture = textures_[request.textureIndex];
                LaunchParams::VirtualTextureMipInfo &mipInfo = texture.mipInfos[request.mip];
                const uint32_t udimCount = static_cast<uint32_t>(texture.activeUdims.size());
                const uint32_t localOffset = request.mip * udimCount + request.localUdim;
                if (localOffset >= texture.udimInfos.size())
                {
                    packed = false;
                    break;
                }
                LaunchParams::VirtualTextureUdimInfo &udimInfo = texture.udimInfos[localOffset];
                mipInfo.level = level;
                udimInfo.basePageX = cursorX;
                udimInfo.basePageY = cursorY;
                udimInfo.pageCountX = request.width;
                udimInfo.pageCountY = request.height;
                cursorX += request.width;
                rowHeight = std::max(rowHeight, request.height);
            }
            if (packed)
            {
                pageTableMipWidths_[level] = width;
                pageTableMipHeights_[level] = std::max(1u, cursorY + rowHeight);
                break;
            }
            width = std::max(width * 2u, 8u);
            if (width > (1u << 16u))
            {
                if (outError)
                {
                    *outError = "VirtualTextureManager: page table level too wide";
                }
                return false;
            }
        }
    }

    uint32_t offset = 0u;
    for (uint32_t level = 0u; level < maxMipCount; ++level)
    {
        pageTableMipOffsets_[level] = offset;
        const uint64_t levelSize =
            static_cast<uint64_t>(pageTableMipWidths_[level]) * static_cast<uint64_t>(pageTableMipHeights_[level]);
        offset += static_cast<uint32_t>(levelSize);
    }
    pageTableEntriesHost_.assign(offset, 0u);

    mipInfosHost_.clear();
    udimInfosHost_.clear();
    for (TextureState &texture : textures_)
    {
        texture.mipInfoOffset = static_cast<uint32_t>(mipInfosHost_.size());
        texture.udimInfoOffset = static_cast<uint32_t>(udimInfosHost_.size());
        for (LaunchParams::VirtualTextureMipInfo mipInfo : texture.mipInfos)
        {
            mipInfo.udimInfoOffset += texture.udimInfoOffset;
            mipInfosHost_.push_back(mipInfo);
        }
        udimInfosHost_.insert(udimInfosHost_.end(), texture.udimInfos.begin(), texture.udimInfos.end());
    }

    uint32_t maxTextureId = 0u;
    for (const TextureState &texture : textures_)
    {
        maxTextureId = std::max(maxTextureId, texture.textureId);
    }
    textureMetaHost_.assign(maxTextureId + 1u, LaunchParams::VirtualTextureTextureMeta{});
    for (LaunchParams::VirtualTextureTextureMeta &meta : textureMetaHost_)
    {
        for (short &v : meta.udimToLocal)
        {
            v = -1;
        }
    }

    for (const TextureState &texture : textures_)
    {
        LaunchParams::VirtualTextureTextureMeta &meta = textureMetaHost_[texture.textureId];
        meta.mipInfoOffset = texture.mipInfoOffset;
        meta.mipCount = texture.mipCount;
        meta.tailFirstMip = texture.tailFirstMip;
        meta.activeUdimCount = static_cast<uint32_t>(texture.activeUdims.size());
        meta.tailPixels = 0ull;
        meta.tailPageCountX = texture.tailPageCountX;
        meta.tailPageCountY = texture.tailPageCountY;
        for (size_t i = 0; i < texture.udimToLocal.size(); ++i)
        {
            meta.udimToLocal[i] = texture.udimToLocal[i];
        }
    }

    for (const TextureState &texture : textures_)
    {
        if (texture.tailPageIndexByLocalUdim.empty())
        {
            continue;
        }
        for (uint32_t mip = texture.tailFirstMip; mip < texture.mipCount; ++mip)
        {
            const LaunchParams::VirtualTextureMipInfo &mipInfo =
                mipInfosHost_[texture.mipInfoOffset + mip];
            const uint32_t level = mipInfo.level;
            for (uint32_t localUdim = 0u; localUdim < texture.tailPageIndexByLocalUdim.size(); ++localUdim)
            {
                const uint32_t tailPageIndex = texture.tailPageIndexByLocalUdim[localUdim];
                if (tailPageIndex == UINT32_MAX)
                {
                    continue;
                }
                if (localUdim >= mipInfo.udimInfoCount)
                {
                    continue;
                }
                const uint32_t udimInfoIndex = mipInfo.udimInfoOffset + localUdim;
                if (udimInfoIndex >= udimInfosHost_.size())
                {
                    continue;
                }
                const LaunchParams::VirtualTextureUdimInfo &udimInfo = udimInfosHost_[udimInfoIndex];
                const uint32_t pageX = tailPageIndex;
                const uint32_t pageY = 0u;
                for (uint32_t ty = 0u; ty < udimInfo.pageCountY; ++ty)
                {
                    for (uint32_t tx = 0u; tx < udimInfo.pageCountX; ++tx)
                    {
                        const uint32_t vaX = udimInfo.basePageX + tx;
                        const uint32_t vaY = udimInfo.basePageY + ty;
                        const uint32_t packed = PackPageTableEntry(pageX, pageY, kPageTypeTail, 1u);
                        if (level >= pageTableMipOffsets_.size() || level >= pageTableMipWidths_.size() ||
                            level >= pageTableMipHeights_.size() || vaX >= pageTableMipWidths_[level] ||
                            vaY >= pageTableMipHeights_[level])
                        {
                            continue;
                        }
                        const uint32_t index = pageTableMipOffsets_[level] +
                                               vaY * pageTableMipWidths_[level] + vaX;
                        if (index < pageTableEntriesHost_.size())
                        {
                            pageTableEntriesHost_[index] = packed;
                        }
                    }
                }
            }
        }
    }
    return true;
}

} // namespace texture
} // namespace ybi
