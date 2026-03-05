#include "texture/virtual_texture/manager.h"

#include "texture/virtual_texture/feedback.h"
#include "texture/virtual_texture/key.h"

#include <algorithm>
#include <cstring>

namespace ybi
{
namespace texture
{

namespace
{
static constexpr uint32_t kPageTypeStream = 1u;
static constexpr uint32_t kUdimMin = 1001u;
static constexpr uint32_t kUdimMax = 1128u;
} // namespace

bool VirtualTextureManager::OpenTileFileIfNeeded(TextureState *texture, std::string *outError)
{
    if (!texture)
    {
        return false;
    }
    if (texture->tileFileOpen)
    {
        return true;
    }
    if (texture->tileFilePath.empty())
    {
        if (outError)
        {
            *outError = "VirtualTextureManager: empty tile file path";
        }
        return false;
    }
    std::string openError;
    if (!OpenVirtualTextureTileFile(texture->tileFilePath, &texture->tileFile, &openError))
    {
        if (outError)
        {
            *outError = openError;
        }
        return false;
    }
    texture->tileFileOpen = true;
    return true;
}

bool VirtualTextureManager::ResolveKey(unsigned long long key, KeyVirtualInfo *outInfo) const
{
    if (!outInfo)
    {
        return false;
    }
    *outInfo = {};
    const uint32_t textureId = FeedbackTextureId(key);
    auto textureIt = textureIdToIndex_.find(textureId);
    if (textureIt == textureIdToIndex_.end())
    {
        return false;
    }
    const TextureState &texture = textures_[textureIt->second];
    const uint32_t mip = FeedbackMip(key);
    if (mip >= texture.mipCount)
    {
        return false;
    }

    const uint32_t udim = FeedbackUdim(key);
    if (udim < kUdimMin || udim > kUdimMax)
    {
        return false;
    }
    const int16_t local = texture.udimToLocal[static_cast<size_t>(udim - kUdimMin)];
    if (local < 0)
    {
        return false;
    }

    const LaunchParams::VirtualTextureMipInfo &mipInfo = texture.mipInfos[mip];
    const uint32_t tileX = FeedbackTileX(key);
    const uint32_t tileY = FeedbackTileY(key);
    if (tileX >= mipInfo.pagesPerUdimX || tileY >= mipInfo.pagesPerUdimY)
    {
        return false;
    }

    outInfo->valid = true;
    outInfo->isTail = (mip >= texture.tailFirstMip);
    outInfo->textureIndex = textureIt->second;
    outInfo->textureId = textureId;
    outInfo->mip = mip;
    outInfo->udim = udim;
    outInfo->udimLocal = static_cast<uint32_t>(local);
    outInfo->tileX = tileX;
    outInfo->tileY = tileY;
    outInfo->vaX = mipInfo.basePageX + outInfo->udimLocal * mipInfo.pagesPerUdimX + tileX;
    outInfo->vaY = mipInfo.basePageY + tileY;
    return true;
}

bool VirtualTextureManager::UploadStreamPage(uint32_t slotIndex,
                                             const std::vector<unsigned char> &rgba8,
                                             uint32_t width,
                                             uint32_t height,
                                             std::string *outError)
{
    if (slotIndex >= streamPageCount_)
    {
        if (outError)
        {
            *outError = "VirtualTextureManager: stream slot out of range";
        }
        return false;
    }
    const uint32_t copyW = std::min(width, config_.pageSize);
    const uint32_t copyH = std::min(height, config_.pageSize);
    const size_t pageBytes =
        static_cast<size_t>(config_.pageSize) * static_cast<size_t>(config_.pageSize) * 4u;
    const size_t pageOffset = static_cast<size_t>(slotIndex) * pageBytes;
    std::memset(streamPixelsHost_.data() + pageOffset, 0, pageBytes);
    for (uint32_t y = 0u; y < copyH; ++y)
    {
        const size_t srcRow = static_cast<size_t>(y) * static_cast<size_t>(width) * 4u;
        const size_t dstRow = pageOffset + static_cast<size_t>(y) * static_cast<size_t>(config_.pageSize) * 4u;
        std::memcpy(streamPixelsHost_.data() + dstRow, rgba8.data() + srcRow, static_cast<size_t>(copyW) * 4u);
    }

    DeviceMemoryView<uint8_t> dst = {streamPixelsDevice_.data() + pageOffset, pageBytes};
    device_->CopyBytesToDevice(dst, streamPixelsHost_.data() + pageOffset, pageBytes);
    return true;
}

bool VirtualTextureManager::LoadStreamPageForKey(const KeyVirtualInfo &info,
                                                 std::vector<unsigned char> *outRgba8,
                                                 uint32_t *outWidth,
                                                 uint32_t *outHeight,
                                                 std::string *outError)
{
    if (!outRgba8 || !outWidth || !outHeight)
    {
        return false;
    }
    TextureState &texture = textures_[info.textureIndex];
    if (!OpenTileFileIfNeeded(&texture, outError))
    {
        return false;
    }

    std::vector<unsigned char> rgba8;
    uint32_t width = 0u;
    uint32_t height = 0u;
    uint64_t sourceBytes = 0u;
    std::string readError;

    uint32_t sourceTileX = info.tileX;
    uint32_t sourceTileY = info.tileY;
    if (info.mip > 0u)
    {
        sourceTileX = info.tileX << info.mip;
        sourceTileY = info.tileY << info.mip;
    }

    if (!ReadVirtualTextureTile(&texture.tileFile,
                                info.udim,
                                sourceTileX,
                                sourceTileY,
                                &rgba8,
                                &width,
                                &height,
                                &sourceBytes,
                                &readError))
    {
        if (info.mip > 0u &&
            ReadVirtualTextureTile(&texture.tileFile,
                                   info.udim,
                                   info.tileX,
                                   info.tileY,
                                   &rgba8,
                                   &width,
                                   &height,
                                   &sourceBytes,
                                   &readError))
        {
            *outRgba8 = std::move(rgba8);
            *outWidth = width;
            *outHeight = height;
            return true;
        }
        if (outError)
        {
            *outError = readError;
        }
        return false;
    }

    *outRgba8 = std::move(rgba8);
    *outWidth = width;
    *outHeight = height;
    return true;
}

bool VirtualTextureManager::TouchResidentKey(unsigned long long key)
{
    auto slotIt = keyToStreamSlot_.find(key);
    if (slotIt == keyToStreamSlot_.end())
    {
        return false;
    }
    const uint32_t slotIndex = slotIt->second;
    if (slotIndex >= streamSlots_.size())
    {
        return false;
    }
    StreamSlotState &slot = streamSlots_[slotIndex];
    if (!slot.occupied)
    {
        return false;
    }
    lruSlots_.splice(lruSlots_.begin(), lruSlots_, slot.lruIt);
    slot.lruIt = lruSlots_.begin();
    return true;
}

bool VirtualTextureManager::EvictOne(uint32_t *outSlotIndex, std::string *outError)
{
    if (lruSlots_.empty())
    {
        if (outError)
        {
            *outError = "VirtualTextureManager: LRU empty and no free stream pages";
        }
        return false;
    }
    const uint32_t slotIndex = lruSlots_.back();
    lruSlots_.pop_back();
    if (slotIndex >= streamSlots_.size())
    {
        return false;
    }
    StreamSlotState &slot = streamSlots_[slotIndex];
    const unsigned long long oldKey = slot.key;
    slot.occupied = false;
    keyToStreamSlot_.erase(oldKey);

    KeyVirtualInfo info = {};
    if (ResolveKey(oldKey, &info) && info.valid && !info.isTail)
    {
        (void)UpdatePageTableTexel(info.mip, info.vaX, info.vaY, 0u, outError);
    }
    *outSlotIndex = slotIndex;
    return true;
}

bool VirtualTextureManager::AllocateStreamSlot(uint32_t *outSlotIndex,
                                               bool *outEvicted,
                                               std::string *outError)
{
    if (!freeSlots_.empty())
    {
        const uint32_t slot = freeSlots_.back();
        freeSlots_.pop_back();
        *outSlotIndex = slot;
        *outEvicted = false;
        return true;
    }
    if (EvictOne(outSlotIndex, outError))
    {
        *outEvicted = true;
        return true;
    }
    return false;
}

bool VirtualTextureManager::ProcessFeedback(const unsigned long long *keys,
                                            uint32_t keyCount,
                                            VirtualTextureUpdateStats *outStats,
                                            std::string *outError)
{
    if (outStats)
    {
        *outStats = {};
    }
    if (!finalized_ || !keys || keyCount == 0u)
    {
        return true;
    }

    std::vector<unsigned long long> keyVec(keys, keys + keyCount);
    std::vector<VirtualTextureFeedbackEntry> histogram;
    BuildFeedbackHistogram(keyVec, keyCount, &histogram);

    VirtualTextureUpdateStats stats = {};
    stats.feedbackCount = keyCount;
    stats.uniqueCount = static_cast<uint32_t>(histogram.size());

    uint32_t uploadsBudget = config_.maxUploadsPerPass;
    for (const VirtualTextureFeedbackEntry &entry : histogram)
    {
        const unsigned long long key = entry.key;
        KeyVirtualInfo info = {};
        if (!ResolveKey(key, &info) || !info.valid)
        {
            stats.failed++;
            continue;
        }
        if (info.isTail)
        {
            stats.hits++;
            continue;
        }
        if (TouchResidentKey(key))
        {
            stats.hits++;
            continue;
        }

        stats.misses++;
        if (uploadsBudget == 0u)
        {
            continue;
        }

        uint32_t slotIndex = 0u;
        bool evicted = false;
        if (!AllocateStreamSlot(&slotIndex, &evicted, outError))
        {
            stats.failed++;
            continue;
        }
        if (evicted)
        {
            stats.evictions++;
        }

        std::vector<unsigned char> rgba8;
        uint32_t width = 0u;
        uint32_t height = 0u;
        if (!LoadStreamPageForKey(info, &rgba8, &width, &height, outError))
        {
            stats.failed++;
            freeSlots_.push_back(slotIndex);
            continue;
        }
        if (!UploadStreamPage(slotIndex, rgba8, width, height, outError))
        {
            stats.failed++;
            freeSlots_.push_back(slotIndex);
            continue;
        }

        const uint32_t pageX = slotIndex % streamPageCountX_;
        const uint32_t pageY = slotIndex / streamPageCountX_;
        const uint32_t packed = PackPageTableEntry(pageX, pageY, kPageTypeStream, 1u);
        if (!UpdatePageTableTexel(info.mip, info.vaX, info.vaY, packed, outError))
        {
            stats.failed++;
            freeSlots_.push_back(slotIndex);
            continue;
        }

        StreamSlotState &slot = streamSlots_[slotIndex];
        slot.occupied = true;
        slot.key = key;
        lruSlots_.push_front(slotIndex);
        slot.lruIt = lruSlots_.begin();
        keyToStreamSlot_[key] = slotIndex;

        stats.uploads++;
        uploadsBudget--;
    }

    if (outStats)
    {
        *outStats = stats;
    }
    return true;
}

} // namespace texture
} // namespace ybi
