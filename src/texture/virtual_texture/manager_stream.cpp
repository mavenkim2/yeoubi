#include "texture/virtual_texture/manager.h"

#include "texture/virtual_texture/feedback.h"
#include "texture/virtual_texture/key.h"

#include <tbb/parallel_for.h>

#include <algorithm>
#include <chrono>
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

using Clock = std::chrono::steady_clock;
static constexpr size_t kParallelLoadBatchSize = 512u;

double ElapsedMs(Clock::time_point start, Clock::time_point end)
{
    return std::chrono::duration<double, std::milli>(end - start).count();
}

void RecordFirstError(std::string *dst, const std::string &src)
{
    if (!dst || dst->empty() == false || src.empty())
    {
        return;
    }
    *dst = src;
}
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
    bool havePixelFormat = false;
    TextureFormat pixelFormat = TextureFormat::RGBA32_FLOAT;
    for (const auto &entry : texture->tileFile.udims)
    {
        if (!havePixelFormat)
        {
            pixelFormat = entry.second.pixelFormat;
            havePixelFormat = true;
            continue;
        }
        if (entry.second.pixelFormat != pixelFormat)
        {
            if (outError)
            {
                *outError = "VirtualTextureManager: tile file has mixed pixel formats: " +
                            texture->tileFilePath;
            }
            return false;
        }
    }
    texture->pixelFormat = pixelFormat;
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
    if (static_cast<uint32_t>(local) >= mipInfo.udimInfoCount)
    {
        return false;
    }
    const uint32_t udimInfoIndex = mipInfo.udimInfoOffset + static_cast<uint32_t>(local);
    if (udimInfoIndex >= texture.udimInfos.size())
    {
        return false;
    }
    const LaunchParams::VirtualTextureUdimInfo &udimInfo = texture.udimInfos[udimInfoIndex];
    const uint32_t tileX = FeedbackTileX(key);
    const uint32_t tileY = FeedbackTileY(key);
    if (tileX >= udimInfo.pageCountX || tileY >= udimInfo.pageCountY)
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
    outInfo->vaX = udimInfo.basePageX + tileX;
    outInfo->vaY = udimInfo.basePageY + tileY;
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
    const uint32_t physicalTextureID = slotIndex / streamPageCountX_;
    const uint32_t page = slotIndex % streamPageCountX_;
    if (physicalTextureID >= streamTextures_.size() || !streamTextures_[physicalTextureID].valid)
    {
        if (outError)
        {
            *outError = "VirtualTextureManager: stream texture out of range";
        }
        return false;
    }

    std::vector<uint8_t> pagePixels(pageBytes, 0u);
    for (uint32_t y = 0u; y < copyH; ++y)
    {
        const size_t srcRow = static_cast<size_t>(y) * static_cast<size_t>(width) * 4u;
        const size_t dstRow = static_cast<size_t>(y) * static_cast<size_t>(config_.pageSize) * 4u;
        std::memcpy(pagePixels.data() + dstRow, rgba8.data() + srcRow, static_cast<size_t>(copyW) * 4u);
    }
    return device_->UpdateTextureRegion(streamTextures_[physicalTextureID],
                                        page * config_.pageSize,
                                        0u,
                                        config_.pageSize,
                                        config_.pageSize,
                                        pagePixels.data(),
                                        pagePixels.size(),
                                        outError);
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

    if (!ReadVirtualTextureTile(&texture.tileFile,
                                info.udim,
                                info.mip,
                                info.tileX,
                                info.tileY,
                                &rgba8,
                                &width,
                                &height,
                                &sourceBytes,
                                &readError))
    {
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
    if (streamTextures_.size() < maxStreamTextureCount_)
    {
        if (!AllocateStreamTexture(outError))
        {
            return false;
        }
        if (!freeSlots_.empty())
        {
            const uint32_t slot = freeSlots_.back();
            freeSlots_.pop_back();
            *outSlotIndex = slot;
            *outEvicted = false;
            return true;
        }
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
    struct PendingStreamLoad
    {
        unsigned long long key = 0ull;
        KeyVirtualInfo info = {};
        uint32_t slotIndex = 0u;
    };

    struct CompletedStreamLoad
    {
        unsigned long long key = 0ull;
        KeyVirtualInfo info = {};
        uint32_t slotIndex = 0u;
        bool success = false;
        uint32_t width = 0u;
        uint32_t height = 0u;
        std::vector<unsigned char> rgba8;
        std::string error;
        double loadMs = 0.0;
    };

    const Clock::time_point processStart = Clock::now();
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
    const Clock::time_point histogramStart = Clock::now();
    BuildFeedbackHistogram(keyVec, keyCount, &histogram);
    const Clock::time_point histogramEnd = Clock::now();

    VirtualTextureUpdateStats stats = {};
    stats.feedbackCount = keyCount;
    stats.uniqueCount = static_cast<uint32_t>(histogram.size());
    stats.histogramMs = ElapsedMs(histogramStart, histogramEnd);

    std::string firstError;
    uint32_t uploadsBudget =
        (config_.maxUploadsPerPass == 0u) ? UINT32_MAX : config_.maxUploadsPerPass;
    size_t histogramIndex = 0u;
    while (histogramIndex < histogram.size() && uploadsBudget > 0u)
    {
        const size_t remainingHistogram = histogram.size() - histogramIndex;
        const size_t remainingBudget =
            uploadsBudget == UINT32_MAX ? kParallelLoadBatchSize : size_t(uploadsBudget);
        const size_t batchLimit = std::min(remainingHistogram, std::min(kParallelLoadBatchSize, remainingBudget));
        if (batchLimit == 0u)
        {
            break;
        }

        std::vector<PendingStreamLoad> pending;
        pending.reserve(batchLimit);
        while (histogramIndex < histogram.size() && pending.size() < batchLimit)
        {
            const VirtualTextureFeedbackEntry &entry = histogram[histogramIndex++];
            const unsigned long long key = entry.key;
            KeyVirtualInfo info = {};
            const Clock::time_point resolveStart = Clock::now();
            const bool resolved = ResolveKey(key, &info);
            stats.resolveKeyMs += ElapsedMs(resolveStart, Clock::now());
            if (!resolved || !info.valid)
            {
                stats.failed++;
                continue;
            }
            if (info.isTail)
            {
                stats.hits++;
                continue;
            }
            const Clock::time_point touchStart = Clock::now();
            const bool resident = TouchResidentKey(key);
            stats.touchResidentMs += ElapsedMs(touchStart, Clock::now());
            if (resident)
            {
                stats.hits++;
                continue;
            }

            stats.misses++;

            uint32_t slotIndex = 0u;
            bool evicted = false;
            stats.allocateStreamSlotCalls++;
            const Clock::time_point allocateStart = Clock::now();
            if (!AllocateStreamSlot(&slotIndex, &evicted, outError))
            {
                stats.allocateStreamSlotMs += ElapsedMs(allocateStart, Clock::now());
                stats.failed++;
                RecordFirstError(&firstError, outError ? *outError : std::string());
                continue;
            }
            stats.allocateStreamSlotMs += ElapsedMs(allocateStart, Clock::now());
            if (evicted)
            {
                stats.evictions++;
            }

            TextureState &texture = textures_[info.textureIndex];
            if (!OpenTileFileIfNeeded(&texture, outError))
            {
                stats.failed++;
                freeSlots_.push_back(slotIndex);
                RecordFirstError(&firstError, outError ? *outError : std::string());
                continue;
            }

            PendingStreamLoad load = {};
            load.key = key;
            load.info = info;
            load.slotIndex = slotIndex;
            pending.push_back(std::move(load));
        }

        if (pending.empty())
        {
            continue;
        }

        std::vector<CompletedStreamLoad> completed(pending.size());
        tbb::parallel_for(size_t(0), pending.size(), [&](size_t i) {
            const PendingStreamLoad &load = pending[i];
            CompletedStreamLoad result = {};
            result.key = load.key;
            result.info = load.info;
            result.slotIndex = load.slotIndex;

            const Clock::time_point loadStart = Clock::now();
            result.success = LoadStreamPageForKey(
                load.info, &result.rgba8, &result.width, &result.height, &result.error);
            result.loadMs = ElapsedMs(loadStart, Clock::now());
            completed[i] = std::move(result);
        });

        for (const CompletedStreamLoad &result : completed)
        {
            stats.loadStreamPageCalls++;
            stats.loadStreamPageMs += result.loadMs;
            if (!result.success)
            {
                stats.failed++;
                freeSlots_.push_back(result.slotIndex);
                RecordFirstError(&firstError, result.error);
                continue;
            }

            stats.uploadStreamPageCalls++;
            const Clock::time_point uploadStart = Clock::now();
            if (!UploadStreamPage(result.slotIndex, result.rgba8, result.width, result.height, outError))
            {
                stats.uploadStreamPageMs += ElapsedMs(uploadStart, Clock::now());
                stats.failed++;
                freeSlots_.push_back(result.slotIndex);
                RecordFirstError(&firstError, outError ? *outError : std::string());
                continue;
            }
            stats.uploadStreamPageMs += ElapsedMs(uploadStart, Clock::now());

            const uint32_t page = result.slotIndex % streamPageCountX_;
            const uint32_t physicalTextureID = result.slotIndex / streamPageCountX_;
            const uint32_t packed = PackPageTableEntry(page, physicalTextureID, kPageTypeStream);
            stats.updatePageTableCalls++;
            const Clock::time_point pageTableStart = Clock::now();
            if (!UpdatePageTableTexel(result.info.mip, result.info.vaX, result.info.vaY, packed, outError))
            {
                stats.updatePageTableMs += ElapsedMs(pageTableStart, Clock::now());
                stats.failed++;
                freeSlots_.push_back(result.slotIndex);
                RecordFirstError(&firstError, outError ? *outError : std::string());
                continue;
            }
            stats.updatePageTableMs += ElapsedMs(pageTableStart, Clock::now());

            StreamSlotState &slot = streamSlots_[result.slotIndex];
            slot.occupied = true;
            slot.key = result.key;
            lruSlots_.push_front(result.slotIndex);
            slot.lruIt = lruSlots_.begin();
            keyToStreamSlot_[result.key] = result.slotIndex;

            stats.uploads++;
            if (uploadsBudget != UINT32_MAX)
            {
                uploadsBudget--;
            }
        }
    }

    stats.totalMs = ElapsedMs(processStart, Clock::now());
    if (outStats)
    {
        *outStats = stats;
    }
    if (outError)
    {
        RecordFirstError(outError, firstError);
    }
    return true;
}

} // namespace texture
} // namespace ybi
