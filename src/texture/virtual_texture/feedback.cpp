#include "texture/virtual_texture/feedback.h"

#include "texture/virtual_texture/key.h"

#include <algorithm>
#include <cstdio>

namespace ybi
{
namespace texture
{

unsigned int FeedbackTileX(unsigned long long key)
{
    return ybi::texture::UnpackVirtualTextureKeyTileX(key);
}

unsigned int FeedbackTileY(unsigned long long key)
{
    return ybi::texture::UnpackVirtualTextureKeyTileY(key);
}

unsigned int FeedbackUdim(unsigned long long key)
{
    return ybi::texture::UnpackVirtualTextureKeyUdim(key);
}

unsigned int FeedbackTextureId(unsigned long long key)
{
    return ybi::texture::UnpackVirtualTextureKeyTextureId(key);
}

unsigned int FeedbackMip(unsigned long long key)
{
    return ybi::texture::UnpackVirtualTextureKeyMip(key);
}

void BuildFeedbackHistogram(const std::vector<unsigned long long> &keys,
                            unsigned int count,
                            std::vector<VirtualTextureFeedbackEntry> *out)
{
    if (!out)
    {
        return;
    }
    out->clear();
    if (keys.empty() || count == 0)
    {
        return;
    }

    const unsigned int clampedCount = std::min<unsigned int>(count, static_cast<unsigned int>(keys.size()));
    std::vector<unsigned long long> sorted(keys.begin(), keys.begin() + clampedCount);
    std::sort(sorted.begin(), sorted.end());

    unsigned int runCount = 1;
    for (size_t i = 1; i < sorted.size(); ++i)
    {
        if (sorted[i] == sorted[i - 1])
        {
            ++runCount;
            continue;
        }
        out->push_back({sorted[i - 1], runCount});
        runCount = 1;
    }
    out->push_back({sorted.back(), runCount});
}

bool WriteFeedbackFile(const std::filesystem::path &feedbackPath,
                       int sppIndex,
                       unsigned int sampledCount,
                       unsigned int copyCount,
                       unsigned int overflowCount,
                       const std::vector<unsigned long long> &keys,
                       std::vector<VirtualTextureFeedbackEntry> *outHistogram,
                       std::string *outError)
{
    if (outHistogram)
    {
        BuildFeedbackHistogram(keys, copyCount, outHistogram);
    }

    std::FILE *feedbackFile = std::fopen(feedbackPath.string().c_str(), "w");
    if (!feedbackFile)
    {
        if (outError)
        {
            *outError = "failed to open feedback file for write: " + feedbackPath.string();
        }
        return false;
    }

    const size_t uniqueCount = outHistogram ? outHistogram->size() : 0u;
    std::fprintf(feedbackFile,
                 "spp=%d sampled=%u stored=%u overflow=%u unique=%zu\n",
                 sppIndex,
                 sampledCount,
                 copyCount,
                 overflowCount,
                 uniqueCount);
    std::fprintf(feedbackFile, "textureId udim tileX tileY mip count\n");

    if (outHistogram)
    {
        for (const VirtualTextureFeedbackEntry &entry : *outHistogram)
        {
            std::fprintf(feedbackFile,
                         "%u %u %u %u %u %u\n",
                         FeedbackTextureId(entry.key),
                         FeedbackUdim(entry.key),
                         FeedbackTileX(entry.key),
                         FeedbackTileY(entry.key),
                         FeedbackMip(entry.key),
                         entry.count);
        }
    }

    std::fclose(feedbackFile);
    return true;
}

int BuildVirtualTextureHashTable(
    const std::unordered_map<unsigned long long, HostVirtualTextureTile> &tiles,
    std::vector<LaunchParams::VirtualTextureTileEntry> *outEntries,
    std::vector<unsigned char> *outPixels,
    std::string *outError)
{
    if (!outEntries || !outPixels)
    {
        if (outError)
        {
            *outError = "null output buffers for virtual texture hash table";
        }
        return 0;
    }

    outEntries->clear();
    outPixels->clear();
    if (tiles.empty())
    {
        return 0;
    }

    int capacity = 1;
    const int minCapacity = static_cast<int>(tiles.size()) * 2;
    while (capacity < minCapacity)
    {
        capacity <<= 1;
    }

    outEntries->resize(static_cast<size_t>(capacity));
    for (LaunchParams::VirtualTextureTileEntry &entry : *outEntries)
    {
        entry.key = kVirtualTextureEmptyKey;
        entry.pixelOffset = 0ull;
        entry.width = 0u;
        entry.height = 0u;
    }

    size_t pixelBytes = 0;
    for (const auto &it : tiles)
    {
        pixelBytes += it.second.rgba8.size();
    }
    outPixels->reserve(pixelBytes);

    const unsigned int mask = static_cast<unsigned int>(capacity - 1);
    for (const auto &it : tiles)
    {
        const HostVirtualTextureTile &tile = it.second;
        const unsigned long long key = tile.key;
        if (tile.rgba8.empty() || tile.width == 0u || tile.height == 0u)
        {
            continue;
        }

        const unsigned long long offset = static_cast<unsigned long long>(outPixels->size());
        outPixels->insert(outPixels->end(), tile.rgba8.begin(), tile.rgba8.end());

        unsigned int slot = ybi::texture::HashVirtualTextureKey(key) & mask;
        bool inserted = false;
        for (int probe = 0; probe < capacity; ++probe)
        {
            LaunchParams::VirtualTextureTileEntry &entry =
                (*outEntries)[static_cast<size_t>(slot)];
            if (entry.key == kVirtualTextureEmptyKey)
            {
                entry.key = key;
                entry.pixelOffset = offset;
                entry.width = tile.width;
                entry.height = tile.height;
                inserted = true;
                break;
            }
            slot = (slot + 1u) & mask;
        }

        if (!inserted)
        {
            outEntries->clear();
            outPixels->clear();
            if (outError)
            {
                *outError = "failed to insert virtual texture tile into hash table";
            }
            return 0;
        }
    }

    return capacity;
}

} // namespace texture
} // namespace ybi
