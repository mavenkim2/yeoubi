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

} // namespace texture
} // namespace ybi
