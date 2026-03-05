#pragma once

#include <cstdint>
#include <filesystem>
#include <string>
#include <vector>

namespace ybi
{
namespace texture
{

struct VirtualTextureFeedbackEntry
{
    unsigned long long key = 0ull;
    unsigned int count = 0u;
};

unsigned int FeedbackTileX(unsigned long long key);
unsigned int FeedbackTileY(unsigned long long key);
unsigned int FeedbackUdim(unsigned long long key);
unsigned int FeedbackTextureId(unsigned long long key);
unsigned int FeedbackMip(unsigned long long key);

void BuildFeedbackHistogram(const std::vector<unsigned long long> &keys,
                            unsigned int count,
                            std::vector<VirtualTextureFeedbackEntry> *out);

bool WriteFeedbackFile(const std::filesystem::path &feedbackPath,
                       int sppIndex,
                       unsigned int sampledCount,
                       unsigned int copyCount,
                       unsigned int overflowCount,
                       const std::vector<unsigned long long> &keys,
                       std::vector<VirtualTextureFeedbackEntry> *outHistogram,
                       std::string *outError);

} // namespace texture
} // namespace ybi
