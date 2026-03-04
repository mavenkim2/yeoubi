#pragma once

#include "render/launch_params.h"

#include <cstdint>
#include <filesystem>
#include <string>
#include <unordered_map>
#include <vector>

namespace ybi
{
namespace texture
{

constexpr unsigned long long kVirtualTextureEmptyKey = ~0ull;

struct VirtualTextureFeedbackEntry
{
    unsigned long long key = 0ull;
    unsigned int count = 0u;
};

struct HostVirtualTextureTile
{
    unsigned long long key = 0ull;
    uint32_t width = 0u;
    uint32_t height = 0u;
    uint64_t sourceBytes = 0ull;
    std::vector<unsigned char> rgba8;
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

int BuildVirtualTextureHashTable(
    const std::unordered_map<unsigned long long, HostVirtualTextureTile> &tiles,
    std::vector<LaunchParams::VirtualTextureTileEntry> *outEntries,
    std::vector<unsigned char> *outPixels,
    std::string *outError);

} // namespace texture
} // namespace ybi
