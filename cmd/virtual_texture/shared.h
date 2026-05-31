#pragma once

#include <pxr/usd/usd/stage.h>

#include <filesystem>
#include <string>
#include <unordered_map>
#include <vector>

namespace fs = std::filesystem;

struct Cli
{
    std::string usdPath;
    std::string outDir;
    int maxMaterials = 0;
    std::vector<std::string> purposes = {"default", "render"};
    int tileSize = 128;
    int tileVerifyCount = 0;
    bool tileVerifyPass = true;
    float tileVerifyEps = 0.0f;
};

struct ChannelTexture
{
    std::string texturePath;
    std::string swizzle;
};

struct MaterialChannels
{
    std::string materialPath;
    std::unordered_map<std::string, ChannelTexture> channels;
};

void PrintUsage(const char *exe);
bool ParseCli(int argc, char **argv, Cli &out);

std::vector<MaterialChannels> CollectMaterialChannels(const pxr::UsdStageRefPtr &stage,
                                                      const std::vector<std::string> &purposes);
std::string Sanitize(const std::string &s);

bool PrepareTexturesForStreamingTiles(const std::vector<MaterialChannels> &materials,
                                      const Cli &cli,
                                      std::string *outError);
