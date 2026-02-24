#pragma once

#include <pxr/usd/usd/stage.h>

#include <libntc/ntc.h>

#include <cstdint>
#include <filesystem>
#include <string>
#include <unordered_map>
#include <vector>

namespace fs = std::filesystem;

struct Cli
{
    std::string usdPath;
    std::string outDir;
    std::string outUsdPath;
    float bitsPerPixel = 4.0f;
    int trainingSteps = 1000;
    int stepsPerIter = 100;
    int cudaDevice = 0;
    int maxMaterials = 0;
    std::vector<std::string> purposes = {"default", "render"};
    bool noEncode = false;
    bool prepareTiles = false;
    int tileSize = 128;
    int tileVerifyCount = 4;
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

struct EncodeStats
{
    uint64_t sourceFileBytes = 0;
    uint64_t decodedBytes = 0;
    uint64_t ntcBytes = 0;
};

void PrintUsage(const char *exe);
bool ParseCli(int argc, char **argv, Cli &out);

std::vector<MaterialChannels> CollectMaterialChannels(const pxr::UsdStageRefPtr &stage,
                                                      const std::vector<std::string> &purposes);
bool WriteNtcBindingsToUsd(const pxr::UsdStageRefPtr &stage,
                           const std::unordered_map<std::string, std::string> &materialToNtcFile,
                           const std::string &outUsdPath,
                           std::string *outError);
std::string Sanitize(const std::string &s);

bool EncodeMaterial(ntc::IContext *context,
                    const MaterialChannels &mat,
                    const fs::path &outPath,
                    float bitsPerPixel,
                    int trainingSteps,
                    int stepsPerIter,
                    int materialIndex,
                    int materialCount,
                    EncodeStats &outStats,
                    float &outActualBpp,
                    std::string &errorOut);

bool PrepareTexturesForStreamingTiles(const std::vector<MaterialChannels> &materials,
                                      const Cli &cli,
                                      std::string *outError);
