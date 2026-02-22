#pragma once

#include <pxr/usd/usd/stage.h>

#include <libntc/ntc.h>

#include <filesystem>
#include <string>
#include <unordered_map>
#include <vector>

namespace fs = std::filesystem;

struct Cli
{
    std::string usdPath;
    std::string outDir;
    float bitsPerPixel = 4.0f;
    int trainingSteps = 1000;
    int stepsPerIter = 100;
    int cudaDevice = 0;
    int maxMaterials = 0;
    bool noEncode = false;
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

std::vector<MaterialChannels> CollectMaterialChannels(const pxr::UsdStageRefPtr &stage);
std::string Sanitize(const std::string &s);
bool WriteManifest(const fs::path &path, const MaterialChannels &mat);

bool EncodeMaterial(ntc::IContext *context,
                    const MaterialChannels &mat,
                    const fs::path &outPath,
                    float bitsPerPixel,
                    int trainingSteps,
                    int stepsPerIter,
                    float &outActualBpp,
                    std::string &errorOut);
