#include "shared.h"

#include <libntc/wrappers.h>

#include <algorithm>
#include <array>
#include <cctype>
#include <chrono>
#include <cstdio>
#include <cstdlib>
#include <filesystem>
#include <string>
#include <unordered_set>
#include <vector>

#define STB_IMAGE_IMPLEMENTATION
#include <stb_image.h>

#define TINYEXR_IMPLEMENTATION
#include <tinyexr.h>

namespace
{

struct LoadedChannel
{
    std::string inputName;
    std::string texturePath;
    ntc::ChannelFormat channelFormat = ntc::ChannelFormat::UNORM8;
    bool isSrgb = false;
    int width = 0;
    int height = 0;
    int numChannels = 0;
    uint64_t sourceFileBytes = 0;
    uint64_t decodedBytes = 0;
    std::vector<unsigned char> unorm8;
    std::vector<float> float32;
};

bool IsSrgbInput(const std::string &inputName)
{
    return inputName == "diffuseColor" || inputName == "emissiveColor";
}

int InferChannelCount(const std::string &inputName, const ChannelTexture &texture)
{
    if (!texture.swizzle.empty())
    {
        return static_cast<int>(texture.swizzle.size());
    }
    if (inputName == "roughness" || inputName == "metallic" || inputName == "occlusion" || inputName == "opacity")
    {
        return 1;
    }
    return 3;
}

int SwizzleCharToChannel(char c)
{
    switch (c)
    {
    case 'R':
        return 0;
    case 'G':
        return 1;
    case 'B':
        return 2;
    case 'A':
        return 3;
    default:
        return -1;
    }
}

std::vector<int> BuildSourceChannels(const std::string &inputName, const ChannelTexture &texture, int numChannels)
{
    std::vector<int> channels;
    channels.reserve(numChannels);

    if (!texture.swizzle.empty())
    {
        for (char c : texture.swizzle)
        {
            const int ch = SwizzleCharToChannel(c);
            if (ch >= 0)
            {
                channels.push_back(ch);
            }
        }
    }

    if (channels.empty())
    {
        if (inputName == "roughness" || inputName == "metallic" || inputName == "occlusion")
        {
            channels.push_back(0);
        }
        else if (inputName == "opacity")
        {
            channels.push_back(3);
        }
        else
        {
            channels = {0, 1, 2};
        }
    }

    while (static_cast<int>(channels.size()) < numChannels)
    {
        channels.push_back(channels.back());
    }
    if (static_cast<int>(channels.size()) > numChannels)
    {
        channels.resize(numChannels);
    }

    return channels;
}

std::string LowerExt(const std::string &path)
{
    std::string ext = std::filesystem::path(path).extension().string();
    for (char &c : ext)
    {
        c = static_cast<char>(std::tolower(static_cast<unsigned char>(c)));
    }
    return ext;
}

std::string ResolveUdimTilePath(const std::string &path)
{
    const std::string tokenUpper = "<UDIM>";
    const std::string tokenLower = "<udim>";

    size_t pos = path.find(tokenUpper);
    size_t tokenLen = tokenUpper.size();
    if (pos == std::string::npos)
    {
        pos = path.find(tokenLower);
        tokenLen = tokenLower.size();
    }
    if (pos == std::string::npos)
    {
        return path;
    }

    for (int tile = 1001; tile <= 1199; ++tile)
    {
        std::string candidate = path;
        candidate.replace(pos, tokenLen, std::to_string(tile));
        if (std::filesystem::exists(candidate))
        {
            return candidate;
        }
    }

    std::string fallback = path;
    fallback.replace(pos, tokenLen, "1001");
    return fallback;
}

bool MakeError(ntc::Status status, std::string &out)
{
    out = std::string(ntc::StatusToString(status)) + ": " + ntc::GetLastErrorMessage();
    return false;
}

bool LoadExrRgba(const std::string &path, int &width, int &height, std::vector<float> &rgba, std::string &reason)
{
    float *raw = nullptr;
    const char *err = nullptr;
    width = 0;
    height = 0;

    const int rc = LoadEXR(&raw, &width, &height, path.c_str(), &err);
    if (rc != TINYEXR_SUCCESS || !raw)
    {
        reason = "LoadEXR failed for " + path;
        if (err)
        {
            reason += " (";
            reason += err;
            reason += ")";
            FreeEXRErrorMessage(err);
        }
        return false;
    }

    const size_t pixelCount = static_cast<size_t>(width) * static_cast<size_t>(height);
    rgba.assign(raw, raw + pixelCount * 4);
    std::free(raw);
    return true;
}

bool LoadTexture(const std::string &inputName, const ChannelTexture &texture, LoadedChannel &out, std::string &reason)
{
    out = {};
    out.inputName = inputName;
    out.texturePath = ResolveUdimTilePath(texture.texturePath);

    if (!std::filesystem::exists(out.texturePath))
    {
        reason = "texture file not found: " + out.texturePath;
        return false;
    }
    std::error_code fileSizeEc;
    out.sourceFileBytes = std::filesystem::file_size(out.texturePath, fileSizeEc);
    if (fileSizeEc)
    {
        out.sourceFileBytes = 0;
    }

    const int numChannels = InferChannelCount(inputName, texture);
    const std::vector<int> srcChannels = BuildSourceChannels(inputName, texture, numChannels);

    if (LowerExt(out.texturePath) == ".exr")
    {
        std::vector<float> rgba;
        int width = 0;
        int height = 0;
        if (!LoadExrRgba(out.texturePath, width, height, rgba, reason))
        {
            return false;
        }

        out.channelFormat = ntc::ChannelFormat::FLOAT32;
        out.isSrgb = false;
        out.width = width;
        out.height = height;
        out.numChannels = numChannels;
        out.float32.resize(static_cast<size_t>(width) * static_cast<size_t>(height) * static_cast<size_t>(numChannels));
        out.decodedBytes = static_cast<uint64_t>(out.float32.size()) * sizeof(float);

        const size_t pixelCount = static_cast<size_t>(width) * static_cast<size_t>(height);
        for (size_t i = 0; i < pixelCount; ++i)
        {
            for (int c = 0; c < numChannels; ++c)
            {
                const int src = std::clamp(srcChannels[c], 0, 3);
                out.float32[i * static_cast<size_t>(numChannels) + static_cast<size_t>(c)] = rgba[i * 4 + static_cast<size_t>(src)];
            }
        }

        return true;
    }

    int width = 0;
    int height = 0;
    int channelsInFile = 0;
    unsigned char *pixels = stbi_load(out.texturePath.c_str(), &width, &height, &channelsInFile, 4);
    if (!pixels)
    {
        reason = "stbi_load failed for " + out.texturePath;
        const char *stbReason = stbi_failure_reason();
        if (stbReason)
        {
            reason += " (";
            reason += stbReason;
            reason += ")";
        }
        return false;
    }

    out.channelFormat = ntc::ChannelFormat::UNORM8;
    out.isSrgb = IsSrgbInput(inputName);
    out.width = width;
    out.height = height;
    out.numChannels = numChannels;
    out.unorm8.resize(static_cast<size_t>(width) * static_cast<size_t>(height) * static_cast<size_t>(numChannels));
    out.decodedBytes = static_cast<uint64_t>(out.unorm8.size());

    const size_t pixelCount = static_cast<size_t>(width) * static_cast<size_t>(height);
    for (size_t i = 0; i < pixelCount; ++i)
    {
        for (int c = 0; c < numChannels; ++c)
        {
            const int src = std::clamp(srcChannels[c], 0, 3);
            out.unorm8[i * static_cast<size_t>(numChannels) + static_cast<size_t>(c)] = pixels[i * 4 + static_cast<size_t>(src)];
        }
    }

    stbi_image_free(pixels);
    return true;
}

} // namespace

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
                    std::string &errorOut)
{
    outStats = {};
    const auto t0 = std::chrono::steady_clock::now();
    std::fprintf(stderr,
                 "NTC stage: [%d/%d] load textures: %s\n",
                 materialIndex,
                 materialCount,
                 mat.materialPath.c_str());
    std::fflush(stderr);

    std::vector<std::pair<std::string, ChannelTexture>> sorted(mat.channels.begin(), mat.channels.end());
    std::sort(sorted.begin(), sorted.end(), [](const auto &a, const auto &b) {
        return a.first < b.first;
    });

    std::vector<LoadedChannel> loaded;
    loaded.reserve(sorted.size());

    for (const auto &kv : sorted)
    {
        LoadedChannel chan;
        std::string reason;
        if (!LoadTexture(kv.first, kv.second, chan, reason))
        {
            errorOut = reason;
            return false;
        }
        loaded.push_back(std::move(chan));
    }

    if (loaded.empty())
    {
        errorOut = "no loadable channels";
        return false;
    }

    const int width = loaded[0].width;
    const int height = loaded[0].height;
    int totalChannels = 0;
    std::unordered_set<std::string> uniqueTexturePaths;
    for (const LoadedChannel &ch : loaded)
    {
        if (ch.width != width || ch.height != height)
        {
            errorOut = "texture dimensions mismatch in material";
            return false;
        }
        totalChannels += ch.numChannels;
        outStats.decodedBytes += ch.decodedBytes;
        if (uniqueTexturePaths.insert(ch.texturePath).second)
        {
            outStats.sourceFileBytes += ch.sourceFileBytes;
        }
    }

    if (totalChannels <= 0 || totalChannels > NTC_MAX_CHANNELS)
    {
        errorOut = "total channels out of range (1..16)";
        return false;
    }
    std::fprintf(stderr,
                 "NTC stage: [%d/%d] textures loaded: %zu, dims=%dx%d, channels=%d, source=%.2f MiB, decoded=%.2f MiB\n",
                 materialIndex,
                 materialCount,
                 loaded.size(),
                 width,
                 height,
                 totalChannels,
                 static_cast<double>(outStats.sourceFileBytes) / (1024.0 * 1024.0),
                 static_cast<double>(outStats.decodedBytes) / (1024.0 * 1024.0));
    std::fflush(stderr);

    ntc::TextureSetDesc desc = {};
    desc.width = width;
    desc.height = height;
    desc.channels = totalChannels;
    desc.mips = 1;

    ntc::TextureSetFeatures features = {};
    features.stagingBytesPerPixel = 16; // allow up to float4 transfers through host staging

    ntc::TextureSetWrapper textureSet(context);
    std::fprintf(stderr, "NTC stage: [%d/%d] setup texture set\n", materialIndex, materialCount);
    std::fflush(stderr);
    ntc::Status status = context->CreateTextureSet(desc, features, textureSet.ptr());
    if (status != ntc::Status::Ok)
    {
        return MakeError(status, errorOut);
    }

    ntc::LatentShape latentShape = {};
    status = ntc::PickLatentShape(bitsPerPixel, outActualBpp, latentShape);
    if (status != ntc::Status::Ok)
    {
        return MakeError(status, errorOut);
    }

    status = textureSet->SetLatentShape(latentShape);
    if (status != ntc::Status::Ok)
    {
        return MakeError(status, errorOut);
    }

    int firstChannel = 0;
    std::fprintf(stderr, "NTC stage: [%d/%d] upload channels\n", materialIndex, materialCount);
    std::fflush(stderr);
    for (const LoadedChannel &ch : loaded)
    {
        ntc::ITextureMetadata *meta = textureSet->AddTexture();
        meta->SetName(ch.inputName.c_str());

        status = meta->SetChannels(firstChannel, ch.numChannels);
        if (status != ntc::Status::Ok)
        {
            return MakeError(status, errorOut);
        }

        meta->SetChannelFormat(ch.channelFormat);
        meta->SetBlockCompressedFormat(ntc::BlockCompressedFormat::None);
        meta->SetRgbColorSpace(ch.isSrgb ? ntc::ColorSpace::sRGB : ntc::ColorSpace::Linear);
        meta->SetAlphaColorSpace(ntc::ColorSpace::Linear);

        std::array<ntc::ColorSpace, 4> colorSpaces = {
            ch.isSrgb ? ntc::ColorSpace::sRGB : ntc::ColorSpace::Linear,
            ch.isSrgb ? ntc::ColorSpace::sRGB : ntc::ColorSpace::Linear,
            ch.isSrgb ? ntc::ColorSpace::sRGB : ntc::ColorSpace::Linear,
            ntc::ColorSpace::Linear,
        };

        ntc::WriteChannelsParameters params = {};
        params.mipLevel = 0;
        params.firstChannel = firstChannel;
        params.numChannels = ch.numChannels;
        params.addressSpace = ntc::AddressSpace::Host;
        params.width = ch.width;
        params.height = ch.height;
        params.channelFormat = ch.channelFormat;
        params.srcColorSpaces = colorSpaces.data();
        params.dstColorSpaces = colorSpaces.data();

        if (ch.channelFormat == ntc::ChannelFormat::FLOAT32)
        {
            params.pData = reinterpret_cast<const unsigned char *>(ch.float32.data());
            params.pixelStride = sizeof(float) * static_cast<size_t>(ch.numChannels);
            params.rowPitch = sizeof(float) * static_cast<size_t>(ch.width) * static_cast<size_t>(ch.numChannels);
        }
        else
        {
            params.pData = ch.unorm8.data();
            params.pixelStride = static_cast<size_t>(ch.numChannels);
            params.rowPitch = static_cast<size_t>(ch.width) * static_cast<size_t>(ch.numChannels);
        }

        status = textureSet->WriteChannels(params);
        if (status != ntc::Status::Ok)
        {
            return MakeError(status, errorOut);
        }

        firstChannel += ch.numChannels;
    }

    ntc::CompressionSettings settings = {};
    settings.trainingSteps = trainingSteps;
    settings.stepsPerIteration = std::max(1, std::min(stepsPerIter, trainingSteps));

    status = textureSet->BeginCompression(settings);
    if (status != ntc::Status::Ok)
    {
        return MakeError(status, errorOut);
    }
    std::fprintf(stderr,
                 "NTC stage: [%d/%d] start compression (steps=%d, stepPerIter=%d)\n",
                 materialIndex,
                 materialCount,
                 settings.trainingSteps,
                 settings.stepsPerIteration);
    std::fflush(stderr);

    ntc::CompressionStats stats = {};
    int nextPercentPrint = 0;
    do
    {
        status = textureSet->RunCompressionSteps(&stats);

        if (settings.trainingSteps > 0)
        {
            const int done = std::max(0, std::min(stats.currentStep, settings.trainingSteps));
            const int pct = (done * 100) / settings.trainingSteps;
            if (pct >= nextPercentPrint || status == ntc::Status::Ok)
            {
                const int barWidth = 24;
                const int filled = (done * barWidth) / settings.trainingSteps;
                std::fprintf(stderr, "  [");
                for (int i = 0; i < barWidth; ++i)
                {
                    std::fprintf(stderr, "%c", i < filled ? '=' : ' ');
                }
                std::fprintf(stderr,
                             "] set %d/%d (%d remaining) step %d/%d\n",
                             materialIndex,
                             materialCount,
                             std::max(0, materialCount - materialIndex),
                             done,
                             settings.trainingSteps);
                std::fflush(stderr);
                nextPercentPrint += 10;
            }
        }
    } while (status == ntc::Status::Incomplete);

    if (status != ntc::Status::Ok)
    {
        return MakeError(status, errorOut);
    }

    status = textureSet->FinalizeCompression();
    if (status != ntc::Status::Ok)
    {
        return MakeError(status, errorOut);
    }
    std::fprintf(stderr, "NTC stage: [%d/%d] finalize compression\n", materialIndex, materialCount);
    std::fflush(stderr);

    status = textureSet->SaveToFile(outPath.string().c_str());
    if (status != ntc::Status::Ok)
    {
        return MakeError(status, errorOut);
    }
    std::error_code ntcSizeEc;
    outStats.ntcBytes = fs::file_size(outPath, ntcSizeEc);
    if (ntcSizeEc)
    {
        outStats.ntcBytes = 0;
    }
    const auto t1 = std::chrono::steady_clock::now();
    const auto elapsedMs = std::chrono::duration_cast<std::chrono::milliseconds>(t1 - t0).count();
    std::fprintf(stderr,
                 "NTC stage: [%d/%d] done %s in %lld ms, ntc=%.2f MiB\n",
                 materialIndex,
                 materialCount,
                 outPath.string().c_str(),
                 static_cast<long long>(elapsedMs),
                 static_cast<double>(outStats.ntcBytes) / (1024.0 * 1024.0));
    std::fflush(stderr);

    return true;
}
