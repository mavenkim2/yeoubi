#include "optix_ntc_runtime.h"

#include <libntc/ntc.h>
#include <libntc/wrappers.h>

#include <algorithm>
#include <cctype>
#include <cstdio>
#include <filesystem>
#include <string>

namespace fs = std::filesystem;

namespace testbvh
{
namespace
{

std::string Sanitize(const std::string &s)
{
    std::string out = s;
    for (char &c : out)
    {
        if (!(std::isalnum((unsigned char)c) || c == '_' || c == '-' || c == '.'))
        {
            c = '_';
        }
    }
    if (out.empty())
    {
        return "material";
    }
    return out;
}

std::string NtcStatusString(ntc::Status status)
{
    std::string out = ntc::StatusToString(status);
    const char *last = ntc::GetLastErrorMessage();
    if (last && last[0] != '\0')
    {
        out += ": ";
        out += last;
    }
    return out;
}

int FindTextureByName(ntc::ITextureSet *textureSet, const char *name)
{
    const int textureCount = textureSet->GetTextureCount();
    for (int i = 0; i < textureCount; ++i)
    {
        const ntc::ITextureMetadata *meta = textureSet->GetTexture(i);
        if (!meta)
        {
            continue;
        }
        const char *textureName = meta->GetName();
        if (textureName && std::string(textureName) == name)
        {
            return i;
        }
    }
    return -1;
}

int FindBestAlbedoTexture(ntc::ITextureSet *textureSet)
{
    const int diffuseIndex = FindTextureByName(textureSet, "diffuseColor");
    if (diffuseIndex >= 0)
    {
        return diffuseIndex;
    }

    int fallback = -1;
    const int textureCount = textureSet->GetTextureCount();
    for (int i = 0; i < textureCount; ++i)
    {
        ntc::ITextureMetadata *meta = textureSet->GetTexture(i);
        if (!meta)
        {
            continue;
        }
        int firstChannel = 0;
        int numChannels = 0;
        meta->GetChannels(firstChannel, numChannels);
        if (numChannels >= 3)
        {
            return i;
        }
        if (fallback < 0 && numChannels > 0)
        {
            fallback = i;
        }
    }
    return fallback;
}

void ExpandToRgba(const std::vector<unsigned char> &src,
                  int width,
                  int height,
                  int srcChannels,
                  std::vector<unsigned char> *dst)
{
    dst->assign(static_cast<size_t>(width) * static_cast<size_t>(height) * 4u, 255u);
    if (srcChannels <= 0)
    {
        return;
    }

    const size_t pixelCount = static_cast<size_t>(width) * static_cast<size_t>(height);
    for (size_t i = 0; i < pixelCount; ++i)
    {
        const unsigned char *s = src.data() + i * static_cast<size_t>(srcChannels);
        unsigned char *d = dst->data() + i * 4u;
        if (srcChannels == 1)
        {
            d[0] = s[0];
            d[1] = s[0];
            d[2] = s[0];
            d[3] = 255u;
            continue;
        }
        if (srcChannels == 2)
        {
            d[0] = s[0];
            d[1] = s[1];
            d[2] = s[0];
            d[3] = 255u;
            continue;
        }

        d[0] = s[0];
        d[1] = s[1];
        d[2] = s[2];
        d[3] = (srcChannels >= 4) ? s[3] : 255u;
    }
}

bool DecodeOneMaterial(ntc::IContext *context,
                       const fs::path &ntcPath,
                       DecodedMaterialTexture &outTexture,
                       std::string &outError)
{
    ntc::TextureSetFeatures features = {};
    features.enableCompression = false;
    features.stagingBytesPerPixel = 4;
    features.separateRefOutData = false;

    ntc::TextureSetWrapper textureSet(context);
    ntc::Status status = context->CreateCompressedTextureSetFromFile(
        ntcPath.string().c_str(), features, textureSet.ptr());
    if (status != ntc::Status::Ok)
    {
        outError = "CreateCompressedTextureSetFromFile failed: " + NtcStatusString(status);
        return false;
    }

    ntc::DecompressionStats decompressStats = {};
    status = textureSet->Decompress(&decompressStats, false);
    if (status != ntc::Status::Ok)
    {
        outError = "Decompress failed: " + NtcStatusString(status);
        return false;
    }

    const int textureIndex = FindBestAlbedoTexture(textureSet);
    if (textureIndex < 0)
    {
        outError = "No texture found in set";
        return false;
    }

    ntc::ITextureMetadata *meta = textureSet->GetTexture(textureIndex);
    if (!meta)
    {
        outError = "Texture metadata lookup failed";
        return false;
    }

    int firstChannel = 0;
    int numChannels = 0;
    meta->GetChannels(firstChannel, numChannels);
    if (numChannels <= 0)
    {
        outError = "Texture has no channels";
        return false;
    }

    const ntc::TextureSetDesc desc = textureSet->GetDesc();
    const int readChannels = std::max(1, std::min(numChannels, 4));
    std::vector<unsigned char> packed(static_cast<size_t>(desc.width) * static_cast<size_t>(desc.height) *
                                      static_cast<size_t>(readChannels));

    ntc::ReadChannelsParameters readParams = {};
    readParams.page = ntc::TextureDataPage::Output;
    readParams.mipLevel = 0;
    readParams.firstChannel = firstChannel;
    readParams.numChannels = readChannels;
    readParams.pOutData = packed.data();
    readParams.addressSpace = ntc::AddressSpace::Host;
    readParams.width = desc.width;
    readParams.height = desc.height;
    readParams.pixelStride = static_cast<size_t>(readChannels);
    readParams.rowPitch = static_cast<size_t>(desc.width) * static_cast<size_t>(readChannels);
    readParams.channelFormat = ntc::ChannelFormat::UNORM8;
    readParams.dstColorSpaces = nullptr;
    readParams.useDithering = false;

    status = textureSet->ReadChannels(readParams);
    if (status != ntc::Status::Ok)
    {
        outError = "ReadChannels failed: " + NtcStatusString(status);
        return false;
    }

    outTexture.valid = true;
    outTexture.width = desc.width;
    outTexture.height = desc.height;
    outTexture.ntcPath = ntcPath.string();
    const char *textureName = meta->GetName();
    outTexture.textureName = textureName ? textureName : "";
    ExpandToRgba(packed, desc.width, desc.height, readChannels, &outTexture.rgba8);
    return true;
}

} // namespace

bool DecodeNtcDiffuseTextures(const std::vector<ybi::ScenePool::MaterialInfo> &materials,
                              const std::string &ntcDir,
                              std::vector<DecodedMaterialTexture> *outTextures,
                              std::string *outError)
{
    if (!outTextures)
    {
        if (outError)
        {
            *outError = "DecodeNtcDiffuseTextures: null outTextures";
        }
        return false;
    }

    outTextures->clear();
    outTextures->resize(materials.size());
    if (ntcDir.empty())
    {
        return true;
    }

    ntc::ContextParameters contextParams = {};
    contextParams.cudaDevice = 0;
    contextParams.graphicsApi = ntc::GraphicsAPI::None;

    ntc::ContextWrapper context;
    ntc::Status status = ntc::CreateContext(context.ptr(), contextParams);
    if (status != ntc::Status::Ok)
    {
        if (outError)
        {
            *outError = "CreateContext failed: " + NtcStatusString(status);
        }
        return false;
    }

    int loadedCount = 0;
    for (size_t materialIndex = 0; materialIndex < materials.size(); ++materialIndex)
    {
        const ybi::ScenePool::MaterialInfo &material = materials[materialIndex];
        const fs::path ntcPath = fs::path(ntcDir) / (Sanitize(material.materialPath) + ".ntc");
        if (!fs::exists(ntcPath))
        {
            continue;
        }

        std::string decodeError;
        if (DecodeOneMaterial(context, ntcPath, (*outTextures)[materialIndex], decodeError))
        {
            loadedCount++;
            std::printf("NTC runtime: material %zu -> %s (%dx%d, texture=%s)\n",
                        materialIndex,
                        ntcPath.string().c_str(),
                        (*outTextures)[materialIndex].width,
                        (*outTextures)[materialIndex].height,
                        (*outTextures)[materialIndex].textureName.c_str());
            continue;
        }

        std::printf("NTC runtime: material %zu decode failed (%s): %s\n",
                    materialIndex,
                    ntcPath.string().c_str(),
                    decodeError.c_str());
    }

    std::printf("NTC runtime: decoded %d/%zu material textures from %s\n",
                loadedCount,
                materials.size(),
                ntcDir.c_str());
    return true;
}

} // namespace testbvh
