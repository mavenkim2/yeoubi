#include "optix_ntc_runtime.h"

#include <libntc/ntc.h>
#include <libntc/wrappers.h>

#include <algorithm>
#include <cstdio>
#include <filesystem>
#include <string>
#include <unordered_map>

namespace fs = std::filesystem;

namespace testbvh
{
namespace
{

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

bool CheckCuda(cudaError_t result, std::string &outError, const char *callName)
{
    if (result == cudaSuccess)
    {
        return true;
    }
    outError = std::string(callName) + " failed: " + cudaGetErrorString(result);
    return false;
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

std::string FindDiffuseTexturePath(const ybi::ScenePool::MaterialInfo &material)
{
    for (const ybi::ScenePool::MaterialTextureInput &input : material.textureInputs)
    {
        if (input.inputName == "diffuseColor")
        {
            return input.texturePath;
        }
    }
    return "";
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
    std::unordered_map<std::string, size_t> decodedByDiffusePath;
    int withNtcPath = 0;
    for (size_t materialIndex = 0; materialIndex < materials.size(); ++materialIndex)
    {
        const ybi::ScenePool::MaterialInfo &material = materials[materialIndex];
        if (material.ntcDiffuseFile.empty())
        {
            continue;
        }
        withNtcPath++;
        const fs::path ntcPath(material.ntcDiffuseFile);
        if (!fs::exists(ntcPath))
        {
            std::printf("NTC runtime: material %zu missing file from USD (%s)\n",
                        materialIndex,
                        ntcPath.string().c_str());
            continue;
        }

        std::string decodeError;
        if (DecodeOneMaterial(context, ntcPath, (*outTextures)[materialIndex], decodeError))
        {
            loadedCount++;
            const std::string diffusePath = FindDiffuseTexturePath(material);
            if (!diffusePath.empty())
            {
                decodedByDiffusePath.emplace(diffusePath, materialIndex);
            }
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

    for (size_t materialIndex = 0; materialIndex < materials.size(); ++materialIndex)
    {
        if ((*outTextures)[materialIndex].valid)
        {
            continue;
        }
        const std::string diffusePath = FindDiffuseTexturePath(materials[materialIndex]);
        if (diffusePath.empty())
        {
            continue;
        }
        auto found = decodedByDiffusePath.find(diffusePath);
        if (found == decodedByDiffusePath.end())
        {
            continue;
        }
        const size_t sourceIndex = found->second;
        if (sourceIndex >= outTextures->size() || !(*outTextures)[sourceIndex].valid)
        {
            continue;
        }
        (*outTextures)[materialIndex] = (*outTextures)[sourceIndex];
        loadedCount++;
        std::printf("NTC runtime: material %zu aliased to material %zu (diffuse path match)\n",
                    materialIndex,
                    sourceIndex);
    }

    std::printf("NTC runtime: decoded %d/%zu material textures from %s\n",
                loadedCount,
                materials.size(),
                "USD material attributes");
    std::printf("NTC runtime: materials with ntc path: %d\n", withNtcPath);
    return true;
}

bool UploadDecodedTexturesToCuda(const std::vector<DecodedMaterialTexture> &decodedTextures,
                                 UploadedMaterialTextures *outTextures,
                                 std::string *outError)
{
    if (!outTextures)
    {
        if (outError)
        {
            *outError = "UploadDecodedTexturesToCuda: null outTextures";
        }
        return false;
    }

    DestroyUploadedTextures(outTextures);
    outTextures->refs.resize(decodedTextures.size());

    for (size_t i = 0; i < decodedTextures.size(); ++i)
    {
        const DecodedMaterialTexture &src = decodedTextures[i];
        if (!src.valid || src.width <= 0 || src.height <= 0 || src.rgba8.empty())
        {
            continue;
        }

        cudaChannelFormatDesc channelDesc =
            cudaCreateChannelDesc(8, 8, 8, 8, cudaChannelFormatKindUnsigned);
        cudaArray_t array = nullptr;
        std::string error;
        if (!CheckCuda(cudaMallocArray(&array, &channelDesc, src.width, src.height), error, "cudaMallocArray"))
        {
            if (outError)
            {
                *outError = error;
            }
            DestroyUploadedTextures(outTextures);
            return false;
        }

        if (!CheckCuda(cudaMemcpy2DToArray(array,
                                           0,
                                           0,
                                           src.rgba8.data(),
                                           static_cast<size_t>(src.width) * 4u,
                                           static_cast<size_t>(src.width) * 4u,
                                           src.height,
                                           cudaMemcpyHostToDevice),
                       error,
                       "cudaMemcpy2DToArray"))
        {
            cudaFreeArray(array);
            if (outError)
            {
                *outError = error;
            }
            DestroyUploadedTextures(outTextures);
            return false;
        }

        cudaResourceDesc resourceDesc = {};
        resourceDesc.resType = cudaResourceTypeArray;
        resourceDesc.res.array.array = array;

        cudaTextureDesc textureDesc = {};
        textureDesc.addressMode[0] = cudaAddressModeWrap;
        textureDesc.addressMode[1] = cudaAddressModeWrap;
        textureDesc.filterMode = cudaFilterModePoint;
        textureDesc.readMode = cudaReadModeNormalizedFloat;
        textureDesc.normalizedCoords = 1;

        cudaTextureObject_t textureObject = 0;
        if (!CheckCuda(cudaCreateTextureObject(&textureObject, &resourceDesc, &textureDesc, nullptr),
                       error,
                       "cudaCreateTextureObject"))
        {
            cudaFreeArray(array);
            if (outError)
            {
                *outError = error;
            }
            DestroyUploadedTextures(outTextures);
            return false;
        }

        outTextures->arrays.push_back(array);
        outTextures->textureObjects.push_back(textureObject);
        MaterialTextureRef &dstRef = outTextures->refs[i];
        dstRef.textureObject = static_cast<unsigned long long>(textureObject);
        dstRef.width = src.width;
        dstRef.height = src.height;
        dstRef.valid = 1;
    }

    return true;
}

void DestroyUploadedTextures(UploadedMaterialTextures *textures)
{
    if (!textures)
    {
        return;
    }
    for (cudaTextureObject_t textureObject : textures->textureObjects)
    {
        if (textureObject != 0)
        {
            cudaDestroyTextureObject(textureObject);
        }
    }
    for (cudaArray_t array : textures->arrays)
    {
        if (array)
        {
            cudaFreeArray(array);
        }
    }
    textures->textureObjects.clear();
    textures->arrays.clear();
    textures->refs.clear();
}

} // namespace testbvh
