#include "texture/virtual_texture/material_metadata.h"

#include "texture/virtual_texture/tile_file.h"

#include <algorithm>
#include <sstream>

namespace ybi
{
namespace texture
{
namespace
{

static constexpr uint32_t kUdimMin = 1001u;

} // namespace

bool BuildVirtualTextureMaterialMetadata(
    size_t materialCount,
    int semanticIndex,
    int materialTextureRefSemanticCount,
    int materialTextureRefStride,
    const std::vector<VirtualTextureMaterialSource> &sources,
    const std::unordered_map<unsigned int, std::string> &tileFiles,
    VirtualTextureMaterialBuildResult *out,
    std::string *outError)
{
    if (!out)
    {
        if (outError)
        {
            *outError = "BuildVirtualTextureMaterialMetadata: output is null";
        }
        return false;
    }

    *out = {};
    if (materialCount == 0u || materialTextureRefSemanticCount <= 0 || materialTextureRefStride <= 0)
    {
        return true;
    }
    if (semanticIndex < 0 || semanticIndex >= materialTextureRefSemanticCount)
    {
        if (outError)
        {
            *outError = "BuildVirtualTextureMaterialMetadata: semantic index out of range";
        }
        return false;
    }

    const size_t refCount = materialCount * static_cast<size_t>(materialTextureRefSemanticCount) *
                            static_cast<size_t>(materialTextureRefStride);
    out->materialTextureRefs.assign(refCount, LaunchParams::MaterialTextureRef{});
    out->registrations.reserve(sources.size());

    for (const VirtualTextureMaterialSource &source : sources)
    {
        if (source.materialIndex >= materialCount || source.texturePath.empty())
        {
            continue;
        }

        const auto tileFileIt = tileFiles.find(source.materialIndex);
        if (tileFileIt == tileFiles.end() || tileFileIt->second.empty())
        {
            if (outError)
            {
                std::ostringstream oss;
                oss << "virtual-texture metadata: missing tile bin for material "
                    << source.materialIndex << " path " << source.texturePath;
                *outError = oss.str();
            }
            return false;
        }

        VirtualTextureTileFile tileFile = {};
        std::string openError;
        if (!OpenVirtualTextureTileFile(tileFileIt->second, &tileFile, &openError))
        {
            if (outError)
            {
                std::ostringstream oss;
                oss << "virtual-texture metadata: failed opening tile bin for material "
                    << source.materialIndex << ": " << openError;
                *outError = oss.str();
            }
            return false;
        }
        if (tileFile.udims.empty())
        {
            if (outError)
            {
                std::ostringstream oss;
                oss << "virtual-texture metadata: tile bin has no UDIMs for material "
                    << source.materialIndex << ": " << tileFileIt->second;
                *outError = oss.str();
            }
            return false;
        }

        std::vector<uint32_t> udims;
        udims.reserve(tileFile.udims.size());
        for (const auto &entry : tileFile.udims)
        {
            udims.push_back(entry.first);
        }
        std::sort(udims.begin(), udims.end());

        VirtualTextureRegisterInput reg = {};
        reg.textureId = source.materialIndex;
        reg.tileFilePath = tileFileIt->second;
        reg.activeUdims.reserve(udims.size());
        reg.udimExtents.reserve(udims.size());

        const size_t base = (static_cast<size_t>(source.materialIndex) *
                             static_cast<size_t>(materialTextureRefSemanticCount) +
                             static_cast<size_t>(semanticIndex)) *
                            static_cast<size_t>(materialTextureRefStride);

        for (uint32_t udim : udims)
        {
            if (udim < kUdimMin)
            {
                continue;
            }
            const uint32_t udimSlot = udim - kUdimMin;
            if (udimSlot >= static_cast<uint32_t>(materialTextureRefStride))
            {
                if (outError)
                {
                    std::ostringstream oss;
                    oss << "virtual-texture metadata: udim " << udim
                        << " exceeds material ref stride for material " << source.materialIndex;
                    *outError = oss.str();
                }
                return false;
            }

            const VirtualTextureUdimTable &udimTable = tileFile.udims.at(udim);
            if (udimTable.imageWidth == 0u || udimTable.imageHeight == 0u)
            {
                if (outError)
                {
                    std::ostringstream oss;
                    oss << "virtual-texture metadata: invalid UDIM extent for material "
                        << source.materialIndex << " udim " << udim;
                    *outError = oss.str();
                }
                return false;
            }

            LaunchParams::MaterialTextureRef &ref =
                out->materialTextureRefs[base + static_cast<size_t>(udimSlot)];
            ref.textureObject = 0ull;
            ref.width = static_cast<int>(udimTable.imageWidth);
            ref.height = static_cast<int>(udimTable.imageHeight);
            ref.valid = 1;
            ref.wrapS = static_cast<int>(source.wrapS);
            ref.wrapT = static_cast<int>(source.wrapT);
            ref._padding0 = 0;
            ref._padding1 = 0;

            reg.activeUdims.push_back(udim);
            VirtualTextureUdimExtent extent = {};
            extent.udim = udim;
            extent.width = udimTable.imageWidth;
            extent.height = udimTable.imageHeight;
            reg.udimExtents.push_back(extent);
            reg.width = std::max(reg.width, udimTable.imageWidth);
            reg.height = std::max(reg.height, udimTable.imageHeight);
            out->activeUdimCount++;
        }

        if (reg.activeUdims.empty() || reg.width == 0u || reg.height == 0u)
        {
            if (outError)
            {
                std::ostringstream oss;
                oss << "virtual-texture metadata: tile bin has no usable UDIM slots for material "
                    << source.materialIndex << ": " << tileFileIt->second;
                *outError = oss.str();
            }
            return false;
        }

        out->registrations.push_back(std::move(reg));
        out->totalVirtualTextureBytes += tileFile.totalTextureBytes;
        out->mappedMaterialCount++;
    }

    return true;
}

} // namespace texture
} // namespace ybi
