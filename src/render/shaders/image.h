#pragma once

#include "core/assert.h"
#include "device/cuda_compat.h"
#include "device/kernel_params.h"
#include "render/ray_cone_utils.h"
#include "render/shaders/shading.h"
#include "render/shaders/shading_helpers.h"
#include "render/virtual_texture_interop.h"
#include "render/virtual_texture_sampling.h"

namespace ybi
{

YBI_CPU_GPU MaterialTextureHandle LoadTextureHandle(const uint32_t *nodes, uint32_t &offset)
{
    MaterialTextureHandle result;
    result.udimBase = nodes[offset++];
    result.udimTableOffset = nodes[offset++];
    return result;
}

YBI_CPU_GPU Vec4f VirtualTextureSampleImage(const KernelParams &params,
                                            MaterialTextureHandle handle,
                                            const ShadingInfo &shadingInfo)
{
    Vec2f flooredUv = Floor(shadingInfo.uv);

    uint32_t udimTableIndex = 0u;
    if (handle.udimBase != ~0u)
    {
        const uint32_t udim =
            1001 + static_cast<uint32_t>(flooredUv.x) + 10 * static_cast<uint32_t>(flooredUv.y);
        udimTableIndex = handle.udimTableOffset + udim - handle.udimBase;
    }
    else
    {
        udimTableIndex = handle.udimTableOffset;
    }

    DeviceTextureDescriptor textureDescriptor = params.textureDescriptors[udimTableIndex];
    const uint32_t pageShift = 7;

    float lod = ComputeTextureLOD(shadingInfo.positions[0],
                                  shadingInfo.positions[1],
                                  shadingInfo.positions[2],
                                  shadingInfo.uvs[0],
                                  shadingInfo.uvs[1],
                                  shadingInfo.uvs[2],
                                  shadingInfo.rayDirection,
                                  shadingInfo.shadingNormal,
                                  shadingInfo.rayCone.width,
                                  textureDescriptor.width,
                                  textureDescriptor.height);

    const uint32_t maxMip = CalculateMaxVirtualTextureMip(textureDescriptor, pageShift);
    lod = Clamp(lod, 0.f, float(maxMip));

    uint32_t baseMip = uint32_t(Floor(lod));
    uint32_t nextMip = Min(baseMip + 1, maxMip);

    // TODO: don't hardcode
    UnpackedFeedbackRequest lowR, highR;
    Vec4f lowerValue = SamplePointVirtualTexture(params,
                                                 textureDescriptor,
                                                 pageShift,
                                                 shadingInfo.uv,
                                                 baseMip,
                                                 &lowR,
                                                 udimTableIndex,
                                                 handle.udimBase,
                                                 handle.udimTableOffset);
    Vec4f higherValue = SamplePointVirtualTexture(params,
                                                  textureDescriptor,
                                                  pageShift,
                                                  shadingInfo.uv,
                                                  nextMip,
                                                  &highR,
                                                  udimTableIndex,
                                                  handle.udimBase,
                                                  handle.udimTableOffset);

    WriteUniqueFeedback(params, lowR);
    WriteUniqueFeedback(params, highR);

    float t = lod - float(baseMip);
    Vec4f value = Lerp(lowerValue, higherValue, t);
    return value;
}

YBI_CPU_GPU Vec4f VirtualTextureSampleImage(const KernelParams &params,
                                            uint32_t &nodeOffset,
                                            const ShadingInfo &shadingInfo)
{
    const MaterialTextureHandle handle = LoadTextureHandle(params.scene.shaderNodes, nodeOffset);
    return VirtualTextureSampleImage(params, handle, shadingInfo);
}

YBI_CPU_GPU void ShaderImageVector3(const KernelParams &params,
                                    const ShadingInfo &info,
                                    float *registers,
                                    uint32_t &nodeOffset)
{
    Vec4f result = VirtualTextureSampleImage(params, nodeOffset, info);
    uint32_t dst = params.scene.shaderNodes[nodeOffset++];
    StoreFloat3Register(registers, result.Xyz(), dst);
}

YBI_CPU_GPU void ShaderImageFloat(const KernelParams &params,
                                  const ShadingInfo &info,
                                  float *registers,
                                  uint32_t &nodeOffset)
{
    YBI_ASSERT(0);
    VirtualTextureSampleImage(params, nodeOffset, info);
}

} // namespace ybi
