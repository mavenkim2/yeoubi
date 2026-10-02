#pragma once

#include "core/assert.h"
#include "device/cuda_compat.h"
#include "kernel/bsdf/closure_types.h"
#include "kernel/integrator/ray_cone.h"
#include "kernel/launch/shading_context.h"
#include "kernel/shading/shading.h"
#include "kernel/shading/shading_helpers.h"
#include "kernel/texture/types.h"
#include "kernel/texture/virtual_sampling.h"

namespace ybi
{

YBI_CPU_GPU MaterialTextureHandle LoadTextureHandle(const uint32_t *nodes, uint32_t &offset)
{
    MaterialTextureHandle result;
    result.udimBase = nodes[offset++];
    result.udimTableOffset = nodes[offset++];
    return result;
}

YBI_CPU_GPU uint32_t GetUdimTableIndex(MaterialTextureHandle handle, const Vec2f &uv)
{
    Vec2f flooredUv = Floor(uv);

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
    return udimTableIndex;
}

YBI_CPU_GPU Vec4f VirtualTextureSampleImage(
    const ShadingContext &context,
    MaterialTextureHandle handle,
    uint32_t udimTableIndex,
    const Vec2f &uv,
    float triangleArea,
    float uvArea,
    const Vec3f &rayDirection,
    const Vec3f &shadingNormal,
    float rayConeWidth)
{
    DeviceTextureDescriptor textureDescriptor = context.textureDescriptors[udimTableIndex];
    const uint32_t pageShift = 7;

    float lod = ComputeTextureLOD(triangleArea,
                                  uvArea,
                                  rayDirection,
                                  shadingNormal,
                                  rayConeWidth,
                                  textureDescriptor.width,
                                  textureDescriptor.height);

    const uint32_t maxMip = CalculateMaxVirtualTextureMip(textureDescriptor, pageShift);
    lod = Clamp(lod, 0.f, float(maxMip));

    const uint32_t baseMip = uint32_t(Floor(lod));
    const uint32_t nextMip = Min(baseMip + 1, maxMip);

    UnpackedFeedbackRequest lowRequest, highRequest;
    const Vec4f lowerValue = SamplePointVirtualTexture(context,
                                                        textureDescriptor,
                                                        pageShift,
                                                        uv,
                                                        baseMip,
                                                        &lowRequest,
                                                        udimTableIndex,
                                                        handle.udimBase,
                                                        handle.udimTableOffset);
    const Vec4f higherValue = SamplePointVirtualTexture(context,
                                                         textureDescriptor,
                                                         pageShift,
                                                         uv,
                                                         nextMip,
                                                         &highRequest,
                                                         udimTableIndex,
                                                         handle.udimBase,
                                                         handle.udimTableOffset);

    WriteUniqueFeedback(context, lowRequest);
    WriteUniqueFeedback(context, highRequest);
    return Lerp(lowerValue, higherValue, lod - float(baseMip));
}

YBI_CPU_GPU Vec4f VirtualTextureSampleImage(const ShadingContext &context,
                                            MaterialTextureHandle handle,
                                            const ClosureInitializationInfo &info)
{
    return VirtualTextureSampleImage(context,
                                     handle,
                                     GetUdimTableIndex(handle, info.uv),
                                     info.uv,
                                     info.triangleArea,
                                     info.uvArea,
                                     info.rayDirection,
                                     info.shadingNormal,
                                     info.rayConeWidth);
}

YBI_CPU_GPU Vec4f
VirtualTextureSampleTiledImage(const ShadingContext &context,
                               MaterialTextureHandle handle,
                               const ClosureInitializationInfo &info,
                               Vec2f uvTiling,
                               Vec2f uvOffset,
                               Vec2f realWorldImageSize,
                               Vec2f realWorldTileSize)
{
    const uint32_t udimTableIndex = GetUdimTableIndex(handle, info.uv);
    const Vec2f scale = uvTiling / realWorldImageSize * realWorldTileSize;
    const Vec2f offset = uvOffset / realWorldImageSize * realWorldTileSize;
    const Vec2f uv = info.uv * scale - offset;
    const float uvArea = info.uvArea * Abs(scale.x * scale.y);
    return VirtualTextureSampleImage(context,
                                     handle,
                                     udimTableIndex,
                                     uv,
                                     info.triangleArea,
                                     uvArea,
                                     info.rayDirection,
                                     info.shadingNormal,
                                     info.rayConeWidth);
}

YBI_CPU_GPU Vec4f VirtualTextureSampleImage(const ShadingContext &context,
                                            MaterialTextureHandle handle,
                                            const ShadingInfo &shadingInfo,
                                            uint32_t udimTableIndex,
                                            const Vec2f &uv,
                                            const Vec2f (&uvs)[3])
{
    const float triangleArea =
        Length(Cross(shadingInfo.positions[1] - shadingInfo.positions[0],
                     shadingInfo.positions[2] - shadingInfo.positions[0]));
    const float uvArea =
        Abs((uvs[1].x - uvs[0].x) * (uvs[2].y - uvs[0].y) -
            (uvs[1].y - uvs[0].y) * (uvs[2].x - uvs[0].x));
    return VirtualTextureSampleImage(context,
                                     handle,
                                     udimTableIndex,
                                     uv,
                                     triangleArea,
                                     uvArea,
                                     shadingInfo.rayDirection,
                                     shadingInfo.shadingNormal,
                                     shadingInfo.rayCone.width);
}

YBI_CPU_GPU Vec4f VirtualTextureSampleImage(const ShadingContext &context,
                                            MaterialTextureHandle handle,
                                            const ShadingInfo &shadingInfo)
{
    return VirtualTextureSampleImage(context,
                                     handle,
                                     shadingInfo,
                                     GetUdimTableIndex(handle, shadingInfo.uv),
                                     shadingInfo.uv,
                                     shadingInfo.uvs);
}

YBI_CPU_GPU Vec4f
VirtualTextureSampleTiledImage(const ShadingContext &context,
                               MaterialTextureHandle handle,
                               const ShadingInfo &shadingInfo,
                               Vec2f uvTiling,
                               Vec2f uvOffset,
                               Vec2f realWorldImageSize,
                               Vec2f realWorldTileSize)
{
    // The texture file is selected from the original UV tile. Tiling only changes
    // coordinates within that image; keep vertex UVs unwrapped for the LOD area.
    const uint32_t udimTableIndex = GetUdimTableIndex(handle, shadingInfo.uv);
    const Vec2f scale = uvTiling / realWorldImageSize * realWorldTileSize;
    const Vec2f offset = uvOffset / realWorldImageSize * realWorldTileSize;
    const Vec2f uv = shadingInfo.uv * scale - offset;
    const Vec2f uvs[3] = {shadingInfo.uvs[0] * scale - offset,
                          shadingInfo.uvs[1] * scale - offset,
                          shadingInfo.uvs[2] * scale - offset};
    return VirtualTextureSampleImage(context, handle, shadingInfo, udimTableIndex, uv, uvs);
}

YBI_CPU_GPU Vec4f VirtualTextureSampleImage(const ShadingContext &context,
                                            uint32_t &nodeOffset,
                                            const ShadingInfo &shadingInfo)
{
    const MaterialTextureHandle handle = LoadTextureHandle(context.scene.shaderNodes, nodeOffset);
    return VirtualTextureSampleImage(context, handle, shadingInfo);
}

YBI_CPU_GPU void ShaderImageVector3(const ShadingContext &context,
                                    const ShadingInfo &info,
                                    float *registers,
                                    uint32_t &nodeOffset)
{
    Vec4f result = VirtualTextureSampleImage(context, nodeOffset, info);
    uint32_t dst = context.scene.shaderNodes[nodeOffset++];
    StoreFloat3Register(registers, result.Xyz(), dst);
}

YBI_CPU_GPU void ShaderImageVector2(const ShadingContext &context,
                                    const ShadingInfo &info,
                                    float *registers,
                                    uint32_t &nodeOffset)
{
    Vec4f result = VirtualTextureSampleImage(context, nodeOffset, info);
    uint32_t dst = context.scene.shaderNodes[nodeOffset++];
    StoreFloat2Register(registers, result.Xy(), dst);
}

YBI_CPU_GPU void ShaderImageVector4(const ShadingContext &context,
                                    const ShadingInfo &info,
                                    float *registers,
                                    uint32_t &nodeOffset)
{
    Vec4f result = VirtualTextureSampleImage(context, nodeOffset, info);
    uint32_t dst = context.scene.shaderNodes[nodeOffset++];
    StoreFloat4Register(registers, result, dst);
}

YBI_CPU_GPU void ShaderImageFloat(const ShadingContext &context,
                                  const ShadingInfo &info,
                                  float *registers,
                                  uint32_t &nodeOffset)
{
    Vec4f result = VirtualTextureSampleImage(context, nodeOffset, info);
    uint32_t dst = context.scene.shaderNodes[nodeOffset++];
    StoreFloatRegister(registers, result.x, dst);
}

} // namespace ybi
