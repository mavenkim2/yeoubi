#pragma once

#include "core/math.h"
#include "core/vec3f.h"
#include "device/cuda_compat.h"
#include "kernel/launch/kernel_params.h"
#include "kernel/shading/shading.h"
#include "kernel/shading/shading_helpers.h"

namespace ybi
{

YBI_GPU void ShaderHeightToNormalVector3(const KernelParams &params,
                                         const ShadingInfo &,
                                         float *registers,
                                         uint32_t &nodeOffset)
{
    const uint32_t *nodes = params.scene.shaderNodes;
    (void)LoadFloatShaderOperand(nodes, nodeOffset, registers);
    (void)LoadFloatShaderOperand(nodes, nodeOffset, registers);
    const uint32_t dst = nodes[nodeOffset++];

    // A generic height graph needs ray differentials (MaterialX Dx/Dy), which
    // this path tracer does not carry. Keep the graph executable and return the
    // MaterialX flat-normal encoding until geometry-normal graphs are enabled.
    StoreFloat3Register(registers, Vec3f(0.5f, 0.5f, 1.f), dst);
}

YBI_GPU void ShaderInvertFloat(const KernelParams &params,
                               const ShadingInfo &,
                               float *registers,
                               uint32_t &nodeOffset)
{
    const uint32_t *nodes = params.scene.shaderNodes;
    const float value = LoadFloatShaderOperand(nodes, nodeOffset, registers);
    const float amount = LoadFloatShaderOperand(nodes, nodeOffset, registers);
    const uint32_t dst = nodes[nodeOffset++];
    StoreFloatRegister(registers, amount - value, dst);
}

YBI_GPU void ShaderCombine3Vector3(const KernelParams &params,
                                   const ShadingInfo &,
                                   float *registers,
                                   uint32_t &nodeOffset)
{
    const uint32_t *nodes = params.scene.shaderNodes;
    const float x = LoadFloatShaderOperand(nodes, nodeOffset, registers);
    const float y = LoadFloatShaderOperand(nodes, nodeOffset, registers);
    const float z = LoadFloatShaderOperand(nodes, nodeOffset, registers);
    const uint32_t dst = nodes[nodeOffset++];
    StoreFloat3Register(registers, Vec3f(x, y, z), dst);
}

YBI_GPU void ShaderRamp4Float(const KernelParams &params,
                              const ShadingInfo &,
                              float *registers,
                              uint32_t &nodeOffset)
{
    const uint32_t *nodes = params.scene.shaderNodes;
    const float valueTopLeft = LoadFloatShaderOperand(nodes, nodeOffset, registers);
    const float valueTopRight = LoadFloatShaderOperand(nodes, nodeOffset, registers);
    const float valueBottomLeft = LoadFloatShaderOperand(nodes, nodeOffset, registers);
    const float valueBottomRight = LoadFloatShaderOperand(nodes, nodeOffset, registers);
    Vec2f texcoord = LoadFloat2ShaderOperand(nodes, nodeOffset, registers);
    const uint32_t dst = nodes[nodeOffset++];

    texcoord.x = Clamp(texcoord.x, 0.f, 1.f);
    texcoord.y = Clamp(texcoord.y, 0.f, 1.f);
    const float top = Lerp(valueTopLeft, valueTopRight, texcoord.x);
    const float bottom = Lerp(valueBottomLeft, valueBottomRight, texcoord.x);
    StoreFloatRegister(registers, Lerp(top, bottom, texcoord.y), dst);
}

} // namespace ybi
