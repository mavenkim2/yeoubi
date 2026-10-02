#pragma once

#include "device/cuda_compat.h"
#include "kernel/launch/kernel_params.h"
#include "kernel/shading/shading.h"
#include "kernel/shading/shading_helpers.h"

namespace ybi
{

YBI_GPU void ShaderRemapColor3(const KernelParams &params,
                               const ShadingInfo &,
                               float *registers,
                               uint32_t &nodeOffset)
{
    const uint32_t *nodes = params.scene.shaderNodes;
    const Vec3f value = LoadFloat3ShaderOperand(nodes, nodeOffset, registers);
    const Vec3f inLow = LoadFloat3ShaderOperand(nodes, nodeOffset, registers);
    const Vec3f inHigh = LoadFloat3ShaderOperand(nodes, nodeOffset, registers);
    const Vec3f outLow = LoadFloat3ShaderOperand(nodes, nodeOffset, registers);
    const Vec3f outHigh = LoadFloat3ShaderOperand(nodes, nodeOffset, registers);
    const uint32_t dst = nodes[nodeOffset++];

    // ND_remap_color3 deliberately does not clamp.
    const Vec3f result = outLow + (value - inLow) * (outHigh - outLow) / (inHigh - inLow);
    StoreFloat3Register(registers, result, dst);
}

} // namespace ybi
