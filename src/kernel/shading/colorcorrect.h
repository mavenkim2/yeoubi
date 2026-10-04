#pragma once

#include "core/assert.h"
#include "device/cuda_compat.h"
#include "kernel/launch/kernel_params.h"
#include "kernel/shading/shading.h"
#include "kernel/shading/shading_helpers.h"

namespace ybi
{

YBI_GPU void ShaderColorCorrectColor3(const KernelParams &params,
                                      const ShadingInfo &info,
                                      float *registers,
                                      uint32_t &nodeOffset)
{
    uint32_t src = params.scene.shaderNodes[nodeOffset++];
    bool isNode = params.scene.shaderNodes[nodeOffset++];
    float gamma;

    if (!isNode)
    {
        gamma = UintAsFloat(params.scene.shaderNodes[nodeOffset++]);
    }
    else
    {
        // uint32_t gammaRegister = params.scene.shaderNodes[nodeOffset++];
        // gamma = LoadFloatRegister(...);
        YBI_ASSERT(0);
    }
    uint32_t dst = params.scene.shaderNodes[nodeOffset++];

    Vec3f v = LoadFloat3Register(registers, src);

    float ex = 1.f / gamma;
    v.x = powf(v.x, ex);
    v.y = powf(v.y, ex);
    v.z = powf(v.z, ex);

    StoreFloat3Register(registers, v, dst);
}

// Materials

} // namespace ybi
