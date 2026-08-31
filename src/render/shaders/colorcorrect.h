#pragma once

#include "core/assert.h"
#include "device/cuda_compat.h"
#include "device/kernel_params.h"
#include "render/shaders/shading.h"
#include "render/shaders/shading_helpers.h"
#include <cmath>

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

    v.x = powf(v.x, gamma);
    v.y = powf(v.y, gamma);
    v.z = powf(v.z, gamma);

    StoreFloat3Register(registers, v, dst);
}

// Materials

} // namespace ybi
