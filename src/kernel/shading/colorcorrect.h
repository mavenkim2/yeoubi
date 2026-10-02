#pragma once

#include "core/assert.h"
#include "device/cuda_compat.h"
#include "kernel/launch/shading_context.h"
#include "kernel/shading/shading.h"
#include "kernel/shading/shading_helpers.h"

namespace ybi
{

YBI_GPU void ShaderColorCorrectColor3(const ShadingContext &context,
                                      const ShadingInfo &info,
                                      float *registers,
                                      uint32_t &nodeOffset)
{
    uint32_t src = context.scene.shaderNodes[nodeOffset++];
    bool isNode = context.scene.shaderNodes[nodeOffset++];
    float gamma;

    if (!isNode)
    {
        gamma = UintAsFloat(context.scene.shaderNodes[nodeOffset++]);
    }
    else
    {
        // uint32_t gammaRegister = params.scene.shaderNodes[nodeOffset++];
        // gamma = LoadFloatRegister(...);
        YBI_ASSERT(0);
    }
    uint32_t dst = context.scene.shaderNodes[nodeOffset++];

    Vec3f v = LoadFloat3Register(registers, src);

    float ex = 1.f / gamma;
    v.x = powf(v.x, ex);
    v.y = powf(v.y, ex);
    v.z = powf(v.z, ex);

    StoreFloat3Register(registers, v, dst);
}

// Materials

} // namespace ybi
