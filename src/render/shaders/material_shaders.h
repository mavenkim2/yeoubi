#pragma once

#include "core/assert.h"
#include "device/cuda_compat.h"
#include "device/kernel_params.h"
#include "render/material_opcodes.h"
#include "render/material_types.h"
#include "render/shaders/attribute.h"
#include "render/shaders/colorcorrect.h"
#include "render/shaders/image.h"
#include "render/shaders/shading_helpers.h"

namespace ybi
{

template <typename T>
YBI_CPU_GPU T EvaluateShader(const KernelParams &params,
                             const ShadingInfo &info,
                             const DeviceMaterialExpression &shader)
{
    float registers[SHADER_NUM_REGISTERS];
    uint32_t offset = shader.bytecodeOffset;

    for (;;)
    {
        MaterialOpcode opcode = static_cast<MaterialOpcode>(params.scene.shaderNodes[offset++]);
        switch (opcode)
        {
            case MaterialOpcode::ConstantColor:
            {
                const Vec3f value(UintAsFloat(params.scene.shaderNodes[offset++]),
                                  UintAsFloat(params.scene.shaderNodes[offset++]),
                                  UintAsFloat(params.scene.shaderNodes[offset++]));
                const uint32_t dst = params.scene.shaderNodes[offset++];
                StoreFloat3Register(registers, value, dst);
                break;
            }
            case MaterialOpcode::ImageVector3:
            {
                ShaderImageVector3(params, info, registers, offset);
                break;
            }
            case MaterialOpcode::ImageFloat:
            {
                ShaderImageFloat(params, info, registers, offset);
                break;
            }
            case MaterialOpcode::AttributeValueVector2:
            {
                ShaderAttributeValueVector2(params, info, registers, offset);
                break;
            }
            case MaterialOpcode::ColorCorrectColor3:
            {
                ShaderColorCorrectColor3(params, info, registers, offset);
                break;
            }
            case MaterialOpcode::Count:
            {
                return LoadTemplateRegister<T>(registers, shader.outputRegister);
            }
            default:
                YBI_ASSERT(0);
        }
    }
}

} // namespace ybi
