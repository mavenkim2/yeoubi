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

YBI_GPU void EvaluateShader(const KernelParams &params,
                            const ShadingInfo &info,
                            const DeviceMaterialExpression &shader,
                            float *registers)
{
    uint32_t offset = shader.bytecodeOffset;

    for (;;)
    {
        MaterialOpcode opcode = static_cast<MaterialOpcode>(params.scene.shaderNodes[offset++]);
        switch (opcode)
        {
            case MaterialOpcode::ConstantBool:
            {
                const bool value = params.scene.shaderNodes[offset++] != 0u;
                const uint32_t dst = params.scene.shaderNodes[offset++];
                StoreBoolRegister(registers, value, dst);
                break;
            }
            case MaterialOpcode::ConstantFloat:
            {
                float value = UintAsFloat(params.scene.shaderNodes[offset++]);
                const uint32_t dst = params.scene.shaderNodes[offset++];
                StoreFloatRegister(registers, value, dst);
                break;
            }
            case MaterialOpcode::ConstantColor:
            {
                const Vec3f value(UintAsFloat(params.scene.shaderNodes[offset++]),
                                  UintAsFloat(params.scene.shaderNodes[offset++]),
                                  UintAsFloat(params.scene.shaderNodes[offset++]));
                const uint32_t dst = params.scene.shaderNodes[offset++];
                StoreFloat3Register(registers, value, dst);
                break;
            }
            case MaterialOpcode::ImageFloat:
            {
                ShaderImageFloat(params, info, registers, offset);
                break;
            }
            case MaterialOpcode::ImageVector3:
            {
                ShaderImageVector3(params, info, registers, offset);
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
            case MaterialOpcode::ExtractVector3:
            {
                const uint32_t src = params.scene.shaderNodes[offset++];
                const uint32_t index = params.scene.shaderNodes[offset++];
                const uint32_t dst = params.scene.shaderNodes[offset++];

                const Vec3f value = LoadFloat3Register(registers, src);
                StoreFloatRegister(registers, value[index], dst);
                break;
            }
            case MaterialOpcode::Count:
            {
                return;
            }
            default:
                YBI_ASSERT(0);
        }
    }
}

YBI_GPU bool EvaluateShaderBool(const KernelParams &params,
                                const ShadingInfo &info,
                                const DeviceMaterialExpression &shader)
{
    float registers[SHADER_NUM_REGISTERS];
    EvaluateShader(params, info, shader, registers);
    return LoadBoolRegister(registers, shader.outputRegister);
}

YBI_GPU Vec2f EvaluateShaderVec2(const KernelParams &params,
                                 const ShadingInfo &info,
                                 const DeviceMaterialExpression &shader)
{
    float registers[SHADER_NUM_REGISTERS];
    EvaluateShader(params, info, shader, registers);
    return LoadFloat2Register(registers, shader.outputRegister);
}

YBI_GPU Vec3f EvaluateShaderVec3(const KernelParams &params,
                                 const ShadingInfo &info,
                                 const DeviceMaterialExpression &shader)
{
    float registers[SHADER_NUM_REGISTERS];
    EvaluateShader(params, info, shader, registers);
    return LoadFloat3Register(registers, shader.outputRegister);
}

YBI_GPU Vec4f EvaluateShaderVec4(const KernelParams &params,
                                 const ShadingInfo &info,
                                 const DeviceMaterialExpression &shader)
{
    float registers[SHADER_NUM_REGISTERS];
    EvaluateShader(params, info, shader, registers);
    return LoadFloat4Register(registers, shader.outputRegister);
}

YBI_GPU float EvaluateShaderFloat(const KernelParams &params,
                                  const ShadingInfo &info,
                                  const DeviceMaterialExpression &shader)
{
    float registers[SHADER_NUM_REGISTERS];
    EvaluateShader(params, info, shader, registers);
    return LoadFloatRegister(registers, shader.outputRegister);
}

} // namespace ybi
