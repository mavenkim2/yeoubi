#include "core/assert.h"
#include "device/cuda_compat.h"
#include "device/kernel_params.h"
#include "render/material_opcodes.h"
#include "render/material_types.h"
#include "render/shaders/attribute.h"
#include "render/shaders/colorcorrect.h"
#include "render/shaders/image.h"

namespace ybi
{

struct Something
{
    MaterialOpcode opcode;
    int dst : 8;
    int src0 : 8;
    int src1 : 8;
};

#if 0
YBI_CPU_GPU void ShaderExtractVec3(float *registers, int dst, int src, int extractSrc)
{
    int index = LoadIntRegister(registers, srcRegister);
    float value = LoadFloatRegister(registers, src + index);
    StoreFloatRegister(registers, value, dst);
}
#endif

#if 0
X(ImageVector3, Vec3) \
X(ImageFloat, Float) \
X(AttributeValueVector2, Vec2) \
X(HeightToNormalVector3, Vec3) \
X(InvertFloat, Float) \
X(ExtractVector3, Float) \
X(Combine3Vector3, Vec3) \
X(NormalMapVector2, Vec3) \
X(ColorCorrectColor3, Vec3) \
X(RemapColor3, Vec3) \
X(Ramp4Float, Float) \
X(Separate3Color3, Float)
#endif

YBI_CPU_GPU Vec3f EvaluateShader(const KernelParams &params,
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
                return LoadFloat3Register(registers, shader.outputRegister);
            }
            default:
                YBI_ASSERT(0);
        }
    }
}

} // namespace ybi
