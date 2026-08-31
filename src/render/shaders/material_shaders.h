#include "core/assert.h"
#include "core/bit.h"
#include "device/cuda_compat.h"
#include "device/kernel_params.h"
#include "render/material_opcodes.h"
#include "render/material_types.h"
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

YBI_CPU_GPU int LoadIntRegister(float *registers, int srcRegister)
{
    YBI_ASSERT(srcRegister < SHADER_NUM_REGISTERS);
    int result = FloatAsInt(registers[srcRegister]);
    return result;
}

YBI_CPU_GPU float LoadFloatRegister(float *registers, int srcRegister)
{
    YBI_ASSERT(srcRegister < SHADER_NUM_REGISTERS);
    float result = registers[srcRegister];
    return result;
}

YBI_CPU_GPU void StoreFloatRegister(float *registers, float value, int dstRegister)
{
    YBI_ASSERT(dstRegister < SHADER_NUM_REGISTERS);
    registers[dstRegister] = value;
}

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

YBI_CPU_GPU void EvaluateShader(const KernelParams &params, int offset)
{
    float registers[SHADER_NUM_REGISTERS];

    // TODO
    MaterialOpcode opcode;

    switch (opcode)
    {
        case MaterialOpcode::ExtractVector3:
            // TODO
            // ShaderExtractVec3(registers, 0, 0);
            break;
        case MaterialOpcode::ImageVector3:
            break;
    }
}

} // namespace ybi
