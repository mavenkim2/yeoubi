#include "core/assert.h"
#include "core/bit.h"
#include "device/cuda_compat.h"
#include "render/material_opcodes.h"
#include "render/material_types.h"

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

YBI_CPU_GPU int LoadFloatRegister(float *registers, int srcRegister)
{
    YBI_ASSERT(srcRegister < SHADER_NUM_REGISTERS);
    int result = registers[srcRegister]);
    return result;
}

YBI_CPU_GPU void StoreFloatRegister(float *registers, float value, int dstRegister)
{
    YBI_ASSERT(dstRegister < SHADER_NUM_REGISTERS);
    registers[dstRegister] = value;
}

YBI_CPU_GPU void ShaderExtractVec3(float *registers, int dst, int src, int extractSrc)
{
    int index = LoadIntRegister(registers, srcRegister);
    float value = LoadFloatRegister(registers, src + index);
    StoreFloatRegister(registers, value, dst);
}

YBI_CPU_GPU void EvaluateShader(KernelParams &params, int offset)
{
    float registers[SHADER_NUM_REGISTERS];

    // TODO
    MaterialOpcode opcode;

    switch (opcode)
    {
        case MaterialOpcode::ExtractVector3:
            // TODO
            ShaderExtractVec3(registers, 0, 0);
            break;
        case MaterialOpcode::Image:
            break;
    }
}

} // namespace ybi
