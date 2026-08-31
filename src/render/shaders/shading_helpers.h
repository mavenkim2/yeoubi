#pragma once

#include "core/assert.h"
#include "core/bit.h"
#include "device/cuda_compat.h"

namespace ybi
{

YBI_GPU int LoadIntRegister(float *registers, int srcRegister)
{
    YBI_ASSERT(srcRegister < SHADER_NUM_REGISTERS);
    int result = FloatAsInt(registers[srcRegister]);
    return result;
}

YBI_GPU float LoadFloatRegister(float *registers, int srcRegister)
{
    YBI_ASSERT(srcRegister < SHADER_NUM_REGISTERS);
    float result = registers[srcRegister];
    return result;
}

YBI_GPU Vec3f LoadFloat3Register(float *registers, int srcRegister)
{
    YBI_ASSERT(srcRegister + 2 < SHADER_NUM_REGISTERS);
    Vec3f v;
    v[0] = registers[srcRegister];
    v[1] = registers[srcRegister + 1];
    v[2] = registers[srcRegister + 2];
    return v;
}

YBI_GPU void StoreFloatRegister(float *registers, float value, int dstRegister)
{
    YBI_ASSERT(dstRegister < SHADER_NUM_REGISTERS);
    registers[dstRegister] = value;
}

YBI_GPU void StoreFloat2Register(float *registers, const Vec2f &value, int dstRegister)
{
    YBI_ASSERT(dstRegister + 1 < SHADER_NUM_REGISTERS);
    registers[dstRegister] = value[0];
    registers[dstRegister + 1] = value[1];
}

YBI_GPU void StoreFloat3Register(float *registers, const Vec3f &value, int dstRegister)
{
    YBI_ASSERT(dstRegister + 2 < SHADER_NUM_REGISTERS);
    registers[dstRegister] = value[0];
    registers[dstRegister + 1] = value[1];
    registers[dstRegister + 2] = value[2];
}

} // namespace ybi
