#pragma once

#include "core/assert.h"
#include "core/bit.h"
#include "device/cuda_compat.h"
#include <type_traits>

namespace ybi
{

// Booleans occupy one float register, with canonical values 0.0f and 1.0f.
YBI_GPU bool LoadBoolRegister(float *registers, int srcRegister)
{
    YBI_ASSERT(srcRegister < SHADER_NUM_REGISTERS);
    return registers[srcRegister] != 0.f;
}

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

YBI_GPU Vec2f LoadFloat2Register(float *registers, int srcRegister)
{
    YBI_ASSERT(srcRegister + 1 < SHADER_NUM_REGISTERS);
    Vec2f v;
    v[0] = registers[srcRegister];
    v[1] = registers[srcRegister + 1];
    return v;
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

YBI_GPU Vec4f LoadFloat4Register(float *registers, int srcRegister)
{
    YBI_ASSERT(srcRegister + 3 < SHADER_NUM_REGISTERS);
    Vec4f v;
    v[0] = registers[srcRegister];
    v[1] = registers[srcRegister + 1];
    v[2] = registers[srcRegister + 2];
    v[3] = registers[srcRegister + 3];
    return v;
}

YBI_GPU void StoreFloatRegister(float *registers, float value, int dstRegister)
{
    YBI_ASSERT(dstRegister < SHADER_NUM_REGISTERS);
    registers[dstRegister] = value;
}

YBI_GPU void StoreBoolRegister(float *registers, bool value, int dstRegister)
{
    YBI_ASSERT(dstRegister < SHADER_NUM_REGISTERS);
    registers[dstRegister] = value ? 1.f : 0.f;
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

YBI_GPU void StoreFloat4Register(float *registers, const Vec4f &value, int dstRegister)
{
    YBI_ASSERT(dstRegister + 3 < SHADER_NUM_REGISTERS);
    registers[dstRegister] = value[0];
    registers[dstRegister + 1] = value[1];
    registers[dstRegister + 2] = value[2];
    registers[dstRegister + 3] = value[3];
}

} // namespace ybi
