#pragma once

#include "core/vec2f.h"
#include "core/vec2u.h"
#include "core/vec3f.h"
#include "core/vec3u.h"
#include "core/vec4f.h"
#include "core/vec4u.h"
#include "device/cuda_compat.h"
#include <cstdint>

namespace ybi
{

struct RNG
{
    static YBI_CPU_GPU uint32_t PCG(uint32_t x)
    {
        uint32_t state = x * 747796405u + 2891336453u;
        uint32_t word = ((state >> ((state >> 28u) + 4u)) ^ state) * 277803737u;
        return (word >> 22u) ^ word;
    }

    // Ref: M. Jarzynski and M. Olano, "Hash Functions for GPU Rendering," Journal of Computer Graphics Techniques, 2020.
    static YBI_CPU_GPU Vec3u PCG3d(Vec3u v)
    {
        v = v * 1664525u + Vec3u(1013904223u);
        v.x += v.y * v.z;
        v.y += v.z * v.x;
        v.z += v.x * v.y;
        v ^= v >> 16u;
        v.x += v.y * v.z;
        v.y += v.z * v.x;
        v.z += v.x * v.y;
        return v;
    }

    // Ref: M. Jarzynski and M. Olano, "Hash Functions for GPU Rendering," Journal of Computer Graphics Techniques, 2020.
    static YBI_CPU_GPU Vec4u PCG4d(Vec4u v)
    {
        v = v * 1664525u + Vec4u(1013904223u);
        v.x += v.y * v.w;
        v.y += v.z * v.x;
        v.z += v.x * v.y;
        v.w += v.y * v.z;
        v ^= v >> 16u;
        v.x += v.y * v.w;
        v.y += v.z * v.x;
        v.z += v.x * v.y;
        v.w += v.y * v.z;
        return v;
    }

    static YBI_CPU_GPU RNG Init(Vec2u pixel, uint32_t frame)
    {
        RNG rng;
        rng.State = RNG::PCG3d(Vec3u(pixel.x, pixel.y, frame)).x;

        return rng;
    }

    static YBI_CPU_GPU RNG Init(Vec2u pixel, uint32_t frame, uint32_t idx)
    {
        RNG rng;
        rng.State = RNG::PCG4d(Vec4u(pixel.x, pixel.y, frame, idx)).x;

        return rng;
    }

    static YBI_CPU_GPU RNG Init(uint32_t idx, uint32_t frame)
    {
        RNG rng;
        rng.State = PCG(idx + PCG(frame));

        return rng;
    }

    static YBI_CPU_GPU RNG Init(uint32_t seed)
    {
        RNG rng;
        rng.State = seed;

        return rng;
    }

    YBI_CPU_GPU uint32_t UniformUint()
    {
        State = State * 747796405u + 2891336453u;
        uint32_t word = ((State >> ((State >> 28u) + 4u)) ^ State) * 277803737u;

        return (word >> 22u) ^ word;
    }

    YBI_CPU_GPU float Uniform()
    {
        // For 32-bit floats, any integer in [0, 2^24] can be represented exactly and
        // there may be rounding errors for anything larger, e.g. 2^24 + 1 is rounded
        // down to 2^24.
        // Given random integers, we can right shift by 8 bits to get integers in
        // [0, 2^24 - 1]. After multiplication by 2^-24, we have uniform numbers in [0, 1).
        // Ref: https://prng.di.unimi.it/
        return float(UniformUint() >> 8) * 0x1p-24f;
    }

    // Returns samples in [0, bound). Requires bound > 0.
    YBI_CPU_GPU uint32_t UniformUintBounded(uint32_t bound)
    {
        YBI_ASSERT(bound > 0);
        uint32_t threshold = (~bound + 1u) % bound;

        for (;;)
        {
            uint32_t r = UniformUint();

            if (r >= threshold)
                return r % bound;
        }
    }

    // Returns samples in [0, bound). Requires bound > 0. Biased but faster than UniformUintBounded():
    // https://en.wikipedia.org/wiki/Fisher%E2%80%93Yates_shuffle
    YBI_CPU_GPU uint32_t UniformUintBounded_Faster(uint32_t bound)
    {
        YBI_ASSERT(bound > 0);
        return uint32_t((uint64_t(UniformUint()) * bound) >> 32u);
    }

    YBI_CPU_GPU Vec2f Uniform2D()
    {
        float u0 = Uniform();
        float u1 = Uniform();

        return Vec2f(u0, u1);
    }

    YBI_CPU_GPU Vec3f Uniform3D()
    {
        float u0 = Uniform();
        float u1 = Uniform();
        float u2 = Uniform();

        return Vec3f(u0, u1, u2);
    }

    YBI_CPU_GPU Vec4f Uniform4D()
    {
        float u0 = Uniform();
        float u1 = Uniform();
        float u2 = Uniform();
        float u3 = Uniform();

        return Vec4f(u0, u1, u2, u3);
    }

    uint32_t State = 0;
};

} // namespace ybi
