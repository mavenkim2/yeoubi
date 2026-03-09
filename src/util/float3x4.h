#pragma once

#include "util/vec4.h"

YBI_NAMESPACE_BEGIN

// NOTE: row major, multiplication order right to left
struct Float3x4
{
    float m[3][4];

    YBI_INTEGRATOR_HD Float3x4()
        : m{{0.0f, 0.0f, 0.0f, 0.0f}, {0.0f, 0.0f, 0.0f, 0.0f}, {0.0f, 0.0f, 0.0f, 0.0f}}
    {
    }
    YBI_INTEGRATOR_HD Float3x4(float m00,
                               float m01,
                               float m02,
                               float m03,
                               float m10,
                               float m11,
                               float m12,
                               float m13,
                               float m20,
                               float m21,
                               float m22,
                               float m23)
    {
        m[0][0] = m00;
        m[0][1] = m01;
        m[0][2] = m02;
        m[0][3] = m03;

        m[1][0] = m10;
        m[1][1] = m11;
        m[1][2] = m12;
        m[1][3] = m13;

        m[2][0] = m20;
        m[2][1] = m21;
        m[2][2] = m22;
        m[2][3] = m23;
    }

    YBI_INTEGRATOR_HD Float3x4(Vec4 r0, Vec4 r1, Vec4 r2)
    {
        for (int c = 0; c < 4; ++c)
        {
            m[0][c] = r0[c];
            m[1][c] = r1[c];
            m[2][c] = r2[c];
        }
    }
};

using float3x4 = Float3x4;

YBI_INTEGRATOR_HD Vec3 operator*(const Float3x4 &m, const Vec4 &p)
{
    return Vec3(Dot(Vec4(m.m[0][0], m.m[0][1], m.m[0][2], m.m[0][3]), p),
                Dot(Vec4(m.m[1][0], m.m[1][1], m.m[1][2], m.m[1][3]), p),
                Dot(Vec4(m.m[2][0], m.m[2][1], m.m[2][2], m.m[2][3]), p));
}

YBI_NAMESPACE_END
