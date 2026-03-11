#pragma once

#include "util/float3x4.h"

namespace ybi
{

// Row-major 4x4. Matrix multiplication is right-to-left.
struct Float4x4
{
    float m[4][4];

    YBI_INTEGRATOR_HD Float4x4()
        : m{{0.0f, 0.0f, 0.0f, 0.0f},
            {0.0f, 0.0f, 0.0f, 0.0f},
            {0.0f, 0.0f, 0.0f, 0.0f},
            {0.0f, 0.0f, 0.0f, 0.0f}}
    {
    }
    YBI_INTEGRATOR_HD Float4x4(float m00,
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
                               float m23,
                               float m30,
                               float m31,
                               float m32,
                               float m33)
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

        m[3][0] = m30;
        m[3][1] = m31;
        m[3][2] = m32;
        m[3][3] = m33;
    }

    YBI_INTEGRATOR_HD Float4x4(Vec4 r0, Vec4 r1, Vec4 r2, Vec4 r3)
    {
        for (int c = 0; c < 4; ++c)
        {
            m[0][c] = r0[c];
            m[1][c] = r1[c];
            m[2][c] = r2[c];
            m[3][c] = r3[c];
        }
    }

    YBI_INTEGRATOR_HD static Float4x4 Identity()
    {
        return Float4x4(1.0f,
                        0.0f,
                        0.0f,
                        0.0f,
                        0.0f,
                        1.0f,
                        0.0f,
                        0.0f,
                        0.0f,
                        0.0f,
                        1.0f,
                        0.0f,
                        0.0f,
                        0.0f,
                        0.0f,
                        1.0f);
    }
};

using float4x4 = Float4x4;

YBI_INTEGRATOR_HD Float4x4 operator*(const Float4x4 &a, const Float4x4 &b)
{
    Float4x4 c;
    for (int i = 0; i < 4; ++i)
    {
        for (int j = 0; j < 4; ++j)
        {
            c.m[i][j] = a.m[i][0] * b.m[0][j] + a.m[i][1] * b.m[1][j] + a.m[i][2] * b.m[2][j] +
                        a.m[i][3] * b.m[3][j];
        }
    }
    return c;
}

YBI_INTEGRATOR_HD Vec4 operator*(const Float4x4 &mat, const Vec4 &p)
{
    return Vec4(Dot(Vec4(mat.m[0][0], mat.m[0][1], mat.m[0][2], mat.m[0][3]), p),
                Dot(Vec4(mat.m[1][0], mat.m[1][1], mat.m[1][2], mat.m[1][3]), p),
                Dot(Vec4(mat.m[2][0], mat.m[2][1], mat.m[2][2], mat.m[2][3]), p),
                Dot(Vec4(mat.m[3][0], mat.m[3][1], mat.m[3][2], mat.m[3][3]), p));
}

} // namespace ybi
