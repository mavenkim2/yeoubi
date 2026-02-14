#pragma once

#include "util/float4.h"

YBI_NAMESPACE_BEGIN

// Row-major 4x4. Matrix multiplication is right-to-left (see .cursor/rules).
struct float4x4
{
    union
    {
        float m[4][4];
        float4 rows[4];
    };

    float4x4() = default;
    float4x4(float m00,
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
    float4x4(float4 r0, float4 r1, float4 r2, float4 r3)
    {
        for (int c = 0; c < 4; c++)
            rows[0][c] = r0[c];
        for (int c = 0; c < 4; c++)
            rows[1][c] = r1[c];
        for (int c = 0; c < 4; c++)
            rows[2][c] = r2[c];
        for (int c = 0; c < 4; c++)
            rows[3][c] = r3[c];
    }

    static float4x4 Identity()
    {
        return float4x4(1, 0, 0, 0, 0, 1, 0, 0, 0, 0, 1, 0, 0, 0, 0, 1);
    }
};

__forceinline float4x4 mul(const float4x4 &a, const float4x4 &b)
{
    float4x4 c;
    for (int i = 0; i < 4; i++)
    {
        for (int j = 0; j < 4; j++)
        {
            c.m[i][j] = a.m[i][0] * b.m[0][j] + a.m[i][1] * b.m[1][j] + a.m[i][2] * b.m[2][j] +
                        a.m[i][3] * b.m[3][j];
        }
    }
    return c;
}

__forceinline float4 mul(const float4x4 &mat, const float4 &p)
{
    return make_float4(
        dot(mat.rows[0], p), dot(mat.rows[1], p), dot(mat.rows[2], p), dot(mat.rows[3], p));
}

YBI_NAMESPACE_END
