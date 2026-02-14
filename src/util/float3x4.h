#pragma once

#include "util/float4.h"

YBI_NAMESPACE_BEGIN

// NOTE: row major, multiplication order right to left
struct float3x4
{
    union
    {
        float m[3][4];
        float4 rows[3];
    };

    float3x4() = default;
    float3x4(float m00,
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
    float3x4(float4 r0, float4 r1, float4 r2)
    {
        for (int c = 0; c < 4; c++)
        {
            m[0][c] = r0[c];
        }
        for (int c = 0; c < 4; c++)
        {
            m[1][c] = r1[c];
        }
        for (int c = 0; c < 4; c++)
        {
            m[2][c] = r2[c];
        }
    }
};

__forceinline float3 mul(const float3x4 &m, const float4 &p)
{
    return make_float3(dot(m.rows[0], p), dot(m.rows[1], p), dot(m.rows[2], p));
}

YBI_NAMESPACE_END
