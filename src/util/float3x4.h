#pragma once

#include "util/vec3.h"
#include "util/vec4.h"

namespace ybi
{

// NOTE: row major, multiplication order right to left
struct Float3x4
{
    float m[3][4];

    YBI_DEVICE Float3x4()
        : m{{0.0f, 0.0f, 0.0f, 0.0f}, {0.0f, 0.0f, 0.0f, 0.0f}, {0.0f, 0.0f, 0.0f, 0.0f}}
    {
    }
    YBI_DEVICE Float3x4(float m00,
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

    YBI_DEVICE Float3x4(Vec4 r0, Vec4 r1, Vec4 r2)
    {
        for (int c = 0; c < 4; ++c)
        {
            m[0][c] = r0[c];
            m[1][c] = r1[c];
            m[2][c] = r2[c];
        }
    }

    YBI_DEVICE static Float3x4 Identity()
    {
        return Float3x4(1.0f, 0.0f, 0.0f, 0.0f, 0.0f, 1.0f, 0.0f, 0.0f, 0.0f, 0.0f, 1.0f, 0.0f);
    }
};

using float3x4 = Float3x4;

YBI_DEVICE Vec3 operator*(const Float3x4 &m, const Vec4 &p)
{
    return Vec3(Dot(Vec4(m.m[0][0], m.m[0][1], m.m[0][2], m.m[0][3]), p),
                Dot(Vec4(m.m[1][0], m.m[1][1], m.m[1][2], m.m[1][3]), p),
                Dot(Vec4(m.m[2][0], m.m[2][1], m.m[2][2], m.m[2][3]), p));
}

YBI_DEVICE Vec3 TransformPointAffine(const Float3x4 &mat, const Vec3 &p)
{
    return mat * Vec4(p, 1.0f);
}

YBI_DEVICE Vec3 TransformVectorAffine(const Float3x4 &mat, const Vec3 &v)
{
    return Vec3(mat.m[0][0] * v.x + mat.m[0][1] * v.y + mat.m[0][2] * v.z,
                mat.m[1][0] * v.x + mat.m[1][1] * v.y + mat.m[1][2] * v.z,
                mat.m[2][0] * v.x + mat.m[2][1] * v.y + mat.m[2][2] * v.z);
}

} // namespace ybi
