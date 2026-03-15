#pragma once

#include "util/math_constants.h"
#include "util/float3x4.h"
#include <cmath>

namespace ybi
{

// Row-major 4x4. Matrix multiplication is right-to-left.
struct Float4x4
{
    float m[4][4];

    Float4x4() = default;
    YBI_DEVICE constexpr Float4x4(float m00,
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
        : m{{m00, m01, m02, m03}, {m10, m11, m12, m13}, {m20, m21, m22, m23}, {m30, m31, m32, m33}}
    {
    }

    YBI_DEVICE Float4x4(Vec4 r0, Vec4 r1, Vec4 r2, Vec4 r3)
    {
        for (int c = 0; c < 4; ++c)
        {
            m[0][c] = r0[c];
            m[1][c] = r1[c];
            m[2][c] = r2[c];
            m[3][c] = r3[c];
        }
    }

    YBI_DEVICE static constexpr Float4x4 Identity()
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

YBI_DEVICE Float4x4 operator*(const Float4x4 &a, const Float4x4 &b)
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

YBI_DEVICE Vec4 operator*(const Float4x4 &mat, const Vec4 &p)
{
    return Vec4(Dot(Vec4(mat.m[0][0], mat.m[0][1], mat.m[0][2], mat.m[0][3]), p),
                Dot(Vec4(mat.m[1][0], mat.m[1][1], mat.m[1][2], mat.m[1][3]), p),
                Dot(Vec4(mat.m[2][0], mat.m[2][1], mat.m[2][2], mat.m[2][3]), p),
                Dot(Vec4(mat.m[3][0], mat.m[3][1], mat.m[3][2], mat.m[3][3]), p));
}

YBI_DEVICE Vec3 TransformPointAffine(const Float4x4 &mat, const Vec3 &p)
{
    const Vec4 hp = mat * Vec4(p, 1.0f);
    return Vec3(hp.x, hp.y, hp.z);
}

YBI_DEVICE Vec3 TransformVectorAffine(const Float4x4 &mat, const Vec3 &v)
{
    return Vec3(mat.m[0][0] * v.x + mat.m[0][1] * v.y + mat.m[0][2] * v.z,
                mat.m[1][0] * v.x + mat.m[1][1] * v.y + mat.m[1][2] * v.z,
                mat.m[2][0] * v.x + mat.m[2][1] * v.y + mat.m[2][2] * v.z);
}

YBI_DEVICE Vec3 TransformPointPerspective(const Float4x4 &mat, const Vec3 &p)
{
    const Vec4 hp = mat * Vec4(p, 1.0f);
    if (fabsf(hp.w) <= 1.0e-20f)
    {
        return Vec3(NAN, NAN, NAN);
    }

    const float invW = 1.0f / hp.w;
    return Vec3(hp.x * invW, hp.y * invW, hp.z * invW);
}

YBI_DEVICE bool Invert(const Float4x4 &mat, Float4x4 *out)
{
    if (!out)
    {
        return false;
    }

    float aug[4][8] = {};
    for (int row = 0; row < 4; ++row)
    {
        for (int col = 0; col < 4; ++col)
        {
            aug[row][col] = mat.m[row][col];
            aug[row][4 + col] = row == col ? 1.0f : 0.0f;
        }
    }

    for (int pivotCol = 0; pivotCol < 4; ++pivotCol)
    {
        int pivotRow = pivotCol;
        float pivotAbs = fabsf(aug[pivotRow][pivotCol]);
        for (int row = pivotCol + 1; row < 4; ++row)
        {
            const float candidateAbs = fabsf(aug[row][pivotCol]);
            if (candidateAbs > pivotAbs)
            {
                pivotAbs = candidateAbs;
                pivotRow = row;
            }
        }

        if (pivotAbs <= 1.0e-20f)
        {
            return false;
        }

        if (pivotRow != pivotCol)
        {
            for (int col = 0; col < 8; ++col)
            {
                const float tmp = aug[pivotCol][col];
                aug[pivotCol][col] = aug[pivotRow][col];
                aug[pivotRow][col] = tmp;
            }
        }

        const float pivotInv = 1.0f / aug[pivotCol][pivotCol];
        for (int col = 0; col < 8; ++col)
        {
            aug[pivotCol][col] *= pivotInv;
        }

        for (int row = 0; row < 4; ++row)
        {
            if (row == pivotCol)
            {
                continue;
            }

            const float factor = aug[row][pivotCol];
            if (factor == 0.0f)
            {
                continue;
            }

            for (int col = 0; col < 8; ++col)
            {
                aug[row][col] -= factor * aug[pivotCol][col];
            }
        }
    }

    for (int row = 0; row < 4; ++row)
    {
        for (int col = 0; col < 4; ++col)
        {
            out->m[row][col] = aug[row][4 + col];
        }
    }
    return true;
}

YBI_DEVICE Float4x4 Transpose(const Float4x4 &mat)
{
    return Float4x4(mat.m[0][0],
                    mat.m[1][0],
                    mat.m[2][0],
                    mat.m[3][0],
                    mat.m[0][1],
                    mat.m[1][1],
                    mat.m[2][1],
                    mat.m[3][1],
                    mat.m[0][2],
                    mat.m[1][2],
                    mat.m[2][2],
                    mat.m[3][2],
                    mat.m[0][3],
                    mat.m[1][3],
                    mat.m[2][3],
                    mat.m[3][3]);
}

YBI_DEVICE Float4x4 BuildCameraFromWorld(const Vec3 &eye, const Vec3 &lookAt)
{
    Vec3 forward = Normalize(lookAt - eye);
    if (Length(forward) <= 1.0e-8f)
    {
        forward = Vec3(0.0f, 0.0f, 1.0f);
    }
    Vec3 worldUp = Vec3(0.0f, 0.0f, 1.0f);
    if (fabsf(Dot(forward, worldUp)) > 0.999f)
    {
        worldUp = Vec3(0.0f, 1.0f, 0.0f);
    }
    const Vec3 right = Normalize(Cross(forward, worldUp));
    const Vec3 up = Normalize(Cross(right, forward));
    return Float4x4(right.x,
                    right.y,
                    right.z,
                    -Dot(right, eye),
                    up.x,
                    up.y,
                    up.z,
                    -Dot(up, eye),
                    forward.x,
                    forward.y,
                    forward.z,
                    -Dot(forward, eye),
                    0.0f,
                    0.0f,
                    0.0f,
                    1.0f);
}

YBI_DEVICE Float4x4 BuildPerspectiveClipFromCamera(float verticalFovDegrees,
                                                          int viewportWidth,
                                                          int viewportHeight,
                                                          float nearPlane = 1.0f,
                                                          float farPlane = 1.0e6f)
{
    const float fovY = verticalFovDegrees * kDegToRad;
    const float tanHalfFovY = fmaxf(1.0e-8f, tanf(0.5f * fovY));
    const float safeWidth = float(viewportWidth > 0 ? viewportWidth : 1);
    const float safeHeight = float(viewportHeight > 0 ? viewportHeight : 1);
    const float aspect = fmaxf(1.0e-8f, safeWidth / safeHeight);
    const float m00 = 1.0f / (tanHalfFovY * aspect);
    const float m11 = 1.0f / tanHalfFovY;
    const float m22 = (farPlane + nearPlane) / (farPlane - nearPlane);
    const float m23 = (-2.0f * farPlane * nearPlane) / (farPlane - nearPlane);
    return Float4x4(m00,
                    0.0f,
                    0.0f,
                    0.0f,
                    0.0f,
                    m11,
                    0.0f,
                    0.0f,
                    0.0f,
                    0.0f,
                    m22,
                    m23,
                    0.0f,
                    0.0f,
                    1.0f,
                    0.0f);
}

} // namespace ybi
