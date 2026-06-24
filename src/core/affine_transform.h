#pragma once

#include "core/assert.h"
#include "core/vec3f.h"
#include "device/cuda_compat.h"
#include <cmath>
#include <cstddef>

namespace ybi
{

class AffineTransform
{
public:
    float m[3][4];

    AffineTransform() = default;

    YBI_CPU_GPU constexpr AffineTransform(float m00,
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
        : m{{m00, m01, m02, m03}, {m10, m11, m12, m13}, {m20, m21, m22, m23}}
    {
    }

    YBI_CPU_GPU explicit AffineTransform(const float values[3][4])
    {
        for (size_t row = 0; row < 3; ++row)
        {
            for (size_t column = 0; column < 4; ++column)
            {
                m[row][column] = values[row][column];
            }
        }
    }

    YBI_CPU_GPU static constexpr AffineTransform Identity()
    {
        return AffineTransform(1.f, 0.f, 0.f, 0.f, 0.f, 1.f, 0.f, 0.f, 0.f, 0.f, 1.f, 0.f);
    }

    YBI_CPU_GPU static constexpr AffineTransform Zero()
    {
        return AffineTransform(
            0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f);
    }

    YBI_CPU_GPU static constexpr AffineTransform Translation(float x, float y, float z)
    {
        return AffineTransform(1.0f, 0.0f, 0.0f, x, 0.0f, 1.0f, 0.0f, y, 0.0f, 0.0f, 1.0f, z);
    }

    YBI_CPU_GPU static constexpr AffineTransform Translation(Vec3f translation)
    {
        return Translation(translation.x, translation.y, translation.z);
    }

    YBI_CPU_GPU static constexpr AffineTransform Scale(float x, float y, float z)
    {
        return AffineTransform(x, 0.0f, 0.0f, 0.0f, 0.0f, y, 0.0f, 0.0f, 0.0f, 0.0f, z, 0.0f);
    }

    YBI_CPU_GPU static constexpr AffineTransform Scale(float value)
    {
        return Scale(value, value, value);
    }

    YBI_CPU_GPU float *Data()
    {
        return &m[0][0];
    }

    YBI_CPU_GPU const float *Data() const
    {
        return &m[0][0];
    }

    YBI_CPU_GPU float *operator[](size_t row)
    {
        YBI_ASSERT(row < 3);
        return m[row];
    }

    YBI_CPU_GPU const float *operator[](size_t row) const
    {
        YBI_ASSERT(row < 3);
        return m[row];
    }

    YBI_CPU_GPU float &operator()(size_t row, size_t column)
    {
        YBI_ASSERT(row < 3);
        YBI_ASSERT(column < 4);
        return m[row][column];
    }

    YBI_CPU_GPU const float &operator()(size_t row, size_t column) const
    {
        YBI_ASSERT(row < 3);
        YBI_ASSERT(column < 4);
        return m[row][column];
    }

    YBI_CPU_GPU bool operator==(const AffineTransform &other) const
    {
        for (int r = 0; r < 3; r++)
        {
            for (int c = 0; c < 4; c++)
            {
                if (m[r][c] != other.m[r][c])
                    return false;
            }
        }
        return true;
    }

    YBI_CPU_GPU constexpr Vec3f GetTranslation() const
    {
        return Vec3f(m[0][3], m[1][3], m[2][3]);
    }

    YBI_CPU_GPU void SetTranslation(Vec3f translation)
    {
        m[0][3] = translation.x;
        m[1][3] = translation.y;
        m[2][3] = translation.z;
    }

    YBI_CPU_GPU constexpr Vec3f TransformPoint(Vec3f point) const
    {
        return Vec3f(m[0][0] * point.x + m[0][1] * point.y + m[0][2] * point.z + m[0][3],
                     m[1][0] * point.x + m[1][1] * point.y + m[1][2] * point.z + m[1][3],
                     m[2][0] * point.x + m[2][1] * point.y + m[2][2] * point.z + m[2][3]);
    }

    YBI_CPU_GPU constexpr Vec3f TransformVector(Vec3f vector) const
    {
        return Vec3f(m[0][0] * vector.x + m[0][1] * vector.y + m[0][2] * vector.z,
                     m[1][0] * vector.x + m[1][1] * vector.y + m[1][2] * vector.z,
                     m[2][0] * vector.x + m[2][1] * vector.y + m[2][2] * vector.z);
    }

    YBI_CPU_GPU bool TryInverse(AffineTransform *out, float epsilon = 1.0e-8f) const
    {
        YBI_ASSERT(out);

        const float a = m[0][0];
        const float b = m[0][1];
        const float c = m[0][2];
        const float d = m[1][0];
        const float e = m[1][1];
        const float f = m[1][2];
        const float g = m[2][0];
        const float h = m[2][1];
        const float i = m[2][2];

        const float det = a * (e * i - f * h) - b * (d * i - f * g) + c * (d * h - e * g);
        if (std::fabs(det) <= epsilon)
        {
            return false;
        }

        const float invDet = 1.0f / det;

        AffineTransform inverse = AffineTransform::Identity();
        inverse.m[0][0] = (e * i - f * h) * invDet;
        inverse.m[0][1] = (c * h - b * i) * invDet;
        inverse.m[0][2] = (b * f - c * e) * invDet;
        inverse.m[1][0] = (f * g - d * i) * invDet;
        inverse.m[1][1] = (a * i - c * g) * invDet;
        inverse.m[1][2] = (c * d - a * f) * invDet;
        inverse.m[2][0] = (d * h - e * g) * invDet;
        inverse.m[2][1] = (b * g - a * h) * invDet;
        inverse.m[2][2] = (a * e - b * d) * invDet;

        const Vec3f translation = GetTranslation();
        const Vec3f inverseTranslation =
            inverse.TransformVector(Vec3f(-translation.x, -translation.y, -translation.z));
        inverse.SetTranslation(inverseTranslation);

        *out = inverse;
        return true;
    }

    YBI_CPU_GPU AffineTransform Inverse(float epsilon = 1.0e-8f) const
    {
        AffineTransform inverse;
        const bool success = TryInverse(&inverse, epsilon);
        YBI_ASSERT(success);
        return inverse;
    }
};

YBI_CPU_GPU constexpr AffineTransform operator*(const AffineTransform &a, const AffineTransform &b)
{
    return AffineTransform(
        a.m[0][0] * b.m[0][0] + a.m[0][1] * b.m[1][0] + a.m[0][2] * b.m[2][0],
        a.m[0][0] * b.m[0][1] + a.m[0][1] * b.m[1][1] + a.m[0][2] * b.m[2][1],
        a.m[0][0] * b.m[0][2] + a.m[0][1] * b.m[1][2] + a.m[0][2] * b.m[2][2],
        a.m[0][0] * b.m[0][3] + a.m[0][1] * b.m[1][3] + a.m[0][2] * b.m[2][3] + a.m[0][3],
        a.m[1][0] * b.m[0][0] + a.m[1][1] * b.m[1][0] + a.m[1][2] * b.m[2][0],
        a.m[1][0] * b.m[0][1] + a.m[1][1] * b.m[1][1] + a.m[1][2] * b.m[2][1],
        a.m[1][0] * b.m[0][2] + a.m[1][1] * b.m[1][2] + a.m[1][2] * b.m[2][2],
        a.m[1][0] * b.m[0][3] + a.m[1][1] * b.m[1][3] + a.m[1][2] * b.m[2][3] + a.m[1][3],
        a.m[2][0] * b.m[0][0] + a.m[2][1] * b.m[1][0] + a.m[2][2] * b.m[2][0],
        a.m[2][0] * b.m[0][1] + a.m[2][1] * b.m[1][1] + a.m[2][2] * b.m[2][1],
        a.m[2][0] * b.m[0][2] + a.m[2][1] * b.m[1][2] + a.m[2][2] * b.m[2][2],
        a.m[2][0] * b.m[0][3] + a.m[2][1] * b.m[1][3] + a.m[2][2] * b.m[2][3] + a.m[2][3]);
}

} // namespace ybi
