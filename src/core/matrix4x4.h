#pragma once

#include "core/assert.h"
#include "core/vec3f.h"
#include "core/vec4f.h"
#include "device/cuda_compat.h"
#include <cmath>
#include <cstddef>

namespace ybi
{

class Matrix4x4
{
public:
    float m[4][4];

    Matrix4x4() = default;

    YBI_CPU_GPU constexpr Matrix4x4(float m00,
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

    YBI_CPU_GPU explicit Matrix4x4(const float values[4][4])
    {
        for (size_t row = 0; row < 4; ++row)
        {
            for (size_t column = 0; column < 4; ++column)
            {
                m[row][column] = values[row][column];
            }
        }
    }

    YBI_CPU_GPU static constexpr Matrix4x4 Identity()
    {
        return Matrix4x4(
            1.f, 0.f, 0.f, 0.f, 0.f, 1.f, 0.f, 0.f, 0.f, 0.f, 1.f, 0.f, 0.f, 0.f, 0.f, 1.f);
    }

    YBI_CPU_GPU static constexpr Matrix4x4 Zero()
    {
        return Matrix4x4(0.0f,
                         0.0f,
                         0.0f,
                         0.0f,
                         0.0f,
                         0.0f,
                         0.0f,
                         0.0f,
                         0.0f,
                         0.0f,
                         0.0f,
                         0.0f,
                         0.0f,
                         0.0f,
                         0.0f,
                         0.0f);
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
        YBI_ASSERT(row < 4);
        return m[row];
    }

    YBI_CPU_GPU const float *operator[](size_t row) const
    {
        YBI_ASSERT(row < 4);
        return m[row];
    }

    YBI_CPU_GPU float &operator()(size_t row, size_t column)
    {
        YBI_ASSERT(row < 4);
        YBI_ASSERT(column < 4);
        return m[row][column];
    }

    YBI_CPU_GPU const float &operator()(size_t row, size_t column) const
    {
        YBI_ASSERT(row < 4);
        YBI_ASSERT(column < 4);
        return m[row][column];
    }

    YBI_CPU_GPU bool operator==(const Matrix4x4 &other) const
    {
        for (int row = 0; row < 4; ++row)
        {
            for (int column = 0; column < 4; ++column)
            {
                if (m[row][column] != other.m[row][column])
                {
                    return false;
                }
            }
        }
        return true;
    }

    YBI_CPU_GPU constexpr Vec4f Transform(Vec4f vector) const
    {
        return Vec4f(
            m[0][0] * vector.x + m[0][1] * vector.y + m[0][2] * vector.z + m[0][3] * vector.w,
            m[1][0] * vector.x + m[1][1] * vector.y + m[1][2] * vector.z + m[1][3] * vector.w,
            m[2][0] * vector.x + m[2][1] * vector.y + m[2][2] * vector.z + m[2][3] * vector.w,
            m[3][0] * vector.x + m[3][1] * vector.y + m[3][2] * vector.z + m[3][3] * vector.w);
    }

    YBI_CPU_GPU constexpr Matrix4x4 Transpose() const
    {
        return Matrix4x4(m[0][0],
                         m[1][0],
                         m[2][0],
                         m[3][0],
                         m[0][1],
                         m[1][1],
                         m[2][1],
                         m[3][1],
                         m[0][2],
                         m[1][2],
                         m[2][2],
                         m[3][2],
                         m[0][3],
                         m[1][3],
                         m[2][3],
                         m[3][3]);
    }

    YBI_CPU_GPU bool TryInverse(Matrix4x4 *out, float epsilon = 1.0e-8f) const
    {
        YBI_ASSERT(out);

        float augmented[4][8];
        for (size_t row = 0; row < 4; ++row)
        {
            for (size_t column = 0; column < 4; ++column)
            {
                augmented[row][column] = m[row][column];
                augmented[row][column + 4] = row == column ? 1.0f : 0.0f;
            }
        }

        for (size_t column = 0; column < 4; ++column)
        {
            size_t pivotRow = column;
            float pivotMagnitude = std::fabs(augmented[pivotRow][column]);
            for (size_t row = column + 1; row < 4; ++row)
            {
                const float magnitude = std::fabs(augmented[row][column]);
                if (magnitude > pivotMagnitude)
                {
                    pivotRow = row;
                    pivotMagnitude = magnitude;
                }
            }

            if (pivotMagnitude <= epsilon)
            {
                return false;
            }

            if (pivotRow != column)
            {
                for (size_t index = 0; index < 8; ++index)
                {
                    const float value = augmented[column][index];
                    augmented[column][index] = augmented[pivotRow][index];
                    augmented[pivotRow][index] = value;
                }
            }

            const float inversePivot = 1.0f / augmented[column][column];
            for (size_t index = 0; index < 8; ++index)
            {
                augmented[column][index] *= inversePivot;
            }

            for (size_t row = 0; row < 4; ++row)
            {
                if (row == column)
                {
                    continue;
                }

                const float factor = augmented[row][column];
                for (size_t index = 0; index < 8; ++index)
                {
                    augmented[row][index] -= factor * augmented[column][index];
                }
            }
        }

        Matrix4x4 inverse = Matrix4x4::Zero();
        for (size_t row = 0; row < 4; ++row)
        {
            for (size_t column = 0; column < 4; ++column)
            {
                inverse.m[row][column] = augmented[row][column + 4];
            }
        }

        *out = inverse;
        return true;
    }

    YBI_CPU_GPU Matrix4x4 Inverse(float epsilon = 1.0e-8f) const
    {
        Matrix4x4 inverse;
        const bool success = TryInverse(&inverse, epsilon);
        YBI_ASSERT(success);
        return inverse;
    }

    static YBI_CPU_GPU Matrix4x4 LookAt(const Vec3f &to,
                                        const Vec3f &worldUp,
                                        const Vec3f &cameraP)
    {
        Vec3f forward = Normalize(cameraP - to);
        Vec3f right = Normalize(Cross(worldUp, forward));
        Vec3f up = Normalize(Cross(forward, right));

        Matrix4x4 result(right.x,
                         right.y,
                         right.z,
                         -Dot(right, cameraP),
                         up.x,
                         up.y,
                         up.z,
                         -Dot(up, cameraP),
                         forward.x,
                         forward.y,
                         forward.z,
                         -Dot(forward, cameraP),
                         0.f,
                         0.f,
                         0.f,
                         1.f);
        return result;
    }
};

YBI_CPU_GPU constexpr Matrix4x4 operator*(const Matrix4x4 &a, const Matrix4x4 &b)
{
    Matrix4x4 result = Matrix4x4::Zero();
    for (size_t row = 0; row < 4; ++row)
    {
        for (size_t column = 0; column < 4; ++column)
        {
            for (size_t index = 0; index < 4; ++index)
            {
                result.m[row][column] += a.m[row][index] * b.m[index][column];
            }
        }
    }
    return result;
}

YBI_CPU_GPU constexpr Vec4f operator*(const Matrix4x4 &matrix, Vec4f vector)
{
    return matrix.Transform(vector);
}

YBI_CPU_GPU Matrix4x4 Perspective(float fov, float aspectRatio, float n = 1e-2f, float f = 1000.f)
{
    float divTanHalf = 1.f / tanf(fov / 2.f);
    Matrix4x4 result(divTanHalf / aspectRatio,
                     0.f,
                     0.f,
                     0.f,
                     0.f,
                     divTanHalf,
                     0.f,
                     0.f,
                     0.f,
                     0.f,
                     f / (n - f),
                     n * f / (n - f),
                     0.f,
                     0.f,
                     -1.f,
                     0.f);
    return result;
}

} // namespace ybi
