#pragma once

#include "core/math.h"
#include "core/vec2f.h"
#include "core/vec3f.h"
#include "device/cuda_compat.h"
#include "kernel/light/types.h"
#include "kernel/sampling/cdf.h"
#include "kernel/texture/sampling.h"

namespace ybi
{

YBI_CPU_GPU Vec2f EnvmapDirectionToUV(const Vec3f &direction)
{
    float u = atan2f(direction.x, -direction.z) / (2.f * YBI_PI);
    u = Min(u - Floor(u), FLOAT_BEFORE_ONE);
    const float sinTheta = Sqrt(direction.x * direction.x + direction.z * direction.z);
    const float v = atan2f(sinTheta, direction.y) / YBI_PI;
    return Vec2f(u, Clamp(v, 0.f, 1.f));
}

YBI_CPU_GPU Vec3f EnvmapUVToDirection(const Vec2f &uv)
{
    const float theta = YBI_PI * uv.y;
    const float phi = 2.f * YBI_PI * uv.x;
    const float sinTheta = Sin(theta);
    return Vec3f(sinTheta * Sin(phi), Cos(theta), -sinTheta * Cos(phi));
}

YBI_GPU Vec3f LookupEnvmapRadiance(const DeviceEnvmapLight &light, const Vec2f &uv)
{
    return SampleDeviceTexture(light.textureObject, light.textureFormat, uv.x, uv.y, 0).Xyz();
}

YBI_GPU Vec3f EvaluateEnvmapLight(const DeviceEnvmapLight &light, const Vec3f &worldDirection)
{
    const Vec3f localDirection = light.localFromWorld.TransformVector(worldDirection);
    return LookupEnvmapRadiance(light, EnvmapDirectionToUV(localDirection));
}

YBI_CPU_GPU float PDFEnvmapLight(const DeviceEnvmapLight &light, const Vec3f &worldDirection)
{
    Vec3f localDirection = light.localFromWorld.TransformVector(worldDirection);
    const float lengthSquared = Dot(localDirection, localDirection);
    localDirection /= Sqrt(lengthSquared);

    const float sinTheta =
        Sqrt(localDirection.x * localDirection.x + localDirection.z * localDirection.z);
    if (sinTheta == 0.f)
        return 0.f;

    Vec2f uv = EnvmapDirectionToUV(localDirection);
    uv.y = Min(uv.y, FLOAT_BEFORE_ONE);
    const float pdfUV =
        PDFPiecewise2DCDF(light.conditionalCDF, light.marginalCDF, light.width, light.height, uv);
    return pdfUV / (2.f * YBI_PI * YBI_PI * sinTheta);
}

YBI_GPU bool SampleEnvmapLight(const DeviceEnvmapLight &light,
                               const Vec2f &u,
                               Vec3f &lightDirection,
                               Vec3f &radiance,
                               float &pdf)
{
    lightDirection = Vec3f(0.f);
    radiance = Vec3f(0.f);
    pdf = 0.f;
    Vec2f uv;
    float pdfUV;
    if (!SamplePiecewise2DCDF(
            light.conditionalCDF, light.marginalCDF, light.width, light.height, u, uv, pdfUV))
        return false;

    const float sinTheta = Sin(YBI_PI * uv.y);
    if (sinTheta == 0.f)
        return false;

    const Vec3f worldDirection = light.worldFromLocal.TransformVector(EnvmapUVToDirection(uv));
    const float lengthSquared = Dot(worldDirection, worldDirection);
    lightDirection = worldDirection / Sqrt(lengthSquared);
    radiance = LookupEnvmapRadiance(light, uv);
    pdf = pdfUV / (2.f * YBI_PI * YBI_PI * sinTheta);
    return pdf > 0.f;
}

} // namespace ybi
