#pragma once

#include "kernel/scene.h"

namespace ybi
{

YBI_GPU bool IntersectRectLight(const Vec3f &rayDirection,
                                const Vec3f &rayOrigin,
                                float tMin,
                                float tMax,
                                const DeviceLight &light,
                                float &hitT)
{
    Vec3f localRayOrigin = light.localFromWorld.TransformPoint(rayOrigin);
    Vec3f localRayDirection = light.localFromWorld.TransformVector(rayDirection);

    // ray = o + td
    // planed = 0
    // dot(n, ray) = 0
    // dot(n, o) + t * dot(n, d) = 0

    // Ray parallel to light, or traveling in same direction as plane normal (-Z)
    if (Abs(localRayDirection.z) < 1e-7f || localRayDirection.z <= 0.f)
    {
        return false;
    }

    float t = -localRayOrigin.z / localRayDirection.z;

    if (t < tMin || t > tMax)
    {
        return false;
    }

    const Vec3f p = localRayOrigin + t * localRayDirection;

    if (Abs(p.x) > light.width * 0.5f || Abs(p.y) > light.height * 0.5f)
    {
        return false;
    }

    hitT = t;
    return true;
}

YBI_GPU float
PDFRectLight(const DeviceLight &light, const Vec3f &hitPosition, const Vec3f &lightPosition)
{
    const Vec3f toLight = lightPosition - hitPosition;
    const float distanceSquared = Dot(toLight, toLight);
    if (distanceSquared <= 0.f)
        return 0.f;

    const Vec3f lightDirection = toLight / sqrtf(distanceSquared);

    const Vec3f worldNormal = Normalize(-Vec3f(
        light.localFromWorld.m[2][0], light.localFromWorld.m[2][1], light.localFromWorld.m[2][2]));

    const float cosLight = Dot(worldNormal, -lightDirection);
    if (cosLight <= 0.f)
        return 0.f;

    const Vec3f worldEdgeX = light.worldFromLocal.TransformVector(Vec3f(light.width, 0.f, 0.f));
    const Vec3f worldEdgeY = light.worldFromLocal.TransformVector(Vec3f(0.f, light.height, 0.f));
    const float area = Length(Cross(worldEdgeX, worldEdgeY));
    if (area <= 0.f)
        return 0.f;

    return distanceSquared / (cosLight * area);
}

YBI_GPU Vec3f SampleRectLight(const DeviceLight &light,
                              const Vec3f &hitPosition,
                              const Vec2f &u,
                              Vec3f &lightPosition,
                              float &pdf)
{
    const Vec3f localPosition((u.x - 0.5f) * light.width, (u.y - 0.5f) * light.height, 0.f);
    lightPosition = light.worldFromLocal.TransformPoint(localPosition);

    const Vec3f toLight = lightPosition - hitPosition;
    const float distanceSquared = Dot(toLight, toLight);
    if (distanceSquared <= 0.f)
    {
        pdf = 0.f;
        return Vec3f(0.f);
    }

    pdf = PDFRectLight(light, hitPosition, lightPosition);
    return toLight / sqrtf(distanceSquared);
}

} // namespace ybi
