#pragma once

#include "scene/device_scene.h"

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

// Carlos Ureña, Marcos Fajardo, and Alan King
// "An Area-Preserving Parametrization for Spherical Rectangles", Eurographics Symposium on
// Rendering, 2013.
//
#if 0
YBI_GPU void SampleRectLight(const DeviceLight &light, const Vec2f &u)
{
    squad.o = o;
    float exl = length(ex), eyl = length(ey);
    // compute local reference system ’R’
    squad.x = ex / exl;
    squad.y = ey / eyl;
    squad.z = cross(squad.x, squad.y);
    // compute rectangle coords in local reference system
    vec3 d = s - o;
    squad.z0 = dot(d, squad.z);
    // flip ’z’ to make it point against ’Q’
    if (squad.z0 > 0)
    {
        squad.z *= -1;
        squad.z0 *= -1;
    }
    squad.z0sq = squad.z0 * squad.z0;
    squad.x0 = dot(d, squad.x);
    squad.y0 = dot(d, squad.y);
    squad.x1 = squad.x0 + exl;
    squad.y1 = squad.y0 + eyl;
    squad.y0sq = squad.y0 * squad.y0;
    squad.y1sq = squad.y1 * squad.y1;
    // create vectors to four vertices
    vec3 v00 = {squad.x0, squad.y0, squad.z0};
    vec3 v01 = {squad.x0, squad.y1, squad.z0};
    vec3 v10 = {squad.x1, squad.y0, squad.z0};
    vec3 v11 = {squad.x1, squad.y1, squad.z0};
    // compute normals to edges
    vec3 n0 = normalize(cross(v00, v10));
    vec3 n1 = normalize(cross(v10, v11));
    vec3 n2 = normalize(cross(v11, v01));
    vec3 n3 = normalize(cross(v01, v00));
    // compute internal angles (gamma_i)
    float g0 = acos(-dot(n0, n1));
    float g1 = acos(-dot(n1, n2));
    float g2 = acos(-dot(n2, n3));
    float g3 = acos(-dot(n3, n0));
    // compute predefined constants
    squad.b0 = n0.z;
    squad.b1 = n2.z;
    squad.b0sq = squad.b0 * squad.b0;
    squad.k = 2 * PI - g2 - g3;
    // compute solid angle from internal angles
    squad.S = g0 + g1 - squad.k;

    float au = u * squad.S + squad.k;
    float fu = (cos(au) * squad.b0 - squad.b1) / sin(au);
    float cu = 1 / sqrt(fu * fu + squad.b0sq) * (fu > 0 ? +1 : -1);
    cu = clamp(cu, -1, 1);

    float xu = -(cu * squad.z0) / sqrt(1 - cu * cu);
    xu = clamp(xu, squad.x0, squad.x1); // avoid Infs
    // 3. compute ’yv’
    float d = sqrt(xu * xu + squad.z0sq);
    float h0 = squad.y0 / sqrt(d * d + squad.y0sq);
    float h1 = squad.y1 / sqrt(d * d + squad.y1sq);
    float hv = h0 + v * (h1 - h0), hv2 = hv * hv;
    float yv = (hv2 < 1 - eps) ? (hv * d) / sqrt(1 - hv2) : squad.y1;
    // 4. transform (xu,yv,z0) to world coords
    return (squad.o + xu * squad.x + yv * squad.y + z0 * squad.z);
}
#endif

} // namespace ybi
