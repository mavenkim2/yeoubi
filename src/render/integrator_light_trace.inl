YBI_DEVICE bool IntersectRectLight(const PackedLight &light,
                                   const Vec3 &rayOrigin,
                                   const Vec3 &rayDir,
                                   float tMin,
                                   float tMax,
                                   float pickPdf,
                                   LightRayHit *outHit)
{
    Vec3 tangent = {};
    Vec3 bitangent = {};
    Vec3 normal = {};
    GetLightFrame(light, &tangent, &bitangent, &normal);
    normal = -normal;

    const float denom = Dot(rayDir, normal);
    if (fabsf(denom) <= 1.0e-8f)
    {
        return false;
    }

    const Vec3 center = GetLightPosition(light);
    const float t = Dot(center - rayOrigin, normal) / denom;
    if (t <= tMin || t >= tMax)
    {
        return false;
    }

    const Vec3 hitPoint = rayOrigin + rayDir * t;
    const Vec3 rel = hitPoint - center;
    const float x = Dot(rel, tangent);
    const float y = Dot(rel, bitangent);
    if (fabsf(x) > light.width * 0.5f || fabsf(y) > light.height * 0.5f)
    {
        return false;
    }

    const float cosLight = Dot(normal, -rayDir);
    const float pdf = SolidAnglePdfFromArea(pickPdf, light.areaScale, t * t, cosLight);
    if (!outHit || pdf <= 0.0f)
    {
        return pdf > 0.0f;
    }

    outHit->radiance = LightEmission(light);
    outHit->distance = t;
    outHit->pdf = pdf;
    outHit->isDeltaLight = false;
    return true;
}

YBI_DEVICE bool IntersectDiskLight(const PackedLight &light,
                                   const Vec3 &rayOrigin,
                                   const Vec3 &rayDir,
                                   float tMin,
                                   float tMax,
                                   float pickPdf,
                                   LightRayHit *outHit)
{
    Vec3 tangent = {};
    Vec3 bitangent = {};
    Vec3 normal = {};
    GetLightFrame(light, &tangent, &bitangent, &normal);
    normal = -normal;

    const float denom = Dot(rayDir, normal);
    if (fabsf(denom) <= 1.0e-8f)
    {
        return false;
    }

    const Vec3 center = GetLightPosition(light);
    const float t = Dot(center - rayOrigin, normal) / denom;
    if (t <= tMin || t >= tMax)
    {
        return false;
    }

    const Vec3 hitPoint = rayOrigin + rayDir * t;
    const Vec3 rel = hitPoint - center;
    const float x = Dot(rel, tangent);
    const float y = Dot(rel, bitangent);
    if (x * x + y * y > light.radius * light.radius)
    {
        return false;
    }

    const float cosLight = Dot(normal, -rayDir);
    const float pdf = SolidAnglePdfFromArea(pickPdf, light.areaScale, t * t, cosLight);
    if (!outHit || pdf <= 0.0f)
    {
        return pdf > 0.0f;
    }

    outHit->radiance = LightEmission(light);
    outHit->distance = t;
    outHit->pdf = pdf;
    outHit->isDeltaLight = false;
    return true;
}

YBI_DEVICE bool IntersectSphereLight(const PackedLight &light,
                                     const Vec3 &rayOrigin,
                                     const Vec3 &rayDir,
                                     float tMin,
                                     float tMax,
                                     float pickPdf,
                                     LightRayHit *outHit)
{
    const Vec3 center = GetLightPosition(light);
    const Vec3 oc = rayOrigin - center;
    const float a = Dot(rayDir, rayDir);
    const float b = 2.0f * Dot(oc, rayDir);
    const float c = Dot(oc, oc) - light.radius * light.radius;
    const float discriminant = b * b - 4.0f * a * c;
    if (discriminant < 0.0f)
    {
        return false;
    }

    const float sqrtDiscriminant = sqrtf(discriminant);
    const float invDenom = 0.5f / MaxF(a, 1.0e-8f);
    float t = (-b - sqrtDiscriminant) * invDenom;
    if (t <= tMin || t >= tMax)
    {
        t = (-b + sqrtDiscriminant) * invDenom;
    }
    if (t <= tMin || t >= tMax)
    {
        return false;
    }

    const Vec3 hitPoint = rayOrigin + rayDir * t;
    const Vec3 normal = Normalize(hitPoint - center);
    const float cosLight = Dot(normal, -rayDir);
    const float pdf = SolidAnglePdfFromArea(pickPdf, light.areaScale, t * t, cosLight);
    if (!outHit || pdf <= 0.0f)
    {
        return pdf > 0.0f;
    }

    outHit->radiance = LightEmission(light);
    outHit->distance = t;
    outHit->pdf = pdf;
    outHit->isDeltaLight = false;
    return true;
}

YBI_DEVICE bool IntersectCylinderLight(const PackedLight &light,
                                       const Vec3 &rayOrigin,
                                       const Vec3 &rayDir,
                                       float tMin,
                                       float tMax,
                                       float pickPdf,
                                       LightRayHit *outHit)
{
    Vec3 tangent = {};
    Vec3 bitangent = {};
    Vec3 axis = {};
    GetLightFrame(light, &tangent, &bitangent, &axis);

    const Vec3 localOrigin = rayOrigin - GetLightPosition(light);
    const float ox = Dot(localOrigin, tangent);
    const float oy = Dot(localOrigin, bitangent);
    const float oz = Dot(localOrigin, axis);
    const float dx = Dot(rayDir, tangent);
    const float dy = Dot(rayDir, bitangent);
    const float dz = Dot(rayDir, axis);

    const float a = dx * dx + dy * dy;
    if (a <= 1.0e-8f)
    {
        return false;
    }

    const float b = 2.0f * (ox * dx + oy * dy);
    const float c = ox * ox + oy * oy - light.radius * light.radius;
    const float discriminant = b * b - 4.0f * a * c;
    if (discriminant < 0.0f)
    {
        return false;
    }

    const float sqrtDiscriminant = sqrtf(discriminant);
    const float invDenom = 0.5f / a;
    float t = (-b - sqrtDiscriminant) * invDenom;
    for (int candidate = 0; candidate < 2; ++candidate)
    {
        if (t > tMin && t < tMax)
        {
            const float z = oz + t * dz;
            if (fabsf(z) <= light.length * 0.5f)
            {
                const float x = ox + t * dx;
                const float y = oy + t * dy;
                const Vec3 normal = Normalize(tangent * x + bitangent * y);
                const float cosLight = Dot(normal, -rayDir);
                const float pdf = SolidAnglePdfFromArea(pickPdf, light.areaScale, t * t, cosLight);
                if (pdf <= 0.0f)
                {
                    return false;
                }
                if (outHit)
                {
                    outHit->radiance = LightEmission(light);
                    outHit->distance = t;
                    outHit->pdf = pdf;
                    outHit->isDeltaLight = false;
                }
                return true;
            }
        }
        t = (-b + sqrtDiscriminant) * invDenom;
    }

    return false;
}

YBI_DEVICE bool TraceAnalyticLight(const LaunchParams &params,
                                   const Vec3 &rayOrigin,
                                   const Vec3 &rayDir,
                                   float tMin,
                                   float tMax,
                                   LightRayHit *outHit)
{
    if (!outHit || params.lights == 0ull || params.lightCount <= 0)
    {
        return false;
    }

    const PackedLight *lights = reinterpret_cast<const PackedLight *>(params.lights);
    bool found = false;
    float bestDistance = tMax;

    for (int i = 0; i < params.lightCount; ++i)
    {
        const PackedLight &light = lights[i];
        if (light.type == static_cast<unsigned int>(LightType::Dome))
        {
            continue;
        }
        const float pickPdf = DirectLightPickPdf(params, i);
        if (pickPdf <= 0.0f)
        {
            continue;
        }

        LightRayHit candidate = {};
        candidate.lightIndex = i;
        bool hit = false;
        switch (light.type)
        {
            case static_cast<unsigned int>(LightType::Rect):
                hit = IntersectRectLight(
                    light, rayOrigin, rayDir, tMin, bestDistance, pickPdf, &candidate);
                break;
            case static_cast<unsigned int>(LightType::Disk):
                hit = IntersectDiskLight(
                    light, rayOrigin, rayDir, tMin, bestDistance, pickPdf, &candidate);
                break;
            case static_cast<unsigned int>(LightType::Sphere):
                hit = IntersectSphereLight(
                    light, rayOrigin, rayDir, tMin, bestDistance, pickPdf, &candidate);
                break;
            case static_cast<unsigned int>(LightType::Cylinder):
                hit = IntersectCylinderLight(
                    light, rayOrigin, rayDir, tMin, bestDistance, pickPdf, &candidate);
                break;
            default:
                break;
        }

        if (!hit)
        {
            continue;
        }

        found = true;
        bestDistance = candidate.distance;
        *outHit = candidate;
    }

    return found;
}

YBI_DEVICE bool TraceAnalyticLight(const LaunchParams &params,
                                   const Vec3 &rayOrigin,
                                   const Vec3 &rayDir,
                                   float tMax,
                                   LightRayHit *outHit)
{
    return TraceAnalyticLight(params, rayOrigin, rayDir, 1.0e-4f, tMax, outHit);
}
