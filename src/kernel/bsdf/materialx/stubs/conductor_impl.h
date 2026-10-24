YBI_GPU bool YBI_MATERIAL_SELF_SAMPLE(const KernelParams &params,
                                      const ClosureEvaluationInfo &info,
                                      const Vec3f &viewDirection,
                                      BsdfSample &sample,
                                      Vec3f u)
{
    const OrthonormalBasis frame =
        MaterialXBsdfFrame(info, YBI_MATERIAL_NORMAL(info), viewDirection);
    const Vec3f localViewDirection = frame.WorldToLocal(viewDirection);
    const Vec2f roughness = YBI_MATERIAL_ROUGHNESS(info);
    if (IsSmoothMicrofacet(roughness))
    {
        if (localViewDirection.z <= 0.f)
            return false;
        sample.direction = frame.LocalToWorld(
            Vec3f(-localViewDirection.x, -localViewDirection.y, localViewDirection.z));
        sample.weight = YBI_MATERIAL_WEIGHT(info) *
                        FresnelConductor(Min(localViewDirection.z, 1.f),
                                         YBI_MATERIAL_IOR(info),
                                         YBI_MATERIAL_EXTINCTION(info));
        sample.pdf = 1.f;
        sample.eta = 1.f;
        sample.refractionDistanceScale = 1.f;
        sample.flags = BSDF_REFLECTION | BSDF_DELTA;
        sample.sampledRoughness = 0.f;
        return true;
    }
    Vec3f localLightDirection;
    if (!SampleConductorBsdf(localViewDirection,
                             u.Xy(),
                             roughness,
                             localLightDirection))
        return false;
    sample.direction = frame.LocalToWorld(localLightDirection);
    sample.eta = 1.f;
    sample.refractionDistanceScale = 1.f;
    sample.flags = BSDF_REFLECTION | BSDF_GLOSSY;
    sample.sampledRoughness = Min(roughness.x, roughness.y);
    return true;
}

YBI_GPU Vec3f YBI_MATERIAL_SELF_EVALUATE(const KernelParams &params,
                                         const ClosureEvaluationInfo &info,
                                         const Vec3f &viewDirection,
                                         const Vec3f &lightDirection,
                                         BsdfEvaluationProperties *properties = nullptr)
{
    const Vec2f roughness = YBI_MATERIAL_ROUGHNESS(info);
    if (IsSmoothMicrofacet(roughness))
        return Vec3f(0.f);
    const OrthonormalBasis frame =
        MaterialXBsdfFrame(info, YBI_MATERIAL_NORMAL(info), viewDirection);
    const Vec3f value = EvaluateConductorBsdf(frame.WorldToLocal(viewDirection),
                                 frame.WorldToLocal(lightDirection),
                                 YBI_MATERIAL_WEIGHT(info),
                                 YBI_MATERIAL_IOR(info),
                                 YBI_MATERIAL_EXTINCTION(info),
                                 roughness);
    if (properties && value.MaxComponent() > 0.f)
    {
        properties->flags |= BSDF_REFLECTION | BSDF_GLOSSY;
        properties->maxRoughness =
            Max(properties->maxRoughness, Max(roughness.x, roughness.y));
    }
    return value;
}

YBI_GPU float YBI_MATERIAL_SELF_PDF(const KernelParams &params,
                                    const ClosureEvaluationInfo &info,
                                    const Vec3f &viewDirection,
                                    const Vec3f &lightDirection)
{
    if (IsSmoothMicrofacet(YBI_MATERIAL_ROUGHNESS(info)))
        return 0.f;
    const OrthonormalBasis frame =
        MaterialXBsdfFrame(info, YBI_MATERIAL_NORMAL(info), viewDirection);
    return PDFConductorBsdf(frame.WorldToLocal(viewDirection),
                            frame.WorldToLocal(lightDirection),
                            YBI_MATERIAL_ROUGHNESS(info));
}

YBI_GPU float YBI_MATERIAL_SELF_ALBEDO(const KernelParams &params,
                                       const ClosureEvaluationInfo &info,
                                       const Vec3f &viewDirection)
{
    const OrthonormalBasis frame =
        MaterialXBsdfFrame(info, YBI_MATERIAL_NORMAL(info), viewDirection);
    ConductorBsdf bsdf;
    bsdf.weight = YBI_MATERIAL_WEIGHT(info);
    bsdf.ior = YBI_MATERIAL_IOR(info);
    bsdf.extinction = YBI_MATERIAL_EXTINCTION(info);
    bsdf.alpha = YBI_MATERIAL_ROUGHNESS(info);
    if (IsSmoothMicrofacet(bsdf.alpha))
        return (bsdf.weight * FresnelConductor(Min(frame.WorldToLocal(viewDirection).z, 1.f),
                                              bsdf.ior,
                                              bsdf.extinction)).MaxComponent();
    return ConductorDirectionalAlbedo(bsdf, frame.WorldToLocal(viewDirection).z).MaxComponent();
}

YBI_GPU float YBI_MATERIAL_SELF_TRANSMITTANCE(const KernelParams &params,
                                              const ClosureEvaluationInfo &info,
                                              const Vec3f &viewDirection)
{
    return 1.f - YBI_MATERIAL_WEIGHT(info);
}
