YBI_GPU bool YBI_MATERIAL_SELF_SAMPLE(const KernelParams &params,
                                      const ClosureEvaluationInfo &info,
                                      const Vec3f &viewDirection,
                                      BsdfSample &sample,
                                      Vec3f u)
{
    const OrthonormalBasis frame =
        MaterialXBsdfFrame(info, YBI_MATERIAL_NORMAL(info), viewDirection);
    const float roughness = YBI_MATERIAL_ROUGHNESS(info);
    Vec3f localLightDirection;
    if (!SampleContyKullaSheen(frame.WorldToLocal(viewDirection),
                               u.Xy(),
                               roughness,
                               localLightDirection))
        return false;
    sample.direction = frame.LocalToWorld(localLightDirection);
    sample.eta = 1.f;
    sample.refractionDistanceScale = 1.f;
    sample.flags = BSDF_REFLECTION | BSDF_GLOSSY;
    sample.sampledRoughness = roughness;
    return true;
}

YBI_GPU Vec3f YBI_MATERIAL_SELF_EVALUATE(const KernelParams &params,
                                         const ClosureEvaluationInfo &info,
                                         const Vec3f &viewDirection,
                                         const Vec3f &lightDirection,
                                         BsdfEvaluationProperties *properties = nullptr)
{
    const OrthonormalBasis frame =
        MaterialXBsdfFrame(info, YBI_MATERIAL_NORMAL(info), viewDirection);
    const float roughness = YBI_MATERIAL_ROUGHNESS(info);
    const Vec3f value = EvaluateContyKullaSheen(frame.WorldToLocal(viewDirection),
                                   frame.WorldToLocal(lightDirection),
                                   YBI_MATERIAL_COLOR(info),
                                   YBI_MATERIAL_WEIGHT(info),
                                   roughness);
    if (properties && value.MaxComponent() > 0.f)
    {
        properties->flags |= BSDF_REFLECTION | BSDF_GLOSSY;
        properties->maxRoughness = Max(properties->maxRoughness, roughness);
    }
    return value;
}

YBI_GPU float YBI_MATERIAL_SELF_PDF(const KernelParams &params,
                                    const ClosureEvaluationInfo &info,
                                    const Vec3f &viewDirection,
                                    const Vec3f &lightDirection)
{
    const OrthonormalBasis frame =
        MaterialXBsdfFrame(info, YBI_MATERIAL_NORMAL(info), viewDirection);
    return PDFContyKullaSheen(frame.WorldToLocal(viewDirection),
                              frame.WorldToLocal(lightDirection),
                              YBI_MATERIAL_ROUGHNESS(info));
}

YBI_GPU float YBI_MATERIAL_SELF_ALBEDO(const KernelParams &params,
                                       const ClosureEvaluationInfo &info,
                                       const Vec3f &viewDirection)
{
    const OrthonormalBasis frame =
        MaterialXBsdfFrame(info, YBI_MATERIAL_NORMAL(info), viewDirection);
    return YBI_MATERIAL_WEIGHT(info) *
           SheenAlbedo(YBI_MATERIAL_ROUGHNESS(info), frame.WorldToLocal(viewDirection).z);
}

YBI_GPU float YBI_MATERIAL_SELF_TRANSMITTANCE(const KernelParams &params,
                                              const ClosureEvaluationInfo &info,
                                              const Vec3f &viewDirection)
{
    const float albedo = YBI_MATERIAL_SELF_ALBEDO(params, info, viewDirection);
    const float transmittance = 1.f - albedo;
    YBI_MATERIAL_DEBUG_PRINT_TRANSMITTANCE(YBI_MATERIAL_WEIGHT(info), albedo, transmittance);
    return transmittance;
}
