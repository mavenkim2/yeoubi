YBI_GPU bool YBI_MATERIAL_SELF_SAMPLE(const KernelParams &params,
                                      const ClosureEvaluationInfo &info,
                                      const Vec3f &viewDirection,
                                      BsdfSample &sample,
                                      Vec3f u)
{
    const OrthonormalBasis frame =
        MaterialXBsdfFrame(info, YBI_MATERIAL_NORMAL(info), viewDirection);
    Vec3f localLightDirection;
    if (!SampleOrenNayarBsdf(frame.WorldToLocal(viewDirection), u.Xy(), localLightDirection))
        return false;
    sample.direction = frame.LocalToWorld(localLightDirection);
    sample.eta = 1.f;
    sample.refractionDistanceScale = 1.f;
    sample.flags = BSDF_REFLECTION | BSDF_DIFFUSE;
    sample.sampledRoughness = 1.f;
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
    const Vec3f value = EvaluateOrenNayarBsdf(frame.WorldToLocal(viewDirection),
                                 frame.WorldToLocal(lightDirection),
                                 YBI_MATERIAL_WEIGHT(info),
                                 YBI_MATERIAL_COLOR(info),
                                 YBI_MATERIAL_ROUGHNESS(info));
    if (properties && value.MaxComponent() > 0.f)
    {
        properties->flags |= BSDF_REFLECTION | BSDF_DIFFUSE;
        properties->maxRoughness = Max(properties->maxRoughness, 1.f);
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
    return PDFOrenNayarBsdf(frame.WorldToLocal(viewDirection),
                            frame.WorldToLocal(lightDirection));
}

YBI_GPU float YBI_MATERIAL_SELF_ALBEDO(const KernelParams &params,
                                       const ClosureEvaluationInfo &info,
                                       const Vec3f &viewDirection)
{
    const OrthonormalBasis frame =
        MaterialXBsdfFrame(info, YBI_MATERIAL_NORMAL(info), viewDirection);
    return YBI_MATERIAL_WEIGHT(info) *
           OrenNayarDirectionalAlbedo(frame.WorldToLocal(viewDirection).z,
                                      YBI_MATERIAL_ROUGHNESS(info));
}

YBI_GPU float YBI_MATERIAL_SELF_TRANSMITTANCE(const KernelParams &params,
                                              const ClosureEvaluationInfo &info,
                                              const Vec3f &viewDirection)
{
    const float weight = YBI_MATERIAL_WEIGHT(info);
    const float transmittance = 1.f - weight;
    YBI_MATERIAL_DEBUG_PRINT_TRANSMITTANCE(weight, transmittance);
    return transmittance;
}
