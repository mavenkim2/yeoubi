YBI_GPU bool YBI_MATERIAL_SELF_SAMPLE(const KernelParams &params,
                                      const ClosureEvaluationInfo &info,
                                      const Vec3f &viewDirection,
                                      BsdfSample &sample,
                                      Vec3f u)
{
    const OrthonormalBasis frame =
        MaterialXBsdfFrame(info, YBI_MATERIAL_NORMAL(info), viewDirection);
    Vec3f localLightDirection;
    if (!SampleContyKullaSheen(frame.WorldToLocal(viewDirection),
                               u.Xy(),
                               YBI_MATERIAL_ROUGHNESS(info),
                               localLightDirection))
        return false;
    sample.direction = frame.LocalToWorld(localLightDirection);
    sample.eta = 1.f;
    sample.flags = BSDF_REFLECTION;
    return true;
}

YBI_GPU Vec3f YBI_MATERIAL_SELF_EVALUATE(const KernelParams &params,
                                         const ClosureEvaluationInfo &info,
                                         const Vec3f &viewDirection,
                                         const Vec3f &lightDirection)
{
    const OrthonormalBasis frame =
        MaterialXBsdfFrame(info, YBI_MATERIAL_NORMAL(info), viewDirection);
    return EvaluateContyKullaSheen(frame.WorldToLocal(viewDirection),
                                   frame.WorldToLocal(lightDirection),
                                   YBI_MATERIAL_COLOR(info),
                                   YBI_MATERIAL_WEIGHT(info),
                                   YBI_MATERIAL_ROUGHNESS(info));
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
