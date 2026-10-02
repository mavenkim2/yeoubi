YBI_GPU bool YBI_MATERIAL_SELF_SAMPLE(const KernelParams &params,
                                      const ClosureEvaluationInfo &info,
                                      const Vec3f &viewDirection,
                                      Vec3f &lightDirection,
                                      Vec3f u)
{
    const OrthonormalBasis frame =
        MaterialXBsdfFrame(info, YBI_MATERIAL_NORMAL(info), viewDirection);
    Vec3f localLightDirection;
    if (!SampleConductorBsdf(frame.WorldToLocal(viewDirection),
                             u.Xy(),
                             YBI_MATERIAL_ROUGHNESS(info),
                             localLightDirection))
        return false;
    lightDirection = frame.LocalToWorld(localLightDirection);
    return true;
}

YBI_GPU Vec3f YBI_MATERIAL_SELF_EVALUATE(const KernelParams &params,
                                         const ClosureEvaluationInfo &info,
                                         const Vec3f &viewDirection,
                                         const Vec3f &lightDirection)
{
    const OrthonormalBasis frame =
        MaterialXBsdfFrame(info, YBI_MATERIAL_NORMAL(info), viewDirection);
    return EvaluateConductorBsdf(frame.WorldToLocal(viewDirection),
                                 frame.WorldToLocal(lightDirection),
                                 YBI_MATERIAL_WEIGHT(info),
                                 YBI_MATERIAL_IOR(info),
                                 YBI_MATERIAL_EXTINCTION(info),
                                 YBI_MATERIAL_ROUGHNESS(info));
}

YBI_GPU float YBI_MATERIAL_SELF_PDF(const KernelParams &params,
                                    const ClosureEvaluationInfo &info,
                                    const Vec3f &viewDirection,
                                    const Vec3f &lightDirection)
{
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
    return ConductorDirectionalAlbedo(bsdf, frame.WorldToLocal(viewDirection).z).MaxComponent();
}

YBI_GPU float YBI_MATERIAL_SELF_TRANSMITTANCE(const KernelParams &params,
                                              const ClosureEvaluationInfo &info,
                                              const Vec3f &viewDirection)
{
    return 1.f - YBI_MATERIAL_WEIGHT(info);
}
