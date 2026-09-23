YBI_GPU bool YBI_MATERIAL_SELF_SAMPLE(const KernelParams &params,
                                      const ShadingInfo &info,
                                      const Vec3f &viewDirection,
                                      Vec3f &lightDirection,
                                      Vec3f u)
{
    const OrthonormalBasis frame =
        MaterialXBsdfFrame(info, YBI_MATERIAL_NORMAL(info), viewDirection);
    Vec3f localLightDirection;
    if (!SampleOrenNayarBsdf(frame.WorldToLocal(viewDirection), u.Xy(), localLightDirection))
        return false;
    lightDirection = frame.LocalToWorld(localLightDirection);
    return true;
}

YBI_GPU Vec3f YBI_MATERIAL_SELF_EVALUATE(const KernelParams &params,
                                         const ShadingInfo &info,
                                         const Vec3f &viewDirection,
                                         const Vec3f &lightDirection)
{
    const OrthonormalBasis frame =
        MaterialXBsdfFrame(info, YBI_MATERIAL_NORMAL(info), viewDirection);
    return EvaluateOrenNayarBsdf(frame.WorldToLocal(viewDirection),
                                 frame.WorldToLocal(lightDirection),
                                 YBI_MATERIAL_WEIGHT(info),
                                 YBI_MATERIAL_COLOR(info),
                                 YBI_MATERIAL_ROUGHNESS(info));
}

YBI_GPU float YBI_MATERIAL_SELF_PDF(const KernelParams &params,
                                    const ShadingInfo &info,
                                    const Vec3f &viewDirection,
                                    const Vec3f &lightDirection)
{
    const OrthonormalBasis frame =
        MaterialXBsdfFrame(info, YBI_MATERIAL_NORMAL(info), viewDirection);
    return PDFOrenNayarBsdf(frame.WorldToLocal(viewDirection),
                            frame.WorldToLocal(lightDirection));
}

YBI_GPU float YBI_MATERIAL_SELF_ALBEDO(const KernelParams &params,
                                       const ShadingInfo &info,
                                       const Vec3f &viewDirection)
{
    const OrthonormalBasis frame =
        MaterialXBsdfFrame(info, YBI_MATERIAL_NORMAL(info), viewDirection);
    return YBI_MATERIAL_WEIGHT(info) *
           OrenNayarDirectionalAlbedo(frame.WorldToLocal(viewDirection).z,
                                      YBI_MATERIAL_ROUGHNESS(info));
}

YBI_GPU float YBI_MATERIAL_SELF_TRANSMITTANCE(const KernelParams &params,
                                              const ShadingInfo &info,
                                              const Vec3f &viewDirection)
{
    return 1.f - YBI_MATERIAL_WEIGHT(info);
}
