YBI_GPU bool YBI_MATERIAL_SELF_SAMPLE(const KernelParams &params,
                                      const ShadingInfo &info,
                                      const Vec3f &viewDirection,
                                      Vec3f &lightDirection,
                                      Vec3f u)
{
    return SampleDielectricBsdf(
        viewDirection, u.Xy(), YBI_MATERIAL_ROUGHNESS(params, info), lightDirection);
}

YBI_GPU Vec3f YBI_MATERIAL_SELF_EVALUATE(const KernelParams &params,
                                         const ShadingInfo &info,
                                         const Vec3f &viewDirection,
                                         const Vec3f &lightDirection)
{
    return EvaluateDielectricBsdf(viewDirection,
                                  lightDirection,
                                  YBI_MATERIAL_WEIGHT(params, info),
                                  YBI_MATERIAL_TINT(params, info),
                                  YBI_MATERIAL_IOR(params, info),
                                  YBI_MATERIAL_ROUGHNESS(params, info));
}

YBI_GPU float YBI_MATERIAL_SELF_PDF(const KernelParams &params,
                                    const ShadingInfo &info,
                                    const Vec3f &viewDirection,
                                    const Vec3f &lightDirection)
{
    return PDFDielectricBsdf(
        viewDirection, lightDirection, YBI_MATERIAL_ROUGHNESS(params, info));
}

YBI_GPU float YBI_MATERIAL_SELF_ALBEDO(const KernelParams &params,
                                       const ShadingInfo &info,
                                       const Vec3f &viewDirection)
{
    return DielectricAlbedo(viewDirection.z,
                            YBI_MATERIAL_WEIGHT(params, info),
                            YBI_MATERIAL_IOR(params, info),
                            YBI_MATERIAL_ROUGHNESS(params, info));
}

YBI_GPU float YBI_MATERIAL_SELF_TRANSMITTANCE(const KernelParams &params,
                                              const ShadingInfo &info,
                                              const Vec3f &viewDirection)
{
    return 1.f - YBI_MATERIAL_SELF_ALBEDO(params, info, viewDirection);
}
