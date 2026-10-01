YBI_GPU bool YBI_MATERIAL_SELF_SAMPLE(const KernelParams &params,
                                      const ShadingInfo &info,
                                      const Vec3f &viewDirection,
                                      Vec3f &lightDirection,
                                      Vec3f u)
{
    const float mix = Clamp(YBI_MATERIAL_MIX(info), 0.f, 1.f);
    const float bgProbability = 1.f - mix;
    if (u.x < bgProbability)
    {
        u.x /= bgProbability;
        YBI_MATERIAL_DEBUG_PRINT_0(bgProbability, mix, bgProbability, mix);
        return YBI_MATERIAL_SAMPLE_BG(params, info, viewDirection, lightDirection, u);
    }

    u.x = (u.x - bgProbability) / mix;
    YBI_MATERIAL_DEBUG_PRINT_1(mix, bgProbability, bgProbability, mix);
    return YBI_MATERIAL_SAMPLE_FG(params, info, viewDirection, lightDirection, u);
}

YBI_GPU Vec3f YBI_MATERIAL_SELF_EVALUATE(const KernelParams &params,
                                         const ShadingInfo &info,
                                         const Vec3f &viewDirection,
                                         const Vec3f &lightDirection)
{
    const float mix = Clamp(YBI_MATERIAL_MIX(info), 0.f, 1.f);
    return (1.f - mix) *
               YBI_MATERIAL_EVALUATE_BG(params, info, viewDirection, lightDirection) +
           mix * YBI_MATERIAL_EVALUATE_FG(params, info, viewDirection, lightDirection);
}

YBI_GPU float YBI_MATERIAL_SELF_PDF(const KernelParams &params,
                                    const ShadingInfo &info,
                                    const Vec3f &viewDirection,
                                    const Vec3f &lightDirection)
{
    const float mix = Clamp(YBI_MATERIAL_MIX(info), 0.f, 1.f);
    return (1.f - mix) * YBI_MATERIAL_PDF_BG(params, info, viewDirection, lightDirection) +
           mix * YBI_MATERIAL_PDF_FG(params, info, viewDirection, lightDirection);
}

YBI_GPU float YBI_MATERIAL_SELF_ALBEDO(const KernelParams &params,
                                       const ShadingInfo &info,
                                       const Vec3f &viewDirection)
{
    const float mix = Clamp(YBI_MATERIAL_MIX(info), 0.f, 1.f);
    return Lerp(YBI_MATERIAL_ALBEDO_BG(params, info, viewDirection),
                YBI_MATERIAL_ALBEDO_FG(params, info, viewDirection),
                mix);
}

YBI_GPU float YBI_MATERIAL_SELF_TRANSMITTANCE(const KernelParams &params,
                                              const ShadingInfo &info,
                                              const Vec3f &viewDirection)
{
    const float mix = Clamp(YBI_MATERIAL_MIX(info), 0.f, 1.f);
    const float bgTransmittance =
        YBI_MATERIAL_TRANSMITTANCE_BG(params, info, viewDirection);
    const float fgTransmittance =
        YBI_MATERIAL_TRANSMITTANCE_FG(params, info, viewDirection);
    const float transmittance = Lerp(bgTransmittance, fgTransmittance, mix);
    YBI_MATERIAL_DEBUG_PRINT_TRANSMITTANCE(
        mix, bgTransmittance, fgTransmittance, transmittance);
    return transmittance;
}
