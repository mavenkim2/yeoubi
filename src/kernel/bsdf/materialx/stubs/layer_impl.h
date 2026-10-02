YBI_GPU bool YBI_MATERIAL_SELF_SAMPLE(const KernelParams &params,
                                      const ClosureEvaluationInfo &info,
                                      const Vec3f &viewDirection,
                                      Vec3f &lightDirection,
                                      Vec3f u)
{
    float albedoTop = YBI_MATERIAL_ALBEDO_TOP(params, info, viewDirection);
    float transmittanceTop = YBI_MATERIAL_TRANSMITTANCE_TOP(params, info, viewDirection);
    float albedoBottom = YBI_MATERIAL_ALBEDO_BOTTOM(params, info, viewDirection);
    float weightedBottomAlbedo = transmittanceTop * albedoBottom;

    float totalAlbedo = albedoTop + weightedBottomAlbedo;
    float prob = totalAlbedo > 0.f ? albedoTop / totalAlbedo : 0.5f;

    if (u.x < prob)
    {
        u.x = u.x / prob;
        YBI_MATERIAL_DEBUG_PRINT_0(
            prob, 1.f - prob, albedoTop, transmittanceTop, albedoBottom, weightedBottomAlbedo);
        return YBI_MATERIAL_SAMPLE_TOP(params, info, viewDirection, lightDirection, u);
    }
    else
    {
        u.x = (u.x - prob) / (1.f - prob);
        YBI_MATERIAL_DEBUG_PRINT_1(
            1.f - prob, prob, albedoTop, transmittanceTop, albedoBottom, weightedBottomAlbedo);
        return YBI_MATERIAL_SAMPLE_BOTTOM(params, info, viewDirection, lightDirection, u);
    }
}

YBI_GPU Vec3f YBI_MATERIAL_SELF_EVALUATE(const KernelParams &params,
                                         const ClosureEvaluationInfo &info,
                                         const Vec3f &viewDirection,
                                         const Vec3f &lightDirection)
{
    float transmittanceTop = YBI_MATERIAL_TRANSMITTANCE_TOP(params, info, viewDirection);
    return YBI_MATERIAL_EVALUATE_TOP(params, info, viewDirection, lightDirection) +
           transmittanceTop *
               YBI_MATERIAL_EVALUATE_BOTTOM(params, info, viewDirection, lightDirection);
}

YBI_GPU float YBI_MATERIAL_SELF_PDF(const KernelParams &params,
                                    const ClosureEvaluationInfo &info,
                                    const Vec3f &viewDirection,
                                    const Vec3f &lightDirection)
{
    float albedoTop = YBI_MATERIAL_ALBEDO_TOP(params, info, viewDirection);
    float transmittanceTop = YBI_MATERIAL_TRANSMITTANCE_TOP(params, info, viewDirection);
    float albedoBottom =
        transmittanceTop * YBI_MATERIAL_ALBEDO_BOTTOM(params, info, viewDirection);
    float totalAlbedo = albedoTop + albedoBottom;
    float prob = totalAlbedo > 0.f ? albedoTop / totalAlbedo : 0.5f;
    return prob * YBI_MATERIAL_PDF_TOP(params, info, viewDirection, lightDirection) +
           (1.f - prob) *
               YBI_MATERIAL_PDF_BOTTOM(params, info, viewDirection, lightDirection);
}

YBI_GPU float YBI_MATERIAL_SELF_TRANSMITTANCE(const KernelParams &params,
                                              const ClosureEvaluationInfo &info,
                                              const Vec3f &viewDirection)
{
    const float topTransmittance =
        YBI_MATERIAL_TRANSMITTANCE_TOP(params, info, viewDirection);
    const float bottomTransmittance =
        YBI_MATERIAL_TRANSMITTANCE_BOTTOM(params, info, viewDirection);
    const float transmittance = topTransmittance * bottomTransmittance;
    YBI_MATERIAL_DEBUG_PRINT_TRANSMITTANCE(
        topTransmittance, bottomTransmittance, transmittance);
    return transmittance;
}

YBI_GPU float YBI_MATERIAL_SELF_ALBEDO(const KernelParams &params,
                                       const ClosureEvaluationInfo &info,
                                       const Vec3f &viewDirection)
{
    float transmittanceTop = YBI_MATERIAL_TRANSMITTANCE_TOP(params, info, viewDirection);
    return YBI_MATERIAL_ALBEDO_TOP(params, info, viewDirection) +
           transmittanceTop * YBI_MATERIAL_ALBEDO_BOTTOM(params, info, viewDirection);
}
