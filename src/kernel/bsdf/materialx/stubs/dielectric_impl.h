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
    const float ior = YBI_MATERIAL_IOR(info);
    const float eta = Dot(viewDirection, info.geometricNormal) > 0.f ? ior : 1.f / ior;
    if (IsSmoothMicrofacet(roughness))
    {
        if (!SampleSmoothDielectricBsdf(localViewDirection,
                                       u.x,
                                       YBI_MATERIAL_WEIGHT(info),
                                       YBI_MATERIAL_TINT(info),
                                       eta,
                                       YBI_MATERIAL_SCATTER_FLAGS,
                                       sample))
            return false;
        sample.direction = frame.LocalToWorld(sample.direction);
        return true;
    }

    YBI_ASSERT(YBI_MATERIAL_SCATTER_FLAGS == BSDF_REFLECTION &&
               "Rough dielectric transmission is not supported");
    if (YBI_MATERIAL_SCATTER_FLAGS != BSDF_REFLECTION)
        return false;
    Vec3f localLightDirection;
    if (!SampleDielectricBsdf(localViewDirection,
                              u.Xy(),
                              roughness,
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
    if (IsSmoothMicrofacet(YBI_MATERIAL_ROUGHNESS(info)) ||
        YBI_MATERIAL_SCATTER_FLAGS != BSDF_REFLECTION)
        return Vec3f(0.f);
    const OrthonormalBasis frame =
        MaterialXBsdfFrame(info, YBI_MATERIAL_NORMAL(info), viewDirection);
    const float ior = YBI_MATERIAL_IOR(info);
    const float eta = Dot(viewDirection, info.geometricNormal) > 0.f ? ior : 1.f / ior;
    return EvaluateDielectricBsdf(frame.WorldToLocal(viewDirection),
                                  frame.WorldToLocal(lightDirection),
                                  YBI_MATERIAL_WEIGHT(info),
                                  YBI_MATERIAL_TINT(info),
                                  eta,
                                  YBI_MATERIAL_ROUGHNESS(info));
}

YBI_GPU float YBI_MATERIAL_SELF_PDF(const KernelParams &params,
                                    const ClosureEvaluationInfo &info,
                                    const Vec3f &viewDirection,
                                    const Vec3f &lightDirection)
{
    if (IsSmoothMicrofacet(YBI_MATERIAL_ROUGHNESS(info)) ||
        YBI_MATERIAL_SCATTER_FLAGS != BSDF_REFLECTION)
        return 0.f;
    const OrthonormalBasis frame =
        MaterialXBsdfFrame(info, YBI_MATERIAL_NORMAL(info), viewDirection);
    return PDFDielectricBsdf(frame.WorldToLocal(viewDirection),
                             frame.WorldToLocal(lightDirection),
                             YBI_MATERIAL_ROUGHNESS(info));
}

YBI_GPU float YBI_MATERIAL_SELF_ALBEDO(const KernelParams &params,
                                       const ClosureEvaluationInfo &info,
                                       const Vec3f &viewDirection)
{
    const OrthonormalBasis frame =
        MaterialXBsdfFrame(info, YBI_MATERIAL_NORMAL(info), viewDirection);
    const float ior = YBI_MATERIAL_IOR(info);
    const float eta = Dot(viewDirection, info.geometricNormal) > 0.f ? ior : 1.f / ior;
    if (IsSmoothMicrofacet(YBI_MATERIAL_ROUGHNESS(info)))
    {
        const float fresnel = FresnelDielectric(frame.WorldToLocal(viewDirection).z, eta);
        const float reflection = (YBI_MATERIAL_SCATTER_FLAGS & BSDF_REFLECTION) ? fresnel : 0.f;
        const float transmission =
            (YBI_MATERIAL_SCATTER_FLAGS & BSDF_TRANSMISSION) ? 1.f - fresnel : 0.f;
        return YBI_MATERIAL_WEIGHT(info) * (reflection + transmission) *
               Max(YBI_MATERIAL_TINT(info), 0.f).MaxComponent();
    }
    return DielectricAlbedo(frame.WorldToLocal(viewDirection).z,
                            YBI_MATERIAL_WEIGHT(info),
                            eta,
                            YBI_MATERIAL_ROUGHNESS(info));
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
