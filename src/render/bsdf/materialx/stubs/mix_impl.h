bool YBI_MATERIAL_SELF_SAMPLE(const Vec3f &viewDirection,
                              const Vec3f &lightDirection,
                              Vec3f u)
{
    float albedoBg = YBI_MATERIAL_ALBEDO_BG();
    float albedoFg = YBI_MATERIAL_ALBEDO_FG();

    float prob = albedoBg / (albedoBg + albedoFg);

    if (u.x < prob)
    {
        u.x = u.x / prob;
        YBI_MATERIAL_SAMPLE_BG(viewDirection, lightDirection, u);
    }
    else 
    {
        u.x = (u.x - prob) / (1.f - prob);
        YBI_MATERIAL_SAMPLE_FG(viewDirection, lightDirection, u);
    }
}

float YBI_MATERIAL_SELF_TRANSMITTANCE(const YBI_MATERIAL_SELF_STRUCT &bsdf)
{
    const float transmittanceBg = YBI_MATERIAL_TRANSMITTANCE_BG();
    const float transmittanceFg = YBI_MATERIAL_TRANSMITTANCE_FG();

    float transmittance = Lerp(transmittanceBg, transmittanceFg, bsdf.weight);
    return transmittance;
}

float YBI_MATERIAL_SELF_ALBEDO()
{
    const float albedoBg = YBI_MATERIAL_ALBEDO_BG();
    const float albedoFg = YBI_MATERIAL_ALBEDO_FG();

    const float albedo = Lerp(albedoBg, albedoFg, bsdf.weight);
    return albedo;
}
