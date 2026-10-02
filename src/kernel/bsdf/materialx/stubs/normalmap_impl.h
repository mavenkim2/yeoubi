YBI_GPU Vec3f YBI_MATERIAL_SELF(const KernelParams &params, const ShadingInfo &info)
{
    Vec3f value = YBI_MATERIAL_INPUT;
    value = value == Vec3f(0.f) ? Vec3f(0.f, 0.f, 1.f) : value * 2.f - Vec3f(1.f);

    const Vec3f normal = Normalize(YBI_MATERIAL_NORMAL);
    const Vec3f tangent = Normalize(YBI_MATERIAL_TANGENT);
    const Vec3f bitangent = Normalize(YBI_MATERIAL_BITANGENT);
    const float scale = YBI_MATERIAL_SCALE;
    return Normalize(tangent * (value.x * scale) + bitangent * (value.y * scale) +
                     normal * value.z);
}
