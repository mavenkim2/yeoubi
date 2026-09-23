YBI_GPU Vec3f YBI_MATERIAL_SELF(const KernelParams &params, const ShadingInfo &info)
{
    const Vec3f reflectivity = Clamp(YBI_MATERIAL_REFLECTIVITY, 0.f, 0.99f);
    const Vec3f edgeColor = YBI_MATERIAL_EDGE_COLOR;
    const Vec3f reflectivitySqrt = Sqrt(reflectivity);
    const Vec3f nMin = (1.f - reflectivity) / (1.f + reflectivity);
    const Vec3f nMax = (1.f + reflectivitySqrt) / (1.f - reflectivitySqrt);
    return nMax + (nMin - nMax) * edgeColor;
}
