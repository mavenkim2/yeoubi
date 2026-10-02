YBI_GPU Vec3f YBI_MATERIAL_SELF(const KernelParams &params,
                                const ClosureInitializationInfo &info)
{
    const Vec3f inputReflectivity = YBI_MATERIAL_REFLECTIVITY;
    const Vec3f reflectivity = Clamp(inputReflectivity, 0.f, 0.99f);
    const Vec3f edgeColor = YBI_MATERIAL_EDGE_COLOR;
    const Vec3f reflectivitySqrt = Sqrt(reflectivity);
    const Vec3f nMin = (1.f - reflectivity) / (1.f + reflectivity);
    const Vec3f nMax = (1.f + reflectivitySqrt) / (1.f - reflectivitySqrt);
    const Vec3f output = nMax + (nMin - nMax) * edgeColor;
    YBI_MATERIAL_DEBUG_PRINT(inputReflectivity, edgeColor, output);
    return output;
}
