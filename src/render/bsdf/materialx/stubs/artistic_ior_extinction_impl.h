YBI_GPU Vec3f YBI_MATERIAL_SELF(const KernelParams &params, const ShadingInfo &info)
{
    const Vec3f inputReflectivity = YBI_MATERIAL_REFLECTIVITY;
    const Vec3f reflectivity = Clamp(inputReflectivity, 0.f, 0.99f);
    const Vec3f edgeColor = YBI_MATERIAL_EDGE_COLOR;
    const Vec3f reflectivitySqrt = Sqrt(reflectivity);
    const Vec3f nMin = (1.f - reflectivity) / (1.f + reflectivity);
    const Vec3f nMax = (1.f + reflectivitySqrt) / (1.f - reflectivitySqrt);
    const Vec3f ior = nMax + (nMin - nMax) * edgeColor;

    const Vec3f nPlusOne = ior + 1.f;
    const Vec3f nMinusOne = ior - 1.f;
    const Vec3f extinctionSquared =
        Max((nPlusOne * nPlusOne * reflectivity - nMinusOne * nMinusOne) /
                (1.f - reflectivity),
            0.f);
    const Vec3f output = Sqrt(extinctionSquared);
    YBI_MATERIAL_DEBUG_PRINT(inputReflectivity, edgeColor, output);
    return output;
}
