YBI_GPU YBI_MATERIAL_RETURN
YBI_MATERIAL_SELF(const KernelParams &params, const ShadingInfo &info)
{
    constexpr MaterialTextureHandle handle = YBI_MATERIAL_HANDLE;
    const Vec4f value = VirtualTextureSampleImage(params, handle, info);
    return YBI_MATERIAL_RESULT;
}
