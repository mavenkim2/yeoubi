YBI_GPU YBI_MATERIAL_RETURN
YBI_MATERIAL_SELF(const KernelParams &params, const ShadingInfo &info)
{
    constexpr MaterialTextureHandle handle = YBI_MATERIAL_HANDLE;
    VirtualTextureMipSelection mipSelection;
    const Vec4f value = VirtualTextureSampleTiledImage(params,
                                                       handle,
                                                       info,
                                                       YBI_MATERIAL_UV_TILING,
                                                       YBI_MATERIAL_UV_OFFSET,
                                                       YBI_MATERIAL_REALWORLD_IMAGE_SIZE,
                                                       YBI_MATERIAL_REALWORLD_TILE_SIZE,
                                                       YBI_MATERIAL_MIP_SELECTION);
    YBI_MATERIAL_DEBUG_PRINT_MIP(mipSelection);
    return YBI_MATERIAL_RESULT;
}
