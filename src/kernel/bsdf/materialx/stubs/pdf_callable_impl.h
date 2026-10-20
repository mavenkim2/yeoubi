YBI_MATERIAL_CALLABLE_SIGNATURE float
YBI_MATERIAL_CALLABLE_NAME(const KernelParams *params,
                           const ClosureEvaluationInfo *info,
                           Vec3f viewDirection,
                           Vec3f lightDirection)
{
    if (Dot(viewDirection, info->geometricNormal) * Dot(lightDirection, info->geometricNormal) <=
        0.f)
        return 0.f;
    return YBI_MATERIAL_CALLABLE_IMPL(*params, *info, viewDirection, lightDirection);
}
