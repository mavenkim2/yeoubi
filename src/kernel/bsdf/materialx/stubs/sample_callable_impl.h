YBI_MATERIAL_CALLABLE_SIGNATURE bool
YBI_MATERIAL_CALLABLE_NAME(const KernelParams *params,
                           const ClosureEvaluationInfo *info,
                           Vec3f viewDirection,
                           Vec3f *lightDirection,
                           Vec3f u)
{
    return YBI_MATERIAL_CALLABLE_IMPL(
        *params, *info, viewDirection, *lightDirection, u);
}
