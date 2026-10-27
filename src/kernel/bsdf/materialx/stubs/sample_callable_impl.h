YBI_MATERIAL_CALLABLE_SIGNATURE bool
YBI_MATERIAL_CALLABLE_NAME(const KernelParams *params,
                           const ClosureEvaluationInfo *info,
                           Vec3f viewDirection,
                           BsdfSample *sample,
                           Vec3f u)
{
    if (!YBI_MATERIAL_CALLABLE_IMPL(*params, *info, viewDirection, *sample, u))
        return false;

    if (!(sample->flags & BSDF_DELTA))
    {
        // Evaluate the complete graph once for continuous samples.
        sample->pdf = YBI_MATERIAL_CALLABLE_PDF(*params, *info, viewDirection, sample->direction);
        if (sample->pdf <= 0.f)
            return false;
        sample->weight =
            YBI_MATERIAL_CALLABLE_EVALUATE(*params, *info, viewDirection, sample->direction) /
            sample->pdf;
    }
    sample->stochasticPdf = sample->pdf;
    sample->bsdfPdf = sample->pdf;
    return true;
}
