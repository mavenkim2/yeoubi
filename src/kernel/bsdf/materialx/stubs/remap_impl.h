YBI_GPU YBI_VALUE_TYPE YBI_SELF(const KernelParams &params,
                                const ClosureInitializationInfo &info)
{
    const YBI_VALUE_TYPE value = YBI_VALUE;
    const YBI_VALUE_TYPE inLow = YBI_IN_LOW;
    const YBI_VALUE_TYPE inHigh = YBI_IN_HIGH;
    const YBI_VALUE_TYPE outLow = YBI_OUT_LOW;
    const YBI_VALUE_TYPE outHigh = YBI_OUT_HIGH;

    const YBI_VALUE_TYPE result = outLow + (value - inLow) * (outHigh - outLow) / (inHigh - inLow);
    return result;
}
