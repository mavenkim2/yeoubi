YBI_GPU YBI_MATERIAL_VALUE_TYPE
YBI_MATERIAL_SELF_GET_CLOSURE(const ClosureEvaluationInfo &info)
{
    return *(YBI_MATERIAL_VALUE_TYPE *)(info.closureData + YBI_MATERIAL_CLOSURE_OFFSET);
}

YBI_GPU void YBI_MATERIAL_SELF_SET_CLOSURE(const KernelParams &params,
                                           ClosureInitializationInfo &info)
{
    *(YBI_MATERIAL_VALUE_TYPE *)(info.closureData + YBI_MATERIAL_CLOSURE_OFFSET) = YBI_MATERIAL_VALUE(params, info);
}
