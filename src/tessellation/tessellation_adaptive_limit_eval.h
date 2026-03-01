static bool IsSupportedLimitEvalInterpolation(PrimvarInterpolation interpolation)
{
    return interpolation == PrimvarInterpolation::Vertex ||
           interpolation == PrimvarInterpolation::Varying ||
           interpolation == PrimvarInterpolation::FaceVarying;
}

static int GetTotalRefinerValueCount(const Far::TopologyRefiner &refiner,
                                     PrimvarInterpolation interpolation,
                                     int channel)
{
    switch (interpolation)
    {
        case PrimvarInterpolation::Vertex:
        case PrimvarInterpolation::Varying:
            return refiner.GetNumVerticesTotal();
        case PrimvarInterpolation::FaceVarying:
            return refiner.GetNumFVarValuesTotal(channel);
        default:
            return 0;
    }
}

static int GetLevelValueCount(const Far::TopologyRefiner &refiner,
                              int level,
                              PrimvarInterpolation interpolation,
                              int channel)
{
    const Far::TopologyLevel &topologyLevel = refiner.GetLevel(level);
    switch (interpolation)
    {
        case PrimvarInterpolation::Vertex:
        case PrimvarInterpolation::Varying:
            return topologyLevel.GetNumVertices();
        case PrimvarInterpolation::FaceVarying:
            return topologyLevel.GetNumFVarValues(channel);
        default:
            return 0;
    }
}

static int GetLocalPointValueCount(const Far::PatchTable &patchTable,
                                   PrimvarInterpolation interpolation,
                                   int channel)
{
    switch (interpolation)
    {
        case PrimvarInterpolation::Vertex:
            return patchTable.GetNumLocalPoints();
        case PrimvarInterpolation::Varying:
            return patchTable.GetNumLocalPointsVarying();
        case PrimvarInterpolation::FaceVarying:
            return patchTable.GetNumLocalPointsFaceVarying(channel);
        default:
            return 0;
    }
}

template <typename LimitEvalT, typename ValueT>
static std::vector<LimitEvalT> BuildLimitEvalValues(const Far::TopologyRefiner &refiner,
                                                    const Far::PatchTable &patchTable,
                                                    const Array<ValueT> &coarseValues,
                                                    PrimvarInterpolation interpolation,
                                                    int channel = 0)
{
    if (!IsSupportedLimitEvalInterpolation(interpolation))
    {
        return {};
    }
    if (interpolation == PrimvarInterpolation::FaceVarying &&
        patchTable.GetNumFVarChannels() <= channel)
    {
        return {};
    }

    const int numRefinerValues = GetTotalRefinerValueCount(refiner, interpolation, channel);
    const int numLocalPoints = GetLocalPointValueCount(patchTable, interpolation, channel);
    if (numRefinerValues <= 0)
    {
        return {};
    }
    std::vector<LimitEvalT> values(numRefinerValues + numLocalPoints);

    const int numBaseValues = GetLevelValueCount(refiner, 0, interpolation, channel);
    const int copyCount = std::min(numBaseValues, int(coarseValues.size()));
    for (int i = 0; i < copyCount; ++i)
    {
        values[i].value = coarseValues[i];
    }

    Far::PrimvarRefiner primvarRefiner(refiner);
    LimitEvalT *src = values.data();
    for (int level = 1; level < refiner.GetNumLevels(); ++level)
    {
        LimitEvalT *dst = src + GetLevelValueCount(refiner, level - 1, interpolation, channel);
        if (interpolation == PrimvarInterpolation::Vertex)
        {
            primvarRefiner.Interpolate(level, src, dst);
        }
        else if (interpolation == PrimvarInterpolation::Varying)
        {
            primvarRefiner.InterpolateVarying(level, src, dst);
        }
        else
        {
            primvarRefiner.InterpolateFaceVarying(level, src, dst, channel);
        }
        src = dst;
    }

    if (numLocalPoints > 0)
    {
        if (interpolation == PrimvarInterpolation::Vertex)
        {
            patchTable.ComputeLocalPointValues(values.data(), values.data() + numRefinerValues);
        }
        else if (interpolation == PrimvarInterpolation::Varying)
        {
            patchTable.ComputeLocalPointValuesVarying(values.data(), values.data() + numRefinerValues);
        }
        else
        {
            patchTable.ComputeLocalPointValuesFaceVarying(values.data(),
                                                          values.data() + numRefinerValues,
                                                          channel);
        }
    }
    return values;
}

template <typename LimitEvalT, typename ValueT>
static bool EvaluateLimitValue(const Far::PatchMap &patchMap,
                               const Far::PatchTable &patchTable,
                               const std::vector<LimitEvalT> &limitValues,
                               int ptexFaceId,
                               const pxr::GfVec2f &uv,
                               PrimvarInterpolation interpolation,
                               ValueT *outValue,
                               int channel = 0)
{
    if (!outValue || !IsSupportedLimitEvalInterpolation(interpolation))
    {
        return false;
    }
    if (interpolation == PrimvarInterpolation::FaceVarying &&
        patchTable.GetNumFVarChannels() <= channel)
    {
        return false;
    }

    const Far::PatchTable::PatchHandle *handle = patchMap.FindPatch(ptexFaceId, uv[0], uv[1]);
    if (!handle)
    {
        return false;
    }

    float pWeights[20] = {0.0f};
    Far::ConstIndexArray cvs;
    if (interpolation == PrimvarInterpolation::Vertex)
    {
        patchTable.EvaluateBasis(*handle, uv[0], uv[1], pWeights);
        cvs = patchTable.GetPatchVertices(*handle);
    }
    else if (interpolation == PrimvarInterpolation::Varying)
    {
        patchTable.EvaluateBasisVarying(*handle, uv[0], uv[1], pWeights);
        cvs = patchTable.GetPatchVaryingVertices(*handle);
    }
    else
    {
        patchTable.EvaluateBasisFaceVarying(*handle, uv[0], uv[1], pWeights, 0, 0, 0, 0, 0, channel);
        cvs = patchTable.GetPatchFVarValues(*handle, channel);
    }
    if (cvs.size() == 0)
    {
        return false;
    }

    ValueT value = LimitEvalValueTraits<ValueT>::Zero();
    for (int i = 0; i < cvs.size(); ++i)
    {
        if (cvs[i] < 0 || cvs[i] >= int(limitValues.size()))
        {
            return false;
        }
        value += limitValues[cvs[i]].value * pWeights[i];
    }
    *outValue = value;
    return true;
}

static std::vector<LimitEvalVertex> BuildLimitEvalVertices(const Far::TopologyRefiner &refiner,
                                                           const Far::PatchTable &patchTable,
                                                           const Array<float3> &coarsePoints)
{
    return BuildLimitEvalValues<LimitEvalVertex, float3>(
        refiner, patchTable, coarsePoints, PrimvarInterpolation::Vertex);
}

static std::vector<LimitEvalFVar2> BuildLimitEvalFVarValues(
    const Far::TopologyRefiner &refiner,
    const Far::PatchTable &patchTable,
    const Array<float2> &coarseValues,
    PrimvarInterpolation interpolation = PrimvarInterpolation::FaceVarying,
    int channel = 0)
{
    return BuildLimitEvalValues<LimitEvalFVar2, float2>(
        refiner, patchTable, coarseValues, interpolation, channel);
}

static bool EvaluateLimitPosition(const Far::PatchMap &patchMap,
                                  const Far::PatchTable &patchTable,
                                  const std::vector<LimitEvalVertex> &limitValues,
                                  int ptexFaceId,
                                  const pxr::GfVec2f &uv,
                                  float3 *outP)
{
    return EvaluateLimitValue<LimitEvalVertex, float3>(patchMap,
                                                       patchTable,
                                                       limitValues,
                                                       ptexFaceId,
                                                       uv,
                                                       PrimvarInterpolation::Vertex,
                                                       outP);
}

static bool EvaluateLimitFVar2(const Far::PatchMap &patchMap,
                               const Far::PatchTable &patchTable,
                               const std::vector<LimitEvalFVar2> &limitValues,
                               int ptexFaceId,
                               const pxr::GfVec2f &uv,
                               float2 *outUV,
                               PrimvarInterpolation interpolation = PrimvarInterpolation::FaceVarying,
                               int channel = 0)
{
    return EvaluateLimitValue<LimitEvalFVar2, float2>(
        patchMap, patchTable, limitValues, ptexFaceId, uv, interpolation, outUV, channel);
}

static float4x4 BuildFallbackCameraFromWorld(const float3 &eye, const float3 &lookAt)
{
    float3 forward = normalize(lookAt - eye);
    if (length(forward) <= 1e-8f)
    {
        forward = make_float3(0.0f, 0.0f, 1.0f);
    }
    float3 worldUp = make_float3(0.0f, 0.0f, 1.0f);
    if (std::abs(dot(forward, worldUp)) > 0.999f)
    {
        worldUp = make_float3(0.0f, 1.0f, 0.0f);
    }
    const float3 right = normalize(cross(forward, worldUp));
    const float3 up = normalize(cross(right, forward));
    return float4x4(right.x,
                    right.y,
                    right.z,
                    -dot(right, eye),
                    up.x,
                    up.y,
                    up.z,
                    -dot(up, eye),
                    forward.x,
                    forward.y,
                    forward.z,
                    -dot(forward, eye),
                    0.0f,
                    0.0f,
                    0.0f,
                    1.0f);
}

static float4x4 BuildFallbackClipFromCamera(float verticalFovDegrees, int viewportWidth, int viewportHeight)
{
    const float fovY = verticalFovDegrees * 3.14159265358979323846f / 180.0f;
    const float tanHalfFovY = std::max(1e-8f, std::tan(0.5f * fovY));
    const float aspect = std::max(1e-8f, float(viewportWidth) / float(viewportHeight));
    const float nearPlane = 1.0f;
    const float farPlane = 1.0e6f;
    const float m00 = 1.0f / (tanHalfFovY * aspect);
    const float m11 = 1.0f / tanHalfFovY;
    const float m22 = (farPlane + nearPlane) / (farPlane - nearPlane);
    const float m23 = (-2.0f * farPlane * nearPlane) / (farPlane - nearPlane);
    return float4x4(m00,
                    0.0f,
                    0.0f,
                    0.0f,
                    0.0f,
                    m11,
                    0.0f,
                    0.0f,
                    0.0f,
                    0.0f,
                    m22,
                    m23,
                    0.0f,
                    0.0f,
                    1.0f,
                    0.0f);
}

static float EvaluateFrustumPlane(const float4 &p, int planeIndex)
{
    switch (planeIndex)
    {
        case 0:
            return p.x + p.w;
        case 1:
            return -p.x + p.w;
        case 2:
            return p.y + p.w;
        case 3:
            return -p.y + p.w;
        case 4:
            return p.z + p.w;
        default:
            return -p.z + p.w;
    }
}

static bool ClipSegmentToFrustum(float4 *a, float4 *b)
{
    YBI_ASSERT(a);
    YBI_ASSERT(b);
    float t0 = 0.0f;
    float t1 = 1.0f;
    const float4 ab = *b - *a;
    for (int plane = 0; plane < 6; ++plane)
    {
        const float fa = EvaluateFrustumPlane(*a, plane);
        const float fb = EvaluateFrustumPlane(*b, plane);
        if (fa < 0.0f && fb < 0.0f)
        {
            return false;
        }
        if ((fa < 0.0f) != (fb < 0.0f))
        {
            const float denom = fa - fb;
            if (std::abs(denom) <= 1e-20f)
            {
                return false;
            }
            const float t = fa / denom;
            if (fa < 0.0f)
            {
                t0 = std::max(t0, t);
            }
            else
            {
                t1 = std::min(t1, t);
            }
            if (t0 > t1)
            {
                return false;
            }
        }
    }
    *a = *a + ab * t0;
    *b = *a + ab * (t1 - t0);
    return true;
}

static bool ProjectClipPointToScreen(const float4 &clipP,
                                     int viewportWidth,
                                     int viewportHeight,
                                     pxr::GfVec2f *outScreen)
{
    YBI_ASSERT(outScreen);
    if (!std::isfinite(clipP.w) || std::abs(clipP.w) <= 1e-20f)
    {
        return false;
    }
    const float invW = 1.0f / clipP.w;
    const float ndcX = clipP.x * invW;
    const float ndcY = clipP.y * invW;
    if (!std::isfinite(ndcX) || !std::isfinite(ndcY))
    {
        return false;
    }
    const float sx = (ndcX * 0.5f + 0.5f) * float(viewportWidth);
    const float sy = (1.0f - (ndcY * 0.5f + 0.5f)) * float(viewportHeight);
    *outScreen = pxr::GfVec2f(sx, sy);
    return std::isfinite(sx) && std::isfinite(sy);
}

enum DiagSplitNonUniformReason
{
    DIAGSPLIT_NON_UNIFORM_NONE = 0,
    DIAGSPLIT_NON_UNIFORM_EVAL_FAIL = 1,
    DIAGSPLIT_NON_UNIFORM_VARIANCE_THRESHOLD = 2,
};

static int ComputeDiagSplitPatchEdgeFactor(const Far::PatchMap &patchMap,
                                           const Far::PatchTable &patchTable,
                                           const std::vector<LimitEvalVertex> &limitValues,
                                           int ptexFaceId,
                                           const pxr::GfVec2f &uvStart,
                                           const pxr::GfVec2f &uvEnd,
                                           int sampleStepsN,
                                           float targetPixelSpacing,
                                           int splitThreshold,
                                           const float3 &eye,
                                           const float3 &lookAt,
                                           int viewportWidth,
                                           int viewportHeight,
                                           float verticalFovDegrees,
                                           bool useCameraMatrices,
                                           const float4x4 &cameraFromWorld,
                                           const float4x4 &clipFromCamera,
                                           DiagSplitNonUniformReason *nonUniformReasonOut)
{
    if (nonUniformReasonOut)
    {
        *nonUniformReasonOut = DIAGSPLIT_NON_UNIFORM_NONE;
    }
    if (sampleStepsN < 2 || targetPixelSpacing <= 0.0f)
    {
        return 1;
    }

    if (viewportWidth <= 0 || viewportHeight <= 0)
    {
        return 1;
    }
    const float4x4 usedCameraFromWorld =
        useCameraMatrices ? cameraFromWorld : BuildFallbackCameraFromWorld(eye, lookAt);
    const float4x4 usedClipFromCamera =
        useCameraMatrices
            ? clipFromCamera
            : BuildFallbackClipFromCamera(verticalFovDegrees, viewportWidth, viewportHeight);
    const float4x4 clipFromWorld = mul(usedClipFromCamera, usedCameraFromWorld);

    float maxLi = 0.0f;
    float sumLi = 0.0f;
    bool hadVisibleSegment = false;
    float3 p0 = make_float3(0.0f);
    if (!EvaluateLimitPosition(patchMap, patchTable, limitValues, ptexFaceId, uvStart, &p0))
    {
        if (nonUniformReasonOut)
        {
            *nonUniformReasonOut = DIAGSPLIT_NON_UNIFORM_EVAL_FAIL;
        }
        return SUBDIV_EDGE_FACTOR_NON_UNIFORM;
    }
    float4 prevClip = mul(clipFromWorld, make_float4(p0.x, p0.y, p0.z, 1.0f));
    for (int i = 1; i < sampleStepsN; ++i)
    {
        const float t = float(i) / float(sampleStepsN - 1);
        const pxr::GfVec2f uv = uvStart * (1.0f - t) + uvEnd * t;
        float3 p = make_float3(0.0f);
        if (!EvaluateLimitPosition(patchMap, patchTable, limitValues, ptexFaceId, uv, &p))
        {
            if (nonUniformReasonOut)
            {
                *nonUniformReasonOut = DIAGSPLIT_NON_UNIFORM_EVAL_FAIL;
            }
            return SUBDIV_EDGE_FACTOR_NON_UNIFORM;
        }
        const float4 currClip = mul(clipFromWorld, make_float4(p.x, p.y, p.z, 1.0f));
        float4 segA = prevClip;
        float4 segB = currClip;
        if (ClipSegmentToFrustum(&segA, &segB))
        {
            pxr::GfVec2f sa(0.0f), sb(0.0f);
            if (ProjectClipPointToScreen(segA, viewportWidth, viewportHeight, &sa) &&
                ProjectClipPointToScreen(segB, viewportWidth, viewportHeight, &sb))
            {
                hadVisibleSegment = true;
                const pxr::GfVec2f d = sb - sa;
                const float li = std::sqrt(d[0] * d[0] + d[1] * d[1]);
                sumLi += li;
                maxLi = std::max(maxLi, li);
            }
        }
        prevClip = currClip;
    }
    if (!hadVisibleSegment)
    {
        return 1;
    }

    const int tMin = std::max(1, int(std::ceil(sumLi / targetPixelSpacing)));
    const int tMax =
        std::max(1, int(std::ceil(float(sampleStepsN - 1) * maxLi / targetPixelSpacing)));
    if ((tMax - tMin) > splitThreshold)
    {
        if (nonUniformReasonOut)
        {
            *nonUniformReasonOut = DIAGSPLIT_NON_UNIFORM_VARIANCE_THRESHOLD;
        }
        return SUBDIV_EDGE_FACTOR_NON_UNIFORM;
    }
    return tMax;
}

static float ComputeLength(const Far::PatchMap &patchMap,
                           const Far::PatchTable &patchTable,
                           const std::vector<LimitEvalVertex> &limitValues,
                           int ptexFaceId,
                           const pxr::GfVec2f &uvStart,
                           const pxr::GfVec2f &uvEnd,
                           int sampleStepsN)
{
    if (sampleStepsN < 2)
    {
        return 1;
    }

    float sumLi = 0.0f;
    float3 p0 = make_float3(0.0f);
    if (!EvaluateLimitPosition(patchMap, patchTable, limitValues, ptexFaceId, uvStart, &p0))
    {
        return 0.f;
    }
    for (int i = 1; i < sampleStepsN; ++i)
    {
        const float t = float(i) / float(sampleStepsN - 1);
        const pxr::GfVec2f uv = uvStart * (1.0f - t) + uvEnd * t;
        float3 p = make_float3(0.0f);
        if (!EvaluateLimitPosition(patchMap, patchTable, limitValues, ptexFaceId, uv, &p))
        {
            return 0.f;
        }
        sumLi += length(p - p0);
        p0 = p;
    }

    return sumLi;
}
