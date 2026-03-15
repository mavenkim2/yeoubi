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
                                                           const Array<Vec3> &coarsePoints)
{
    return BuildLimitEvalValues<LimitEvalVertex, Vec3>(
        refiner, patchTable, coarsePoints, PrimvarInterpolation::Vertex);
}

static std::vector<LimitEvalFVar2> BuildLimitEvalFVarValues(
    const Far::TopologyRefiner &refiner,
    const Far::PatchTable &patchTable,
    const Array<Vec2> &coarseValues,
    PrimvarInterpolation interpolation = PrimvarInterpolation::FaceVarying,
    int channel = 0)
{
    return BuildLimitEvalValues<LimitEvalFVar2, Vec2>(
        refiner, patchTable, coarseValues, interpolation, channel);
}

static bool EvaluateLimitPosition(const Far::PatchMap &patchMap,
                                  const Far::PatchTable &patchTable,
                                  const std::vector<LimitEvalVertex> &limitValues,
                                  int ptexFaceId,
                                  const pxr::GfVec2f &uv,
                                  Vec3 *outP)
{
    return EvaluateLimitValue<LimitEvalVertex, Vec3>(patchMap,
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
                               Vec2 *outUV,
                               PrimvarInterpolation interpolation = PrimvarInterpolation::FaceVarying,
                               int channel = 0)
{
    return EvaluateLimitValue<LimitEvalFVar2, Vec2>(
        patchMap, patchTable, limitValues, ptexFaceId, uv, interpolation, outUV, channel);
}

static float EvaluateFrustumPlane(const Vec4 &p, int planeIndex)
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

static bool ClipSegmentToFrustum(Vec4 *a, Vec4 *b)
{
    YBI_ASSERT(a);
    YBI_ASSERT(b);
    float t0 = 0.0f;
    float t1 = 1.0f;
    const Vec4 ab = *b - *a;
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

static bool ProjectClipPointToScreen(const Vec4 &clipP,
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
                                           int viewportWidth,
                                           int viewportHeight,
                                           const Float4x4 &cameraFromWorld,
                                           const Float4x4 &clipFromCamera,
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
    const Float4x4 clipFromWorld = clipFromCamera * cameraFromWorld;

    float maxLi = 0.0f;
    float sumLi = 0.0f;
    bool hadVisibleSegment = false;
    Vec3 p0 = Vec3(0.0f);
    if (!EvaluateLimitPosition(patchMap, patchTable, limitValues, ptexFaceId, uvStart, &p0))
    {
        if (nonUniformReasonOut)
        {
            *nonUniformReasonOut = DIAGSPLIT_NON_UNIFORM_EVAL_FAIL;
        }
        return SUBDIV_EDGE_FACTOR_NON_UNIFORM;
    }
    Vec4 prevClip = clipFromWorld * Vec4(p0.x, p0.y, p0.z, 1.0f);
    for (int i = 1; i < sampleStepsN; ++i)
    {
        const float t = float(i) / float(sampleStepsN - 1);
        const pxr::GfVec2f uv = uvStart * (1.0f - t) + uvEnd * t;
        Vec3 p = Vec3(0.0f);
        if (!EvaluateLimitPosition(patchMap, patchTable, limitValues, ptexFaceId, uv, &p))
        {
            if (nonUniformReasonOut)
            {
                *nonUniformReasonOut = DIAGSPLIT_NON_UNIFORM_EVAL_FAIL;
            }
            return SUBDIV_EDGE_FACTOR_NON_UNIFORM;
        }
        const Vec4 currClip = clipFromWorld * Vec4(p.x, p.y, p.z, 1.0f);
        Vec4 segA = prevClip;
        Vec4 segB = currClip;
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
    Vec3 p0 = Vec3(0.0f);
    if (!EvaluateLimitPosition(patchMap, patchTable, limitValues, ptexFaceId, uvStart, &p0))
    {
        return 0.f;
    }
    for (int i = 1; i < sampleStepsN; ++i)
    {
        const float t = float(i) / float(sampleStepsN - 1);
        const pxr::GfVec2f uv = uvStart * (1.0f - t) + uvEnd * t;
        Vec3 p = Vec3(0.0f);
        if (!EvaluateLimitPosition(patchMap, patchTable, limitValues, ptexFaceId, uv, &p))
        {
            return 0.f;
        }
        sumLi += Length(p - p0);
        p0 = p;
    }

    return sumLi;
}
