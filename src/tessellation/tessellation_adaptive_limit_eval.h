static std::vector<LimitEvalVertex> BuildLimitEvalVertices(const Far::TopologyRefiner &refiner,
                                                           const Far::PatchTable &patchTable,
                                                           const Array<float3> &coarsePoints)
{
    const int numRefinerVerts = refiner.GetNumVerticesTotal();
    const int numLocalPoints = patchTable.GetNumLocalPoints();
    std::vector<LimitEvalVertex> values(numRefinerVerts + numLocalPoints);

    const int numCoarseVerts = refiner.GetLevel(0).GetNumVertices();
    const int copyCount = std::min(numCoarseVerts, int(coarsePoints.size()));
    for (int i = 0; i < copyCount; ++i)
    {
        values[i].p = coarsePoints[i];
    }

    Far::PrimvarRefiner primvarRefiner(refiner);
    LimitEvalVertex *src = values.data();
    for (int level = 1; level < refiner.GetNumLevels(); ++level)
    {
        LimitEvalVertex *dst = src + refiner.GetLevel(level - 1).GetNumVertices();
        primvarRefiner.Interpolate(level, src, dst);
        src = dst;
    }

    if (numLocalPoints > 0)
    {
        patchTable.ComputeLocalPointValues(values.data(), values.data() + numRefinerVerts);
    }
    return values;
}

static bool EvaluateLimitPosition(const Far::PatchMap &patchMap,
                                  const Far::PatchTable &patchTable,
                                  const std::vector<LimitEvalVertex> &limitValues,
                                  int ptexFaceId,
                                  const pxr::GfVec2f &uv,
                                  float3 *outP)
{
    if (!outP)
    {
        return false;
    }
    const Far::PatchTable::PatchHandle *handle = patchMap.FindPatch(ptexFaceId, uv[0], uv[1]);
    if (!handle)
    {
        return false;
    }

    float pWeights[20] = {0.0f};
    patchTable.EvaluateBasis(*handle, uv[0], uv[1], pWeights);
    Far::ConstIndexArray cvs = patchTable.GetPatchVertices(*handle);

    float3 p = make_float3(0.0f);
    for (int i = 0; i < cvs.size(); ++i)
    {
        p += limitValues[cvs[i]].p * pWeights[i];
    }
    *outP = p;
    return true;
}

static pxr::GfVec2f ProjectToScreen(const float3 &p,
                                    const float3 &eye,
                                    const float3 &lookAt,
                                    int viewportWidth,
                                    int viewportHeight,
                                    float verticalFovDegrees)
{
    const float3 forward = normalize(lookAt - eye);
    float3 worldUp = make_float3(0.0f, 0.0f, 1.0f);
    if (std::abs(dot(forward, worldUp)) > 0.999f)
    {
        worldUp = make_float3(0.0f, 1.0f, 0.0f);
    }
    const float3 right = normalize(cross(forward, worldUp));
    const float3 up = normalize(cross(right, forward));

    const float3 v = p - eye;
    const float x = dot(v, right);
    const float y = dot(v, up);
    const float z = std::max(1e-6f, dot(v, forward));

    const float fovY = verticalFovDegrees * 3.14159265358979323846f / 180.0f;
    const float tanHalfFovY = std::tan(0.5f * fovY);
    const float aspect = float(viewportWidth) / float(viewportHeight);
    const float ndcX = x / (z * tanHalfFovY * aspect);
    const float ndcY = y / (z * tanHalfFovY);

    const float sx = (ndcX * 0.5f + 0.5f) * float(viewportWidth);
    const float sy = (1.0f - (ndcY * 0.5f + 0.5f)) * float(viewportHeight);
    return pxr::GfVec2f(sx, sy);
}

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
                                           float verticalFovDegrees)
{
    if (sampleStepsN < 2 || targetPixelSpacing <= 0.0f)
    {
        return 1;
    }

    float maxLi = 0.0f;
    float sumLi = 0.0f;
    float3 p0 = make_float3(0.0f);
    if (!EvaluateLimitPosition(patchMap, patchTable, limitValues, ptexFaceId, uvStart, &p0))
    {
        return SUBDIV_EDGE_FACTOR_NON_UNIFORM;
    }
    pxr::GfVec2f prev =
        ProjectToScreen(p0, eye, lookAt, viewportWidth, viewportHeight, verticalFovDegrees);
    for (int i = 1; i < sampleStepsN; ++i)
    {
        const float t = float(i) / float(sampleStepsN - 1);
        const pxr::GfVec2f uv = uvStart * (1.0f - t) + uvEnd * t;
        float3 p = make_float3(0.0f);
        if (!EvaluateLimitPosition(patchMap, patchTable, limitValues, ptexFaceId, uv, &p))
        {
            return SUBDIV_EDGE_FACTOR_NON_UNIFORM;
        }
        const pxr::GfVec2f s =
            ProjectToScreen(p, eye, lookAt, viewportWidth, viewportHeight, verticalFovDegrees);
        const pxr::GfVec2f d = s - prev;
        const float li = std::sqrt(d[0] * d[0] + d[1] * d[1]);
        sumLi += li;
        maxLi = std::max(maxLi, li);
        prev = s;
    }

    const int tMin = std::max(1, int(std::ceil(sumLi / targetPixelSpacing)));
    const int tMax =
        std::max(1, int(std::ceil(float(sampleStepsN - 1) * maxLi / targetPixelSpacing)));
    if ((tMax - tMin) >= splitThreshold)
    {
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
