static float Dot(const pxr::GfVec3f &a, const pxr::GfVec3f &b)
{
    return a[0] * b[0] + a[1] * b[1] + a[2] * b[2];
}

static pxr::GfVec3f Cross(const pxr::GfVec3f &a, const pxr::GfVec3f &b)
{
    return pxr::GfVec3f(
        a[1] * b[2] - a[2] * b[1], a[2] * b[0] - a[0] * b[2], a[0] * b[1] - a[1] * b[0]);
}

static pxr::GfVec3f Normalize(const pxr::GfVec3f &v)
{
    const float lenSq = Dot(v, v);
    if (lenSq <= 1e-12f)
    {
        return pxr::GfVec3f(0.0f, 0.0f, 1.0f);
    }
    return v * (1.0f / std::sqrt(lenSq));
}

static std::vector<LimitEvalVertex> BuildLimitEvalVertices(const Far::TopologyRefiner &refiner,
                                                           const Far::PatchTable &patchTable,
                                                           const pxr::VtVec3fArray &coarsePoints)
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
                                  pxr::GfVec3f *outP)
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

    pxr::GfVec3f p(0.0f);
    for (int i = 0; i < cvs.size(); ++i)
    {
        p += limitValues[cvs[i]].p * pWeights[i];
    }
    *outP = p;
    return true;
}

static pxr::GfVec2f ProjectToScreen(const pxr::GfVec3f &p,
                                    const pxr::GfVec3f &eye,
                                    const pxr::GfVec3f &lookAt,
                                    int viewportWidth,
                                    int viewportHeight,
                                    float verticalFovDegrees)
{
    const pxr::GfVec3f forward = Normalize(lookAt - eye);
    pxr::GfVec3f worldUp(0.0f, 0.0f, 1.0f);
    if (std::abs(Dot(forward, worldUp)) > 0.999f)
    {
        worldUp = pxr::GfVec3f(0.0f, 1.0f, 0.0f);
    }
    const pxr::GfVec3f right = Normalize(Cross(forward, worldUp));
    const pxr::GfVec3f up = Normalize(Cross(right, forward));

    const pxr::GfVec3f v = p - eye;
    const float x = Dot(v, right);
    const float y = Dot(v, up);
    const float z = std::max(1e-6f, Dot(v, forward));

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
                                           const pxr::GfVec3f &eye,
                                           const pxr::GfVec3f &lookAt,
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
    pxr::GfVec3f p0(0.0f);
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
        pxr::GfVec3f p(0.0f);
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
    pxr::GfVec3f p0(0.0f);
    if (!EvaluateLimitPosition(patchMap, patchTable, limitValues, ptexFaceId, uvStart, &p0))
    {
        return 0.f;
    }
    for (int i = 1; i < sampleStepsN; ++i)
    {
        const float t = float(i) / float(sampleStepsN - 1);
        const pxr::GfVec2f uv = uvStart * (1.0f - t) + uvEnd * t;
        pxr::GfVec3f p(0.0f);
        if (!EvaluateLimitPosition(patchMap, patchTable, limitValues, ptexFaceId, uv, &p))
        {
            return 0.f;
        }
        sumLi += (p - p0).GetLength();
        p0 = p;
    }

    return sumLi;
}
