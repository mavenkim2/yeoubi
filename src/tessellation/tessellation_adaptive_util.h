struct CreasePairs
{
    std::vector<int> pairs;
    std::vector<float> weights;
};

static CreasePairs BuildCreasePairs(const SubdivisionMesh &m)
{
    CreasePairs out = {};
    size_t c = 0;
    for (size_t i = 0; i < m.creaseLengths.size(); i++)
    {
        const int len = m.creaseLengths[i];
        if (len < 2 || c + size_t(len) > m.creaseIndices.size())
        {
            c += std::max(0, len);
            continue;
        }
        const float w = (i < m.creaseSharpnesses.size()) ? m.creaseSharpnesses[i] : 0.0f;
        for (int j = 0; j + 1 < len; j++)
        {
            out.pairs.push_back(m.creaseIndices[c + j]);
            out.pairs.push_back(m.creaseIndices[c + j + 1]);
            out.weights.push_back(w);
        }
        c += len;
    }
    return out;
}

static Sdc::SchemeType SchemeFromString(const std::string &s)
{
    if (s == "loop")
    {
        return Sdc::SCHEME_LOOP;
    }
    if (s == "bilinear")
    {
        return Sdc::SCHEME_BILINEAR;
    }
    return Sdc::SCHEME_CATMARK;
}

static Sdc::Options::CreasingMethod ToSdcCreasingMethod(SubdivisionCreasingMethod method)
{
    if (method == SUBDIVISION_CREASING_CHAIKIN)
    {
        return Sdc::Options::CREASE_CHAIKIN;
    }
    return Sdc::Options::CREASE_UNIFORM;
}

static Sdc::Options::TriangleSubdivision ToSdcTriangleSubdivision(TriangleSubdivisionRule rule)
{
    if (rule == TRIANGLE_SUBDIVISION_SMOOTH)
    {
        return Sdc::Options::TRI_SUB_SMOOTH;
    }
    return Sdc::Options::TRI_SUB_CATMARK;
}

static int CountNonManifoldEdges(const SubdivisionEdgeMap &edgeMap)
{
    int count = 0;
    for (const auto &it : edgeMap)
    {
        if (it.second.nonManifold)
        {
            count++;
        }
    }
    return count;
}

static int CountBoundaryEdges(const SubdivisionEdgeMap &edgeMap)
{
    int count = 0;
    for (const auto &it : edgeMap)
    {
        if (it.second.boundary)
        {
            count++;
        }
    }
    return count;
}

static int CountEdgesWithMidpointVertex(const SubdivisionEdgeMap &edgeMap)
{
    int count = 0;
    for (const auto &it : edgeMap)
    {
        if (it.second.midpointVertex >= 0)
        {
            count++;
        }
    }
    return count;
}

static int CountEdgesWithComputedTMax(const SubdivisionEdgeMap &edgeMap)
{
    int count = 0;
    for (const auto &it : edgeMap)
    {
        if (it.second.tmaxEdgeFactor != SUBDIV_EDGE_FACTOR_UNINITIALIZED)
        {
            count++;
        }
    }
    return count;
}

static bool VerifyInitializedUniformEdgesHaveStoredPatchParams(const SubdivisionEdgeMap &edgeMap,
                                                               int *missingStoredParamsOut,
                                                               int *badUniformFactorOut)
{
    int missingStoredParams = 0;
    int badUniformFactor = 0;
    for (const auto &it : edgeMap)
    {
        const SubdivisionEdge &edge = it.second;
        if (!edge.transitionedUninitializedToUniform)
        {
            continue;
        }
        if (edge.tmaxEdgeFactor < 1)
        {
            badUniformFactor++;
        }
        if (!edge.hasStoredPatchParams || edge.storedPtexFaceId < 0 || edge.sampleVStart < 0 ||
            edge.sampleVEnd < 0)
        {
            missingStoredParams++;
        }
    }

    if (missingStoredParamsOut)
    {
        *missingStoredParamsOut = missingStoredParams;
    }
    if (badUniformFactorOut)
    {
        *badUniformFactorOut = badUniformFactor;
    }
    return (missingStoredParams == 0) && (badUniformFactor == 0);
}
