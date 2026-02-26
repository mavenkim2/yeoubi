static SubdivisionEdge &GetOrCreateEdge(SubdivisionEdgeMap &edgeMap,
                                        int v0,
                                        int v1,
                                        int ptexFaceId,
                                        const pxr::GfVec2f &uv0,
                                        const pxr::GfVec2f &uv1)
{
    const uint64_t key = MakeEdgeKey(v0, v1);
    auto it = edgeMap.find(key);
    if (it == edgeMap.end())
    {
        SubdivisionEdge edge = {};
        edge.v0 = v0;
        edge.v1 = v1;
        edge.sampleVStart = v0;
        edge.sampleVEnd = v1;
        edge.storedPtexFaceId = ptexFaceId;
        edge.storedUv0 = uv0;
        edge.storedUv1 = uv1;
        edge.hasStoredPatchParams = (ptexFaceId >= 0);
        it = edgeMap.emplace(key, edge).first;
    }
    return it->second;
}

static SubdivisionEdge &GetEdge(SubdivisionEdgeMap &edgeMap, int v0, int v1)
{
    const uint64_t key = MakeEdgeKey(v0, v1);
    auto it = edgeMap.find(key);
    YBI_ASSERT(it != edgeMap.end());
    return it->second;
}

static const SubdivisionEdge &GetEdge(const SubdivisionEdgeMap &edgeMap, int v0, int v1)
{
    const uint64_t key = MakeEdgeKey(v0, v1);
    const auto it = edgeMap.find(key);
    YBI_ASSERT(it != edgeMap.end());
    return it->second;
}

static void SetEdgeSampleParams(SubdivisionEdge &edge,
                                int vStart,
                                int vEnd,
                                int ptexFaceId,
                                const pxr::GfVec2f &uvStart,
                                const pxr::GfVec2f &uvEnd,
                                bool overwrite)
{
    if (!overwrite && edge.hasStoredPatchParams && edge.sampleVStart >= 0 && edge.sampleVEnd >= 0)
    {
        return;
    }
    edge.sampleVStart = vStart;
    edge.sampleVEnd = vEnd;
    edge.storedPtexFaceId = ptexFaceId;
    edge.storedUv0 = uvStart;
    edge.storedUv1 = uvEnd;
    edge.hasStoredPatchParams = (ptexFaceId >= 0);
}

static int GetEdgeFactor(const SubdivisionEdgeMap &edgeMap, int v0, int v1)
{
    const uint64_t key = MakeEdgeKey(v0, v1);
    const auto it = edgeMap.find(key);
    YBI_ERROR(it != edgeMap.end(), "Missing edge factor lookup for edge (%d, %d)\n", v0, v1);
    return it->second.tmaxEdgeFactor;
}

static int EnsurePatchEdgeFactor(const SelectedSubdivMesh &m,
                                 const SubdivisionPatch &patch,
                                 int edgeIndex,
                                 SubdivisionEdgeMap &edgeMap,
                                 int &nextGeneratedVertexId,
                                 bool allowQuadCoarseUpgrade,
                                 const Far::PatchMap &patchMap,
                                 const Far::PatchTable &patchTable,
                                 const std::vector<LimitEvalVertex> &limitValues,
                                 int sampleStepsN,
                                 float targetPixelSpacing,
                                 int splitThreshold,
                                 const pxr::GfVec3f &eye,
                                 const pxr::GfVec3f &lookAt,
                                 int viewportWidth,
                                 int viewportHeight,
                                 float verticalFovDegrees,
                                 int *computedCountOut)
{
    const int next = (edgeIndex + 1) & 3;
    SubdivisionEdge &edge = GetOrCreateEdge(edgeMap,
                                            patch.verts[edgeIndex],
                                            patch.verts[next],
                                            patch.ptexFaceId,
                                            patch.uv[edgeIndex],
                                            patch.uv[next]);
    if (edge.tmaxEdgeFactor == SUBDIV_EDGE_FACTOR_UNINITIALIZED)
    {
        edge.tmaxEdgeFactor = ComputeDiagSplitPatchEdgeFactor(patchMap,
                                                              patchTable,
                                                              limitValues,
                                                              patch.ptexFaceId,
                                                              patch.uv[edgeIndex],
                                                              patch.uv[next],
                                                              sampleStepsN,
                                                              targetPixelSpacing,
                                                              splitThreshold,
                                                              eye,
                                                              lookAt,
                                                              viewportWidth,
                                                              viewportHeight,
                                                              verticalFovDegrees);
        if (edge.tmaxEdgeFactor != SUBDIV_EDGE_FACTOR_NON_UNIFORM && edge.tmaxEdgeFactor < 1)
        {
            edge.tmaxEdgeFactor = 1;
        }
        if (edge.tmaxEdgeFactor >= 1)
        {
            edge.transitionedUninitializedToUniform = true;
            SetEdgeSampleParams(edge,
                                patch.verts[edgeIndex],
                                patch.verts[next],
                                patch.ptexFaceId,
                                patch.uv[edgeIndex],
                                patch.uv[next],
                                /*overwrite*/ true);
            if (edge.tmaxEdgeFactor >= 2)
            {
                const int splitPosStored = edge.tmaxEdgeFactor / 2;
                edge.edgeVertexIndexStart = nextGeneratedVertexId;
                nextGeneratedVertexId += (edge.tmaxEdgeFactor - 1);
                edge.midpointVertex = edge.edgeVertexIndexStart + (splitPosStored - 1);
            }
            else
            {
                edge.edgeVertexIndexStart = -1;
                edge.midpointVertex = -1;
            }
        }
        else
        {
            edge.edgeVertexIndexStart = -1;
        }
        if (computedCountOut)
        {
            (*computedCountOut)++;
        }
    }
    return edge.tmaxEdgeFactor;
}

static std::vector<SubdivisionPatch>
DiagSplitPatches(const SelectedSubdivMesh &m,
                 const std::vector<SubdivisionPatch> &patches,
                 SubdivisionEdgeMap &edgeMap,
                 int &nextGeneratedVertexId,
                 const Far::PatchMap &patchMap,
                 const Far::PatchTable &patchTable,
                 const std::vector<LimitEvalVertex> &limitValues,
                 int sampleStepsN,
                 float targetPixelSpacing,
                 int splitThreshold,
                 const pxr::GfVec3f &eye,
                 const pxr::GfVec3f &lookAt,
                 int viewportWidth,
                 int viewportHeight,
                 float verticalFovDegrees,
                 int *computedCountOut)
{
    std::vector<SubdivisionPatch> worklist = patches;
    std::vector<SubdivisionPatch> out;
    out.reserve(patches.size());

    size_t patchCursor = 0;
    while (patchCursor < worklist.size())
    {
        const SubdivisionPatch patch = worklist[patchCursor++];

        int edgeFactor[4] = {SUBDIV_EDGE_FACTOR_UNINITIALIZED,
                             SUBDIV_EDGE_FACTOR_UNINITIALIZED,
                             SUBDIV_EDGE_FACTOR_UNINITIALIZED,
                             SUBDIV_EDGE_FACTOR_UNINITIALIZED};
        int splitEdge = -1;
        for (int e = 0; e < 4; ++e)
        {
            edgeFactor[e] = EnsurePatchEdgeFactor(m,
                                                  patch,
                                                  e,
                                                  edgeMap,
                                                  nextGeneratedVertexId,
                                                  /*allowQuadCoarseUpgrade*/ false,
                                                  patchMap,
                                                  patchTable,
                                                  limitValues,
                                                  sampleStepsN,
                                                  targetPixelSpacing,
                                                  splitThreshold,
                                                  eye,
                                                  lookAt,
                                                  viewportWidth,
                                                  viewportHeight,
                                                  verticalFovDegrees,
                                                  computedCountOut);
        }

        bool splitU = edgeFactor[0] == SUBDIV_EDGE_FACTOR_NON_UNIFORM || edgeFactor[2] == SUBDIV_EDGE_FACTOR_NON_UNIFORM;
        bool splitV = edgeFactor[1] == SUBDIV_EDGE_FACTOR_NON_UNIFORM || edgeFactor[3] == SUBDIV_EDGE_FACTOR_NON_UNIFORM;

        if (splitU && splitV)
        {
            float l0 = ComputeLength(patchMap, patchTable, limitValues, patch.ptexFaceId, patch.uv[0], patch.uv[1], sampleStepsN);
            float l1 = ComputeLength(patchMap, patchTable, limitValues, patch.ptexFaceId, patch.uv[1], patch.uv[2], sampleStepsN);
            float l2 = ComputeLength(patchMap, patchTable, limitValues, patch.ptexFaceId, patch.uv[2], patch.uv[3], sampleStepsN);
            float l3 = ComputeLength(patchMap, patchTable, limitValues, patch.ptexFaceId, patch.uv[3], patch.uv[0], sampleStepsN);

            if (l0 + l2 > l1 + l3)
            {
                splitV = false;
            }
            else 
            {
                splitU = false;
            }
        }

        if (splitU)
        {
            splitEdge = edgeFactor[0] == SUBDIV_EDGE_FACTOR_NON_UNIFORM ? 0 : 2;
        }
        else if (splitV)
        {
            splitEdge = edgeFactor[1] == SUBDIV_EDGE_FACTOR_NON_UNIFORM ? 1 : 3;
        }

        if (splitEdge < 0)
        {
            out.push_back(patch);
            continue;
        }

        const int oppositeEdge = (splitEdge + 2) & 3;
        if (edgeFactor[oppositeEdge] == 1)
        {
            const int ov0 = patch.verts[oppositeEdge];
            const int ov1 = patch.verts[(oppositeEdge + 1) & 3];
            const int ovn = (oppositeEdge + 1) & 3;
            SubdivisionEdge &opp = GetOrCreateEdge(
                edgeMap, ov0, ov1, patch.ptexFaceId, patch.uv[oppositeEdge], patch.uv[ovn]);
            opp.tmaxEdgeFactor = 2;
            SetEdgeSampleParams(
                opp, ov0, ov1, patch.ptexFaceId, patch.uv[oppositeEdge], patch.uv[ovn], true);
            if (opp.edgeVertexIndexStart < 0 || opp.midpointVertex < 0)
            {
                opp.edgeVertexIndexStart = nextGeneratedVertexId++;
                opp.midpointVertex = opp.edgeVertexIndexStart;
            }
            edgeFactor[oppositeEdge] = 2;
        }

        const int splitPair[2] = {splitEdge, oppositeEdge};
        int splitVerts[2] = {-1, -1};
        pxr::GfVec2f splitUV[2];
        for (int i = 0; i < 2; ++i)
        {
            const int e = splitPair[i];
            const int n = (e + 1) & 3;
            const int v0 = patch.verts[e];
            const int v1 = patch.verts[n];
            SubdivisionEdge &edge =
                GetOrCreateEdge(edgeMap, v0, v1, patch.ptexFaceId, patch.uv[e], patch.uv[n]);

            int mid = edge.midpointVertex;
            int t = edgeFactor[e];
            float alpha = 0.5f;

            if (t != SUBDIV_EDGE_FACTOR_NON_UNIFORM)
            {
                YBI_ASSERT(t > 1);
                const int splitPosStored = t / 2;
                alpha = float(splitPosStored) / float(t);
                splitUV[i] = patch.uv[e] * (1.0f - alpha) + patch.uv[n] * alpha;
                YBI_ASSERT(mid >= 0);
                YBI_ASSERT(edge.edgeVertexIndexStart != -1);

                if (!edge.split)
                {
                    edge.split = true;

                    // NOTE: these are technically wrong when child edge factor is 1, but the
                    // lookup won't use edgeVertexIndexStart in that case.
                    SubdivisionEdge &edgeA = GetOrCreateEdge(
                        edgeMap, edge.v0, mid, patch.ptexFaceId, patch.uv[e], splitUV[i]);
                    SubdivisionEdge &edgeB = GetOrCreateEdge(
                        edgeMap, mid, edge.v1, patch.ptexFaceId, splitUV[i], patch.uv[n]);
                    edgeA.tmaxEdgeFactor = splitPosStored;
                    edgeA.edgeVertexIndexStart = edge.edgeVertexIndexStart;
                    edgeA.midpointVertex =
                        edgeA.tmaxEdgeFactor > 1
                            ? edgeA.edgeVertexIndexStart + (edgeA.tmaxEdgeFactor / 2) - 1
                            : -1;

                    edgeB.tmaxEdgeFactor = t - splitPosStored;
                    edgeB.edgeVertexIndexStart = edge.edgeVertexIndexStart + splitPosStored;
                    edgeB.midpointVertex =
                        edgeB.tmaxEdgeFactor > 1
                            ? edgeB.edgeVertexIndexStart + (edgeB.tmaxEdgeFactor / 2) - 1
                            : -1;

                    edgeA.boundary = edge.boundary;
                    edgeB.boundary = edge.boundary;
                }
            }
            else
            {
                splitUV[i] = patch.uv[e] * (1.0f - alpha) + patch.uv[n] * alpha;
                if (!edge.split)
                {
                    edge.split = true;
                    mid = nextGeneratedVertexId++;
                    edge.midpointVertex = mid;
                    SubdivisionEdge &edgeA = GetOrCreateEdge(
                        edgeMap, v0, mid, patch.ptexFaceId, patch.uv[e], splitUV[i]);
                    SubdivisionEdge &edgeB = GetOrCreateEdge(
                        edgeMap, mid, v1, patch.ptexFaceId, splitUV[i], patch.uv[n]);
                    edgeA.tmaxEdgeFactor = SUBDIV_EDGE_FACTOR_UNINITIALIZED;
                    edgeB.tmaxEdgeFactor = SUBDIV_EDGE_FACTOR_UNINITIALIZED;
                    edgeA.boundary = edge.boundary;
                    edgeB.boundary = edge.boundary;
                }
            }

            splitVerts[i] = mid;
        }

        SubdivisionEdge &seam = GetOrCreateEdge(
            edgeMap, splitVerts[0], splitVerts[1], patch.ptexFaceId, splitUV[0], splitUV[1]);
        seam.tmaxEdgeFactor = SUBDIV_EDGE_FACTOR_UNINITIALIZED;
        seam.edgeVertexIndexStart = -1;
        seam.boundary = false;

        const int i0 = splitEdge;
        const int i1 = (i0 + 1) & 3;
        const int i2 = (i0 + 2) & 3;
        const int i3 = (i0 + 3) & 3;
        const int midA = splitVerts[0];
        const int midB = splitVerts[1];

        SubdivisionPatch childA = patch;
        childA.verts[0] = patch.verts[i0];
        childA.verts[1] = midA;
        childA.verts[2] = midB;
        childA.verts[3] = patch.verts[i3];
        childA.uv[0] = patch.uv[i0];
        childA.uv[1] = splitUV[0];
        childA.uv[2] = splitUV[1];
        childA.uv[3] = patch.uv[i3];

        SubdivisionPatch childB = patch;
        childB.verts[0] = midA;
        childB.verts[1] = patch.verts[i1];
        childB.verts[2] = patch.verts[i2];
        childB.verts[3] = midB;
        childB.uv[0] = splitUV[0];
        childB.uv[1] = patch.uv[i1];
        childB.uv[2] = patch.uv[i2];
        childB.uv[3] = splitUV[1];

        worklist.push_back(childA);
        worklist.push_back(childB);
    }

    return out;
}
