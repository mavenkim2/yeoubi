#include "util/assert.h"
static SubdivisionEdge &GetOrCreateEdge(SubdivisionEdgeMap &edgeMap,
                                        int v0,
                                        int v1,
                                        int depth,
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
        edge.depth = depth;
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

static int EnsurePatchEdgeFactor(const SubdivisionPatch &patch,
                                 int edgeIndex,
                                 SubdivisionEdgeMap &edgeMap,
                                 int &nextGeneratedVertexId,
                                 int maxDiagSplitDepth,
                                 bool allowQuadCoarseUpgrade,
                                 const Far::PatchMap &patchMap,
                                 const Far::PatchTable &patchTable,
                                 const std::vector<LimitEvalVertex> &limitValues,
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
                                 int *computedCountOut)
{
    auto nonUniformReasonLabel = [](DiagSplitNonUniformReason reason) {
        switch (reason)
        {
            case DIAGSPLIT_NON_UNIFORM_EVAL_FAIL:
                return "eval-fail";
            case DIAGSPLIT_NON_UNIFORM_VARIANCE_THRESHOLD:
                return "variance-threshold";
            default:
                return "unknown";
        }
    };
    const int next = (edgeIndex + 1) & 3;
    SubdivisionEdge &edge = GetOrCreateEdge(edgeMap,
                                            patch.verts[edgeIndex],
                                            patch.verts[next],
                                            /*depth*/ 0,
                                            patch.ptexFaceId,
                                            patch.uv[edgeIndex],
                                            patch.uv[next]);
    const int clampedDepth = std::max(0, edge.depth);
    const int remainingDepth = std::max(0, maxDiagSplitDepth - clampedDepth);
    const int maxFactorAtDepth =
        remainingDepth >= 30 ? std::numeric_limits<int>::max() : (1 << remainingDepth);
    if (edge.tmaxEdgeFactor == SUBDIV_EDGE_FACTOR_UNINITIALIZED)
    {
        DiagSplitNonUniformReason nonUniformReason = DIAGSPLIT_NON_UNIFORM_NONE;
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
                                                              verticalFovDegrees,
                                                              useCameraMatrices,
                                                              cameraFromWorld,
                                                              clipFromCamera,
                                                              &nonUniformReason);
        if (edge.tmaxEdgeFactor == SUBDIV_EDGE_FACTOR_NON_UNIFORM)
        {
            if (nonUniformReason != DIAGSPLIT_NON_UNIFORM_VARIANCE_THRESHOLD)
            {
                std::fprintf(stderr,
                             "DiagSplit non-uniform edge: reason=%s ptex=%d edge=%d v=(%d,%d)"
                             " uv0=(%.9g,%.9g) uv1=(%.9g,%.9g)\n",
                             nonUniformReasonLabel(nonUniformReason),
                             patch.ptexFaceId,
                             edgeIndex,
                             patch.verts[edgeIndex],
                             patch.verts[next],
                             double(patch.uv[edgeIndex][0]),
                             double(patch.uv[edgeIndex][1]),
                             double(patch.uv[next][0]),
                             double(patch.uv[next][1]));
            }
        }
        if (edge.tmaxEdgeFactor == SUBDIV_EDGE_FACTOR_NON_UNIFORM &&
            clampedDepth >= maxDiagSplitDepth)
        {
            edge.tmaxEdgeFactor = 1;
        }
        if (edge.tmaxEdgeFactor != SUBDIV_EDGE_FACTOR_NON_UNIFORM && edge.tmaxEdgeFactor < 1)
        {
            edge.tmaxEdgeFactor = 1;
        }
        if (edge.tmaxEdgeFactor != SUBDIV_EDGE_FACTOR_NON_UNIFORM)
        {
            edge.tmaxEdgeFactor = std::min(edge.tmaxEdgeFactor, maxFactorAtDepth);
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
                YBI_ASSERT(edge.midpointVertex != edge.v0);
                YBI_ASSERT(edge.midpointVertex != edge.v1);
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
    if (edge.tmaxEdgeFactor == SUBDIV_EDGE_FACTOR_NON_UNIFORM && clampedDepth >= maxDiagSplitDepth)
    {
        edge.tmaxEdgeFactor = 1;
        edge.transitionedUninitializedToUniform = true;
        SetEdgeSampleParams(edge,
                            patch.verts[edgeIndex],
                            patch.verts[next],
                            patch.ptexFaceId,
                            patch.uv[edgeIndex],
                            patch.uv[next],
                            /*overwrite*/ true);
        edge.edgeVertexIndexStart = -1;
        edge.midpointVertex = -1;
    }
    // if (edge.tmaxEdgeFactor != SUBDIV_EDGE_FACTOR_NON_UNIFORM)
    // {
    //     edge.tmaxEdgeFactor = std::max(1, std::min(edge.tmaxEdgeFactor, maxFactorAtDepth));
    //     if (edge.tmaxEdgeFactor <= 1)
    //     {
    //         edge.edgeVertexIndexStart = -1;
    //         edge.midpointVertex = -1;
    //     }
    // }
    return edge.tmaxEdgeFactor;
}

static std::vector<SubdivisionPatch>
DiagSplitPatches(const std::vector<SubdivisionPatch> &patches,
                 SubdivisionEdgeMap &edgeMap,
                 int &nextGeneratedVertexId,
                 int maxDiagSplitDepth,
                 const Far::PatchMap &patchMap,
                 const Far::PatchTable &patchTable,
                 const std::vector<LimitEvalVertex> &limitValues,
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
            edgeFactor[e] = EnsurePatchEdgeFactor(patch,
                                                  e,
                                                  edgeMap,
                                                  nextGeneratedVertexId,
                                                  maxDiagSplitDepth,
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
                                                  useCameraMatrices,
                                                  cameraFromWorld,
                                                  clipFromCamera,
                                                  computedCountOut);
        }

        bool splitU = edgeFactor[0] == SUBDIV_EDGE_FACTOR_NON_UNIFORM ||
                      edgeFactor[2] == SUBDIV_EDGE_FACTOR_NON_UNIFORM;
        bool splitV = edgeFactor[1] == SUBDIV_EDGE_FACTOR_NON_UNIFORM ||
                      edgeFactor[3] == SUBDIV_EDGE_FACTOR_NON_UNIFORM;

        if (splitU && splitV)
        {
            float l0 = ComputeLength(patchMap,
                                     patchTable,
                                     limitValues,
                                     patch.ptexFaceId,
                                     patch.uv[0],
                                     patch.uv[1],
                                     sampleStepsN);
            float l1 = ComputeLength(patchMap,
                                     patchTable,
                                     limitValues,
                                     patch.ptexFaceId,
                                     patch.uv[1],
                                     patch.uv[2],
                                     sampleStepsN);
            float l2 = ComputeLength(patchMap,
                                     patchTable,
                                     limitValues,
                                     patch.ptexFaceId,
                                     patch.uv[2],
                                     patch.uv[3],
                                     sampleStepsN);
            float l3 = ComputeLength(patchMap,
                                     patchTable,
                                     limitValues,
                                     patch.ptexFaceId,
                                     patch.uv[3],
                                     patch.uv[0],
                                     sampleStepsN);

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
            SubdivisionEdge &opp = GetEdge(edgeMap, ov0, ov1);
            opp.tmaxEdgeFactor = 2;
            // TODO: unsure if these edges' limit positions are always evaluated. assert should
            // hopefully catch.
            SetEdgeSampleParams(
                opp, ov0, ov1, patch.ptexFaceId, patch.uv[oppositeEdge], patch.uv[ovn], true);
            if (opp.edgeVertexIndexStart < 0 || opp.midpointVertex < 0)
            {
                opp.edgeVertexIndexStart = nextGeneratedVertexId++;
                opp.midpointVertex = opp.edgeVertexIndexStart;
                YBI_ASSERT(opp.midpointVertex != opp.v0);
                YBI_ASSERT(opp.midpointVertex != opp.v1);
            }
            edgeFactor[oppositeEdge] = 2;
        }

        const int splitPair[2] = {splitEdge, oppositeEdge};
        YBI_ASSERT(edgeFactor[splitEdge] == SUBDIV_EDGE_FACTOR_NON_UNIFORM);
        int splitVerts[2] = {-1, -1};
        pxr::GfVec2f splitUV[2];
        for (int i = 0; i < 2; ++i)
        {
            const int e = splitPair[i];
            const int n = (e + 1) & 3;
            const int v0 = patch.verts[e];
            const int v1 = patch.verts[n];
            SubdivisionEdge &edge = GetEdge(edgeMap, v0, v1);

            int mid = edge.midpointVertex;
            YBI_ASSERT(v0 != mid);
            YBI_ASSERT(v1 != mid);

            int t = edgeFactor[e];
            float alpha = 0.5f;
            if (t != SUBDIV_EDGE_FACTOR_NON_UNIFORM)
            {
                YBI_ASSERT(t > 1);
                const int splitPosStored = t / 2;
                const int splitPosPatch = edge.v0 == v0 ? splitPosStored : t - splitPosStored;

                alpha = float(splitPosPatch) / float(t);
                splitUV[i] = patch.uv[e] * (1.0f - alpha) + patch.uv[n] * alpha;
                YBI_ERROR(mid >= 0, "%i %i\n", mid, edge.edgeVertexIndexStart);
                YBI_ASSERT(edge.edgeVertexIndexStart != -1);

                if (!edge.split)
                {
                    edge.split = true;
                    const int parentV0 = edge.v0;
                    const int parentV1 = edge.v1;
                    const int parentEdgeVertexIndexStart = edge.edgeVertexIndexStart;
                    const bool parentBoundary = edge.boundary;
                    YBI_ASSERT(parentEdgeVertexIndexStart >= 0);

                    // NOTE: these are technically wrong when child edge factor is 1, but the
                    // lookup won't use edgeVertexIndexStart in that case.
                    SubdivisionEdge &edgeA =
                        GetOrCreateEdge(edgeMap,
                                        parentV0,
                                        mid,
                                        edge.depth + 1,
                                        patch.ptexFaceId,
                                        parentV0 == v0 ? patch.uv[e] : patch.uv[n],
                                        splitUV[i]);
                    SubdivisionEdge &edgeB =
                        GetOrCreateEdge(edgeMap,
                                        mid,
                                        parentV1,
                                        edge.depth + 1,
                                        patch.ptexFaceId,
                                        splitUV[i],
                                        parentV0 == v0 ? patch.uv[n] : patch.uv[e]);
                    edgeA.tmaxEdgeFactor = splitPosStored;
                    edgeA.edgeVertexIndexStart = parentEdgeVertexIndexStart;
                    edgeA.midpointVertex =
                        edgeA.tmaxEdgeFactor > 1
                            ? edgeA.edgeVertexIndexStart + (edgeA.tmaxEdgeFactor / 2) - 1
                            : -1;
                    if (edgeA.midpointVertex >= 0)
                    {
                        YBI_ASSERT(edgeA.midpointVertex != edgeA.v0);
                        YBI_ASSERT(edgeA.midpointVertex != edgeA.v1);
                    }

                    YBI_ASSERT(edgeA.v0 != edgeA.v1);

                    edgeB.tmaxEdgeFactor = t - splitPosStored;
                    edgeB.edgeVertexIndexStart = parentEdgeVertexIndexStart + splitPosStored;
                    edgeB.midpointVertex =
                        edgeB.tmaxEdgeFactor > 1
                            ? edgeB.edgeVertexIndexStart + (edgeB.tmaxEdgeFactor / 2) - 1
                            : -1;
                    if (edgeB.midpointVertex >= 0)
                    {
                        YBI_ERROR(edgeB.midpointVertex != edgeB.v0 &&
                                      edgeB.midpointVertex != edgeB.v1,
                                  "Invalid edgeB midpoint: ptex=%d splitEdge=%d "
                                  "pairIdx=%d edgeIndex=%d t=%d splitPosStored=%d "
                                  "parent=(v0=%d v1=%d start=%d mid=%d depth=%d) "
                                  "edgeA=(v0=%d v1=%d t=%d start=%d mid=%d depth=%d) "
                                  "edgeB=(v0=%d v1=%d t=%d start=%d mid=%d depth=%d)",
                                  patch.ptexFaceId,
                                  splitEdge,
                                  i,
                                  e,
                                  t,
                                  splitPosStored,
                                  parentV0,
                                  parentV1,
                                  parentEdgeVertexIndexStart,
                                  edge.midpointVertex,
                                  edge.depth,
                                  edgeA.v0,
                                  edgeA.v1,
                                  edgeA.tmaxEdgeFactor,
                                  edgeA.edgeVertexIndexStart,
                                  edgeA.midpointVertex,
                                  edgeA.depth,
                                  edgeB.v0,
                                  edgeB.v1,
                                  edgeB.tmaxEdgeFactor,
                                  edgeB.edgeVertexIndexStart,
                                  edgeB.midpointVertex,
                                  edgeB.depth);
                    }

                    YBI_ASSERT(edgeB.v0 != edgeB.v1);

                    edgeA.boundary = parentBoundary;
                    edgeB.boundary = parentBoundary;
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
                    YBI_ASSERT(edge.midpointVertex != edge.v0);
                    YBI_ASSERT(edge.midpointVertex != edge.v1);
                    SubdivisionEdge &edgeA = GetOrCreateEdge(edgeMap,
                                                             v0,
                                                             mid,
                                                             edge.depth + 1,
                                                             patch.ptexFaceId,
                                                             patch.uv[e],
                                                             splitUV[i]);
                    SubdivisionEdge &edgeB = GetOrCreateEdge(edgeMap,
                                                             mid,
                                                             v1,
                                                             edge.depth + 1,
                                                             patch.ptexFaceId,
                                                             splitUV[i],
                                                             patch.uv[n]);
                    edgeA.tmaxEdgeFactor = SUBDIV_EDGE_FACTOR_UNINITIALIZED;
                    edgeB.tmaxEdgeFactor = SUBDIV_EDGE_FACTOR_UNINITIALIZED;
                    edgeA.boundary = edge.boundary;
                    edgeB.boundary = edge.boundary;
                }
            }

            splitVerts[i] = mid;
        }

        const int splitDepth0 =
            splitPair[0] >= 0
                ? GetEdge(edgeMap, patch.verts[splitPair[0]], patch.verts[(splitPair[0] + 1) & 3])
                      .depth
                : 0;
        const int splitDepth1 =
            splitPair[1] >= 0
                ? GetEdge(edgeMap, patch.verts[splitPair[1]], patch.verts[(splitPair[1] + 1) & 3])
                      .depth
                : 0;
        const int seamDepth = std::max(splitDepth0, splitDepth1) + 1;
        SubdivisionEdge &seam = GetOrCreateEdge(edgeMap,
                                                splitVerts[0],
                                                splitVerts[1],
                                                seamDepth,
                                                patch.ptexFaceId,
                                                splitUV[0],
                                                splitUV[1]);
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
