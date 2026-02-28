static int GetOrAllocateMidpointVertex(SubdivisionEdgeMap &edgeMap,
                                       int v0,
                                       int v1,
                                       int ptexFaceId,
                                       const pxr::GfVec2f &uv0,
                                       const pxr::GfVec2f &uv1,
                                       int &nextVertexId)
{
    SubdivisionEdge &edge = GetOrCreateEdge(edgeMap, v0, v1, ptexFaceId, uv0, uv1);
    if (edge.midpointVertex < 0)
    {
        edge.midpointVertex = nextVertexId++;
        edge.split = true;
    }
    return edge.midpointVertex;
}

static void FinalizeEdgeFlags(SubdivisionEdgeMap &edgeMap)
{
    for (auto &kv : edgeMap)
    {
        SubdivisionEdge &edge = kv.second;
        edge.boundary = (edge.faceCount == 1);
        edge.nonManifold = (edge.faceCount > 2);
    }
}

static std::vector<SubdivisionPatch> BuildSubdivisionPatches(const SubdivisionMesh &m,
                                                             const Far::TopologyRefiner &refiner,
                                                             SubdivisionEdgeMap &edgeMap,
                                                             int &nextVertexId)
{
    std::vector<SubdivisionPatch> patches;
    patches.reserve(m.vertsPerFace.size());

    Far::PtexIndices ptexIndices(refiner);
    std::vector<int> faceCenterVertex(m.vertsPerFace.size(), -1);

    size_t cursor = 0;
    for (size_t f = 0; f < m.vertsPerFace.size(); ++f)
    {
        const int n = m.vertsPerFace[f];
        if (n < 3 || cursor + size_t(n) > m.indices.size())
        {
            cursor += std::max(0, n);
            continue;
        }

        const int basePtexFaceId = ptexIndices.GetFaceId(int(f));
        if (n == 4)
        {
            SubdivisionPatch patch = {};
            patch.verts[0] = m.indices[cursor + 0];
            patch.verts[1] = m.indices[cursor + 1];
            patch.verts[2] = m.indices[cursor + 2];
            patch.verts[3] = m.indices[cursor + 3];
            patch.uv[0] = pxr::GfVec2f(0.0f, 0.0f);
            patch.uv[1] = pxr::GfVec2f(1.0f, 0.0f);
            patch.uv[2] = pxr::GfVec2f(1.0f, 1.0f);
            patch.uv[3] = pxr::GfVec2f(0.0f, 1.0f);
            patch.coarseFace = int(f);
            patch.quadrant = 0;
            patch.ptexFaceId = basePtexFaceId;
            patches.push_back(patch);
        }
        else
        {
            if (faceCenterVertex[f] < 0)
            {
                faceCenterVertex[f] = nextVertexId++;
            }
            const int center = faceCenterVertex[f];

            for (int i = 0; i < n; ++i)
            {
                const int v = m.indices[cursor + i];
                const int vNext = m.indices[cursor + ((i + 1) % n)];
                const int vPrev = m.indices[cursor + ((i + n - 1) % n)];
                const int midNext = GetOrAllocateMidpointVertex(edgeMap,
                                                                v,
                                                                vNext,
                                                                basePtexFaceId + i,
                                                                pxr::GfVec2f(0.0f, 0.0f),
                                                                pxr::GfVec2f(1.0f, 0.0f),
                                                                nextVertexId);
                const int midPrev = GetOrAllocateMidpointVertex(edgeMap,
                                                                vPrev,
                                                                v,
                                                                basePtexFaceId + i,
                                                                pxr::GfVec2f(0.0f, 1.0f),
                                                                pxr::GfVec2f(0.0f, 0.0f),
                                                                nextVertexId);
                SubdivisionEdge &edge = GetEdge(edgeMap, v, vNext);
                edge.tmaxEdgeFactor = SUBDIV_EDGE_FACTOR_NON_UNIFORM;

                SubdivisionPatch patch = {};
                patch.verts[0] = v;
                patch.verts[1] = midNext;
                patch.verts[2] = center;
                patch.verts[3] = midPrev;
                patch.uv[0] = pxr::GfVec2f(0.0f, 0.0f);
                patch.uv[1] = pxr::GfVec2f(1.0f, 0.0f);
                patch.uv[2] = pxr::GfVec2f(1.0f, 1.0f);
                patch.uv[3] = pxr::GfVec2f(0.0f, 1.0f);
                patch.coarseFace = int(f);
                patch.quadrant = i;
                patch.ptexFaceId = basePtexFaceId + i;
                const int patchId = int(patches.size());
                patches.push_back(patch);
            }
        }

        cursor += n;
    }

    FinalizeEdgeFlags(edgeMap);
    return patches;
}
