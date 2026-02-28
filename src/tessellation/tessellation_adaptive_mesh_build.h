#pragma once

static void RecurseChildEdges(
    const SubdivisionEdgeMap &edgeMap, int v0, int v1, std::vector<int> &vertices, int depth = 0)
{
    const SubdivisionEdge &edge = GetEdge(edgeMap, v0, v1);
    if (edge.split)
    {
        YBI_ASSERT(edge.midpointVertex >= 0);
        RecurseChildEdges(edgeMap, v0, edge.midpointVertex, vertices, depth + 1);
        RecurseChildEdges(edgeMap, edge.midpointVertex, v1, vertices, depth + 1);
    }
    else
    {
        YBI_ERROR(vertices.back() == v0 || vertices.back() == v1,
                  "back: %i, edge: %i %i, depth: %i\n",
                  vertices.back(),
                  v0,
                  v1,
                  depth);
        YBI_ASSERT(edge.tmaxEdgeFactor >= 1);
        YBI_ASSERT(edge.edgeVertexIndexStart >= 0 || edge.tmaxEdgeFactor == 1);
        bool forward = vertices.back() == v0;
        YBI_ASSERT(forward);
        bool edgeForward = edge.v0 == v0;

        for (int k = 1; k < edge.tmaxEdgeFactor; k++)
        {
            const int localOffset = edgeForward ? (k - 1) : (edge.tmaxEdgeFactor - 1 - k);
            int vertex = edge.edgeVertexIndexStart + localOffset;
            vertices.push_back(vertex);
        }
        vertices.push_back(forward ? v1 : v0);
    }
}

static bool BuildLeafPatchStitchedMesh(const std::vector<SubdivisionPatch> &leafPatches,
                                       const SubdivisionEdgeMap &edgeMap,
                                       const Far::PatchMap &patchMap,
                                       const Far::PatchTable &patchTable,
                                       const std::vector<LimitEvalVertex> &limitValues,
                                       int nextGeneratedVertexId,
                                       Mesh *outMesh,
                                       std::vector<int> *outTriPatchFaceIds,
                                       std::vector<int> *outTriCoarseFaceIds,
                                       std::vector<int> *outTriPtexFaceIds,
                                       std::vector<int> *outTriQuadrants,
                                       int *outVertexCount,
                                       int *outTriangleCount)
{
    if (!outMesh)
    {
        return false;
    }

    if (outVertexCount)
    {
        *outVertexCount = 0;
    }
    if (outTriangleCount)
    {
        *outTriangleCount = 0;
    }

    int maxVertexId = nextGeneratedVertexId;

    const float inf = std::numeric_limits<float>::infinity();
    std::vector<float3> sharedPos(maxVertexId, make_float3(inf, inf, inf));
    std::vector<uint8_t> sharedInit(maxVertexId, uint8_t(0));
    auto isSharedInit = [&](int id) -> bool {
        return (id >= 0) && (id < int(sharedInit.size())) && (sharedInit[id] != 0);
    };

    for (const auto &it : edgeMap)
    {
        const SubdivisionEdge &edge = it.second;
        if (!edge.transitionedUninitializedToUniform)
        {
            continue;
        }

        YBI_ASSERT(edge.hasStoredPatchParams);
        YBI_ASSERT(edge.storedPtexFaceId >= 0);
        YBI_ASSERT(edge.sampleVStart >= 0 && edge.sampleVStart < int(sharedPos.size()));
        YBI_ASSERT(edge.sampleVEnd >= 0 && edge.sampleVEnd < int(sharedPos.size()));

        int t = edge.tmaxEdgeFactor;
        YBI_ASSERT(t >= 1);
        if (t > 1)
        {
            YBI_ASSERT(edge.edgeVertexIndexStart >= 0);
            YBI_ASSERT(edge.edgeVertexIndexStart < int(sharedPos.size()));
        }

        for (int k = 0; k <= t; ++k)
        {
            const float a = float(k) / float(t);
            const pxr::GfVec2f uv = edge.storedUv0 * (1.0f - a) + edge.storedUv1 * a;
            float3 p = make_float3(0.0f);
            if (!EvaluateLimitPosition(
                    patchMap, patchTable, limitValues, edge.storedPtexFaceId, uv, &p))
            {
                return false;
            }

            int vertexId = -1;
            if (k == 0)
            {
                vertexId = edge.sampleVStart;
            }
            else if (k == t)
            {
                vertexId = edge.sampleVEnd;
            }
            else
            {
                vertexId = edge.edgeVertexIndexStart + (k - 1);
            }
            YBI_ASSERT(vertexId >= 0 && vertexId < int(sharedPos.size()));
            if (!isSharedInit(vertexId))
            {
                sharedPos[vertexId] = p;
                sharedInit[vertexId] = 1;
            }
        }
    }

    std::vector<float3> positions;
    std::vector<int> indices;
    if (outTriPatchFaceIds)
    {
        outTriPatchFaceIds->clear();
    }
    if (outTriCoarseFaceIds)
    {
        outTriCoarseFaceIds->clear();
    }
    if (outTriPtexFaceIds)
    {
        outTriPtexFaceIds->clear();
    }
    if (outTriQuadrants)
    {
        outTriQuadrants->clear();
    }

    std::vector<int> sharedMeshIndex(maxVertexId, -1);
    for (int i = 0; i < int(sharedPos.size()); ++i)
    {
        if (!isSharedInit(i))
        {
            continue;
        }
        const float3 &p = sharedPos[i];
        sharedMeshIndex[i] = int(positions.size());
        positions.push_back(p);
    }

    int currentPatchFaceId = -1;
    int currentCoarseFaceId = -1;
    int currentPtexFaceId = -1;
    int currentQuadrant = -1;

    auto emitTri = [&](int a, int b, int c) {
        YBI_ERROR(a >= 0 && b >= 0 && c >= 0, "%i %i %i\n", a, b, c);
        if (a == b || b == c || a == c)
        {
            return;
        }
        indices.push_back(a);
        indices.push_back(b);
        indices.push_back(c);
        if (outTriPatchFaceIds)
        {
            outTriPatchFaceIds->push_back(currentPatchFaceId);
        }
        if (outTriCoarseFaceIds)
        {
            outTriCoarseFaceIds->push_back(currentCoarseFaceId);
        }
        if (outTriPtexFaceIds)
        {
            outTriPtexFaceIds->push_back(currentPtexFaceId);
        }
        if (outTriQuadrants)
        {
            outTriQuadrants->push_back(currentQuadrant);
        }
    };

    auto stitchStrip =
        [&](const std::vector<int> &outer, const std::vector<int> &inner, int m, int n) {
            YBI_ASSERT(!outer.empty());
            YBI_ASSERT(!inner.empty());

            size_t io = 0;
            size_t ii = 0;
            if (inner.size() == 1)
            {
                while (io + 1 < outer.size())
                {
                    emitTri(outer[io], outer[io + 1], inner[0]);
                    io++;
                }
                return;
            }

            int q = m - 3 * n;
            while (io + 1 < outer.size() || ii + 1 < inner.size())
            {
                const bool canAdvanceOuter = (io + 1 < outer.size());
                const bool canAdvanceInner = (ii + 1 < inner.size());

                if (!canAdvanceInner)
                {
                    emitTri(outer[io], outer[io + 1], inner[ii]);
                    io++;
                    continue;
                }
                if (!canAdvanceOuter)
                {
                    emitTri(outer[io], inner[ii + 1], inner[ii]);
                    ii++;
                    continue;
                }

                if (q >= 0)
                {
                    emitTri(outer[io], inner[ii + 1], inner[ii]);
                    ii++;
                    q -= 2 * n;
                }
                else
                {
                    emitTri(outer[io], outer[io + 1], inner[ii]);
                    io++;
                    q += 2 * m;
                }
            }
        };

    for (size_t patchFaceId = 0; patchFaceId < leafPatches.size(); ++patchFaceId)
    {
        const SubdivisionPatch &patch = leafPatches[patchFaceId];
        currentPatchFaceId = int(patchFaceId);
        currentCoarseFaceId = patch.coarseFace;
        currentPtexFaceId = patch.ptexFaceId;
        currentQuadrant = patch.quadrant;

        auto buildOuter = [&](const SubdivisionEdge &edge, int a, int b, int t) {
            std::vector<int> outer;
            outer.reserve(t + 1);
            outer.push_back(a);
            RecurseChildEdges(edgeMap, a, b, outer);
            YBI_ERROR((!edge.split && outer.size() == t + 1) || (edge.split && outer.size() > t),
                      "split: %i\n",
                      edge.split);

            std::vector<int> out;
            for (int o : outer)
            {
                YBI_ASSERT(o >= 0 && o < int(sharedPos.size()));
                YBI_ASSERT(isSharedInit(o));
                const int meshId = sharedMeshIndex[o];
                YBI_ASSERT(meshId >= 0);
                out.push_back(meshId);
            }
            return out;
        };

        int e0 = GetEdgeFactor(edgeMap, patch.verts[0], patch.verts[1]);
        int e1 = GetEdgeFactor(edgeMap, patch.verts[1], patch.verts[2]);
        int e2 = GetEdgeFactor(edgeMap, patch.verts[2], patch.verts[3]);
        int e3 = GetEdgeFactor(edgeMap, patch.verts[3], patch.verts[0]);

        YBI_ASSERT(e0 != SUBDIV_EDGE_FACTOR_NON_UNIFORM);
        YBI_ASSERT(e1 != SUBDIV_EDGE_FACTOR_NON_UNIFORM);
        YBI_ASSERT(e2 != SUBDIV_EDGE_FACTOR_NON_UNIFORM);
        YBI_ASSERT(e3 != SUBDIV_EDGE_FACTOR_NON_UNIFORM);
        YBI_ASSERT(e0 >= 1 && e1 >= 1 && e2 >= 1 && e3 >= 1);

        const SubdivisionEdge &edge0 = GetEdge(edgeMap, patch.verts[0], patch.verts[1]);
        const SubdivisionEdge &edge1 = GetEdge(edgeMap, patch.verts[1], patch.verts[2]);
        const SubdivisionEdge &edge2 = GetEdge(edgeMap, patch.verts[2], patch.verts[3]);
        const SubdivisionEdge &edge3 = GetEdge(edgeMap, patch.verts[3], patch.verts[0]);

        std::vector<int> outer0 = buildOuter(edge0, patch.verts[0], patch.verts[1], e0);
        std::vector<int> outer1 = buildOuter(edge1, patch.verts[1], patch.verts[2], e1);
        std::vector<int> outer2 = buildOuter(edge2, patch.verts[2], patch.verts[3], e2);
        std::vector<int> outer3 = buildOuter(edge3, patch.verts[3], patch.verts[0], e3);

        e0 = outer0.size() - 1;
        e1 = outer1.size() - 1;
        e2 = outer2.size() - 1;
        e3 = outer3.size() - 1;

        const int nu = std::max(std::max(e0, e2), 2);
        const int nv = std::max(std::max(e1, e3), 2);
        const int cols = nu - 1;
        const int rows = nv - 1;
        std::vector<int> innerGridIndex(cols * rows, -1);
        auto innerAt = [&](int iu, int iv) -> int & {
            YBI_ASSERT(iu >= 1 && iu <= cols);
            YBI_ASSERT(iv >= 1 && iv <= rows);
            return innerGridIndex[(iv - 1) * cols + (iu - 1)];
        };

        for (int iv = 1; iv < nv; ++iv)
        {
            const float sv = float(iv) / float(nv);
            for (int iu = 1; iu < nu; ++iu)
            {
                const float su = float(iu) / float(nu);
                const pxr::GfVec2f uvBottom = patch.uv[0] * (1.0f - su) + patch.uv[1] * su;
                const pxr::GfVec2f uvTop = patch.uv[3] * (1.0f - su) + patch.uv[2] * su;
                const pxr::GfVec2f uv = uvBottom * (1.0f - sv) + uvTop * sv;

                float3 p = make_float3(0.0f);
                if (!EvaluateLimitPosition(
                        patchMap, patchTable, limitValues, patch.ptexFaceId, uv, &p))
                {
                    return false;
                }
                innerAt(iu, iv) = int(positions.size());
                positions.push_back(p);
            }
        }

        if (cols >= 2 && rows >= 2)
        {
            for (int y = 0; y < rows - 1; ++y)
            {
                for (int x = 0; x < cols - 1; ++x)
                {
                    const int i00 = innerAt(x + 1, y + 1);
                    const int i10 = innerAt(x + 2, y + 1);
                    const int i01 = innerAt(x + 1, y + 2);
                    const int i11 = innerAt(x + 2, y + 2);
                    emitTri(i00, i10, i11);
                    emitTri(i00, i11, i01);
                }
            }
        }

        std::vector<int> inner0;
        std::vector<int> inner1;
        std::vector<int> inner2;
        std::vector<int> inner3;
        inner0.reserve(cols);
        inner1.reserve(rows);
        inner2.reserve(cols);
        inner3.reserve(rows);
        for (int iu = 1; iu < nu; ++iu)
        {
            inner0.push_back(innerAt(iu, 1));
        }
        for (int iv = 1; iv < nv; ++iv)
        {
            inner1.push_back(innerAt(cols, iv));
        }
        for (int iu = 1; iu < nu; ++iu)
        {
            inner2.push_back(innerAt(iu, rows));
        }
        for (int iv = 1; iv < nv; ++iv)
        {
            inner3.push_back(innerAt(1, iv));
        }

        std::reverse(outer2.begin(), outer2.end());
        std::reverse(outer3.begin(), outer3.end());
        stitchStrip(outer0, inner0, std::max(e0, e2), e0);
        stitchStrip(outer1, inner1, std::max(e1, e3), e1);
        stitchStrip(outer2, inner2, std::max(e0, e2), e2);
        stitchStrip(outer3, inner3, std::max(e1, e3), e3);
    }

    Array<float3> outPositions(positions);
    Array<int> outIndices(indices);
    YBI_ASSERT(outPositions.size() == positions.size());
    YBI_ASSERT(outIndices.size() == indices.size());
    outMesh->positions = std::move(outPositions);
    outMesh->indices = std::move(outIndices);

    if (outVertexCount)
    {
        *outVertexCount = int(positions.size());
    }
    if (outTriangleCount)
    {
        *outTriangleCount = int(indices.size() / 3);
    }
    return true;
}
