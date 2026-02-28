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

static void RecurseChildEdgeUVs(const SubdivisionEdgeMap &edgeMap,
                                int v0,
                                int v1,
                                const pxr::GfVec2f &uv0,
                                const pxr::GfVec2f &uv1,
                                std::vector<pxr::GfVec2f> &uvs)
{
    const SubdivisionEdge &edge = GetEdge(edgeMap, v0, v1);
    auto GetEdgeVertexUV = [](const SubdivisionEdge &e,
                              int vertexId,
                              const pxr::GfVec2f &fallback) -> pxr::GfVec2f {
        if (!e.hasStoredPatchParams)
        {
            return fallback;
        }
        if (e.sampleVStart == vertexId)
        {
            return e.storedUv0;
        }
        if (e.sampleVEnd == vertexId)
        {
            return e.storedUv1;
        }
        if (e.v0 == vertexId)
        {
            return e.storedUv0;
        }
        if (e.v1 == vertexId)
        {
            return e.storedUv1;
        }
        return fallback;
    };
    if (edge.split)
    {
        YBI_ASSERT(edge.midpointVertex >= 0);
        const SubdivisionEdge &edgeA = GetEdge(edgeMap, v0, edge.midpointVertex);
        const SubdivisionEdge &edgeB = GetEdge(edgeMap, edge.midpointVertex, v1);
        const pxr::GfVec2f uvFallbackMid = uv0 * 0.5f + uv1 * 0.5f;
        const pxr::GfVec2f uvMid = GetEdgeVertexUV(edgeA,
                                                   edge.midpointVertex,
                                                   GetEdgeVertexUV(edgeB,
                                                                   edge.midpointVertex,
                                                                   uvFallbackMid));
        RecurseChildEdgeUVs(edgeMap, v0, edge.midpointVertex, uv0, uvMid, uvs);
        RecurseChildEdgeUVs(edgeMap, edge.midpointVertex, v1, uvMid, uv1, uvs);
    }
    else
    {
        YBI_ASSERT(edge.tmaxEdgeFactor >= 1);
        for (int k = 1; k < edge.tmaxEdgeFactor; ++k)
        {
            const float a = float(k) / float(edge.tmaxEdgeFactor);
            uvs.push_back(uv0 * (1.0f - a) + uv1 * a);
        }
        uvs.push_back(uv1);
    }
}

static bool BuildLeafPatchStitchedMesh(const std::vector<SubdivisionPatch> &leafPatches,
                                       const SubdivisionEdgeMap &edgeMap,
                                       const Far::PatchMap &patchMap,
                                       const Far::PatchTable &patchTable,
                                       const std::vector<LimitEvalVertex> &limitValues,
                                       const std::vector<LimitEvalFVar2> *limitFVarValues,
                                       int nextGeneratedVertexId,
                                       Mesh *outMesh,
                                       std::vector<int> *outTriPatchFaceIds,
                                       std::vector<int> *outTriCoarseFaceIds,
                                       std::vector<int> *outTriPtexFaceIds,
                                       std::vector<int> *outTriQuadrants,
                                       int *outVertexCount,
                                       int *outTriangleCount,
                                       int *outInnerGridTriangleCount,
                                       int *outStitchingTriangleCount)
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
    if (outInnerGridTriangleCount)
    {
        *outInnerGridTriangleCount = 0;
    }
    if (outStitchingTriangleCount)
    {
        *outStitchingTriangleCount = 0;
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

    const bool hasTexcoords = limitFVarValues != nullptr;
    Array<float3> positions;
    Array<int> indices;
    Array<float2> texcoords;
    Array<int> texcoordIndices;
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
        positions.EmplaceBack(p);
    }

    int currentPatchFaceId = -1;
    int currentCoarseFaceId = -1;
    int currentPtexFaceId = -1;
    int currentQuadrant = -1;
    int innerGridTriangleCount = 0;
    int stitchingTriangleCount = 0;
    bool emitFromInnerGrid = false;

    auto emitTri = [&](int a, int b, int c, int ta, int tb, int tc) {
        YBI_ERROR(a >= 0 && b >= 0 && c >= 0, "%i %i %i\n", a, b, c);
        if (a == b || b == c || a == c)
        {
            return;
        }
        indices.EmplaceBack(a);
        indices.EmplaceBack(b);
        indices.EmplaceBack(c);
        if (hasTexcoords)
        {
            YBI_ASSERT(ta >= 0 && tb >= 0 && tc >= 0);
            texcoordIndices.EmplaceBack(ta);
            texcoordIndices.EmplaceBack(tb);
            texcoordIndices.EmplaceBack(tc);
        }
        if (emitFromInnerGrid)
        {
            innerGridTriangleCount++;
        }
        else
        {
            stitchingTriangleCount++;
        }
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

    auto stitchStrip = [&](const std::vector<int> &outer,
                           const std::vector<int> &inner,
                           const std::vector<int> *outerTc,
                           const std::vector<int> *innerTc,
                           int m,
                           int n) {
            YBI_ASSERT(!outer.empty());
            YBI_ASSERT(!inner.empty());
            if (hasTexcoords)
            {
                YBI_ASSERT(outerTc != nullptr);
                YBI_ASSERT(innerTc != nullptr);
                YBI_ASSERT(outerTc->size() == outer.size());
                YBI_ASSERT(innerTc->size() == inner.size());
            }

            size_t io = 0;
            size_t ii = 0;
            if (inner.size() == 1)
            {
                while (io + 1 < outer.size())
                {
                    emitTri(outer[io],
                            outer[io + 1],
                            inner[0],
                            hasTexcoords ? (*outerTc)[io] : -1,
                            hasTexcoords ? (*outerTc)[io + 1] : -1,
                            hasTexcoords ? (*innerTc)[0] : -1);
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
                    emitTri(outer[io],
                            outer[io + 1],
                            inner[ii],
                            hasTexcoords ? (*outerTc)[io] : -1,
                            hasTexcoords ? (*outerTc)[io + 1] : -1,
                            hasTexcoords ? (*innerTc)[ii] : -1);
                    io++;
                    continue;
                }
                if (!canAdvanceOuter)
                {
                    emitTri(outer[io],
                            inner[ii + 1],
                            inner[ii],
                            hasTexcoords ? (*outerTc)[io] : -1,
                            hasTexcoords ? (*innerTc)[ii + 1] : -1,
                            hasTexcoords ? (*innerTc)[ii] : -1);
                    ii++;
                    continue;
                }

                if (q >= 0)
                {
                    emitTri(outer[io],
                            inner[ii + 1],
                            inner[ii],
                            hasTexcoords ? (*outerTc)[io] : -1,
                            hasTexcoords ? (*innerTc)[ii + 1] : -1,
                            hasTexcoords ? (*innerTc)[ii] : -1);
                    ii++;
                    q -= 2 * n;
                }
                else
                {
                    emitTri(outer[io],
                            outer[io + 1],
                            inner[ii],
                            hasTexcoords ? (*outerTc)[io] : -1,
                            hasTexcoords ? (*outerTc)[io + 1] : -1,
                            hasTexcoords ? (*innerTc)[ii] : -1);
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

        auto buildOuter = [&](const SubdivisionEdge &edge,
                              int a,
                              int b,
                              int t,
                              const pxr::GfVec2f &uvA,
                              const pxr::GfVec2f &uvB,
                              std::vector<int> *out,
                              std::vector<int> *outTc) -> bool {
            YBI_ASSERT(out);
            YBI_ASSERT(outTc);

            std::vector<int> outerVertexIds;
            outerVertexIds.reserve(t + 1);
            outerVertexIds.push_back(a);
            RecurseChildEdges(edgeMap, a, b, outerVertexIds);
            YBI_ERROR((!edge.split && outerVertexIds.size() == t + 1) ||
                          (edge.split && outerVertexIds.size() > size_t(t)),
                      "split: %i\n",
                      edge.split);

            std::vector<pxr::GfVec2f> outerPatchUVs;
            if (hasTexcoords)
            {
                outerPatchUVs.reserve(outerVertexIds.size());
                outerPatchUVs.push_back(uvA);
                RecurseChildEdgeUVs(edgeMap, a, b, uvA, uvB, outerPatchUVs);
                YBI_ASSERT(outerPatchUVs.size() == outerVertexIds.size());
            }

            out->clear();
            out->reserve(outerVertexIds.size());
            outTc->clear();
            outTc->reserve(outerVertexIds.size());
            for (size_t i = 0; i < outerVertexIds.size(); ++i)
            {
                const int o = outerVertexIds[i];
                YBI_ASSERT(o >= 0 && o < int(sharedPos.size()));
                YBI_ASSERT(isSharedInit(o));
                const int meshId = sharedMeshIndex[o];
                YBI_ASSERT(meshId >= 0);
                out->push_back(meshId);

                if (hasTexcoords)
                {
                    float2 tc = make_float2(0.0f);
                    if (!EvaluateLimitFVar2(patchMap,
                                            patchTable,
                                            *limitFVarValues,
                                            patch.ptexFaceId,
                                            outerPatchUVs[i],
                                            &tc))
                    {
                        return false;
                    }
                    outTc->push_back(int(texcoords.size()));
                    texcoords.EmplaceBack(tc);
                }
            }
            return true;
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

        std::vector<int> outer0, outer1, outer2, outer3;
        std::vector<int> outer0Tc, outer1Tc, outer2Tc, outer3Tc;
        if (!buildOuter(edge0,
                        patch.verts[0],
                        patch.verts[1],
                        e0,
                        patch.uv[0],
                        patch.uv[1],
                        &outer0,
                        &outer0Tc) ||
            !buildOuter(edge1,
                        patch.verts[1],
                        patch.verts[2],
                        e1,
                        patch.uv[1],
                        patch.uv[2],
                        &outer1,
                        &outer1Tc) ||
            !buildOuter(edge2,
                        patch.verts[2],
                        patch.verts[3],
                        e2,
                        patch.uv[2],
                        patch.uv[3],
                        &outer2,
                        &outer2Tc) ||
            !buildOuter(edge3,
                        patch.verts[3],
                        patch.verts[0],
                        e3,
                        patch.uv[3],
                        patch.uv[0],
                        &outer3,
                        &outer3Tc))
        {
            return false;
        }

        e0 = outer0.size() - 1;
        e1 = outer1.size() - 1;
        e2 = outer2.size() - 1;
        e3 = outer3.size() - 1;

        const int nu = std::max(std::max(e0, e2), 2);
        const int nv = std::max(std::max(e1, e3), 2);
        const int cols = nu - 1;
        const int rows = nv - 1;
        std::vector<int> innerGridIndex(cols * rows, -1);
        std::vector<int> innerGridTcIndex(cols * rows, -1);
        auto innerAt = [&](int iu, int iv) -> int & {
            YBI_ASSERT(iu >= 1 && iu <= cols);
            YBI_ASSERT(iv >= 1 && iv <= rows);
            return innerGridIndex[(iv - 1) * cols + (iu - 1)];
        };
        auto innerTcAt = [&](int iu, int iv) -> int & {
            YBI_ASSERT(iu >= 1 && iu <= cols);
            YBI_ASSERT(iv >= 1 && iv <= rows);
            return innerGridTcIndex[(iv - 1) * cols + (iu - 1)];
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
                positions.EmplaceBack(p);
                if (hasTexcoords)
                {
                    float2 tc = make_float2(0.0f);
                    if (!EvaluateLimitFVar2(
                            patchMap, patchTable, *limitFVarValues, patch.ptexFaceId, uv, &tc))
                    {
                        return false;
                    }
                    innerTcAt(iu, iv) = int(texcoords.size());
                    texcoords.EmplaceBack(tc);
                }
            }
        }

        if (cols >= 2 && rows >= 2)
        {
            emitFromInnerGrid = true;
            for (int y = 0; y < rows - 1; ++y)
            {
                for (int x = 0; x < cols - 1; ++x)
                {
                    const int i00 = innerAt(x + 1, y + 1);
                    const int i10 = innerAt(x + 2, y + 1);
                    const int i01 = innerAt(x + 1, y + 2);
                    const int i11 = innerAt(x + 2, y + 2);
                    const int tc00 = hasTexcoords ? innerTcAt(x + 1, y + 1) : -1;
                    const int tc10 = hasTexcoords ? innerTcAt(x + 2, y + 1) : -1;
                    const int tc01 = hasTexcoords ? innerTcAt(x + 1, y + 2) : -1;
                    const int tc11 = hasTexcoords ? innerTcAt(x + 2, y + 2) : -1;
                    emitTri(i00, i10, i11, tc00, tc10, tc11);
                    emitTri(i00, i11, i01, tc00, tc11, tc01);
                }
            }
            emitFromInnerGrid = false;
        }

        std::vector<int> inner0;
        std::vector<int> inner1;
        std::vector<int> inner2;
        std::vector<int> inner3;
        std::vector<int> inner0Tc;
        std::vector<int> inner1Tc;
        std::vector<int> inner2Tc;
        std::vector<int> inner3Tc;
        inner0.reserve(cols);
        inner1.reserve(rows);
        inner2.reserve(cols);
        inner3.reserve(rows);
        inner0Tc.reserve(cols);
        inner1Tc.reserve(rows);
        inner2Tc.reserve(cols);
        inner3Tc.reserve(rows);
        for (int iu = 1; iu < nu; ++iu)
        {
            inner0.push_back(innerAt(iu, 1));
            if (hasTexcoords)
            {
                inner0Tc.push_back(innerTcAt(iu, 1));
            }
        }
        for (int iv = 1; iv < nv; ++iv)
        {
            inner1.push_back(innerAt(cols, iv));
            if (hasTexcoords)
            {
                inner1Tc.push_back(innerTcAt(cols, iv));
            }
        }
        for (int iu = 1; iu < nu; ++iu)
        {
            inner2.push_back(innerAt(iu, rows));
            if (hasTexcoords)
            {
                inner2Tc.push_back(innerTcAt(iu, rows));
            }
        }
        for (int iv = 1; iv < nv; ++iv)
        {
            inner3.push_back(innerAt(1, iv));
            if (hasTexcoords)
            {
                inner3Tc.push_back(innerTcAt(1, iv));
            }
        }

        std::reverse(outer2.begin(), outer2.end());
        std::reverse(outer3.begin(), outer3.end());
        std::reverse(outer2Tc.begin(), outer2Tc.end());
        std::reverse(outer3Tc.begin(), outer3Tc.end());
        stitchStrip(outer0,
                    inner0,
                    hasTexcoords ? &outer0Tc : nullptr,
                    hasTexcoords ? &inner0Tc : nullptr,
                    std::max(e0, e2),
                    e0);
        stitchStrip(outer1,
                    inner1,
                    hasTexcoords ? &outer1Tc : nullptr,
                    hasTexcoords ? &inner1Tc : nullptr,
                    std::max(e1, e3),
                    e1);
        stitchStrip(outer2,
                    inner2,
                    hasTexcoords ? &outer2Tc : nullptr,
                    hasTexcoords ? &inner2Tc : nullptr,
                    std::max(e0, e2),
                    e2);
        stitchStrip(outer3,
                    inner3,
                    hasTexcoords ? &outer3Tc : nullptr,
                    hasTexcoords ? &inner3Tc : nullptr,
                    std::max(e1, e3),
                    e3);
    }

    outMesh->positions = std::move(positions);
    outMesh->indices = std::move(indices);
    if (hasTexcoords)
    {
        YBI_ASSERT(texcoordIndices.size() == outMesh->indices.size());
        outMesh->texcoords = std::move(texcoords);
        outMesh->texcoordIndices = std::move(texcoordIndices);
    }
    else
    {
        outMesh->texcoords.Resize(0);
        outMesh->texcoordIndices.Resize(0);
    }

    if (outVertexCount)
    {
        *outVertexCount = int(outMesh->positions.size());
    }
    if (outTriangleCount)
    {
        *outTriangleCount = int(outMesh->indices.size() / 3);
    }
    if (outInnerGridTriangleCount)
    {
        *outInnerGridTriangleCount = innerGridTriangleCount;
    }
    if (outStitchingTriangleCount)
    {
        *outStitchingTriangleCount = stitchingTriangleCount;
    }
    return true;
}
