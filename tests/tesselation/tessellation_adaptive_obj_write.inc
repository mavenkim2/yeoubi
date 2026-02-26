static bool WriteLeafPatchCornerQuadsObj(const std::vector<SubdivisionPatch> &leafPatches,
                                         const SubdivisionEdgeMap &edgeMap,
                                         const Far::PatchMap &patchMap,
                                         const Far::PatchTable &patchTable,
                                         const std::vector<LimitEvalVertex> &limitValues,
                                         const std::string &outObjPath,
                                         int *outVertexCount,
                                         int *outQuadCount)
{
    if (outVertexCount)
    {
        *outVertexCount = 0;
    }
    if (outQuadCount)
    {
        *outQuadCount = 0;
    }

    int maxVertexId = -1;
    for (const SubdivisionPatch &patch : leafPatches)
    {
        for (int i = 0; i < 4; ++i)
        {
            maxVertexId = std::max(maxVertexId, patch.verts[i]);
        }
    }
    if (maxVertexId < 0)
    {
        FILE *fEmpty = std::fopen(outObjPath.c_str(), "w");
        if (!fEmpty)
        {
            return false;
        }
        std::fclose(fEmpty);
        return true;
    }

    const float inf = std::numeric_limits<float>::infinity();
    std::vector<pxr::GfVec3f> cornerPos(maxVertexId + 1, pxr::GfVec3f(inf, inf, inf));
    std::vector<uint8_t> cornerInit(maxVertexId + 1, uint8_t(0));

    for (const SubdivisionPatch &patch : leafPatches)
    {
        for (int i = 0; i < 4; ++i)
        {
            const int vid = patch.verts[i];
            YBI_ASSERT(vid >= 0 && vid < int(cornerPos.size()));
            if (cornerInit[vid] != 0)
            {
                continue;
            }

            pxr::GfVec3f p(0.0f);
            if (!EvaluateLimitPosition(
                    patchMap, patchTable, limitValues, patch.ptexFaceId, patch.uv[i], &p))
            {
                return false;
            }
            cornerPos[vid] = p;
            cornerInit[vid] = 1;
        }
    }

    FILE *f = std::fopen(outObjPath.c_str(), "w");
    if (!f)
    {
        return false;
    }

    int vertexCount = 0;
    int quadCount = 0;
    std::vector<int> objIndex(maxVertexId + 1, 0);
    for (int i = 0; i < int(cornerPos.size()); ++i)
    {
        if (cornerInit[i] == 0)
        {
            continue;
        }
        const pxr::GfVec3f &p = cornerPos[i];
        std::fprintf(f, "v %.9g %.9g %.9g\n", double(p[0]), double(p[1]), double(p[2]));
        vertexCount++;
        objIndex[i] = vertexCount;
    }

    int faceId = 0;
    for (const SubdivisionPatch &patch : leafPatches)
    {
        const int i0 = objIndex[patch.verts[0]];
        const int i1 = objIndex[patch.verts[1]];
        const int i2 = objIndex[patch.verts[2]];
        const int i3 = objIndex[patch.verts[3]];
        YBI_ASSERT(i0 > 0 && i1 > 0 && i2 > 0 && i3 > 0);
        const int e0 = GetEdgeFactor(edgeMap, patch.verts[0], patch.verts[1]);
        const int e1 = GetEdgeFactor(edgeMap, patch.verts[1], patch.verts[2]);
        const int e2 = GetEdgeFactor(edgeMap, patch.verts[2], patch.verts[3]);
        const int e3 = GetEdgeFactor(edgeMap, patch.verts[3], patch.verts[0]);
        std::fprintf(f,
                     "# face_id %d corner_vertex_ids %d %d %d %d obj_vertex_ids %d %d %d %d "
                     "edge_factors %d %d %d %d\n",
                     faceId,
                     patch.verts[0],
                     patch.verts[1],
                     patch.verts[2],
                     patch.verts[3],
                     i0,
                     i1,
                     i2,
                     i3,
                     e0,
                     e1,
                     e2,
                     e3);
        std::fprintf(f, "f %d %d %d %d\n", i0, i1, i2, i3);
        quadCount++;
        faceId++;
    }

    std::fclose(f);
    if (outVertexCount)
    {
        *outVertexCount = vertexCount;
    }
    if (outQuadCount)
    {
        *outQuadCount = quadCount;
    }
    return true;
}

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

static bool WriteLeafPatchInnerGridObj(const std::vector<SubdivisionPatch> &leafPatches,
                                       const SubdivisionEdgeMap &edgeMap,
                                       const Far::PatchMap &patchMap,
                                       const Far::PatchTable &patchTable,
                                       const std::vector<LimitEvalVertex> &limitValues,
                                       const std::string &outObjPath,
                                       int *outVertexCount,
                                       int *outTriangleCount)
{
    if (outVertexCount)
    {
        *outVertexCount = 0;
    }
    if (outTriangleCount)
    {
        *outTriangleCount = 0;
    }

    FILE *f = std::fopen(outObjPath.c_str(), "w");
    if (!f)
    {
        return false;
    }

    int maxVertexId = -1;
    for (const SubdivisionPatch &patch : leafPatches)
    {
        for (int i = 0; i < 4; ++i)
        {
            maxVertexId = std::max(maxVertexId, patch.verts[i]);
        }
    }
    for (const auto &it : edgeMap)
    {
        const SubdivisionEdge &edge = it.second;
        maxVertexId = std::max(maxVertexId, edge.v0);
        maxVertexId = std::max(maxVertexId, edge.v1);
        maxVertexId = std::max(maxVertexId, edge.sampleVStart);
        maxVertexId = std::max(maxVertexId, edge.sampleVEnd);
        maxVertexId = std::max(maxVertexId, edge.midpointVertex);
        if (edge.edgeVertexIndexStart >= 0 && edge.tmaxEdgeFactor >= 2)
        {
            maxVertexId =
                std::max(maxVertexId, edge.edgeVertexIndexStart + (edge.tmaxEdgeFactor - 2));
        }
    }
    if (maxVertexId < 0)
    {
        std::fclose(f);
        return true;
    }

    const float inf = std::numeric_limits<float>::infinity();
    std::vector<pxr::GfVec3f> sharedPos(maxVertexId + 1, pxr::GfVec3f(inf, inf, inf));
    std::vector<uint8_t> sharedInit(maxVertexId + 1, uint8_t(0));
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
            pxr::GfVec3f p(0.0f);
            if (!EvaluateLimitPosition(
                    patchMap, patchTable, limitValues, edge.storedPtexFaceId, uv, &p))
            {
                std::fclose(f);
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

    int vertexCount = 0;
    int triangleCount = 0;
    std::vector<int> sharedObjIndex(maxVertexId + 1, 0);
    for (int i = 0; i < int(sharedPos.size()); ++i)
    {
        if (!isSharedInit(i))
        {
            continue;
        }
        const pxr::GfVec3f &p = sharedPos[i];
        std::fprintf(f, "v %.9g %.9g %.9g\n", double(p[0]), double(p[1]), double(p[2]));
        vertexCount++;
        sharedObjIndex[i] = vertexCount;
    }

    auto emitTri = [&](int a, int b, int c) {
        YBI_ERROR(a > 0 && b > 0 && c > 0, "%i %i %i\n", a, b, c);
        if (a == b || b == c || a == c)
        {
            return;
        }
        std::fprintf(f, "f %d %d %d\n", a, b, c);
        triangleCount++;
    };

#if 1
    auto getOrientedEdgeVertexId =
        [&](const SubdivisionEdge &edge, int a, int b, int t, int k) -> int {
        YBI_ASSERT(t >= 1);
        YBI_ASSERT(k >= 0 && k <= t);
        const bool forward = (edge.v0 == a && edge.v1 == b);
        const bool reverse = (edge.v0 == b && edge.v1 == a);
        YBI_ASSERT(forward || reverse);
        if (k == 0)
        {
            return a;
        }
        if (k == t)
        {
            return b;
        }
        YBI_ASSERT(t > 1);
        YBI_ASSERT(edge.edgeVertexIndexStart >= 0);
        const int localOffset = forward ? (k - 1) : (t - 1 - k);
        return edge.edgeVertexIndexStart + localOffset;
    };
#endif

    auto stitchStrip = [&](const std::vector<int> &outer, const std::vector<int> &inner) {
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

            const float uOuter = float(io + 1) / float(outer.size() - 1);
            const float uInner = float(ii + 1) / float(inner.size() - 1);
            if (uOuter <= uInner)
            {
                emitTri(outer[io], outer[io + 1], inner[ii]);
                io++;
            }
            else
            {
                emitTri(outer[io], inner[ii + 1], inner[ii]);
                ii++;
            }
        }
    };

    for (const SubdivisionPatch &patch : leafPatches)
    {
        auto buildOuter = [&](const SubdivisionEdge &edge, int a, int b, int t) {
            std::vector<int> outer;
            outer.reserve(t + 1);
#if 1
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
                const int objId = sharedObjIndex[o];
                YBI_ASSERT(objId > 0);
                out.push_back(objId);
            }
            return out;
#else
            std::vector<int> test;
            for (int k = 0; k <= t; ++k)
            {
                const int vertexId = getOrientedEdgeVertexId(edge, a, b, t, k);
                YBI_ASSERT(vertexId >= 0 && vertexId < int(sharedPos.size()));
                YBI_ASSERT(isSharedInit(vertexId));
                const int objId = sharedObjIndex[vertexId];
                YBI_ASSERT(objId > 0);
                test.push_back(objId);
            }
            return test;
#endif
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
        std::vector<int> innerGridObj(cols * rows, 0);
        auto innerAt = [&](int iu, int iv) -> int & {
            YBI_ASSERT(iu >= 1 && iu <= cols);
            YBI_ASSERT(iv >= 1 && iv <= rows);
            return innerGridObj[(iv - 1) * cols + (iu - 1)];
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

                pxr::GfVec3f p(0.0f);
                if (!EvaluateLimitPosition(
                        patchMap, patchTable, limitValues, patch.ptexFaceId, uv, &p))
                {
                    std::fclose(f);
                    return false;
                }
                std::fprintf(f, "v %.9g %.9g %.9g\n", double(p[0]), double(p[1]), double(p[2]));
                vertexCount++;
                innerAt(iu, iv) = vertexCount;
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
        for (int iu = nu - 1; iu >= 1; --iu)
        {
            inner2.push_back(innerAt(iu, rows));
        }
        for (int iv = nv - 1; iv >= 1; --iv)
        {
            inner3.push_back(innerAt(1, iv));
        }

        stitchStrip(outer0, inner0);
        stitchStrip(outer1, inner1);
        stitchStrip(outer2, inner2);
        stitchStrip(outer3, inner3);
    }

    std::fclose(f);
    if (outVertexCount)
    {
        *outVertexCount = vertexCount;
    }
    if (outTriangleCount)
    {
        *outTriangleCount = triangleCount;
    }
    return true;
}
