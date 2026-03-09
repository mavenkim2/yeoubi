#pragma once

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
    std::vector<Vec3> cornerPos(maxVertexId + 1, Vec3(inf, inf, inf));
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

            Vec3 p = Vec3(0.0f);
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
        const Vec3 &p = cornerPos[i];
        std::fprintf(f, "v %.9g %.9g %.9g\n", double(p.x), double(p.y), double(p.z));
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
