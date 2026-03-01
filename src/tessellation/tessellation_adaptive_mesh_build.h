#pragma once

template <typename T> static int AppendAttributeValue(Array<uint8_t> *dst, const T &value)
{
    YBI_ASSERT(dst);
    const int index = int(dst->size() / sizeof(T));
    const size_t oldBytes = dst->size();
    dst->Resize(oldBytes + sizeof(T));
    memcpy(dst->data() + oldBytes, &value, sizeof(T));
    return index;
}

static bool EvaluateAndAppendAttributeSample(const SubdivisionLimitEvalAttribute &evalAttr,
                                             const Far::PatchMap &patchMap,
                                             const Far::PatchTable &patchTable,
                                             int ptexFaceId,
                                             const pxr::GfVec2f &uv,
                                             Array<uint8_t> *dst,
                                             int *outIndex)
{
    YBI_ASSERT(dst);
    YBI_ASSERT(outIndex);

    if (evalAttr.type == AttributeType::Float)
    {
        float value = 0.0f;
        if (!EvaluateLimitValue<LimitEvalFloat, float>(patchMap,
                                                       patchTable,
                                                       evalAttr.valuesFloat,
                                                       ptexFaceId,
                                                       uv,
                                                       evalAttr.interpolation,
                                                       &value,
                                                       evalAttr.fvarChannel))
        {
            return false;
        }
        *outIndex = AppendAttributeValue<float>(dst, value);
        return true;
    }
    if (evalAttr.type == AttributeType::Float2)
    {
        float2 value = make_float2(0.0f);
        if (!EvaluateLimitValue<LimitEvalFloat2, float2>(patchMap,
                                                         patchTable,
                                                         evalAttr.valuesFloat2,
                                                         ptexFaceId,
                                                         uv,
                                                         evalAttr.interpolation,
                                                         &value,
                                                         evalAttr.fvarChannel))
        {
            return false;
        }
        *outIndex = AppendAttributeValue<float2>(dst, value);
        return true;
    }
    if (evalAttr.type == AttributeType::Float3)
    {
        float3 value = make_float3(0.0f);
        if (!EvaluateLimitValue<LimitEvalFloat3, float3>(patchMap,
                                                         patchTable,
                                                         evalAttr.valuesFloat3,
                                                         ptexFaceId,
                                                         uv,
                                                         evalAttr.interpolation,
                                                         &value,
                                                         evalAttr.fvarChannel))
        {
            return false;
        }
        *outIndex = AppendAttributeValue<float3>(dst, value);
        return true;
    }
    if (evalAttr.type == AttributeType::Float4)
    {
        float4 value = make_float4(0.0f, 0.0f, 0.0f, 0.0f);
        if (!EvaluateLimitValue<LimitEvalFloat4, float4>(patchMap,
                                                         patchTable,
                                                         evalAttr.valuesFloat4,
                                                         ptexFaceId,
                                                         uv,
                                                         evalAttr.interpolation,
                                                         &value,
                                                         evalAttr.fvarChannel))
        {
            return false;
        }
        *outIndex = AppendAttributeValue<float4>(dst, value);
        return true;
    }
    return false;
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

struct TessellationAttributeBuilder
{
    std::string name;
    AttributeType type = AttributeType::Unknown;
    PrimvarInterpolation interpolation = PrimvarInterpolation::Unknown;
    bool isFaceVarying = false;
    Array<uint8_t> values;
    Array<int> indices;
};

static bool BuildLeafPatchStitchedMesh(const std::vector<SubdivisionPatch> &leafPatches,
                                       const SubdivisionEdgeMap &edgeMap,
                                       const Far::PatchMap &patchMap,
                                       const Far::PatchTable &patchTable,
                                       const std::vector<LimitEvalVertex> &limitValues,
                                       const std::vector<SubdivisionLimitEvalAttribute> &limitEvalAttributes,
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
    std::vector<pxr::GfVec2f> sharedUV(maxVertexId, pxr::GfVec2f(0.0f, 0.0f));
    std::vector<int> sharedPtexFace(maxVertexId, -1);
    std::vector<uint8_t> sharedInit(maxVertexId, uint8_t(0));
    std::vector<uint8_t> sharedUVInit(maxVertexId, uint8_t(0));
    auto isSharedInit = [&](int id) -> bool {
        return (id >= 0) && (id < int(sharedInit.size())) && (sharedInit[id] != 0);
    };

    std::vector<TessellationAttributeBuilder> attrBuilders(limitEvalAttributes.size());
    std::vector<int> faceVaryingBuilderIndices;
    std::vector<int> builderToFVarSlot(limitEvalAttributes.size(), -1);
    for (size_t i = 0; i < limitEvalAttributes.size(); ++i)
    {
        attrBuilders[i].name = limitEvalAttributes[i].name;
        attrBuilders[i].type = limitEvalAttributes[i].type;
        attrBuilders[i].interpolation = limitEvalAttributes[i].interpolation;
        attrBuilders[i].isFaceVarying =
            (limitEvalAttributes[i].interpolation == PrimvarInterpolation::FaceVarying);
        if (attrBuilders[i].isFaceVarying)
        {
            builderToFVarSlot[i] = int(faceVaryingBuilderIndices.size());
            faceVaryingBuilderIndices.push_back(int(i));
        }
    }
    const int fvarAttrCount = int(faceVaryingBuilderIndices.size());

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
                sharedUV[vertexId] = uv;
                sharedPtexFace[vertexId] = edge.storedPtexFaceId;
                sharedUVInit[vertexId] = 1;
            }
        }
    }

    Array<float3> positions;
    Array<int> indices;
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
        const int meshIndex = int(positions.size());
        sharedMeshIndex[i] = meshIndex;
        positions.EmplaceBack(sharedPos[i]);
        YBI_ASSERT(sharedUVInit[i] != 0 && sharedPtexFace[i] >= 0);
        for (size_t ai = 0; ai < limitEvalAttributes.size(); ++ai)
        {
            if (attrBuilders[ai].isFaceVarying)
            {
                continue;
            }
            int attrValueIndex = -1;
            if (!EvaluateAndAppendAttributeSample(limitEvalAttributes[ai],
                                                  patchMap,
                                                  patchTable,
                                                  sharedPtexFace[i],
                                                  sharedUV[i],
                                                  &attrBuilders[ai].values,
                                                  &attrValueIndex))
            {
                return false;
            }
            YBI_ASSERT(attrValueIndex == meshIndex);
        }
    }

    int currentPatchFaceId = -1;
    int currentCoarseFaceId = -1;
    int currentPtexFaceId = -1;
    int currentQuadrant = -1;
    int innerGridTriangleCount = 0;
    int stitchingTriangleCount = 0;
    auto emitTriToBuffer = [&](const int refA,
                               const int refB,
                               const int refC,
                               const int outA,
                               const int outB,
                               const int outC,
                               bool trackMetadata,
                               bool fromInnerGrid,
                               Array<int> *outIndexBuffer) -> bool {
        YBI_ASSERT(outIndexBuffer);
        YBI_ERROR(refA >= 0 && refB >= 0 && refC >= 0, "%i %i %i\n", refA, refB, refC);
        if (refA == refB || refB == refC || refA == refC)
        {
            return false;
        }
        outIndexBuffer->EmplaceBack(outA);
        outIndexBuffer->EmplaceBack(outB);
        outIndexBuffer->EmplaceBack(outC);
        if (trackMetadata)
        {
            if (fromInnerGrid)
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
        }
        return true;
    };

    auto stitchStripToBuffer = [&](const std::vector<int> &outerRef,
                                   const std::vector<int> &innerRef,
                                   const std::vector<int> &outerWrite,
                                   const std::vector<int> &innerWrite,
                                   Array<int> *outIndexBuffer,
                                   int m,
                                   int n,
                                   bool trackMetadata) -> int {
            YBI_ASSERT(!outerRef.empty());
            YBI_ASSERT(!innerRef.empty());
            YBI_ASSERT(outerRef.size() == outerWrite.size());
            YBI_ASSERT(innerRef.size() == innerWrite.size());
            YBI_ASSERT(outIndexBuffer);

            const size_t beforeTriCount = outIndexBuffer->size() / 3;
            size_t io = 0;
            size_t ii = 0;
            if (innerRef.size() == 1)
            {
                while (io + 1 < outerRef.size())
                {
                    emitTriToBuffer(outerRef[io],
                                    outerRef[io + 1],
                                    innerRef[0],
                                    outerWrite[io],
                                    outerWrite[io + 1],
                                    innerWrite[0],
                                    trackMetadata,
                                    false,
                                    outIndexBuffer);
                    io++;
                }
                return int(outIndexBuffer->size() / 3 - beforeTriCount);
            }

            int q = m - 3 * n;
            while (io + 1 < outerRef.size() || ii + 1 < innerRef.size())
            {
                const bool canAdvanceOuter = (io + 1 < outerRef.size());
                const bool canAdvanceInner = (ii + 1 < innerRef.size());

                if (!canAdvanceInner)
                {
                    emitTriToBuffer(outerRef[io],
                                    outerRef[io + 1],
                                    innerRef[ii],
                                    outerWrite[io],
                                    outerWrite[io + 1],
                                    innerWrite[ii],
                                    trackMetadata,
                                    false,
                                    outIndexBuffer);
                    io++;
                    continue;
                }
                if (!canAdvanceOuter)
                {
                    emitTriToBuffer(outerRef[io],
                                    innerRef[ii + 1],
                                    innerRef[ii],
                                    outerWrite[io],
                                    innerWrite[ii + 1],
                                    innerWrite[ii],
                                    trackMetadata,
                                    false,
                                    outIndexBuffer);
                    ii++;
                    continue;
                }

                if (q >= 0)
                {
                    emitTriToBuffer(outerRef[io],
                                    innerRef[ii + 1],
                                    innerRef[ii],
                                    outerWrite[io],
                                    innerWrite[ii + 1],
                                    innerWrite[ii],
                                    trackMetadata,
                                    false,
                                    outIndexBuffer);
                    ii++;
                    q -= 2 * n;
                }
                else
                {
                    emitTriToBuffer(outerRef[io],
                                    outerRef[io + 1],
                                    innerRef[ii],
                                    outerWrite[io],
                                    outerWrite[io + 1],
                                    innerWrite[ii],
                                    trackMetadata,
                                    false,
                                    outIndexBuffer);
                    io++;
                    q += 2 * m;
                }
            }
            return int(outIndexBuffer->size() / 3 - beforeTriCount);
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
                              std::vector<std::vector<int>> *outFVar) -> bool {
            YBI_ASSERT(out);
            YBI_ASSERT(outFVar);
            outFVar->assign(size_t(fvarAttrCount), std::vector<int>{});

            std::vector<int> outerVertexIds;
            outerVertexIds.reserve(t + 1);
            outerVertexIds.push_back(a);
            RecurseChildEdges(edgeMap, a, b, outerVertexIds);
            YBI_ERROR((!edge.split && outerVertexIds.size() == size_t(t + 1)) ||
                          (edge.split && outerVertexIds.size() > size_t(t)),
                      "split: %i\n",
                      edge.split);

            const int edgeRate = int(outerVertexIds.size()) - 1;
            YBI_ASSERT(edgeRate >= 1);
            out->clear();
            out->reserve(outerVertexIds.size());
            for (int k = 0; k < fvarAttrCount; ++k)
            {
                (*outFVar)[k].reserve(outerVertexIds.size());
            }

            for (size_t i = 0; i < outerVertexIds.size(); ++i)
            {
                const int o = outerVertexIds[i];
                YBI_ASSERT(o >= 0 && o < int(sharedPos.size()));
                YBI_ASSERT(isSharedInit(o));
                const int meshId = sharedMeshIndex[o];
                YBI_ASSERT(meshId >= 0);
                out->push_back(meshId);

                const float alpha = float(i) / float(edgeRate);
                const pxr::GfVec2f uv = uvA * (1.0f - alpha) + uvB * alpha;
                for (int fvSlot = 0; fvSlot < fvarAttrCount; ++fvSlot)
                {
                    const int builderIndex = faceVaryingBuilderIndices[fvSlot];
                    int valueIndex = -1;
                    if (!EvaluateAndAppendAttributeSample(limitEvalAttributes[builderIndex],
                                                          patchMap,
                                                          patchTable,
                                                          patch.ptexFaceId,
                                                          uv,
                                                          &attrBuilders[builderIndex].values,
                                                          &valueIndex))
                    {
                        return false;
                    }
                    (*outFVar)[fvSlot].push_back(valueIndex);
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
        std::vector<std::vector<int>> outer0FVar, outer1FVar, outer2FVar, outer3FVar;
        if (!buildOuter(edge0,
                        patch.verts[0],
                        patch.verts[1],
                        e0,
                        patch.uv[0],
                        patch.uv[1],
                        &outer0,
                        &outer0FVar) ||
            !buildOuter(edge1,
                        patch.verts[1],
                        patch.verts[2],
                        e1,
                        patch.uv[1],
                        patch.uv[2],
                        &outer1,
                        &outer1FVar) ||
            !buildOuter(edge2,
                        patch.verts[2],
                        patch.verts[3],
                        e2,
                        patch.uv[2],
                        patch.uv[3],
                        &outer2,
                        &outer2FVar) ||
            !buildOuter(edge3,
                        patch.verts[3],
                        patch.verts[0],
                        e3,
                        patch.uv[3],
                        patch.uv[0],
                        &outer3,
                        &outer3FVar))
        {
            return false;
        }

        e0 = int(outer0.size()) - 1;
        e1 = int(outer1.size()) - 1;
        e2 = int(outer2.size()) - 1;
        e3 = int(outer3.size()) - 1;

        const int nu = std::max(std::max(e0, e2), 2);
        const int nv = std::max(std::max(e1, e3), 2);
        const int cols = nu - 1;
        const int rows = nv - 1;
        std::vector<int> innerGridIndex(cols * rows, -1);
        std::vector<std::vector<int>> innerGridFVarIndex(
            size_t(fvarAttrCount), std::vector<int>(cols * rows, -1));
        auto innerAt = [&](int iu, int iv) -> int & {
            YBI_ASSERT(iu >= 1 && iu <= cols);
            YBI_ASSERT(iv >= 1 && iv <= rows);
            return innerGridIndex[(iv - 1) * cols + (iu - 1)];
        };
        auto innerFVarAt = [&](int fvarSlot, int iu, int iv) -> int & {
            YBI_ASSERT(iu >= 1 && iu <= cols);
            YBI_ASSERT(iv >= 1 && iv <= rows);
            return innerGridFVarIndex[size_t(fvarSlot)][(iv - 1) * cols + (iu - 1)];
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
                const int meshIndex = int(positions.size());
                innerAt(iu, iv) = meshIndex;
                positions.EmplaceBack(p);
                for (size_t ai = 0; ai < limitEvalAttributes.size(); ++ai)
                {
                    if (attrBuilders[ai].isFaceVarying)
                    {
                        continue;
                    }
                    int attrValueIndex = -1;
                    if (!EvaluateAndAppendAttributeSample(limitEvalAttributes[ai],
                                                          patchMap,
                                                          patchTable,
                                                          patch.ptexFaceId,
                                                          uv,
                                                          &attrBuilders[ai].values,
                                                          &attrValueIndex))
                    {
                        return false;
                    }
                    YBI_ASSERT(attrValueIndex == meshIndex);
                }
                for (int fvarSlot = 0; fvarSlot < fvarAttrCount; ++fvarSlot)
                {
                    const int builderIndex = faceVaryingBuilderIndices[fvarSlot];
                    int valueIndex = -1;
                    if (!EvaluateAndAppendAttributeSample(limitEvalAttributes[builderIndex],
                                                          patchMap,
                                                          patchTable,
                                                          patch.ptexFaceId,
                                                          uv,
                                                          &attrBuilders[builderIndex].values,
                                                          &valueIndex))
                    {
                        return false;
                    }
                    innerFVarAt(fvarSlot, iu, iv) = valueIndex;
                }
            }
        }

        if (cols >= 2 && rows >= 2)
        {
            const size_t innerGridTriStart = indices.size() / 3;
            std::vector<size_t> innerGridFVarTriStart(size_t(fvarAttrCount), 0);
            for (int fvarSlot = 0; fvarSlot < fvarAttrCount; ++fvarSlot)
            {
                const int builderIndex = faceVaryingBuilderIndices[fvarSlot];
                innerGridFVarTriStart[size_t(fvarSlot)] =
                    attrBuilders[builderIndex].indices.size() / 3;
            }
            for (int y = 0; y < rows - 1; ++y)
            {
                for (int x = 0; x < cols - 1; ++x)
                {
                    const int i00 = innerAt(x + 1, y + 1);
                    const int i10 = innerAt(x + 2, y + 1);
                    const int i01 = innerAt(x + 1, y + 2);
                    const int i11 = innerAt(x + 2, y + 2);
                    const bool emit0 = emitTriToBuffer(
                        i00, i10, i11, i00, i10, i11, true, true, &indices);
                    const bool emit1 = emitTriToBuffer(
                        i00, i11, i01, i00, i11, i01, true, true, &indices);
                    if (emit0)
                    {
                        for (int fvarSlot = 0; fvarSlot < fvarAttrCount; ++fvarSlot)
                        {
                            const int builderIndex = faceVaryingBuilderIndices[fvarSlot];
                            emitTriToBuffer(i00,
                                            i10,
                                            i11,
                                            innerFVarAt(fvarSlot, x + 1, y + 1),
                                            innerFVarAt(fvarSlot, x + 2, y + 1),
                                            innerFVarAt(fvarSlot, x + 2, y + 2),
                                            false,
                                            true,
                                            &attrBuilders[builderIndex].indices);
                        }
                    }
                    if (emit1)
                    {
                        for (int fvarSlot = 0; fvarSlot < fvarAttrCount; ++fvarSlot)
                        {
                            const int builderIndex = faceVaryingBuilderIndices[fvarSlot];
                            emitTriToBuffer(i00,
                                            i11,
                                            i01,
                                            innerFVarAt(fvarSlot, x + 1, y + 1),
                                            innerFVarAt(fvarSlot, x + 2, y + 2),
                                            innerFVarAt(fvarSlot, x + 1, y + 2),
                                            false,
                                            true,
                                            &attrBuilders[builderIndex].indices);
                        }
                    }
                }
            }
            const int emittedInnerGridPosTris = int(indices.size() / 3 - innerGridTriStart);
            for (int fvarSlot = 0; fvarSlot < fvarAttrCount; ++fvarSlot)
            {
                const int builderIndex = faceVaryingBuilderIndices[fvarSlot];
                const int emittedInnerGridFVarTris =
                    int(attrBuilders[builderIndex].indices.size() / 3 -
                        innerGridFVarTriStart[size_t(fvarSlot)]);
                YBI_ASSERT(emittedInnerGridPosTris == emittedInnerGridFVarTris);
            }
        }

        std::vector<int> inner0;
        std::vector<int> inner1;
        std::vector<int> inner2;
        std::vector<int> inner3;
        std::vector<std::vector<int>> inner0FVar{size_t(fvarAttrCount)};
        std::vector<std::vector<int>> inner1FVar{size_t(fvarAttrCount)};
        std::vector<std::vector<int>> inner2FVar{size_t(fvarAttrCount)};
        std::vector<std::vector<int>> inner3FVar{size_t(fvarAttrCount)};
        inner0.reserve(cols);
        inner1.reserve(rows);
        inner2.reserve(cols);
        inner3.reserve(rows);
        for (int f = 0; f < fvarAttrCount; ++f)
        {
            inner0FVar[size_t(f)].reserve(cols);
            inner1FVar[size_t(f)].reserve(rows);
            inner2FVar[size_t(f)].reserve(cols);
            inner3FVar[size_t(f)].reserve(rows);
        }
        for (int iu = 1; iu < nu; ++iu)
        {
            inner0.push_back(innerAt(iu, 1));
            for (int f = 0; f < fvarAttrCount; ++f)
            {
                inner0FVar[size_t(f)].push_back(innerFVarAt(f, iu, 1));
            }
        }
        for (int iv = 1; iv < nv; ++iv)
        {
            inner1.push_back(innerAt(cols, iv));
            for (int f = 0; f < fvarAttrCount; ++f)
            {
                inner1FVar[size_t(f)].push_back(innerFVarAt(f, cols, iv));
            }
        }
        for (int iu = 1; iu < nu; ++iu)
        {
            inner2.push_back(innerAt(iu, rows));
            for (int f = 0; f < fvarAttrCount; ++f)
            {
                inner2FVar[size_t(f)].push_back(innerFVarAt(f, iu, rows));
            }
        }
        for (int iv = 1; iv < nv; ++iv)
        {
            inner3.push_back(innerAt(1, iv));
            for (int f = 0; f < fvarAttrCount; ++f)
            {
                inner3FVar[size_t(f)].push_back(innerFVarAt(f, 1, iv));
            }
        }

        std::reverse(outer2.begin(), outer2.end());
        std::reverse(outer3.begin(), outer3.end());
        for (int f = 0; f < fvarAttrCount; ++f)
        {
            std::reverse(outer2FVar[size_t(f)].begin(), outer2FVar[size_t(f)].end());
            std::reverse(outer3FVar[size_t(f)].begin(), outer3FVar[size_t(f)].end());
        }

        const int mU = std::max(e0, e2);
        const int mV = std::max(e1, e3);
        const size_t stitchingTriStart = indices.size() / 3;
        std::vector<size_t> stitchingFVarTriStart(size_t(fvarAttrCount), 0);
        for (int fvarSlot = 0; fvarSlot < fvarAttrCount; ++fvarSlot)
        {
            const int builderIndex = faceVaryingBuilderIndices[fvarSlot];
            stitchingFVarTriStart[size_t(fvarSlot)] =
                attrBuilders[builderIndex].indices.size() / 3;
        }
        const int emittedStitchPos0 = stitchStripToBuffer(outer0, inner0, outer0, inner0, &indices, mU, e0, true);
        const int emittedStitchPos1 = stitchStripToBuffer(outer1, inner1, outer1, inner1, &indices, mV, e1, true);
        const int emittedStitchPos2 = stitchStripToBuffer(outer2, inner2, outer2, inner2, &indices, mU, e2, true);
        const int emittedStitchPos3 = stitchStripToBuffer(outer3, inner3, outer3, inner3, &indices, mV, e3, true);
        const int emittedStitchPosTotal = emittedStitchPos0 + emittedStitchPos1 + emittedStitchPos2 + emittedStitchPos3;

        for (int fvarSlot = 0; fvarSlot < fvarAttrCount; ++fvarSlot)
        {
            YBI_ASSERT(outer0FVar[size_t(fvarSlot)].size() == outer0.size());
            YBI_ASSERT(inner0FVar[size_t(fvarSlot)].size() == inner0.size());
            YBI_ASSERT(outer1FVar[size_t(fvarSlot)].size() == outer1.size());
            YBI_ASSERT(inner1FVar[size_t(fvarSlot)].size() == inner1.size());
            YBI_ASSERT(outer2FVar[size_t(fvarSlot)].size() == outer2.size());
            YBI_ASSERT(inner2FVar[size_t(fvarSlot)].size() == inner2.size());
            YBI_ASSERT(outer3FVar[size_t(fvarSlot)].size() == outer3.size());
            YBI_ASSERT(inner3FVar[size_t(fvarSlot)].size() == inner3.size());

            const int builderIndex = faceVaryingBuilderIndices[fvarSlot];
            const int emittedFVar0 = stitchStripToBuffer(outer0,
                                                         inner0,
                                                         outer0FVar[size_t(fvarSlot)],
                                                         inner0FVar[size_t(fvarSlot)],
                                                         &attrBuilders[builderIndex].indices,
                                                         mU,
                                                         e0,
                                                         false);
            const int emittedFVar1 = stitchStripToBuffer(outer1,
                                                         inner1,
                                                         outer1FVar[size_t(fvarSlot)],
                                                         inner1FVar[size_t(fvarSlot)],
                                                         &attrBuilders[builderIndex].indices,
                                                         mV,
                                                         e1,
                                                         false);
            const int emittedFVar2 = stitchStripToBuffer(outer2,
                                                         inner2,
                                                         outer2FVar[size_t(fvarSlot)],
                                                         inner2FVar[size_t(fvarSlot)],
                                                         &attrBuilders[builderIndex].indices,
                                                         mU,
                                                         e2,
                                                         false);
            const int emittedFVar3 = stitchStripToBuffer(outer3,
                                                         inner3,
                                                         outer3FVar[size_t(fvarSlot)],
                                                         inner3FVar[size_t(fvarSlot)],
                                                         &attrBuilders[builderIndex].indices,
                                                         mV,
                                                         e3,
                                                         false);
            const int emittedFVarTotal = emittedFVar0 + emittedFVar1 + emittedFVar2 + emittedFVar3;
            YBI_ASSERT(emittedFVarTotal == emittedStitchPosTotal);
            const int emittedFVarBySize =
                int(attrBuilders[builderIndex].indices.size() / 3 -
                    stitchingFVarTriStart[size_t(fvarSlot)]);
            YBI_ASSERT(emittedFVarBySize == emittedStitchPosTotal);
        }
    }

    outMesh->positions = std::move(positions);
    outMesh->indices = std::move(indices);
    outMesh->attributes.clear();
    outMesh->attributes.reserve(attrBuilders.size());
    for (TessellationAttributeBuilder &builder : attrBuilders)
    {
        if (!builder.isFaceVarying)
        {
            builder.indices.Resize(outMesh->indices.size());
            if (outMesh->indices.size() > 0)
            {
                memcpy(builder.indices.data(),
                       outMesh->indices.data(),
                       sizeof(int) * outMesh->indices.size());
            }
        }
        YBI_ASSERT(builder.indices.size() == outMesh->indices.size());
        outMesh->attributes.emplace_back(std::move(builder.values),
                                         std::move(builder.indices),
                                         builder.type,
                                         builder.interpolation,
                                         builder.name);
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
