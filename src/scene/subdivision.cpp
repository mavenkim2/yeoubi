#include "scene/subdivision.h"
#include "scene/scene.h"
#include "util/host_memory_arena.h"
#include <memory>
#include <opensubdiv/far/patchTableFactory.h>
#include <opensubdiv/far/primvarRefiner.h>
#include <opensubdiv/far/topologyDescriptor.h>
#include <opensubdiv/sdc/types.h>
#include <opensubdiv/version.h>
#include <vector>

YBI_NAMESPACE_BEGIN

using namespace OpenSubdiv;

struct OsdVertex
{
    float p[3];
    void Clear()
    {
        p[0] = 0.f;
        p[1] = 0.f;
        p[2] = 0.f;
    }
    void AddWithWeight(const OsdVertex &src, float w)
    {
        p[0] += w * src.p[0];
        p[1] += w * src.p[1];
        p[2] += w * src.p[2];
    }
};

static OpenSubdiv::Sdc::Options::VtxBoundaryInterpolation ToOsdVtxBoundary(BoundaryInterpolation r)
{
    switch (r)
    {
        case BOUNDARY_INTERPOLATION_NONE:
            return Sdc::Options::VTX_BOUNDARY_NONE;
        case BOUNDARY_INTERPOLATION_EDGE:
            return Sdc::Options::VTX_BOUNDARY_EDGE_ONLY;
        case BOUNDARY_INTERPOLATION_EDGE_AND_CORNER:
            return Sdc::Options::VTX_BOUNDARY_EDGE_AND_CORNER;
        default:
            return Sdc::Options::VTX_BOUNDARY_EDGE_AND_CORNER;
    }
}

static OpenSubdiv::Sdc::Options::FVarLinearInterpolation ToOsdFVarLinear(FVarLinearInterpolation f)
{
    switch (f)
    {
        case FVAR_LINEAR_NONE:
            return Sdc::Options::FVAR_LINEAR_NONE;
        case FVAR_LINEAR_CORNERS_ONLY:
            return Sdc::Options::FVAR_LINEAR_CORNERS_ONLY;
        case FVAR_LINEAR_CORNERS_PLUS1:
            return Sdc::Options::FVAR_LINEAR_CORNERS_PLUS1;
        case FVAR_LINEAR_CORNERS_PLUS2:
            return Sdc::Options::FVAR_LINEAR_CORNERS_PLUS2;
        case FVAR_LINEAR_BOUNDARIES:
            return Sdc::Options::FVAR_LINEAR_BOUNDARIES;
        case FVAR_LINEAR_ALL:
            return Sdc::Options::FVAR_LINEAR_ALL;
        default:
            return Sdc::Options::FVAR_LINEAR_CORNERS_ONLY;
    }
}

template <typename T>
static void InterpolateValues(MemoryView<T> &view, const Far::TopologyRefiner *refiner)
{
    const int numLevels = refiner->GetNumLevels();
    T *src = view.data();

    for (int level = 1; level < numLevels; level++)
    {
        T *dst = src + ref->GetLevel(level - 1).GetNumVertices();
        Far::PrimvarRefiner(*refiner).Interpolate(level, src, dst);
        src = dst;
    }
}

void Subdivision(const SubdivisionMesh &mesh, int refineLevel)
{
    using TopologyDescriptor = Far::TopologyDescriptor;

    HostMemoryArena arena;

    const int numFaces = static_cast<int>(mesh.vertsPerFace.size());

    TopologyDescriptor desc;
    desc.numVertices = static_cast<int>(mesh.vertices.size());
    desc.numFaces = numFaces;
    desc.numFVarChannels = 0;
    desc.fvarChannels = nullptr;

    desc.numVertsPerFace = mesh.vertsPerFace.data();
    desc.vertIndicesPerFace = mesh.indices.data();

    Sdc::SchemeType scheme = Sdc::SCHEME_CATMARK;
    Sdc::Options options;
    options.SetVtxBoundaryInterpolation(ToOsdVtxBoundary(mesh.interpolationRule));
    options.SetFVarLinearInterpolation(ToOsdFVarLinear(mesh.fvarLinearInterpolation));

    Far::TopologyRefiner *refiner = Far::TopologyRefinerFactory<TopologyDescriptor>::Create(
        desc, Far::TopologyRefinerFactory<TopologyDescriptor>::Options(scheme, options));
    YBI_ASSERT(refiner);

    Far::TopologyRefiner::AdaptiveOptions adaptiveOptions(refineLevel);
    refiner->RefineAdaptive(adaptiveOptions);

    bool hasCreasesOrCorners = mesh.cornerIndices.size() || mesh.creaseIndices.size();
    Far::PatchTableFactory::Options patchOptions(refineLevel);
    patchOptions.endCapType = Far::PatchTableFactory::Options::ENDCAP_GREGORY_BASIS;
    patchOptions.useInfSharpPatch = true;
    patchOptions.generateFVarLegacyLinearPatches = false;

    const Far::PatchTable *patchTable = Far::PatchTableFactory::Create(*refiner, patchOptions);
    YBI_ASSERT(patchTable);

#if 0

    const int numVerticesTotal = ref->GetNumVerticesTotal();
    const int numLocalPointsTotal = pt->GetNumLocalPoints();

    MemoryView<OsdVertex> verts =
        arena.PushArray<OsdVertex>(numVerticesTotal + numLocalPointsTotal);
    const int numCoarse = ref->GetLevel(0).GetNumVertices();

    for (int i = 0; i < numCoarse && i < static_cast<int>(mesh.vertices.size()); i++)
    {
        verts[i].p[0] = mesh.vertices[i].x;
        verts[i].p[1] = mesh.vertices[i].y;
        verts[i].p[2] = mesh.vertices[i].z;
    }

    const int numLevels = ref->GetNumLevels();
    Far::PrimvarRefinerReal<float> primvarRefiner(*ref);

    OsdVertex *src = verts.data();
    for (int level = 1; level < numLevels; level++)
    {
        OsdVertex *dst = src + ref->GetLevel(level - 1).GetNumVertices();
        primvarRefiner.Interpolate(level, src, dst);
        src = dst;
    }

    pt->ComputeLocalPointValues(src, verts.data() + numVerticesTotal);

    const size_t numPatchesTotal = pt->GetNumPatchesTotal();
    result.limitPositions.Resize(numPatchesTotal);

    const int numArrays = pt->GetNumPatchArrays();
    int globalPatchIndex = 0;
    std::vector<float> wP;

    for (int array = 0; array < numArrays; array++)
    {
        const int numPatchesInArray = pt->GetNumPatches(array);
        const int numCVsPerPatch = pt->GetPatchArrayDescriptor(array).GetNumControlVertices();
        wP.resize(numCVsPerPatch);

        for (int patch = 0; patch < numPatchesInArray; patch++)
        {
            Far::PatchTable::PatchHandle handle;
            handle.arrayIndex = array;
            handle.patchIndex = globalPatchIndex;
            handle.vertIndex = patch * numCVsPerPatch;

            pt->EvaluateBasis(handle, 0.5f, 0.5f, wP.data());

            Far::ConstIndexArray cvs = pt->GetPatchVertices(array, patch);
            float x = 0.f, y = 0.f, z = 0.f;
            for (int cv = 0; cv < cvs.size(); cv++)
            {
                int idx = cvs[cv];
                float w = wP[cv];
                const OsdVertex &v = verts[idx];
                x += w * v.p[0];
                y += w * v.p[1];
                z += w * v.p[2];
            }
            result.limitPositions[globalPatchIndex] = make_float3(x, y, z);
            ++globalPatchIndex;
        }
    }

#endif
}

static void BuildRefinedVertexBuffer(const Far::TopologyRefiner *ref,
                                     const Far::PatchTable *pt,
                                     const SubdivisionMesh &controlMesh,
                                     std::vector<OsdVertex> *vertsOut)
{
    using namespace OpenSubdiv;
    const int nRefinerVerts = ref->GetNumVerticesTotal();
    const int nLocalPoints = pt->GetNumLocalPoints();
    if (nRefinerVerts == 0 || nLocalPoints < 0)
        return;

    std::vector<OsdVertex> &verts = *vertsOut;
    verts.resize(nRefinerVerts + nLocalPoints);

    const int nCoarse = ref->GetLevel(0).GetNumVertices();
    for (int i = 0; i < nCoarse && i < static_cast<int>(controlMesh.vertices.size()); i++)
    {
        verts[i].p[0] = controlMesh.vertices[i].x;
        verts[i].p[1] = controlMesh.vertices[i].y;
        verts[i].p[2] = controlMesh.vertices[i].z;
    }

    const int nLevels = ref->GetNumLevels();
    Far::PrimvarRefinerReal<float> primvarRefiner(*ref);
    OsdVertex *src = verts.data();
    for (int level = 1; level < nLevels; level++)
    {
        OsdVertex *dst = src + ref->GetLevel(level - 1).GetNumVertices();
        primvarRefiner.Interpolate(level, src, dst);
        src = dst;
    }

    if (nLocalPoints > 0 && pt->GetLocalPointStencilTable())
        pt->ComputeLocalPointValues(src, verts.data() + nRefinerVerts);
}

static float3 EvaluateLimitAt(const OpenSubdiv::Far::PatchTable *pt,
                              const std::vector<OsdVertex> &verts,
                              int array,
                              int patch,
                              int globalPatchIndex,
                              int numCVsPerPatch,
                              float u,
                              float v,
                              std::vector<float> *wP)
{
    using namespace OpenSubdiv;
    Far::PatchTable::PatchHandle handle;
    handle.arrayIndex = array;
    handle.patchIndex = globalPatchIndex;
    handle.vertIndex = patch * numCVsPerPatch;

    wP->resize(numCVsPerPatch);
    pt->EvaluateBasis(handle, u, v, wP->data());

    Far::ConstIndexArray cvs = pt->GetPatchVertices(array, patch);
    float x = 0.f, y = 0.f, z = 0.f;
    for (int cv = 0; cv < cvs.size(); cv++)
    {
        const int idx = cvs[cv];
        const float w = (*wP)[cv];
        const OsdVertex &v = verts[idx];
        x += w * v.p[0];
        y += w * v.p[1];
        z += w * v.p[2];
    }
    return make_float3(x, y, z);
}

MicropolygonMesh CreateMicropolygons(const SubdivisionResult &subdivisionResult,
                                     const SubdivisionMesh &controlMesh,
                                     const int edgeRateU,
                                     const int edgeRateV)
{
    using namespace OpenSubdiv;
    MicropolygonMesh result;

    if (!subdivisionResult.refiner || !subdivisionResult.patchTable)
        return result;
    if (edgeRateU < 1 || edgeRateV < 1)
        return result;

    const Far::PatchTable *pt = subdivisionResult.patchTable.get();
    const Far::TopologyRefiner *ref = subdivisionResult.refiner.get();

    std::vector<OsdVertex> verts;
    BuildRefinedVertexBuffer(ref, pt, controlMesh, &verts);
    if (verts.empty())
        return result;

    const int numArrays = pt->GetNumPatchArrays();
    const int nu = edgeRateU;
    const int nv = edgeRateV;
    const int vertsPerPatch = (nu + 1) * (nv + 1);
    const int quadsPerPatch = nu * nv;

    size_t totalPositions = 0;
    size_t totalIndices = 0;
    for (int array = 0; array < numArrays; array++)
    {
        const int numPatchesInArray = pt->GetNumPatches(array);
        totalPositions += static_cast<size_t>(numPatchesInArray) * vertsPerPatch;
        totalIndices += static_cast<size_t>(numPatchesInArray) * quadsPerPatch * 4;
    }

    result.positions.Resize(totalPositions);
    result.indices.Resize(totalIndices);

    std::vector<float> wP;
    size_t posOffset = 0;
    size_t idxOffset = 0;
    int globalPatchIndex = 0;

    for (int array = 0; array < numArrays; array++)
    {
        const int numPatchesInArray = pt->GetNumPatches(array);
        const int numCVsPerPatch = pt->GetPatchArrayDescriptor(array).GetNumControlVertices();

        for (int patch = 0; patch < numPatchesInArray; patch++)
        {
            const size_t patchBase = posOffset;

            for (int j = 0; j <= nv; j++)
            {
                const float v = (nv > 0) ? static_cast<float>(j) / static_cast<float>(nv) : 0.f;
                for (int i = 0; i <= nu; i++)
                {
                    const float u =
                        (nu > 0) ? static_cast<float>(i) / static_cast<float>(nu) : 0.f;
                    const float3 p = EvaluateLimitAt(
                        pt, verts, array, patch, globalPatchIndex, numCVsPerPatch, u, v, &wP);
                    result.positions[patchBase + static_cast<size_t>(j) * (nu + 1) + i] = p;
                }
            }

            for (int j = 0; j < nv; j++)
            {
                for (int i = 0; i < nu; i++)
                {
                    const size_t base = patchBase + static_cast<size_t>(j) * (nu + 1) + i;
                    const size_t i0 = base;
                    const size_t i1 = base + 1;
                    const size_t i2 = base + (nu + 1) + 1;
                    const size_t i3 = base + (nu + 1);
                    result.indices[idxOffset++] = static_cast<int>(i0);
                    result.indices[idxOffset++] = static_cast<int>(i1);
                    result.indices[idxOffset++] = static_cast<int>(i2);
                    result.indices[idxOffset++] = static_cast<int>(i3);
                }
            }

            posOffset += vertsPerPatch;
            ++globalPatchIndex;
        }
    }

    return result;
}

YBI_NAMESPACE_END
