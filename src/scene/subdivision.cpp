#include "scene/subdivision.h"
#include "scene/attributes.h"
#include "scene/scene.h"
#include "util/assert.h"
#include "util/base.h"
#include "util/float2.h"
#include "util/host_memory_arena.h"
#include <algorithm>
#include <opensubdiv/far/patchMap.h>
#include <opensubdiv/far/patchTable.h>
#include <opensubdiv/far/patchTableFactory.h>
#include <opensubdiv/far/primvarRefiner.h>
#include <opensubdiv/far/ptexIndices.h>
#include <opensubdiv/far/topologyDescriptor.h>
#include <opensubdiv/sdc/types.h>
#include <opensubdiv/version.h>
#include <vector>

YBI_NAMESPACE_BEGIN

namespace Far = OpenSubdiv::Far;
namespace Sdc = OpenSubdiv::Sdc;

template <typename T>
struct OsdData
{
    T value;

    void Clear()
    {
        memset(&value, 0, sizeof(T));
    }
    void AddWithWeight(const OsdData &src, float w)
    {
        value += w * src.value;
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
static MemoryView<OsdData> InterpolateVertex(HostMemoryArena &arena,
                                             const T &array,
                                             AttributeType type,
                                             const Far::TopologyRefiner *refiner,
                                             const Far::PatchTable *patchTable)
{

    const int numLevels = refiner->GetNumLevels();
    const int numVertices = refiner->GetNumVerticesTotal();
    const int numLocalPoints = patchTable->GetNumLocalPoints();

    YBI_ASSERT(array.size() == refiner->GetLevel(0).GetNumVertices());

    size_t typeSize = AttributeTypeGetSize(type);
    size_t numFloats = typeSize / sizeof(float);
    MemoryView<OsdData> view = arena.PushArray<OsdData>(numVertices + numLocalPoints);
    MemoryView<float> data = arena.PushArray<float>((numVertices + numLocalPoints) * numFloats);

    memcpy(data.data(), array.data(), array.size() * typeSize);

    for (size_t i = 0; i < numVertices + numLocalPoints; i++)
    {
        view[i].p = &data[i * numFloats];
        view[i].numFloats = numFloats;
    }
    OsdData *src = view.data();

    for (int level = 1; level < numLevels; level++)
    {
        OsdData *dst = src + refiner->GetLevel(level - 1).GetNumVertices();
        Far::PrimvarRefiner(*refiner).Interpolate(level, src, dst);
        src = dst;
    }
    patchTable->ComputeLocalPointValues(view.data(), view.data() + numVertices);
    return view;
}

template <typename T>
static void InterpolateVarying(HostMemoryArena &arena,
                               const T &array,
                               const Far::TopologyRefiner *refiner,
                               const Far::PatchTable *patchTable,
                               int channel)
{
    const int numVerticesTotal = refiner->GetNumFVarValuesTotal(channel);
    const int numLocalPointsTotal = patchTable->GetNumLocalPointsVarying();

    MemoryView<OsdData> view = arena.PushArray<OsdData>(numVerticesTotal + numLocalPointsTotal);

    const int numLevels = refiner->GetNumLevels();
    OsdData *src = view.data();

    for (int level = 1; level < numLevels; level++)
    {
        OsdData *dst = src + refiner->GetLevel(level - 1).GetNumFVarValues(channel);
        Far::PrimvarRefiner(*refiner).InterpolateVarying(level, src, dst);
        src = dst;
    }
    patchTable->ComputeLocalPointValuesVarying(src, view.data() + numVerticesTotal);
}

template <typename T>
static void InterpolateFaceVarying(MemoryView<OsdData> &view,
                                   const T &array,
                                   const Far::TopologyRefiner *refiner,
                                   const Far::PatchTable *patchTable,
                                   int channel)
{
    const int numVerticesTotal = refiner->GetNumFVarValuesTotal(channel);
    const int numLocalPointsTotal = patchTable->GetNumLocalPointsFaceVarying();
    const int numLevels = refiner->GetNumLevels();
    OsdData *src = view.data();

    for (int level = 1; level < numLevels; level++)
    {
        OsdData *dst = src + refiner->GetLevel(level - 1).GetNumFVarValues(channel);
        Far::PrimvarRefiner(*refiner).InterpolateFaceVarying(level, src, dst);
        src = dst;
    }
    patchTable->ComputeLocalPointValuesVarying(src, view.data() + numVerticesTotal);
}

struct LimitSurfaceSample
{
    float2 uv;
    int faceID;
};

static void GenerateLimitSurfaceSamples(const Far::PatchMap *patchMap,
                                        const Far::PatchTable *patchTable)
{
    // for each face
    // get the edge rates of each edge
    // this determines the internal grid
    // determine the edges (if any) that have different edge rates from their opposites

    for (int patch = 0; patch < patchTable->GetNumPatchesTotal(); patch++)
    {
        // const Far::PatchTable::PatchHandle *handle = patchMap->FindPatch(fid);
        // patchTable->GetPatchVertices()
    }
    // patchTable->GetPatchVert
}

namespace DiagSplit
{

static const int NON_UNIFORM = -1;
static const float2 CORNER_TO_UV[4] = {
    make_float2(0, 0),
    make_float2(1, 0),
    make_float2(1, 1),
    make_float2(0, 1),
};
static const float2 UV_DELTA[4] = {
    make_float2(1, 0),
    make_float2(0, 1),
    make_float2(-1, 0),
    make_float2(0, -1),
};

struct DiagSplitParams
{
    const Far::PatchMap *patchMap;
    const Far::PatchTable *patchTable;
    const MemoryView<OsdData<float3>> positions;
    int N;   // number of times to sample edge in T()
    float R; // vertex pixel spacing goal
};

struct SubPatch
{
    float3 corners[4];
};

static void EvaluatePosition() {}

static int T(const DiagSplitParams &params, const float2 &uvStart, const float2 &uvEnd, int edge)
{
    const int N = params.N;
    float R = params.R;
    const int splitThreshold = 0;

    int fid = 0;
    int axis = edge & 1;

    float3 prev;

    YBI_ASSERT(N > 1);

    float maxLi = 0.f;

    for (int i = 1; i < N; i++)
    {
        float2 uv = lerp(uvStart, uvEnd, float(i) / (N - 1));

        const Far::PatchTable::PatchHandle *handle = params.patchMap->FindPatch(fid, uv.x, uv.y);
        auto cvIndices = params.patchTable->GetPatchVertices(*handle);

        float pWeights[20];
        params.patchTable->EvaluateBasis(*handle, uv.x, uv.y, pWeights);
        float3 p = make_float3(0.f);
        for (int cv = 0; cv < cvIndices.size(); cv++)
        {
            p += params.positions[cvIndices[cv]].value * pWeights[cv];
        }

        maxLi = std::max(maxLi, length(p - prev));
        prev = p;
    }

    int tMin = 0.f;
    int tMax = (int)ceilf(N * maxLi / R);

    if (tMax - tMin >= splitThreshold)
    {
        return DiagSplit::NON_UNIFORM;
    }
    return tMax;
}

static void PartitionEdge(const DiagSplitParams &params,
                          float start,
                          float end,
                          int edge,
                          int edgeFactor,
                          int &t0,
                          int &t1)
{
    float2 uvStart;
    float2 uvEnd;

    if (edgeFactor == NON_UNIFORM)
    {
        float2 uv = (uvStart + uvEnd) / 2.f;

        t0 = T(params, uvStart, uv, edge);
        t1 = T(params, uv, uvEnd, edge);
    }
    else
    {
        t0 = (int)floorf(edgeFactor / 2.f);
        t1 = edgeFactor - t0;

        float2 uv = lerp(uvStart, uvEnd, float(t0) / edgeFactor);
    }
}

static void Split() {}

} // namespace DiagSplit

// https://techmatt.github.io/pdfs/diagSplit.pdf

void Subdivision(Scene *scene, const SubdivisionMesh &mesh, int refineLevel)
{
    using TopologyDescriptor = Far::TopologyDescriptor;

    HostMemoryArena arena;

    const int numFaces = static_cast<int>(mesh.vertsPerFace.size());
    MemoryView<TopologyDescriptor::FVarChannel> channels =
        arena.PushArray<TopologyDescriptor::FVarChannel>(mesh.attributeEnd - mesh.attributeStart);
    size_t numChannels = 0;

    for (size_t i = mesh.attributeStart; i < mesh.attributeEnd; i++)
    {
        Attribute &attribute = scene->attributes[i];
        if (AttributeTypeIsFloat(attribute.type))
        {
            size_t numValues = AttributeTypeGetSize(attribute.type) / sizeof(float);
            channels[numChannels].numValues = numValues;
            channels[numChannels].valueIndices = attribute.indices.data();
            numChannels++;
        }
    }

    TopologyDescriptor desc = {};
    desc.numVertices = static_cast<int>(mesh.vertices.size());
    desc.numFaces = numFaces;
    desc.numVertsPerFace = mesh.vertsPerFace.data();
    desc.vertIndicesPerFace = mesh.indices.data();
#if 0
    desc.numFVarChannels = numChannels;
    desc.fvarChannels = channels.data();
#endif
    desc.numCreases = mesh.creaseSharpnesses.size();
    desc.creaseWeights = mesh.creaseSharpnesses.data();
    // TODO: does the spec match?
    desc.creaseVertexIndexPairs = mesh.creaseIndices.data();
    desc.numCorners = mesh.cornerSharpnesses.size();
    desc.cornerWeights = mesh.cornerSharpnesses.data();
    desc.cornerVertexIndices = mesh.cornerIndices.data();
    desc.numHoles = mesh.holeIndices.size();
    desc.holeIndices = mesh.holeIndices.data();

    Sdc::SchemeType scheme = Sdc::SCHEME_CATMARK;
    Sdc::Options options;
    options.SetVtxBoundaryInterpolation(ToOsdVtxBoundary(mesh.interpolationRule));
    options.SetFVarLinearInterpolation(ToOsdFVarLinear(mesh.fvarLinearInterpolation));
    // options.SetCreasingMethod(CreasingMethod c)
    // options.SetTriangleSubdivision(TriangleSubdivision t)

    Far::TopologyRefiner *refiner = Far::TopologyRefinerFactory<TopologyDescriptor>::Create(
        desc, Far::TopologyRefinerFactory<TopologyDescriptor>::Options(scheme, options));
    YBI_ASSERT(refiner);

    Far::TopologyRefiner::AdaptiveOptions adaptiveOptions(1);
    refiner->RefineAdaptive(adaptiveOptions);

    bool hasCreasesOrCorners = bool(desc.numCreases) || bool(desc.numCorners);
    Far::PatchTableFactory::Options patchOptions(refineLevel);
    patchOptions.endCapType = Far::PatchTableFactory::Options::ENDCAP_GREGORY_BASIS;
    patchOptions.useInfSharpPatch = hasCreasesOrCorners;
    patchOptions.generateFVarLegacyLinearPatches = false;

    const Far::PatchTable *patchTable = Far::PatchTableFactory::Create(*refiner, patchOptions);
    Far::PatchMap patchMap(*patchTable);
    YBI_ASSERT(patchTable);

    Far::PtexIndices ptexIndices(*refiner);
    const int numPtexFaces = ptexIndices.GetNumFaces();
    MemoryView<OsdData> positions =
        InterpolateVertex(arena, mesh.vertices, AttributeType::Float3, refiner, patchTable);

#if 0
    const Far::TopologyLevel &level0 = refiner->GetLevel(0);
    const Far::TopologyLevel &level1 = refiner->GetLevel(1);

    // Print vertex positions of quadrilaterals created by quadrangulating non-quad faces.
    // Ordinary quadrilateral faces in the control cage are skipped.
    const int level0VertexOffset = level0.GetNumVertices();
    for (int i = 0; i < numFaces; i++)
    {
        Far::ConstIndexArray faceVerts = level0.GetFaceVertices(i);
        if (faceVerts.size() == 4)
            continue;

        Far::ConstIndexArray childFaces = level0.GetFaceChildFaces(i);
        printf("face %d (%d-gon) -> %d quads:\n",
               i, static_cast<int>(faceVerts.size()), static_cast<int>(childFaces.size()));
        printf("  control verts:");
        for (int j = 0; j < faceVerts.size(); j++)
        {
            const float *p = positions[faceVerts[j]].p;
            printf(" (%.4f, %.4f, %.4f)", p[0], p[1], p[2]);
        }
        printf("\n");

        for (int j = 0; j < childFaces.size(); j++)
        {
            Far::ConstIndexArray quadVerts = level1.GetFaceVertices(childFaces[j]);
            printf("  quad %d:", j);
            for (int k = 0; k < quadVerts.size(); k++)
            {
                const float *p = positions[quadVerts[k] + level0VertexOffset].p;
                printf(" (%.4f, %.4f, %.4f)", p[0], p[1], p[2]);
            }
            printf("\n");
        }
    }

    GenerateLimitSurfaceSamples(&patchMap, patchTable);
#endif

#if 0

    MemoryView<const Attribute *> vertexAttributes =
        arena.PushArray<const Attribute *>(mesh.attributeEnd - mesh.attributeStart);
    size_t numVertexAttributes = 0;
    MemoryView<const Attribute *> varyingAttributes =
        arena.PushArray<const Attribute *>(mesh.attributeEnd - mesh.attributeStart);
    size_t numVaryingAttributes = 0;
    MemoryView<const Attribute *> faceVaryingAttributes =
        arena.PushArray<const Attribute *>(mesh.attributeEnd - mesh.attributeStart);
    size_t numFaceVaryingAttributes = 0;

    auto categorizeAttributes = [&](auto &arr, size_t &count, PrimvarInterpolation interp) {
        for (size_t idx = mesh.attributeStart; idx < mesh.attributeEnd; idx++)
        {
            const Attribute &attr = scene->attributes[idx];
            if (attr.interpolation == interp)
                arr[count++] = &attr;
        }
    };

    categorizeAttributes(vertexAttributes, numVertexAttributes, PrimvarInterpolation::Vertex);
    categorizeAttributes(varyingAttributes, numVaryingAttributes, PrimvarInterpolation::Varying);
    categorizeAttributes(
        faceVaryingAttributes, numFaceVaryingAttributes, PrimvarInterpolation::FaceVarying);

    // CALL_INTERPOLATE(Vertex, InterpolateVertex, vertexAttributes, numVertexAttributes);
    // CALL_INTERPOLATE(Varying, InterpolateVarying, varyingAttributes, numVaryingAttributes);
    // CALL_INTERPOLATE(
    //     FaceVarying, InterpolateFaceVarying, faceVaryingAttributes, numFaceVaryingAttributes);
#endif
}

YBI_NAMESPACE_END
