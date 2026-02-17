#include "bvh/usd_subdiv_select.h"

#include <pxr/usd/usd/stage.h>

#include <opensubdiv/far/patchMap.h>
#include <opensubdiv/far/patchTable.h>
#include <opensubdiv/far/patchTableFactory.h>
#include <opensubdiv/far/primvarRefiner.h>
#include <opensubdiv/far/ptexIndices.h>
#include <opensubdiv/far/topologyDescriptor.h>
#include <opensubdiv/sdc/types.h>

#include <algorithm>
#include <filesystem>
#include <fstream>
#include <string>
#include <vector>

namespace fs = std::filesystem;
namespace Far = OpenSubdiv::Far;
namespace Sdc = OpenSubdiv::Sdc;

template <typename T>
struct OsdData
{
    T value;

    void Clear()
    {
        value = T(0.0f);
    }
    void AddWithWeight(const OsdData &src, float w)
    {
        value = value + src.value * w;
    }
};

struct TessMesh
{
    std::vector<pxr::GfVec3f> positions;
    std::vector<int> indices;
};

static void BuildCreasePairs(const SelectedSubdivMesh &source,
                             std::vector<int> &creaseVertexPairsOut,
                             std::vector<float> &creaseSharpnessOut)
{
    creaseVertexPairsOut.clear();
    creaseSharpnessOut.clear();

    size_t indexOffset = 0;
    const size_t creaseCount = std::min(source.creaseLengths.size(), source.creaseSharpnesses.size());
    for (size_t crease = 0; crease < creaseCount; crease++)
    {
        const int length = source.creaseLengths[crease];
        if (length < 2)
        {
            indexOffset += size_t(std::max(0, length));
            continue;
        }
        if (indexOffset + size_t(length) > source.creaseIndices.size())
        {
            break;
        }

        const float sharpness = source.creaseSharpnesses[crease];
        for (int segment = 0; segment < length - 1; segment++)
        {
            const int i0 = source.creaseIndices[indexOffset + size_t(segment)];
            const int i1 = source.creaseIndices[indexOffset + size_t(segment + 1)];
            creaseVertexPairsOut.push_back(i0);
            creaseVertexPairsOut.push_back(i1);
            creaseSharpnessOut.push_back(sharpness);
        }
        indexOffset += size_t(length);
    }
}

template <typename T>
static std::vector<OsdData<T>> InterpolateVertex(const std::vector<T> &array,
                                                 const Far::TopologyRefiner *refiner,
                                                 const Far::PatchTable *patchTable)
{
    const int numLevels = refiner->GetNumLevels();
    const int numVertices = refiner->GetNumVerticesTotal();
    const int numLocalPoints = patchTable->GetNumLocalPoints();

    std::vector<OsdData<T>> values(size_t(numVertices + numLocalPoints));
    for (size_t i = 0; i < array.size(); i++)
    {
        values[i].value = array[i];
    }

    OsdData<T> *src = values.data();
    for (int level = 1; level < numLevels; level++)
    {
        OsdData<T> *dst = src + refiner->GetLevel(level - 1).GetNumVertices();
        Far::PrimvarRefiner(*refiner).Interpolate(level, src, dst);
        src = dst;
    }
    patchTable->ComputeLocalPointValues(values.data(), values.data() + numVertices);
    return values;
}

static pxr::GfVec3f EvaluatePosition(const Far::PatchMap &patchMap,
                                     const Far::PatchTable &patchTable,
                                     const std::vector<OsdData<pxr::GfVec3f>> &positions,
                                     int ptexFace,
                                     float u,
                                     float v)
{
    const Far::PatchTable::PatchHandle *handle = patchMap.FindPatch(ptexFace, u, v);
    if (!handle)
    {
        return pxr::GfVec3f(0.0f);
    }

    Far::ConstIndexArray cvIndices = patchTable.GetPatchVertices(*handle);
    std::vector<float> pWeights(size_t(cvIndices.size()));
    patchTable.EvaluateBasis(*handle, u, v, pWeights.data());

    pxr::GfVec3f p(0.0f);
    for (int cv = 0; cv < cvIndices.size(); cv++)
    {
        p += positions[size_t(cvIndices[cv])].value * pWeights[size_t(cv)];
    }
    return p;
}

static TessMesh TessellateLimitSurfaceUniform(const Far::PatchMap &patchMap,
                                              const Far::PatchTable &patchTable,
                                              const std::vector<OsdData<pxr::GfVec3f>> &positions,
                                              int numPtexFaces,
                                              int edgeRate)
{
    TessMesh mesh = {};
    for (int face = 0; face < numPtexFaces; face++)
    {
        const int gridStart = int(mesh.positions.size());
        const int width = edgeRate + 1;
        const int height = edgeRate + 1;

        for (int v = 0; v < height; v++)
        {
            for (int u = 0; u < width; u++)
            {
                const float fu = float(u) / float(edgeRate);
                const float fv = float(v) / float(edgeRate);
                mesh.positions.push_back(
                    EvaluatePosition(patchMap, patchTable, positions, face, fu, fv));
            }
        }

        auto GridIndex = [&](int u, int v) -> int { return gridStart + v * width + u; };
        for (int v = 0; v < edgeRate; v++)
        {
            for (int u = 0; u < edgeRate; u++)
            {
                const int i0 = GridIndex(u, v);
                const int i1 = GridIndex(u + 1, v);
                const int i2 = GridIndex(u + 1, v + 1);
                const int i3 = GridIndex(u, v + 1);
                mesh.indices.push_back(i0);
                mesh.indices.push_back(i1);
                mesh.indices.push_back(i2);
                mesh.indices.push_back(i0);
                mesh.indices.push_back(i2);
                mesh.indices.push_back(i3);
            }
        }
    }
    return mesh;
}

static bool WriteObj(const TessMesh &mesh,
                     const fs::path &path,
                     const SelectedSubdivMesh &source,
                     int edgeRate)
{
    std::ofstream out(path, std::ios::out | std::ios::binary);
    if (!out.is_open())
    {
        return false;
    }

    out << "# source_prim " << source.path.GetString() << "\n";
    out << "# scheme " << source.subdivisionScheme << "\n";
    out << "# control_cage_faces " << source.faceVertexCounts.size() << "\n";
    out << "# edge_rate " << edgeRate << "\n";
    out << "# vertices " << mesh.positions.size() << "\n";
    out << "# triangles " << (mesh.indices.size() / 3) << "\n";

    for (const pxr::GfVec3f &p : mesh.positions)
    {
        out << "v " << p[0] << " " << p[1] << " " << p[2] << "\n";
    }
    for (size_t i = 0; i + 2 < mesh.indices.size(); i += 3)
    {
        out << "f " << (mesh.indices[i + 0] + 1) << " " << (mesh.indices[i + 1] + 1) << " "
            << (mesh.indices[i + 2] + 1) << "\n";
    }
    return out.good();
}

static bool WriteControlCageObj(const SelectedSubdivMesh &source, const fs::path &path)
{
    std::ofstream out(path, std::ios::out | std::ios::binary);
    if (!out.is_open())
    {
        return false;
    }

    out << "# source_prim " << source.path.GetString() << "\n";
    out << "# scheme " << source.subdivisionScheme << "\n";
    out << "# control_cage_faces " << source.faceVertexCounts.size() << "\n";
    out << "# control_cage_vertices " << source.points.size() << "\n";

    for (const pxr::GfVec3f &p : source.points)
    {
        out << "v " << p[0] << " " << p[1] << " " << p[2] << "\n";
    }

    int indexOffset = 0;
    for (int faceVertexCount : source.faceVertexCounts)
    {
        if (faceVertexCount < 3)
        {
            indexOffset += faceVertexCount;
            continue;
        }

        const int i0 = source.faceVertexIndices[indexOffset];
        for (int i = 1; i < faceVertexCount - 1; i++)
        {
            const int i1 = source.faceVertexIndices[indexOffset + i];
            const int i2 = source.faceVertexIndices[indexOffset + i + 1];
            out << "f " << (i0 + 1) << " " << (i1 + 1) << " " << (i2 + 1) << "\n";
        }
        indexOffset += faceVertexCount;
    }

    return out.good();
}

static bool ProcessSelectedMesh(const SelectedSubdivMesh &selected,
                                const fs::path &outputDir,
                                const char *controlCageName,
                                const char *rate2Name,
                                const char *rate8Name)
{
    Far::TopologyDescriptor desc = {};
    desc.numVertices = int(selected.points.size());
    desc.numFaces = int(selected.faceVertexCounts.size());
    desc.numVertsPerFace = selected.faceVertexCounts.cdata();
    desc.vertIndicesPerFace = selected.faceVertexIndices.cdata();
    std::vector<int> creaseVertexPairs;
    std::vector<float> creaseSharpness;
    BuildCreasePairs(selected, creaseVertexPairs, creaseSharpness);
    desc.numCreases = int(creaseSharpness.size());
    desc.creaseVertexIndexPairs = creaseVertexPairs.empty() ? nullptr : creaseVertexPairs.data();
    desc.creaseWeights = creaseSharpness.empty() ? nullptr : creaseSharpness.data();
    desc.numCorners = int(std::min(selected.cornerIndices.size(), selected.cornerSharpnesses.size()));
    desc.cornerVertexIndices = selected.cornerIndices.cdata();
    desc.cornerWeights = selected.cornerSharpnesses.cdata();
    desc.numHoles = int(selected.holeIndices.size());
    desc.holeIndices = selected.holeIndices.cdata();

    const Sdc::SchemeType scheme = Sdc::SCHEME_CATMARK;
    Sdc::Options options;
    Far::TopologyRefiner *refiner = Far::TopologyRefinerFactory<Far::TopologyDescriptor>::Create(
        desc, Far::TopologyRefinerFactory<Far::TopologyDescriptor>::Options(scheme, options));
    if (!refiner)
    {
        fprintf(stderr, "Failed to create topology refiner.\n");
        return false;
    }

    Far::TopologyRefiner::AdaptiveOptions adaptiveOptions(1);
    refiner->RefineAdaptive(adaptiveOptions);

    Far::PatchTableFactory::Options patchOptions(1);
    patchOptions.endCapType = Far::PatchTableFactory::Options::ENDCAP_GREGORY_BASIS;
    const Far::PatchTable *patchTable = Far::PatchTableFactory::Create(*refiner, patchOptions);
    if (!patchTable)
    {
        fprintf(stderr, "Failed to create patch table.\n");
        delete refiner;
        return false;
    }

    Far::PatchMap patchMap(*patchTable);
    Far::PtexIndices ptexIndices(*refiner);
    const int numPtexFaces = ptexIndices.GetNumFaces();

    std::vector<pxr::GfVec3f> cagePositions;
    cagePositions.reserve(selected.points.size());
    for (const pxr::GfVec3f &p : selected.points)
    {
        cagePositions.push_back(p);
    }
    const std::vector<OsdData<pxr::GfVec3f>> refinedPositions =
        InterpolateVertex(cagePositions, refiner, patchTable);

    const fs::path controlCagePath = outputDir / controlCageName;
    if (!WriteControlCageObj(selected, controlCagePath))
    {
        fprintf(stderr, "Failed to write control cage OBJ: %s\n", controlCagePath.string().c_str());
        delete patchTable;
        delete refiner;
        return false;
    }
    printf("Wrote %s (verts=%zu faces=%zu)\n",
           controlCagePath.string().c_str(),
           selected.points.size(),
           selected.faceVertexCounts.size());
    printf("Control-cage features: creaseSegments=%zu corners=%d holes=%d\n",
           creaseSharpness.size(),
           desc.numCorners,
           desc.numHoles);

    const int rates[] = {2, 8};
    const char *names[] = {rate2Name, rate8Name};
    for (int i = 0; i < 2; i++)
    {
        const int rate = rates[i];
        const TessMesh tessMesh =
            TessellateLimitSurfaceUniform(patchMap, *patchTable, refinedPositions, numPtexFaces, rate);
        const fs::path outPath = outputDir / names[i];
        if (!WriteObj(tessMesh, outPath, selected, rate))
        {
            fprintf(stderr, "Failed to write OBJ: %s\n", outPath.string().c_str());
            delete patchTable;
            delete refiner;
            return false;
        }
        printf("Wrote %s (verts=%zu tris=%zu)\n",
               outPath.string().c_str(),
               tessMesh.positions.size(),
               tessMesh.indices.size() / 3);
    }

    delete patchTable;
    delete refiner;
    return true;
}

int main()
{
    const std::string usdPath = "C:/Users/maven/Downloads/ALab-2.2.0/ALab/entry.usda";
    pxr::UsdStageRefPtr stage = pxr::UsdStage::Open(usdPath);
    if (!stage)
    {
        fprintf(stderr, "Failed to open USD stage: %s\n", usdPath.c_str());
        return 1;
    }

    const fs::path outputDir = fs::path("tests") / "bvh" / "out";
    fs::create_directories(outputDir);

    SelectedSubdivMesh selected = {};
    if (!SelectLargestCatmullClarkMesh(stage, selected))
    {
        fprintf(stderr, "No Catmull-Clark mesh found in stage.\n");
        return 1;
    }
    if (!ProcessSelectedMesh(selected,
                             outputDir,
                             "subdiv_control_cage.obj",
                             "subdiv_limit_uniform_r2.obj",
                             "subdiv_limit_uniform_r8.obj"))
    {
        return 1;
    }

    SelectedSubdivMesh selectedWithCreases = {};
    if (SelectLargestCatmullClarkMeshWithCreases(stage, selectedWithCreases))
    {
        if (!ProcessSelectedMesh(selectedWithCreases,
                                 outputDir,
                                 "subdiv_control_cage_creased.obj",
                                 "subdiv_limit_uniform_r2_creased.obj",
                                 "subdiv_limit_uniform_r8_creased.obj"))
        {
            return 1;
        }
    }
    else
    {
        printf("No Catmull-Clark mesh with creases found in stage.\n");
    }

    return 0;
}
