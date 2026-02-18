#include <pxr/usd/usd/primRange.h>
#include <pxr/usd/usd/stage.h>
#include <pxr/usd/usdGeom/basisCurves.h>
#include <pxr/usd/usdGeom/mesh.h>
#include <pxr/usd/usdGeom/metrics.h>

#include <filesystem>
#include <fstream>
#include <cctype>
#include <limits>
#include <string>
#include <vector>

#include "usd_subdiv_select.h"
#include "usd_camera_utils.h"

namespace fs = std::filesystem;

namespace
{
struct MeshCandidate
{
    pxr::SdfPath path;
    pxr::VtVec3fArray points;
    pxr::VtIntArray faceVertexCounts;
    pxr::VtIntArray faceVertexIndices;
};

struct CurveCandidate
{
    pxr::SdfPath path;
    pxr::VtIntArray curveVertexCounts;
    pxr::VtVec3fArray points;
    pxr::VtFloatArray widths;
    pxr::TfToken basis;
    pxr::TfToken type;
    pxr::TfToken wrap;
};

struct StoatMeshResult
{
    MeshCandidate mesh;
    int score = std::numeric_limits<int>::min();
};

struct StoatCandidateInfo
{
    std::string path;
    size_t vertices = 0;
    int score = 0;
};

static pxr::GfVec3f RotateToZUpIfNeeded(const pxr::GfVec3f &point, bool rotateYUpToZUp)
{
    if (!rotateYUpToZUp)
    {
        return point;
    }
    // +90 degrees around X: (x, y, z) -> (x, -z, y)
    return pxr::GfVec3f(point[0], -point[2], point[1]);
}

static std::string ToLower(std::string s)
{
    for (char &c : s)
    {
        c = (char)std::tolower((unsigned char)c);
    }
    return s;
}

static bool IsLikelyLightweightMesh(const pxr::VtVec3fArray &points,
                                    const pxr::VtIntArray &faceVertexCounts)
{
    return points.size() > 0 && points.size() <= 50000 && faceVertexCounts.size() <= 50000;
}

static bool IsLikelyLightweightCurve(const pxr::VtVec3fArray &points, const pxr::VtIntArray &curveCounts)
{
    return points.size() > 0 && points.size() <= 50000 && curveCounts.size() <= 20000;
}

static bool WriteMeshObj(const MeshCandidate &mesh,
                         const fs::path &outputPath,
                         bool rotateYUpToZUp,
                         const UsdCameraInfo *usdCameraInfo = nullptr)
{
    std::ofstream out(outputPath, std::ios::out | std::ios::binary);
    if (!out.is_open())
    {
        return false;
    }

    out << "# source_prim " << mesh.path.GetString() << "\n";
    out << "# vertices " << mesh.points.size() << "\n";
    out << "# faces " << mesh.faceVertexCounts.size() << "\n";
    out << "# exported_coordinate_system right_handed_z_up\n";
    if (usdCameraInfo && usdCameraInfo->found)
    {
        out << "# usd_camera_path " << usdCameraInfo->path.GetString() << "\n";
        out << "# usd_camera_distance_to_mesh_center " << usdCameraInfo->distanceToMeshCenter << "\n";
    }

    for (const pxr::GfVec3f &point : mesh.points)
    {
        const pxr::GfVec3f p = RotateToZUpIfNeeded(point, rotateYUpToZUp);
        out << "v " << p[0] << " " << p[1] << " " << p[2] << "\n";
    }

    int indexOffset = 0;
    for (int faceVertexCount : mesh.faceVertexCounts)
    {
        if (faceVertexCount < 3)
        {
            indexOffset += faceVertexCount;
            continue;
        }

        out << "f";
        for (int i = 0; i < faceVertexCount; i++)
        {
            const int index = mesh.faceVertexIndices[indexOffset + i];
            out << " " << (index + 1);
        }
        out << "\n";
        indexOffset += faceVertexCount;
    }

    return out.good();
}

static bool WriteSelectedSubdivJson(const SelectedSubdivMesh &mesh,
                                    const fs::path &outputPath,
                                    bool rotateYUpToZUp,
                                    const UsdCameraInfo *usdCameraInfo = nullptr)
{
    std::ofstream out(outputPath, std::ios::out | std::ios::binary);
    if (!out.is_open())
    {
        return false;
    }

    auto writeIntArray = [&](const char *key, const pxr::VtIntArray &values, bool trailingComma) {
        out << "  \"" << key << "\": [";
        for (size_t i = 0; i < values.size(); i++)
        {
            if (i > 0)
            {
                out << ", ";
            }
            out << values[i];
        }
        out << "]";
        out << (trailingComma ? ",\n" : "\n");
    };

    auto writeFloatArray = [&](const char *key, const pxr::VtFloatArray &values, bool trailingComma) {
        out << "  \"" << key << "\": [";
        for (size_t i = 0; i < values.size(); i++)
        {
            if (i > 0)
            {
                out << ", ";
            }
            out << values[i];
        }
        out << "]";
        out << (trailingComma ? ",\n" : "\n");
    };

    out << "{\n";
    out << "  \"source_prim\": \"" << mesh.path.GetString() << "\",\n";
    out << "  \"scheme\": \"" << mesh.subdivisionScheme << "\",\n";
    out << "  \"points\": [\n";
    for (size_t i = 0; i < mesh.points.size(); i++)
    {
        const pxr::GfVec3f point = RotateToZUpIfNeeded(mesh.points[i], rotateYUpToZUp);
        out << "    [" << point[0] << ", " << point[1] << ", " << point[2] << "]";
        out << (i + 1 < mesh.points.size() ? ",\n" : "\n");
    }
    out << "  ],\n";
    writeIntArray("face_vertex_counts", mesh.faceVertexCounts, true);
    writeIntArray("face_vertex_indices", mesh.faceVertexIndices, true);
    writeIntArray("corner_indices", mesh.cornerIndices, true);
    writeFloatArray("corner_sharpnesses", mesh.cornerSharpnesses, true);
    writeIntArray("crease_indices", mesh.creaseIndices, true);
    writeIntArray("crease_lengths", mesh.creaseLengths, true);
    writeFloatArray("crease_sharpnesses", mesh.creaseSharpnesses, true);
    writeIntArray("hole_indices", mesh.holeIndices, usdCameraInfo && usdCameraInfo->found);
    if (usdCameraInfo && usdCameraInfo->found)
    {
        out << "  \"usd_camera_path\": \"" << usdCameraInfo->path.GetString() << "\",\n";
        out << "  \"usd_camera_distance_to_mesh_center\": " << usdCameraInfo->distanceToMeshCenter << "\n";
    }
    out << "}\n";

    return out.good();
}

static bool WriteCurveJson(const CurveCandidate &curve, const fs::path &outputPath)
{
    std::ofstream out(outputPath, std::ios::out | std::ios::binary);
    if (!out.is_open())
    {
        return false;
    }

    out << "{\n";
    out << "  \"source_prim\": \"" << curve.path.GetString() << "\",\n";
    out << "  \"basis\": \"" << curve.basis.GetString() << "\",\n";
    out << "  \"type\": \"" << curve.type.GetString() << "\",\n";
    out << "  \"wrap\": \"" << curve.wrap.GetString() << "\",\n";
    out << "  \"curve_vertex_counts\": [";
    for (size_t i = 0; i < curve.curveVertexCounts.size(); i++)
    {
        if (i > 0)
        {
            out << ", ";
        }
        out << curve.curveVertexCounts[i];
    }
    out << "],\n";

    out << "  \"points\": [\n";
    for (size_t i = 0; i < curve.points.size(); i++)
    {
        const pxr::GfVec3f &p = curve.points[i];
        out << "    [" << p[0] << ", " << p[1] << ", " << p[2] << "]";
        if (i + 1 < curve.points.size())
        {
            out << ",";
        }
        out << "\n";
    }
    out << "  ],\n";

    out << "  \"widths\": [";
    for (size_t i = 0; i < curve.widths.size(); i++)
    {
        if (i > 0)
        {
            out << ", ";
        }
        out << curve.widths[i];
    }
    out << "]\n";
    out << "}\n";
    return out.good();
}

static int ScoreStoatMeshPath(const std::string &pathLower)
{
    if (pathLower.find("/root/stoat/") != 0)
    {
        return std::numeric_limits<int>::min();
    }

    const char *rejectedTerms[] = {
        "fur",
        "hair",
        "whisker",
        "outfit",
        "cloth",
        "shirt",
        "sweater",
        "pant",
        "shoe",
        "groom",
        "curve",
        "procedural",
    };
    for (const char *term : rejectedTerms)
    {
        if (pathLower.find(term) != std::string::npos)
        {
            return std::numeric_limits<int>::min();
        }
    }

    int score = 0;
    if (pathLower.find("body") != std::string::npos)
    {
        score += 1000;
    }
    if (pathLower.find("body_m_geo") != std::string::npos)
    {
        score += 800;
    }
    if (pathLower.find("body_m_volume_geo") != std::string::npos)
    {
        score -= 300;
    }
    if (pathLower.find("/geo/") != std::string::npos)
    {
        score += 100;
    }
    if (pathLower.find("stoat") != std::string::npos)
    {
        score += 50;
    }
    return score;
}

static void WriteStoatCandidates(const std::vector<StoatCandidateInfo> &candidates,
                                 const fs::path &outputPath)
{
    std::ofstream out(outputPath, std::ios::out | std::ios::binary);
    if (!out.is_open())
    {
        return;
    }
    for (const StoatCandidateInfo &entry : candidates)
    {
        out << entry.path << " vertices=" << entry.vertices << " score=" << entry.score << "\n";
    }
}

static void WriteStoatPrimTypeReport(const std::vector<std::string> &lines, const fs::path &outputPath)
{
    std::ofstream out(outputPath, std::ios::out | std::ios::binary);
    if (!out.is_open())
    {
        return;
    }
    for (const std::string &line : lines)
    {
        out << line << "\n";
    }
}
} // namespace

int main()
{
    const std::string usdPath = "C:/Users/maven/Downloads/ALab-2.2.0/ALab/entry.usda";
    pxr::UsdStageRefPtr stage = pxr::UsdStage::Open(usdPath);
    if (!stage)
    {
        printf("Failed to open USD stage: %s\n", usdPath.c_str());
        return 1;
    }

    const pxr::TfToken stageUpAxis = pxr::UsdGeomGetStageUpAxis(stage);
    const bool rotateYUpToZUp = stageUpAxis == pxr::UsdGeomTokens->y;
    printf("USD stage up-axis: %s (rotate_y_up_to_z_up=%d)\n",
           stageUpAxis.GetString().c_str(),
           rotateYUpToZUp ? 1 : 0);

    const fs::path outputDir = fs::path("tests") / "bvh" / "out";
    fs::create_directories(outputDir);

    SelectedSubdivMesh selectedCatmullClark = {};
    if (SelectLargestCatmullClarkMesh(stage, selectedCatmullClark))
    {
        UsdCameraInfo usdCameraInfo = {};
        GetClosestUsdCameraInfo(stage, selectedCatmullClark.points, usdCameraInfo);

        MeshCandidate selectedMesh = {};
        selectedMesh.path = selectedCatmullClark.path;
        selectedMesh.points = selectedCatmullClark.points;
        selectedMesh.faceVertexCounts = selectedCatmullClark.faceVertexCounts;
        selectedMesh.faceVertexIndices = selectedCatmullClark.faceVertexIndices;

        const fs::path selectedObjPath = outputDir / "selected_catclark_control_cage.obj";
        if (!WriteMeshObj(selectedMesh, selectedObjPath, rotateYUpToZUp, &usdCameraInfo))
        {
            printf("Failed to write Catmull-Clark control cage OBJ: %s\n",
                   selectedObjPath.string().c_str());
            return 1;
        }
        const fs::path selectedJsonPath = outputDir / "selected_catclark_control_cage.json";
        if (!WriteSelectedSubdivJson(
                selectedCatmullClark, selectedJsonPath, rotateYUpToZUp, &usdCameraInfo))
        {
            printf("Failed to write Catmull-Clark control cage JSON: %s\n",
                   selectedJsonPath.string().c_str());
            return 1;
        }
        printf("Catmull-Clark control cage exported: %s\n", selectedObjPath.string().c_str());
        printf("Catmull-Clark control cage data exported: %s\n", selectedJsonPath.string().c_str());
        printf("  prim=%s scheme=%s vertices=%zu faces=%zu\n",
               selectedCatmullClark.path.GetString().c_str(),
               selectedCatmullClark.subdivisionScheme.c_str(),
               selectedCatmullClark.points.size(),
               selectedCatmullClark.faceVertexCounts.size());
        if (usdCameraInfo.found)
        {
            printf("  usd_camera=%s distance_to_mesh_center=%f\n",
                   usdCameraInfo.path.GetString().c_str(),
                   usdCameraInfo.distanceToMeshCenter);
        }
    }
    else
    {
        printf("No Catmull-Clark mesh found in stage.\n");
    }
    SelectedSubdivMesh selectedCatmullClarkWithCreases = {};
    if (SelectLargestCatmullClarkMeshWithCreases(stage, selectedCatmullClarkWithCreases))
    {
        UsdCameraInfo usdCameraInfo = {};
        GetClosestUsdCameraInfo(stage, selectedCatmullClarkWithCreases.points, usdCameraInfo);

        MeshCandidate selectedCreasedMesh = {};
        selectedCreasedMesh.path = selectedCatmullClarkWithCreases.path;
        selectedCreasedMesh.points = selectedCatmullClarkWithCreases.points;
        selectedCreasedMesh.faceVertexCounts = selectedCatmullClarkWithCreases.faceVertexCounts;
        selectedCreasedMesh.faceVertexIndices = selectedCatmullClarkWithCreases.faceVertexIndices;

        const fs::path selectedObjPath = outputDir / "selected_catclark_control_cage_creased.obj";
        if (!WriteMeshObj(selectedCreasedMesh, selectedObjPath, rotateYUpToZUp, &usdCameraInfo))
        {
            printf("Failed to write creased Catmull-Clark control cage OBJ: %s\n",
                   selectedObjPath.string().c_str());
            return 1;
        }
        const fs::path selectedJsonPath = outputDir / "selected_catclark_control_cage_creased.json";
        if (!WriteSelectedSubdivJson(
                selectedCatmullClarkWithCreases, selectedJsonPath, rotateYUpToZUp, &usdCameraInfo))
        {
            printf("Failed to write creased Catmull-Clark control cage JSON: %s\n",
                   selectedJsonPath.string().c_str());
            return 1;
        }
        printf("Creased Catmull-Clark control cage exported: %s\n", selectedObjPath.string().c_str());
        printf("Creased Catmull-Clark control cage data exported: %s\n",
               selectedJsonPath.string().c_str());
        printf("  prim=%s scheme=%s vertices=%zu faces=%zu creases=%zu\n",
               selectedCatmullClarkWithCreases.path.GetString().c_str(),
               selectedCatmullClarkWithCreases.subdivisionScheme.c_str(),
               selectedCatmullClarkWithCreases.points.size(),
               selectedCatmullClarkWithCreases.faceVertexCounts.size(),
               selectedCatmullClarkWithCreases.creaseSharpnesses.size());
        if (usdCameraInfo.found)
        {
            printf("  usd_camera=%s distance_to_mesh_center=%f\n",
                   usdCameraInfo.path.GetString().c_str(),
                   usdCameraInfo.distanceToMeshCenter);
        }
    }
    else
    {
        printf("No Catmull-Clark mesh with creases found in stage.\n");
    }

    const pxr::Usd_PrimFlagsConjunction flags =
        pxr::UsdPrimIsActive && pxr::UsdPrimIsLoaded && !pxr::UsdPrimIsAbstract;
    const pxr::Usd_PrimFlagsPredicate predicate(flags);

    MeshCandidate bestMesh;
    CurveCandidate bestCurve;
    StoatMeshResult bestStoatMesh;
    size_t bestMeshVertexCount = std::numeric_limits<size_t>::max();
    size_t bestCurvePointCount = std::numeric_limits<size_t>::max();
    std::vector<StoatCandidateInfo> stoatCandidateReport;
    std::vector<std::string> stoatPrimTypeLines;

    size_t scannedMeshes = 0;
    size_t scannedCurves = 0;
    const size_t maxMeshesToScan = 48;
    const size_t maxCurvesToScan = 48;

    for (pxr::UsdPrim prim : stage->Traverse(predicate))
    {
        const std::string primPath = prim.GetPath().GetString();
        const std::string pathLower = ToLower(primPath);
        if (pathLower.find("/root/stoat/") == 0)
        {
            const std::string typeName = prim.GetTypeName().GetString();
            std::string line = primPath + " | type=" + typeName;
            if (prim.IsA<pxr::UsdGeomMesh>())
            {
                line += " | IsA:Mesh";
            }
            if (prim.IsA<pxr::UsdGeomBasisCurves>())
            {
                line += " | IsA:BasisCurves";
            }
            stoatPrimTypeLines.push_back(std::move(line));
        }

        if (prim.IsA<pxr::UsdGeomMesh>())
        {
            pxr::UsdGeomMesh mesh(prim);
            const bool isStoatMesh = pathLower.find("/root/stoat/") == 0;

            if (!isStoatMesh && scannedMeshes >= maxMeshesToScan)
            {
                // keep lightweight global sampling bounded, but still allow full stoat search.
            }
            else
            {
                pxr::VtVec3fArray points;
                pxr::VtIntArray faceCounts;
                pxr::VtIntArray faceIndices;
                if (mesh.GetPointsAttr().Get(&points) &&
                    mesh.GetFaceVertexCountsAttr().Get(&faceCounts) &&
                    mesh.GetFaceVertexIndicesAttr().Get(&faceIndices))
                {
                    if (!isStoatMesh)
                    {
                        scannedMeshes++;
                        if (IsLikelyLightweightMesh(points, faceCounts) &&
                            points.size() < bestMeshVertexCount)
                        {
                            bestMesh.path = prim.GetPath();
                            bestMesh.points = points;
                            bestMesh.faceVertexCounts = faceCounts;
                            bestMesh.faceVertexIndices = faceIndices;
                            bestMeshVertexCount = bestMesh.points.size();
                        }
                    }

                    const int baseScore = ScoreStoatMeshPath(pathLower);
                    if (baseScore != std::numeric_limits<int>::min() && points.size() >= 50 &&
                        points.size() <= 500000)
                    {
                        const int score =
                            baseScore + (int)std::min<size_t>(points.size(), 300000) / 100;
                        stoatCandidateReport.push_back({primPath, points.size(), score});
                        if (score > bestStoatMesh.score)
                        {
                            bestStoatMesh.mesh.path = prim.GetPath();
                            bestStoatMesh.mesh.points = std::move(points);
                            bestStoatMesh.mesh.faceVertexCounts = std::move(faceCounts);
                            bestStoatMesh.mesh.faceVertexIndices = std::move(faceIndices);
                            bestStoatMesh.score = score;
                        }
                    }
                }
            }
        }

        if (scannedCurves < maxCurvesToScan && prim.IsA<pxr::UsdGeomBasisCurves>())
        {
            pxr::UsdGeomBasisCurves curves(prim);
            pxr::VtVec3fArray points;
            pxr::VtIntArray curveCounts;
            pxr::VtFloatArray widths;
            pxr::TfToken basis;
            pxr::TfToken type;
            pxr::TfToken wrap;
            if (curves.GetPointsAttr().Get(&points) && curves.GetCurveVertexCountsAttr().Get(&curveCounts))
            {
                curves.GetWidthsAttr().Get(&widths);
                curves.GetBasisAttr().Get(&basis);
                curves.GetTypeAttr().Get(&type);
                curves.GetWrapAttr().Get(&wrap);
                scannedCurves++;
                if (IsLikelyLightweightCurve(points, curveCounts) && points.size() < bestCurvePointCount)
                {
                    bestCurve.path = prim.GetPath();
                    bestCurve.points = std::move(points);
                    bestCurve.curveVertexCounts = std::move(curveCounts);
                    bestCurve.widths = std::move(widths);
                    bestCurve.basis = basis;
                    bestCurve.type = type;
                    bestCurve.wrap = wrap;
                    bestCurvePointCount = bestCurve.points.size();
                }
            }
        }

        if (scannedMeshes >= maxMeshesToScan && scannedCurves >= maxCurvesToScan &&
            bestMeshVertexCount != std::numeric_limits<size_t>::max() &&
            bestCurvePointCount != std::numeric_limits<size_t>::max())
        {
            break;
        }
    }

    WriteStoatCandidates(stoatCandidateReport, outputDir / "stoat_mesh_candidates.txt");
    WriteStoatPrimTypeReport(stoatPrimTypeLines, outputDir / "stoat_prim_types.txt");

    if (bestStoatMesh.score != std::numeric_limits<int>::min())
    {
        const fs::path stoatPath = outputDir / "stoat_body_selected.obj";
        if (!WriteMeshObj(bestStoatMesh.mesh, stoatPath, rotateYUpToZUp))
        {
            printf("Failed to write stoat body OBJ: %s\n", stoatPath.string().c_str());
            return 1;
        }
        printf("Stoat body candidate exported: %s\n", stoatPath.string().c_str());
        printf("  prim=%s vertices=%zu faces=%zu score=%d\n",
               bestStoatMesh.mesh.path.GetString().c_str(),
               bestStoatMesh.mesh.points.size(),
               bestStoatMesh.mesh.faceVertexCounts.size(),
               bestStoatMesh.score);
    }
    else
    {
        printf("No stoat body candidate mesh found.\n");
    }

    if (bestMeshVertexCount != std::numeric_limits<size_t>::max())
    {
        const fs::path meshPath = outputDir / "selected_mesh.obj";
        if (!WriteMeshObj(bestMesh, meshPath, rotateYUpToZUp))
        {
            printf("Failed to write mesh OBJ: %s\n", meshPath.string().c_str());
            return 1;
        }
        printf("Mesh exported: %s\n", meshPath.string().c_str());
        printf("  prim=%s vertices=%zu faces=%zu\n",
               bestMesh.path.GetString().c_str(),
               bestMesh.points.size(),
               bestMesh.faceVertexCounts.size());
    }
    else
    {
        printf("No suitable lightweight UsdGeomMesh found in scanned set.\n");
    }

    if (bestCurvePointCount != std::numeric_limits<size_t>::max())
    {
        const fs::path curvePath = outputDir / "selected_curve.json";
        if (!WriteCurveJson(bestCurve, curvePath))
        {
            printf("Failed to write curve JSON: %s\n", curvePath.string().c_str());
            return 1;
        }
        printf("Curve exported: %s\n", curvePath.string().c_str());
        printf("  prim=%s curves=%zu points=%zu widths=%zu basis=%s type=%s wrap=%s\n",
               bestCurve.path.GetString().c_str(),
               bestCurve.curveVertexCounts.size(),
               bestCurve.points.size(),
               bestCurve.widths.size(),
               bestCurve.basis.GetString().c_str(),
               bestCurve.type.GetString().c_str(),
               bestCurve.wrap.GetString().c_str());
    }
    else
    {
        printf("No suitable lightweight UsdGeomBasisCurves found in scanned set.\n");
    }

    return 0;
}
