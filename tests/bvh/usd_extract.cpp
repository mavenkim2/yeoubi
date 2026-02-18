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
#include "io/usd_mesh_io.h"

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

int main(int argc, char **argv)
{
    if (argc != 2)
    {
        printf("Usage: %s <path-to-entry.usda>\n", argv[0]);
        return 1;
    }

    const std::string usdPath = argv[1];
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
        if (!ybi::testio::WriteMeshObj(selectedMesh.path,
                                       selectedMesh.points,
                                       selectedMesh.faceVertexCounts,
                                       selectedMesh.faceVertexIndices,
                                       selectedObjPath,
                                       rotateYUpToZUp,
                                       &usdCameraInfo))
        {
            printf("Failed to write Catmull-Clark control cage OBJ: %s\n",
                   selectedObjPath.string().c_str());
            return 1;
        }
        const fs::path selectedJsonPath = outputDir / "selected_catclark_control_cage.json";
        if (!ybi::testio::WriteSelectedSubdivJson(
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
        if (!ybi::testio::WriteMeshObj(selectedCreasedMesh.path,
                                       selectedCreasedMesh.points,
                                       selectedCreasedMesh.faceVertexCounts,
                                       selectedCreasedMesh.faceVertexIndices,
                                       selectedObjPath,
                                       rotateYUpToZUp,
                                       &usdCameraInfo))
        {
            printf("Failed to write creased Catmull-Clark control cage OBJ: %s\n",
                   selectedObjPath.string().c_str());
            return 1;
        }
        const fs::path selectedJsonPath = outputDir / "selected_catclark_control_cage_creased.json";
        if (!ybi::testio::WriteSelectedSubdivJson(
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
        if (!ybi::testio::WriteMeshObj(bestStoatMesh.mesh.path,
                                       bestStoatMesh.mesh.points,
                                       bestStoatMesh.mesh.faceVertexCounts,
                                       bestStoatMesh.mesh.faceVertexIndices,
                                       stoatPath,
                                       rotateYUpToZUp))
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
        if (!ybi::testio::WriteMeshObj(bestMesh.path,
                                       bestMesh.points,
                                       bestMesh.faceVertexCounts,
                                       bestMesh.faceVertexIndices,
                                       meshPath,
                                       rotateYUpToZUp))
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
        if (!ybi::testio::WriteCurveJson(bestCurve.path,
                                         bestCurve.curveVertexCounts,
                                         bestCurve.points,
                                         bestCurve.widths,
                                         bestCurve.basis,
                                         bestCurve.type,
                                         bestCurve.wrap,
                                         curvePath))
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
