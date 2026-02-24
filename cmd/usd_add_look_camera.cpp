#include <pxr/base/gf/matrix4d.h>
#include <pxr/base/gf/range3d.h>
#include <pxr/base/gf/vec2f.h>
#include <pxr/base/gf/vec3d.h>
#include <pxr/usd/usd/primRange.h>
#include <pxr/usd/usd/stage.h>
#include <pxr/usd/usdGeom/bboxCache.h>
#include <pxr/usd/usdGeom/camera.h>
#include <pxr/usd/usdGeom/imageable.h>
#include <pxr/usd/usdGeom/metrics.h>
#include <pxr/usd/usdGeom/tokens.h>
#include <pxr/usd/usdGeom/xform.h>
#include <pxr/usd/usdGeom/xformable.h>

#include <algorithm>
#include <cctype>
#include <cstdio>
#include <cstdlib>
#include <string>
#include <vector>

namespace
{

void PrintUsage(const char *exe)
{
    std::fprintf(stderr,
                 "Usage: %s <input.usda> <target-prim-path-or-name> <output.usda> [camera_distance_scale]\n",
                 exe);
    std::fprintf(stderr,
                 "       %s --resolve <input.usda> <target-prim-path-or-name>\n",
                 exe);
}

pxr::GfVec3d ChooseOffsetDirection(const pxr::TfToken &upAxis)
{
    if (upAxis == pxr::UsdGeomTokens->y)
    {
        return pxr::GfVec3d(0.0, 0.25, 1.0).GetNormalized();
    }
    return pxr::GfVec3d(0.0, -1.0, 0.25).GetNormalized();
}

pxr::GfVec3d ChooseUpDirection(const pxr::TfToken &upAxis)
{
    if (upAxis == pxr::UsdGeomTokens->y)
    {
        return pxr::GfVec3d(0.0, 1.0, 0.0);
    }
    return pxr::GfVec3d(0.0, 0.0, 1.0);
}

std::string ToLower(std::string s)
{
    for (char &c : s)
    {
        c = (char)std::tolower((unsigned char)c);
    }
    return s;
}

std::string ExtractHintToken(const std::string &targetPrimPathOrName)
{
    const size_t slash = targetPrimPathOrName.find_last_of('/');
    if (slash == std::string::npos)
    {
        return targetPrimPathOrName;
    }
    if (slash + 1 >= targetPrimPathOrName.size())
    {
        return std::string();
    }
    return targetPrimPathOrName.substr(slash + 1);
}

pxr::SdfPath ResolveTargetPath(const pxr::UsdStageRefPtr &stage, const std::string &targetPrimPathOrName)
{
    if (!stage)
    {
        return pxr::SdfPath();
    }

    pxr::SdfPath asPath(targetPrimPathOrName);
    if (asPath.IsAbsolutePath())
    {
        const pxr::UsdPrim prim = stage->GetPrimAtPath(asPath);
        if (prim)
        {
            return asPath;
        }
    }

    const std::string hint = ExtractHintToken(targetPrimPathOrName);
    if (hint.empty())
    {
        return pxr::SdfPath();
    }
    const std::string hintLower = ToLower(hint);

    std::vector<pxr::UsdPrim> candidates;
    candidates.reserve(32);
    for (const pxr::UsdPrim &prim : stage->Traverse())
    {
        if (!prim.IsA<pxr::UsdGeomImageable>())
        {
            continue;
        }
        const std::string pathLower = ToLower(prim.GetPath().GetString());
        const std::string nameLower = ToLower(prim.GetName().GetString());
        if (nameLower == hintLower || pathLower.find(hintLower) != std::string::npos)
        {
            candidates.push_back(prim);
        }
    }

    if (candidates.empty())
    {
        return pxr::SdfPath();
    }

    const pxr::TfTokenVector includedPurposes = {
        pxr::UsdGeomTokens->default_,
        pxr::UsdGeomTokens->render,
        pxr::UsdGeomTokens->proxy,
        pxr::UsdGeomTokens->guide,
    };
    pxr::UsdGeomBBoxCache bboxCache(pxr::UsdTimeCode::Default(), includedPurposes);

    double bestDiagonal = -1.0;
    pxr::SdfPath bestPath;
    for (const pxr::UsdPrim &candidate : candidates)
    {
        const pxr::GfBBox3d worldBounds = bboxCache.ComputeWorldBound(candidate);
        const pxr::GfRange3d range = worldBounds.ComputeAlignedRange();
        if (range.IsEmpty())
        {
            continue;
        }
        const double diagonal = range.GetSize().GetLength();
        if (diagonal > bestDiagonal)
        {
            bestDiagonal = diagonal;
            bestPath = candidate.GetPath();
        }
    }
    if (!bestPath.IsEmpty())
    {
        return bestPath;
    }

    return candidates.front().GetPath();
}

bool AddCameraToStage(const pxr::UsdStageRefPtr &stage,
                      const std::string &targetPrimPathOrName,
                      double distanceScale)
{
    if (!stage)
    {
        return false;
    }

    const pxr::SdfPath targetPath = ResolveTargetPath(stage, targetPrimPathOrName);
    if (targetPath.IsEmpty())
    {
        std::fprintf(stderr,
                     "Target prim not found by path/name: %s\n",
                     targetPrimPathOrName.c_str());
        return false;
    }

    const pxr::UsdPrim targetPrim = stage->GetPrimAtPath(targetPath);
    if (!targetPrim)
    {
        std::fprintf(stderr, "Target prim not found: %s\n", targetPath.GetText());
        return false;
    }

    const pxr::TfTokenVector includedPurposes = {
        pxr::UsdGeomTokens->default_,
        pxr::UsdGeomTokens->render,
        pxr::UsdGeomTokens->proxy,
        pxr::UsdGeomTokens->guide,
    };
    pxr::UsdGeomBBoxCache bboxCache(pxr::UsdTimeCode::Default(), includedPurposes);
    const pxr::GfBBox3d worldBounds = bboxCache.ComputeWorldBound(targetPrim);
    const pxr::GfRange3d range = worldBounds.ComputeAlignedRange();
    if (range.IsEmpty())
    {
        std::fprintf(stderr, "Target prim has empty world bounds: %s\n", targetPath.GetText());
        return false;
    }

    const pxr::TfToken upAxis = pxr::UsdGeomGetStageUpAxis(stage);
    const pxr::GfVec3d center = range.GetMidpoint();
    const double diagonal = std::max(1e-3, range.GetSize().GetLength());
    const double cameraDistance = std::max(0.1, diagonal * std::max(0.1, distanceScale));
    const pxr::GfVec3d eye = center + ChooseOffsetDirection(upAxis) * cameraDistance;
    const pxr::GfVec3d up = ChooseUpDirection(upAxis);

    const pxr::SdfPath renderPath("/Render");
    const pxr::SdfPath cameraPath("/Render/Camera_yeoubi");
    pxr::UsdGeomXform::Define(stage, renderPath);
    pxr::UsdGeomCamera camera = pxr::UsdGeomCamera::Define(stage, cameraPath);

    pxr::UsdGeomXformable xformable(camera.GetPrim());
    xformable.ClearXformOpOrder();
    pxr::GfMatrix4d view(1.0);
    view.SetLookAt(eye, center, up);
    const pxr::GfMatrix4d cameraLocalToWorld = view.GetInverse();
    xformable.AddTransformOp().Set(cameraLocalToWorld);

    const float nearClip = std::max(0.01f, static_cast<float>(diagonal * 0.001));
    const float farClip = std::max(nearClip + 1.0f, static_cast<float>(cameraDistance * 10.0));
    camera.CreateClippingRangeAttr().Set(pxr::GfVec2f(nearClip, farClip));
    camera.CreateFocalLengthAttr().Set(35.0f);

    std::fprintf(stdout,
                 "Camera authored at %s eye=(%.3f, %.3f, %.3f) target=(%.3f, %.3f, %.3f)\n",
                 cameraPath.GetText(),
                 eye[0],
                 eye[1],
                 eye[2],
                 center[0],
                 center[1],
                 center[2]);
    return true;
}

} // namespace

int main(int argc, char **argv)
{
    if (argc == 4 && std::string(argv[1]) == "--resolve")
    {
        const std::string inputPath = argv[2];
        const std::string targetPrimPathOrName = argv[3];
        pxr::UsdStageRefPtr stage = pxr::UsdStage::Open(inputPath);
        if (!stage)
        {
            std::fprintf(stderr, "Failed to open input stage: %s\n", inputPath.c_str());
            return 1;
        }
        const pxr::SdfPath targetPath = ResolveTargetPath(stage, targetPrimPathOrName);
        if (targetPath.IsEmpty())
        {
            std::fprintf(stderr,
                         "Target prim not found by path/name: %s\n",
                         targetPrimPathOrName.c_str());
            return 1;
        }
        std::printf("%s\n", targetPath.GetText());
        return 0;
    }

    if (argc < 4 || argc > 5)
    {
        PrintUsage(argv[0]);
        return 1;
    }

    const std::string inputPath = argv[1];
    const std::string targetPrimPathOrName = argv[2];
    const std::string outputPath = argv[3];
    const double distanceScale = argc == 5 ? std::atof(argv[4]) : 1.5;

    pxr::UsdStageRefPtr stage = pxr::UsdStage::Open(inputPath);
    if (!stage)
    {
        std::fprintf(stderr, "Failed to open input stage: %s\n", inputPath.c_str());
        return 1;
    }

    if (!AddCameraToStage(stage, targetPrimPathOrName, distanceScale))
    {
        return 1;
    }

    if (!stage->Export(outputPath))
    {
        std::fprintf(stderr, "Failed to export output stage: %s\n", outputPath.c_str());
        return 1;
    }

    std::fprintf(stdout, "Wrote stage with camera: %s\n", outputPath.c_str());
    return 0;
}
