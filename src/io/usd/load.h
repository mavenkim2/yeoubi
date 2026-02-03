#pragma once

#include "../../scene/scene.h"
#include "pxr/base/vt/types.h"
#include <pxr/usd/usd/prim.h>
#include <pxr/usd/usd/primRange.h>
#include <pxr/usd/usd/stage.h>
#include <pxr/usd/usdGeom/basisCurves.h>
#include <pxr/usd/usdGeom/mesh.h>
#include <pxr/usd/usdGeom/metrics.h>
#include <pxr/usd/usdGeom/xform.h>
#include <vector>

YBI_NAMESPACE_BEGIN

template <typename T, typename V>
using vector = std::vector<T, V>;

static void Test()
{
    std::string filePath = "C:/Users/maven/workspace/ALab-2.2.0/ALab/entry.usda";

    pxr::UsdStageRefPtr stage = pxr::UsdStage::Open(filePath.c_str());

    if (!stage)
    {
        printf("error opening usd stage at %s\n", filePath.c_str());
        return;
    }

    if (pxr::UsdGeomGetStageUpAxis(stage) != pxr::UsdGeomTokens->z)
    {
    }

    double startTimeCode     = stage->GetStartTimeCode();
    double endTimeCode       = stage->GetEndTimeCode();
    double fps               = stage->GetFramesPerSecond();
    double tcps              = stage->GetTimeCodesPerSecond();
    double timeCodesPerFrame = tcps / fps;

    printf("start: %f, end: %f, fps: %f, tcps: %f\n", startTimeCode, endTimeCode, fps, tcps);

    std::vector<pxr::UsdGeomMesh> meshes;

    for (pxr::UsdPrim prim : stage->Traverse())
    {
        if (prim.IsA<pxr::UsdGeomMesh>())
        {
            meshes.push_back(pxr::UsdGeomMesh(prim));
#if 0
            std::vector<double> timeSamples;
            pxr::UsdAttribute pointsAttr = mesh.GetPointsAttr();
            if (mesh.GetVelocitiesAttr().IsValid())
            {
                pxr::VtVec3fArray velocities;
                if (mesh.GetVelocitiesAttr().Get(&velocities))
                {
                }
                else
                {
                    // printf("no v?\n");
                }
            }
            if (pointsAttr.ValueMightBeTimeVarying())
            {
                if (pointsAttr.GetTimeSamples(&timeSamples))
                {
                    for (double time : timeSamples)
                    {
                        pxr::VtVec3fArray positions;
                        pointsAttr.Get(&positions, time);
                        printf("time: %f, p: %f %f %f\n", time, positions[0][0],
                               positions[0][1], positions[0][2]);
                    }
                }
                break;
            }
            // if (mesh.GetFaceVertexCountsAttr().ValueMightBeTimeVarying())
            {
            }
#endif
        }
        else if (prim.IsA<pxr::UsdGeomXform>())
        {
            // printf("xform\n");
        }
        else if (prim.IsA<pxr::UsdGeomBasisCurves>())
        {
            // printf("curves\n");
        }
        else
        {
            // printf("type: %s\n", prim.GetTypeName().GetString().c_str());
        }
    }

    Scene scene;
    scene.meshes.reserve(meshes.size());
    for (pxr::UsdGeomMesh &mesh : meshes)
    {
        pxr::VtVec3fArray positions;
        pxr::VtIntArray faceIndices;
        pxr::VtIntArray faceCounts;
        pxr::TfToken scheme;

        mesh.GetSubdivisionSchemeAttr().Get(&scheme);
        if (scheme == pxr::UsdGeomTokens->catmullClark)
        {
            // printf("cat\n");
        }
        else if (scheme == pxr::UsdGeomTokens->none)
        {
            // printf("poly\n");
        }
        else
        {
            printf("%s\n", scheme.GetText());
        }

        bool success = mesh.GetPointsAttr().Get(&positions);
        assert(success);
        success = mesh.GetFaceVertexIndicesAttr().Get(&faceIndices);
        assert(success);
        success = mesh.GetFaceVertexCountsAttr().Get(&faceCounts);
        assert(success);

        int firstFaceCount     = faceCounts[0];
        bool constantFaceCount = true;

        // TODO: cat clark :)
        for (int faceCount : faceCounts)
        {
            if (faceCount != firstFaceCount)
            {
                // printf("no: %i %i %s\n", faceCount, firstFaceCount, scheme.GetText());
                constantFaceCount = false;
                break;
            }
        }
        // break;
    }
}

YBI_NAMESPACE_END
