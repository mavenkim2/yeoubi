#include "load.h"
#include "../../scene/scene.h"
#include "../../util/float3.h"
#include "pxr/usd/usdGeom/subset.h"
#include "pxr/usd/usdShade/shader.h"
#include <pxr/base/vt/types.h>
#include <pxr/usd/usd/prim.h>
#include <pxr/usd/usd/primRange.h>
#include <pxr/usd/usd/relationship.h>
#include <pxr/usd/usd/stage.h>
#include <pxr/usd/usd/variantSets.h>
#include <pxr/usd/usdGeom/basisCurves.h>
#include <pxr/usd/usdGeom/camera.h>
#include <pxr/usd/usdGeom/curves.h>
#include <pxr/usd/usdGeom/hermiteCurves.h>
#include <pxr/usd/usdGeom/mesh.h>
#include <pxr/usd/usdGeom/metrics.h>
#include <pxr/usd/usdGeom/nurbsCurves.h>
#include <pxr/usd/usdGeom/pointInstancer.h>
#include <pxr/usd/usdGeom/scope.h>
#include <pxr/usd/usdGeom/xform.h>
#include <pxr/usd/usdRender/settings.h>
#include <pxr/usd/usdShade/material.h>
#include <vector>

YBI_NAMESPACE_BEGIN

struct USDTraversalState
{
};

#define USD_ASSERT(expr)                                                                      \
    {                                                                                         \
        bool result = expr;                                                                   \
        if (!result)                                                                          \
        {                                                                                     \
            printf("USD Error in %s (%s:%d)\n", #expr, __FILE__, __LINE__);                   \
            assert(false);                                                                    \
        }                                                                                     \
    }

static void ProcessUSDBasisCurve(pxr::UsdGeomBasisCurves &curve)
{
    size_t numCurves = curve.GetCurveCount(0.0);

    pxr::VtIntArray curveVertexCounts;
    pxr::VtVec3fArray points;
    pxr::VtFloatArray widths;
    pxr::TfToken basisToken;
    pxr::TfToken typeToken;
    pxr::TfToken wrapToken;

    USD_ASSERT(curve.GetCurveVertexCountsAttr().Get(&curveVertexCounts, 0.0));
    USD_ASSERT(curve.GetPointsAttr().Get(&points, 0.0));
    USD_ASSERT(curve.GetWidthsAttr().Get(&widths, 0.0));
    USD_ASSERT(curve.GetBasisAttr().Get(&basisToken, 0.0));
    USD_ASSERT(curve.GetTypeAttr().Get(&typeToken, 0.0));
    USD_ASSERT(curve.GetWrapAttr().Get(&wrapToken, 0.0));

    std::vector<pxr::UsdGeomSubset> subsets = pxr::UsdGeomSubset::GetAllGeomSubsets(curve);
    if (subsets.size())
    {
        printf("has subsets\n");
    }

    printf("curve: %s %s %s\n",
           basisToken.GetText(),
           typeToken.GetText(),
           wrapToken.GetText());

    int vStep = 1;

    if (basisToken.GetString() == "bezier")
    {
        vStep = 3;
    }

    int totalNumVertices = 0;
    int *curveOffsets    = (int *)malloc(numCurves + 1);
    int curveOffsetIndex = 0;

    for (int i = 0; i < 10; i++)
    {
        printf("points: %f %f %f\n", points[i][0], points[i][1], points[i][2]);
    }

    for (int curveVertexCount : curveVertexCounts)
    {
        curveOffsets[curveOffsetIndex++] = totalNumVertices;
        totalNumVertices += curveVertexCount;
    }

    curveOffsets[curveOffsetIndex] = totalNumVertices;
    assert(totalNumVertices == points.size());

    int curveFlags = 0;

    if (typeToken == "cubic")
    {
        curveFlags |= CurveFlags::CURVE_FLAGS_CUBIC;
    }
    else
    {
        assert(typeToken == "linear");
        curveFlags |= CurveFlags::CURVE_FLAGS_LINEAR;
    }

    if (curve.GetNormalsAttr().HasValue())
    {
        pxr::VtVec3fArray normals;

        USD_ASSERT(curve.GetNormalsAttr().Get(&normals, 0.0));
        printf("normal %f %f %f\n", normals[0][0], normals[0][1], normals[0][2]);
        curveFlags |= CurveFlags::CURVE_FLAGS_RIBBON;
    }
    else
    {
        curveFlags |= CurveFlags::CURVE_FLAGS_TUBE;
    }

    printf("basis: %zi %zi %s\n", numCurves, points.size(), basisToken.GetText());
}

void Test(Scene *scene)
{
    std::string filePath = "C:/Users/maven/Downloads/ALab-2.2.0/ALab/entry.usda";

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
    std::vector<pxr::UsdGeomBasisCurves> basisCurves;

    pxr::Usd_PrimFlagsConjunction filterFlags =
        pxr::UsdPrimIsActive && pxr::UsdPrimIsLoaded && !pxr::UsdPrimIsAbstract;

    // if (defined_prims_only)
    // {
    //     filter_flags &= pxr::UsdPrimIsDefined;
    // }

    pxr::Usd_PrimFlagsPredicate filterPredicate(filterFlags);

    // if (!params_.support_scene_instancing)
    // {
    // filterPredicate = pxr::UsdTraverseInstanceProxiesfilterPredicate);
    // }

    // pxr::UsdPrim root = stage->GetPrototypes();

    for (pxr::UsdPrim prim : stage->Traverse(filterPredicate))
    {
        pxr::UsdVariantSets variantSets = prim.GetVariantSets();
        for (const std::string &setNames : variantSets.GetNames())
        {
            printf("variants: %s\n", setNames.c_str());
            printf("type: %s\n", prim.GetTypeName().GetString().c_str());
        }

        if (prim.IsInstance())
        {
            pxr::UsdPrim proto = prim.GetPrototype();
        }
        else if (prim.IsA<pxr::UsdGeomMesh>())
        {
            meshes.push_back(pxr::UsdGeomMesh(prim));
#if 0
            std::vector<double> timeSamples;
            pxr::UsdAttribute pointsAttr = mesh.GetPointsAttr();
            pointsAttr
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
            pxr::UsdGeomXform xform(prim);
            // printf("xform\n");
        }
        else if (prim.IsA<pxr::UsdGeomBasisCurves>())
        {
            basisCurves.push_back(pxr::UsdGeomBasisCurves(prim));
        }
        else if (prim.IsA<pxr::UsdGeomNurbsCurves>())
        {
            pxr::UsdGeomNurbsCurves nurbsCurves(prim);
            size_t numCurves = nurbsCurves.GetCurveCount(0.0);
            // printf("nurbs: %zi\n", numCurves);

            // pxr::VtVec3fArray points;
            // nurbsCurves.GetPointsAttr().Get(&points, 0.0);

            // for (auto point : points)
            // {
            //     printf("%f %f %f\n", point[0], point[1], point[2]);
            // }
            // nurbsCurves.GetOrderAttr();
            // nurbsCurves.GetKnotsAttr();
            // nurbsCurves.GetRangesAttr();
            // nurbsCurves.GetPointWeightsAttr();
        }
        else if (prim.IsA<pxr::UsdGeomHermiteCurves>())
        {
            printf("hermite\n");
        }
        else if (prim.IsA<pxr::UsdGeomCurves>())
        {
            printf("curve\n");
        }
        else if (prim.IsA<pxr::UsdGeomCamera>())
        {
            pxr::UsdGeomCamera camera(prim);
            double shutterOpen, shutterClose;

            // if (camera.GetShutterOpenAttr().ValueMightBeTimeVarying() ||
            //     camera.GetShutterCloseAttr().ValueMightBeTimeVarying())
            {
                camera.GetShutterOpenAttr().Get(&shutterOpen, 0.0);
                camera.GetShutterCloseAttr().Get(&shutterClose, 0.0);
                printf("open: %f close: %f\n", shutterOpen, shutterClose);
            }
        }
        else if (prim.IsA<pxr::UsdRenderSettings>())
        {
            pxr::UsdRenderSettings settings(prim);
        }
        else if (prim.IsA<pxr::UsdShadeMaterial>())
        {
            pxr::UsdShadeMaterial material(prim);
        }
        else if (prim.IsA<pxr::UsdShadeShader>())
        {
            pxr::UsdShadeShader shader(prim);
        }
        else if (prim.IsA<pxr::UsdGeomPointInstancer>())
        {
#if 0
            pxr::SdfPathVector targets;
            pxr::UsdGeomPointInstancer pointInstancer(prim);
            pxr::UsdRelationship rel = pointInstancer.GetPrototypesRel();
            USD_SUCCESS(rel.GetTargets(&targets));

            // for (pxr::SdfPath path : targets)
            // {
            //     printf("path: %s\n", path.GetText());
            // }

            pxr::UsdAttribute protoIndicesAttr = pointInstancer.GetProtoIndicesAttr();

            pxr::VtIntArray protoIndices;
            USD_SUCCESS(protoIndicesAttr.Get(&protoIndices));

            for (int i : protoIndices)
            {
            }

            // pointInstancer.GetPrototypesRel();
            // pointInstancer.GetPositionsAttr
#endif
        }
        else if (prim.IsA<pxr::UsdGeomScope>())
        {
            pxr::UsdGeomScope scope(prim);

            if (prim.GetName().GetString() == "alfro_body_main")
            {
                printf("alfro\n");
            }
        }
        else
        {
            printf("type: %s\n", prim.GetTypeName().GetString().c_str());
        }
    }

    for (pxr::UsdGeomBasisCurves &curve : basisCurves)
    {
        ProcessUSDBasisCurve(curve);
    }

    scene->meshes.reserve(meshes.size());
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

        bool success = mesh.GetPointsAttr().Get(&positions, 0.0);
        assert(success);

        // if (mesh.GetPointsAttr().ValueMightBeTimeVarying())
        // {
        //     std::vector<double> timeSamples;
        //
        //     if (mesh.GetPointsAttr.GetTimeSamplesInInterval(&timeSamples))
        //     {
        //         for (double time : timeSamples)
        //         {
        //             pxr::VtVec3fArray p;
        //             mesh.GetPointsAttr().Get(&p, time);
        //             printf("time: %f, p: %f %f %f %f %f %f\n",
        //                    time,
        //                    positions[0][0],
        //                    positions[0][1],
        //                    positions[0][2],
        //                    p[0][0],
        //                    p[0][1],
        //                    p[0][2]);
        //         }
        //     }
        // }
        // printf("\n");

        success = mesh.GetFaceVertexIndicesAttr().Get(&faceIndices, 0.0);
        assert(success);
        success = mesh.GetFaceVertexCountsAttr().Get(&faceCounts, 0.0);
        assert(success);

        bool constantFaceCount = true;
        int numTriangles       = 0;

        for (int faceCount : faceCounts)
        {
            if (faceCount != 3)
            {
                if (faceCount == 4)
                {
                    numTriangles += 2;
                }
                else
                {
                    printf("n-gon found\n");
                    assert(0);
                }
            }
            else
            {
                numTriangles++;
            }
        }
        size_t positionSize    = sizeof(float3) * positions.size();
        float3 *finalPositions = (float3 *)malloc(positionSize);
        memcpy(finalPositions, positions.data(), positionSize);
        int *finalIndices = (int *)malloc(sizeof(int) * 3 * numTriangles);

        int inputOffset = 0;
        int finalOffset = 0;
        for (int faceCount : faceCounts)
        {
            if (faceCount == 3)
            {
                for (int i = 0; i < 3; i++)
                {
                    int index = faceIndices[inputOffset++];
                    assert(index < positions.size());
                    finalIndices[finalOffset++] = index;
                }
            }
            else if (faceCount == 4)
            {
                // TODO: handle concavity/non-planarity?
                int tempIndices[4];
                for (int i = 0; i < 4; i++)
                {
                    int index = faceIndices[inputOffset++];
                    assert(index < positions.size());
                    tempIndices[i] = index;
                }
                finalIndices[finalOffset++] = tempIndices[0];
                finalIndices[finalOffset++] = tempIndices[1];
                finalIndices[finalOffset++] = tempIndices[2];

                finalIndices[finalOffset++] = tempIndices[0];
                finalIndices[finalOffset++] = tempIndices[2];
                finalIndices[finalOffset++] = tempIndices[3];
            }
            else
            {
                assert(0);
            }
        }

        scene->meshes.push_back(
            Mesh(finalPositions, finalIndices, positions.size(), 3 * numTriangles));
    }
}

YBI_NAMESPACE_END
