#include "load.h"
#include "scene/scene.h"
#include "util/float3.h"
#include "vector_functions.h"
#include <algorithm>
#include <cmath>
#include <limits>
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
#include <pxr/usd/usdGeom/subset.h>
#include <pxr/usd/usdGeom/xform.h>
#include <pxr/usd/usdRender/settings.h>
#include <pxr/usd/usdShade/material.h>
#include <pxr/usd/usdShade/shader.h>
#include <utility>

YBI_NAMESPACE_BEGIN

struct USDTraversalState
{
};

#define USD_ASSERT(expr) \
    { \
        bool result = expr; \
        if (!result) \
        { \
            printf("USD Error in %s (%s:%d)\n", #expr, __FILE__, __LINE__); \
            assert(false); \
        } \
    }

static void ProcessUSDBasisCurve(pxr::UsdGeomBasisCurves &curve, Scene *scene)
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

    uint32_t offset = 0;
    uint32_t curveIndex = 0;

#if 0
    for (auto count : curveVertexCounts)
    {
        printf("curve: %u\n", curveIndex++);
        float minF[3];
        float maxF[3];
        uint32_t expMin[3];
        uint32_t expMax[3];
        for (int i = 0; i < 3; i++)
        {
            minF[i] = std::numeric_limits<float>::infinity();
            maxF[i] = -std::numeric_limits<float>::infinity();

            expMin[i] = 512;
            expMax[i] = 0u;
        }
        for (uint32_t i = offset; i < offset + count; i++)
        {
            auto point = points[i];

            for (int j = 0; j < 3; j++)
            {
                minF[j] = std::min(minF[j], point[j]);
                maxF[j] = std::max(maxF[j], point[j]);

                uint32_t test;
                memcpy(&test, &point[j], sizeof(float));

                uint32_t exponent = (test >> 23) & 0x7f;
                if (exponent == 0)
                {
                    printf("0 exp: %f %f %f\n", point[0], point[1], point[2]);
                }
                expMin[j] = std::min(expMin[j], exponent);
                expMax[j] = std::max(expMax[j], exponent);
            }
        }

        for (uint32_t i = offset; i < offset + count; i++)
        {
            auto point = points[i];

            printf("point: %f %f %f\n", point[0], point[1], point[2]);
        }

        for (int j = 0; j < 3; j++)
        {
            float delta = maxF[j] - minF[j];
            int expDelta = expMax[j] - expMin[j];

            int test;
            memcpy(&test, &delta, sizeof(float));

            printf("exponents: %i %i %i\n", expMin[j], expMax[j], expDelta);
            printf("test delta: floats %u %f %f\n", test, maxF[j], delta);
        }
        offset += count;
    }
#endif

    std::vector<float> testX;
    std::vector<float> testY;
    std::vector<float> testZ;

    std::vector<int> indicesX;
    std::vector<int> indicesY;
    std::vector<int> indicesZ;

    int index = 0;

    for (auto &point : points)
    {
        testX.push_back(point[0]);
        testY.push_back(point[1]);
        testZ.push_back(point[2]);

        indicesX.push_back(index);
        indicesY.push_back(index);
        indicesZ.push_back(index);
        index++;
    }

    std::sort(
        indicesX.begin(), indicesX.end(), [&](size_t i, size_t j) { return testX[i] < testX[j]; });
    std::sort(
        indicesY.begin(), indicesY.end(), [&](size_t i, size_t j) { return testY[i] < testY[j]; });
    std::sort(
        indicesZ.begin(), indicesZ.end(), [&](size_t i, size_t j) { return testZ[i] < testZ[j]; });

    std::sort(testX.begin(), testX.end());
    std::sort(testY.begin(), testY.end());
    std::sort(testZ.begin(), testZ.end());

    int numUnique = 0;
    std::vector<int> mappingX(testX.size());
    std::vector<int> mappingY(testY.size());
    std::vector<int> mappingZ(testZ.size());

    for (int j = 0; j < 3; j++)
    {
        auto &test = j == 0 ? testX : (j == 1 ? testY : testZ);
        auto &indices = j == 0 ? indicesX : (j == 1 ? indicesY : indicesZ);
        auto &mapping = j == 0 ? mappingX : (j == 1 ? mappingY : mappingZ);

        int numDuplicated = 0;

        for (int floatIndex = 0; floatIndex < test.size(); floatIndex++)
        {
            mapping[indices[floatIndex]] = floatIndex;
        }

        int temp = -1;
        float prevFloat = std::nanf("0");
        for (int floatIndex = 0; floatIndex < test.size(); floatIndex++)
        {
            float currentFloat = test[floatIndex];
            if (prevFloat != currentFloat)
            {
                prevFloat = currentFloat;
                temp++;
            }
            else
            {
                numDuplicated++;
            }

            mapping[indices[floatIndex]] = temp;
            test[temp] = currentFloat;

            // int val = indices[mapping[floatIndex]];
            // if (val < numDuplicated)
            // {
            //     printf("abort: %i %i\n", numDuplicated, val);
            // }
            // test[temp] = currentFloat;
            // indices[mapping[floatIndex]] -= numDuplicated;
            // indices[mapping[floatIndex]] = temp;
            // printf("test : %i %i %i, check: %f %f\n",
            //        indices[mapping[floatIndex]],
            //        mapping[floatIndex],
            //        val,
            //        test[indices[mapping[floatIndex]]],
            //        currentFloat);
            // indices.push_back(index);
        }
        // test[temp] = prevFloat;
        numUnique += temp;
    }

    int step = 128;
    int numBits = 0;
    for (int axis = 0; axis < 3; axis++)
    {
        auto &indices = axis == 0 ? indicesX : (axis == 1 ? indicesY : indicesZ);
        auto &mapping = axis == 0 ? mappingX : (axis == 1 ? mappingY : mappingZ);
        auto &test = axis == 0 ? testX : (axis == 1 ? testY : testZ);

        int curveOffset = 0;
        int curveIndex = 0;

        for (int curveCount : curveVertexCounts)
        {
            // printf("curve: %i\n", curveIndex++);
            int minIndex = std::numeric_limits<int>::infinity();
            int maxIndex = 0;

            float minF = std::numeric_limits<float>::infinity();
            float maxF = -std::numeric_limits<float>::infinity();

            for (int j = curveOffset; j < curveOffset + curveCount; j++)
            {
                // if (indices[j] == 0)
                {
                    if (test[mapping[j]] != points[j][axis])
                    {
                        printf("aborot\n");
                    }
                }
                minIndex = mapping[j] < minIndex ? mapping[j] : minIndex;
                maxIndex = mapping[j] > maxIndex ? mapping[j] : maxIndex;

                minF = std::min(points[j][axis], minF);
                maxF = std::max(points[j][axis], maxF);
            }
            {
                int delta = maxIndex - minIndex;
                // printf("min: %f, max: %f, check: %f %f\n",
                //        test[minIndex],
                //        test[maxIndex],
                //        minF,
                //        maxF);
                // printf("axis: %i, delta %i\n", axis, delta);
                int r = 0;
                int x = delta;

                // We check ranges and shift x to narrow down the highest bit
                if (x >= 0x10000)
                {
                    x >>= 16;
                    r += 16;
                }
                if (x >= 0x100)
                {
                    x >>= 8;
                    r += 8;
                }
                if (x >= 0x10)
                {
                    x >>= 4;
                    r += 4;
                }
                if (x >= 0x4)
                {
                    x >>= 2;
                    r += 2;
                }
                if (x >= 0x2)
                {
                    r += 1;
                }

                r++;
                numBits += curveCount * r;
            }
            curveOffset += curveCount;
        }
        // printf("test delta %i r %i, %i %i\n", delta, r, minIndex, maxIndex);
        // float delta = test[i + step - 1] - test[i];

        // uint32_t u;
        // memcpy(&u, &delta, sizeof(float));
        // printf("delta: %u, %f\n", u, delta);
    }
    printf("num bytes: %i, uncompressed: %i\n", numUnique * 4 + numBits / 8, 12 * points.size());
    printf("num unique: %i, total: %i\n", numUnique, points.size());

    std::vector<pxr::UsdGeomSubset> subsets = pxr::UsdGeomSubset::GetAllGeomSubsets(curve);
    if (subsets.size())
    {
        printf("has subsets\n");
    }

    printf("curve: %s %s %s\n", basisToken.GetText(), typeToken.GetText(), wrapToken.GetText());

    int vStep = 1;

    if (basisToken.GetString() == "bezier")
    {
        vStep = 3;
    }

    int totalNumVertices = 0;

    // TODO: really need to get rid of this
    Array<float3> positions(points);
    Array<int> curveOffsets(numCurves);
    int curveOffsetIndex = 0;

    for (int curveVertexCount : curveVertexCounts)
    {
        curveOffsets[curveOffsetIndex++] = totalNumVertices;
        totalNumVertices += curveVertexCount;
    }

    assert(totalNumVertices == points.size());
    scene->curves.emplace_back(std::move(positions), std::move(curveOffsets));

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

    double startTimeCode = stage->GetStartTimeCode();
    double endTimeCode = stage->GetEndTimeCode();
    double fps = stage->GetFramesPerSecond();
    double tcps = stage->GetTimeCodesPerSecond();
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

    scene->curves.reserve(basisCurves.size());
    int curveIndex = 0;
    for (pxr::UsdGeomBasisCurves &curve : basisCurves)
    {
        ProcessUSDBasisCurve(curve, scene);
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
        int numTriangles = 0;

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

        Array<float3> finalPositions(positions);
        Array<int> finalIndices(3 * numTriangles);

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

        scene->meshes.emplace_back(std::move(finalPositions), std::move(finalIndices));
    }
}

YBI_NAMESPACE_END
