#include "io/usd/load.h"
#include "scene/scene.h"
#include "util/float3.h"
#include "util/float3x4.h"
#include "util/float4.h"

#include <algorithm>
#include <cmath>
#include <limits>
#include <pxr/base/gf/matrix4f.h>
#include <pxr/base/gf/vec4f.h>
#include <pxr/base/vt/types.h>
#include <pxr/usd/sdf/path.h>
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
#include <pxr/usd/usdGeom/primvarsAPI.h>
#include <pxr/usd/usdGeom/scope.h>
#include <pxr/usd/usdGeom/subset.h>
#include <pxr/usd/usdGeom/xform.h>
#include <pxr/usd/usdRender/settings.h>
#include <pxr/usd/usdShade/material.h>
#include <pxr/usd/usdShade/materialBindingAPI.h>
#include <pxr/usd/usdShade/shader.h>
#include <pxr/usd/usdVol/volume.h>
#include <unordered_map>
#include <utility>

#include <immintrin.h>

YBI_NAMESPACE_BEGIN

#define USD_ASSERT(expr) \
    { \
        bool result = expr; \
        if (!result) \
        { \
            printf("USD Error in %s (%s:%d)\n", #expr, __FILE__, __LINE__); \
            assert(false); \
        } \
    }

struct USDTraversalState
{
    std::vector<pxr::UsdGeomMesh> meshes;
    std::vector<pxr::UsdGeomBasisCurves> basisCurves;
    std::vector<pxr::UsdGeomPointInstancer> pointInstancers;
    std::vector<pxr::UsdPrim> instances;
};

static void TraversePrim(pxr::UsdPrim &root,
                         USDTraversalState &state,
                         pxr::Usd_PrimFlagsPredicate filterPredicate)
{
    const uint32_t stackMax = 1024;
    uint32_t stackTop = 0;
    pxr::UsdPrim stack[stackMax];
    stack[stackTop++] = root;

    while (stackTop)
    {
        uint32_t stackIndex = --stackTop;
        pxr::UsdPrim prim = stack[stackIndex];

        bool pushChildren = true;
        if (prim.IsInstance())
        {
            state.instances.push_back(prim);
        }
        else if (prim.IsA<pxr::UsdGeomPointInstancer>())
        {
            pxr::UsdGeomPointInstancer pointInstancer(prim);
            state.pointInstancers.push_back(pointInstancer);
            pushChildren = false;
        }
        else if (prim.IsA<pxr::UsdGeomMesh>())
        {
            state.meshes.push_back(pxr::UsdGeomMesh(prim));
        }
        else if (prim.IsA<pxr::UsdGeomXform>())
        {
            pxr::UsdGeomXform xform(prim);
            // printf("xform\n");
        }
        else if (prim.IsA<pxr::UsdGeomBasisCurves>())
        {
            state.basisCurves.push_back(pxr::UsdGeomBasisCurves(prim));
        }
        else if (prim.IsA<pxr::UsdGeomNurbsCurves>())
        {
            pxr::UsdGeomNurbsCurves nurbsCurves(prim);
            size_t numCurves = nurbsCurves.GetCurveCount(0.0);
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
        else if (prim.IsA<pxr::UsdShadeShader>())
        {
            pxr::UsdShadeShader shader(prim);
            std::vector<pxr::UsdShadeInput> inputs = shader.GetInputs();
            for (pxr::UsdShadeInput &input : inputs)
            {
                // printf("help: %s %s\n",
                //        input.GetFullName().GetText(),
                //        input.GetTypeName().GetAsToken().GetText());
            }

            pxr::TfToken shaderId;
            if (shader.GetShaderId(&shaderId))
            {
                // printf("help: %s\n", shaderId.GetText());
            }
            else
            {
            }
        }
        else if (prim.IsA<pxr::UsdGeomScope>())
        {
            pxr::UsdGeomScope scope(prim);
        }
        else if (prim.IsA<pxr::UsdVolVolume>())
        {
            pxr::UsdVolVolume volume(prim);
            pxr::TfToken fieldPath;
            volume.GetFieldPath(fieldPath);
            printf("volume: %s\n", fieldPath.GetText());
        }
        else
        {
            // printf("type: %s\n", prim.GetTypeName().GetString().c_str());
        }

        if (pushChildren)
        {
            pxr::UsdPrimSiblingRange range = prim.GetFilteredChildren(filterPredicate);
            for (pxr::UsdPrim &prim : range)
            {
                if (stackTop >= stackMax)
                {
                    printf("help: %i %i\n", stackTop, stackMax);
                }
                assert(stackTop < stackMax);
                stack[stackTop++] = prim;
            }
        }
    }
}

static pxr::UsdShadeMaterial GetPrimMaterial(pxr::UsdPrim &prim,
                                             pxr::TfToken token = pxr::UsdShadeTokens->full)
{
    pxr::UsdShadeMaterialBindingAPI bindingApi(prim);
    // TODO: set this as an option
    pxr::UsdShadeMaterial material = bindingApi.ComputeBoundMaterial(token);

    if (material)
    {
        // printf("Prim %s is bound to material: %s\n",
        //        prim.GetPath().GetText(),
        //        material.GetPath().GetText());
    }
    else
    {
        // printf("Prim %s has no material bound.\n", prim.GetPath().GetText());
    }

    return material;
}

static int AddMaterialToMap(std::unordered_map<std::string, int> &materialMap,
                            std::vector<pxr::UsdShadeMaterial> &materials,
                            const pxr::UsdShadeMaterial &material)
{
    int materialIndex = -1;
    if (material)
    {
        pxr::SdfPath sdfPath = material.GetPath();
        if (!sdfPath.IsEmpty())
        {
            std::string path = sdfPath.GetString();

            auto found = materialMap.find(path);
            if (found == materialMap.end())
            {
                materialIndex = static_cast<int>(materials.size());
                materialMap.emplace(path, materialIndex);
                materials.push_back(material);
            }
            else
            {
                materialIndex = found->second;
            }
        }
    }
    return materialIndex;
}

// NOTE: For shade inputs that most likely have a single UsdUVTexture
static pxr::UsdShadeShader HandleExpectedUVTexture(const pxr::UsdShadeInput &input)
{
    pxr::TfToken token;

    auto sources = input.GetConnectedSources();
    pxr::UsdShadeShader shader;

    if (sources.size() == 1)
    {
        pxr::UsdPrim sourcePrim = sources[0].source.GetPrim();
        if (sourcePrim.IsA<pxr::UsdShadeShader>())
        {
            pxr::UsdShadeShader connectedShader(sourcePrim);
            pxr::TfToken newTok;
            connectedShader.GetShaderId(&newTok);
            if (newTok == pxr::TfToken("UsdUVTexture"))
            {
                shader = connectedShader;
            }
            else
            {
                printf("Expected UsdUvTexture, other cases not handled yet\n");
                assert(0);
            }
        }
        else
        {
            printf("Expected UsdShadeShader, other cases not handled yet. Type name is: %s\n",
                   sourcePrim.GetTypeName().GetText());
            assert(0);
        }
    }
    else
    {
        printf("Expected one source, other cases not handled yet\n");
        assert(0);
    }
    return shader;
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
    uint32_t numBits = 0;
    for (auto count : curveVertexCounts)
    {
        numBits += 32 * 3;

        float minF[3];
        float maxF[3];
        uint32_t expMin[3];
        uint32_t expMax[3];
        uint32_t mantissaMin[3];
        uint32_t mantissaMax[3];
        for (int i = 0; i < 3; i++)
        {
            minF[i] = std::numeric_limits<float>::infinity();
            maxF[i] = -std::numeric_limits<float>::infinity();

            expMin[i] = 512;
            expMax[i] = 0u;

            mantissaMin[i] = ~0u;
            mantissaMax[i] = 0;
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

                uint32_t exponent = (test >> 23) & 0xff;
                uint32_t mantissa = test & 0x7fffff;
                if (exponent == 0)
                {
                    // printf("0 exp: %f %f %f, %u\n", point[0], point[1], point[2], test);
                }
                expMin[j] = std::min(expMin[j], exponent);
                expMax[j] = std::max(expMax[j], exponent);

                mantissaMin[j] = std::min(mantissaMin[j], mantissa);
                mantissaMax[j] = std::max(mantissaMax[j], mantissa);
            }
        }

        if (numCurves == 3689911)
        {
            printf("curve: %u\n", curveIndex++);
            for (uint32_t i = offset; i < offset + count; i++)
            {
                auto point = points[i];
                printf("point: %f %f %f\n", point[0], point[1], point[2]);
            }
        }

        numBits += 8 * 3 + 23 * 3;
        for (int j = 0; j < 3; j++)
        {
            float delta = maxF[j] - minF[j];
            int expDelta = expMax[j] - expMin[j];
            int mantissaDelta = mantissaMax[j] - mantissaMin[j];

            int numExpBits, numMantissaBits;
            if (expDelta == 0)
            {
                numExpBits = 0;
            }
            else
            {
                _BitScanReverse((unsigned long *)&numExpBits, expDelta);
                numExpBits++;
            }

            if (mantissaDelta == 0)
            {
                numMantissaBits = 0;
            }
            else
            {
                _BitScanReverse((unsigned long *)&numMantissaBits, mantissaDelta);
                numMantissaBits++;
            }
            numBits += count * (numExpBits + numMantissaBits);

            int test;
            memcpy(&test, &delta, sizeof(float));

            // printf("exponents: %i %i %i\n", expMin[j], expMax[j], expDelta);
            // printf("mant: %i %i %i\n", mantissaMin[j], mantissaMax[j], mantissaDelta);
            // printf("test delta: floats %u %f %f\n", test, maxF[j], delta);
        }
        offset += count;
    }
    printf("uncompressed: %i, compressed: %i\n", points.size() * 12, numBits / 8);
#endif

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

    Array<float3> positions(points);
    Array<float> curveWidths(widths);
    Array<int> curveOffsets(numCurves);
    int curveOffsetIndex = 0;

    for (int curveVertexCount : curveVertexCounts)
    {
        curveOffsets[curveOffsetIndex++] = totalNumVertices;
        totalNumVertices += curveVertexCount;
    }

    assert(totalNumVertices == points.size());
    scene->curves.emplace_back(
        std::move(positions), std::move(curveWidths), std::move(curveOffsets));

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

static void ProcessUSDPointInstancer(pxr::UsdGeomPointInstancer &pointInstancer,
                                     Scene *scene,
                                     std::unordered_map<std::string, int> &instanceMap,
                                     pxr::UsdTimeCode time = 0.0)
{
    pxr::VtIntArray protoIndices;
    pxr::VtVec3fArray positions;
    pxr::VtQuatfArray orientations;
    pxr::VtVec3fArray scales;
    pxr::SdfPathVector prototypePaths;

    USD_ASSERT(pointInstancer.GetProtoIndicesAttr().Get(&protoIndices, time));
    USD_ASSERT(pointInstancer.GetPositionsAttr().Get(&positions, time));
    USD_ASSERT(pointInstancer.GetPrototypesRel().GetTargets(&prototypePaths));

    pxr::UsdAttribute orientationsAttr = pointInstancer.GetOrientationsfAttr();
    pxr::UsdAttribute scaleAttr = pointInstancer.GetScalesAttr();

    if (orientationsAttr.HasValue())
    {
        printf("has orientations\n");
        orientationsAttr.Get(&orientations, time);
    }
    if (scaleAttr.HasValue())
    {
        scaleAttr.Get(&scales, time);
    }

    Array<float3x4> affineTransforms(protoIndices.size());
    Array<int> objectIDs(protoIndices.size());

    for (pxr::SdfPath &path : prototypePaths)
    {
        printf("%s\n", path.GetText());
    }
    // so whwat's the idea?

    // 1. get all prototype paths, for all types of instances
    // 2. load the prototypes, with whatever data they have
    // 3.

    for (size_t index = 0; index < protoIndices.size(); index++)
    {
        int objectID = 0;
        int instanceIndex = protoIndices[index];

        pxr::GfMatrix4f transform;
        transform.SetIdentity();

        if (!scales.empty())
        {
            transform.SetScale(scales[index]);
        }
        if (!orientations.empty())
        {
            transform.SetRotate(orientations[index]);
        }
        transform.SetTranslateOnly(positions[index]);

        // pxr::SdfPath sdfPath = protoIndices[instanceIndex];
        // auto found = instanceMap.find(sdfPath.GetString());
        //
        // if (found == instanceMap.end())
        // {
        //     printf("not found\n");
        // }
        // assert(found != instanceMap.end());

        // NOTE: OpenUSD uses left to right multiplication. Tranpose to get right
        // to left.
        pxr::GfMatrix4f transposeMatrix = transform.GetTranspose();
        pxr::GfVec4f row0 = transposeMatrix.GetRow(0);
        pxr::GfVec4f row1 = transposeMatrix.GetRow(1);
        pxr::GfVec4f row2 = transposeMatrix.GetRow(2);

        float4 r0 = make_float4(row0[0], row0[1], row0[2], row0[3]);
        float4 r1 = make_float4(row1[0], row1[1], row1[2], row1[3]);
        float4 r2 = make_float4(row2[0], row2[1], row2[2], row2[3]);

        float3x4 affineMatrix(r0, r1, r2);

        affineTransforms[index] = affineMatrix;
        objectIDs[index] = objectID;
    }

    scene->instancesArray.emplace_back(std::move(affineTransforms), std::move(objectIDs));
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

    pxr::Usd_PrimFlagsConjunction filterFlags =
        pxr::UsdPrimIsActive && pxr::UsdPrimIsLoaded && !pxr::UsdPrimIsAbstract;

    pxr::Usd_PrimFlagsPredicate filterPredicate(filterFlags);

    pxr::UsdRenderSettings settings = pxr::UsdRenderSettings::GetStageRenderSettings(stage);
    if (settings)
    {
        printf("found\n");
        float shutterOpen, shutterClose;
        if (settings.GetPrim().GetAttribute(pxr::TfToken("shutter:open")).Get(&shutterOpen))
        {
            printf("Global Render Shutter Open: %f\n", shutterOpen);
        }
    }

    USDTraversalState state;
    pxr::UsdPrim root = stage->GetPseudoRoot();
    TraversePrim(root, state, filterPredicate);

    scene->curves.reserve(state.basisCurves.size());

    std::unordered_map<std::string, int> pathGeomIndexMap;
    std::vector<pxr::UsdPrim> prototypes;

    // Assumes no more than 32 levels of nested instancing. Used to prevent circular references.
    int depth = 0;
    const int maxDepth = 32;
    size_t instanceStart = 0;
    size_t pointInstancerStart = 0;
    while (depth++ < maxDepth && (instanceStart < state.instances.size() ||
                                  pointInstancerStart < state.pointInstancers.size()))
    {
        size_t prototypeStart = prototypes.size();

        // Collect all prototypes
        for (size_t instanceIndex = instanceStart; instanceIndex < state.instances.size();
             instanceIndex++)
        {
            pxr::UsdPrim &instancePrim = state.instances[instanceIndex];
            pxr::UsdPrim proto = instancePrim.GetPrototype();
            std::string pathString = proto.GetPath().GetString();
            auto found = pathGeomIndexMap.find(pathString);
            if (found == pathGeomIndexMap.end())
            {
                int index = static_cast<int>(prototypes.size());
                pathGeomIndexMap.emplace(pathString, index);
                prototypes.push_back(proto);
            }
        }
        for (size_t instanceIndex = pointInstancerStart;
             instanceIndex < state.pointInstancers.size();
             instanceIndex++)
        {
            pxr::UsdGeomPointInstancer &pointInstancer = state.pointInstancers[instanceIndex];
            pxr::SdfPathVector prototypePaths;
            USD_ASSERT(pointInstancer.GetPrototypesRel().GetTargets(&prototypePaths));
            for (pxr::SdfPath &path : prototypePaths)
            {
                std::string pathString = path.GetString();
                auto found = pathGeomIndexMap.find(pathString);
                if (found == pathGeomIndexMap.end())
                {
                    int index = static_cast<int>(prototypes.size());
                    pathGeomIndexMap.emplace(pathString, index);

                    pxr::UsdPrim proto = stage->GetPrimAtPath(path);
                    prototypes.push_back(proto);
                }
            }
        }

        // Prepare for next pass
        instanceStart = state.instances.size();
        pointInstancerStart = state.pointInstancers.size();

        for (size_t protoIndex = prototypeStart; protoIndex < prototypes.size(); protoIndex++)
        {
            pxr::UsdPrim &proto = prototypes[protoIndex];
            TraversePrim(proto, state, filterPredicate);
        }
    }

    for (pxr::UsdGeomPointInstancer &pointInstancer : state.pointInstancers)
    {
        ProcessUSDPointInstancer(pointInstancer, scene, pathGeomIndexMap);
    }

    std::unordered_map<std::string, int> materialMap;
    std::vector<pxr::UsdShadeMaterial> materials;

    Array<int> curveMaterialIndices(state.basisCurves.size());

    // Handle materials
    int curveIndex = 0;
    for (pxr::UsdGeomBasisCurves &curve : state.basisCurves)
    {
        pxr::UsdShadeMaterial material = GetPrimMaterial(curve.GetPrim());
        int materialIndex = AddMaterialToMap(materialMap, materials, material);
        curveMaterialIndices[curveIndex++] = materialIndex;
    }

    for (pxr::UsdGeomMesh &mesh : state.meshes)
    {
        const std::vector<pxr::UsdGeomSubset> subsets =
            pxr::UsdGeomSubset::GetAllGeomSubsets(mesh);

        if (subsets.size())
        {
            printf("subsets not handled yet\n");
        }
        else
        {
            pxr::UsdShadeMaterial material = GetPrimMaterial(mesh.GetPrim());
            int materialIndex = AddMaterialToMap(materialMap, materials, material);
        }
    }

    printf("num materials: %zi\n", materials.size());

    for (pxr::UsdShadeMaterial &material : materials)
    {
        if (pxr::UsdShadeShader shader = material.ComputeSurfaceSource())
        {
            pxr::TfToken token;
            if (shader.GetShaderId(&token) && token == pxr::TfToken("UsdPreviewSurface"))
            {
                if (pxr::UsdShadeInput diffuseInput =
                        shader.GetInput(pxr::TfToken("diffuseColor")))
                {
                    // pxr::UsdShadeShader uvTextureShader = HandleExpectedUVTexture(diffuseInput);
                }
                if (pxr::UsdShadeInput normalInput = shader.GetInput(pxr::TfToken("normal")))
                {
                    pxr::UsdShadeShader uvTextureShader = HandleExpectedUVTexture(normalInput);

                    if (uvTextureShader)
                    {
                    }
                }
            }
            else
            {
                printf("material is not a usdpreviewsurface\n");
                assert(0);
            }
        }
    }

    // Process geometry
    for (pxr::UsdGeomBasisCurves &curve : state.basisCurves)
    {
        ProcessUSDBasisCurve(curve, scene);
    }

    scene->meshes.reserve(state.meshes.size());
    for (pxr::UsdGeomMesh &mesh : state.meshes)
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
