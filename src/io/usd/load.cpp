#include "io/usd/load.h"
#include "io/usd/asset_path_resolve.h"
#include "io/usd/instance_dag_build.h"
#include "pxr/base/gf/vec2f.h"
#include "scene/attributes.h"
#include "scene/scene.h"
#include "util/assert.h"
#include "util/vec2.h"
#include "util/vec3.h"
#include "util/float3x4.h"
#include "util/vec4.h"
#include "util/float4x4.h"

#include <algorithm>
#include <cctype>
#include <cmath>
#include <limits>
#include <pxr/base/gf/frustum.h>
#include <pxr/base/gf/matrix4d.h>
#include <pxr/base/gf/matrix4f.h>
#include <pxr/base/gf/vec4f.h>
#include <pxr/base/vt/types.h>
#include <pxr/base/vt/value.h>
#include <pxr/usd/sdf/path.h>
#include <pxr/usd/sdf/assetPath.h>
#include <pxr/usd/sdf/valueTypeName.h>
#include <pxr/usd/usd/prim.h>
#include <pxr/usd/usd/primRange.h>
#include <pxr/usd/usd/relationship.h>
#include <pxr/usd/usd/stage.h>
#include <pxr/usd/usd/variantSets.h>
#include <pxr/usd/usdGeom/basisCurves.h>
#include <pxr/usd/usdGeom/camera.h>
#include <pxr/usd/usdGeom/curves.h>
#include <pxr/usd/usdGeom/hermiteCurves.h>
#include <pxr/usd/usdGeom/imageable.h>
#include <pxr/usd/usdGeom/mesh.h>
#include <pxr/usd/usdGeom/metrics.h>
#include <pxr/usd/usdGeom/nurbsCurves.h>
#include <pxr/usd/usdGeom/pointInstancer.h>
#include <pxr/usd/usdGeom/primvarsAPI.h>
#include <pxr/usd/usdGeom/scope.h>
#include <pxr/usd/usdGeom/subset.h>
#include <pxr/usd/usdGeom/xformCache.h>
#include <pxr/usd/usdGeom/xform.h>
#include <pxr/usd/usdGeom/xformable.h>
#include <pxr/usd/usdLux/cylinderLight.h>
#include <pxr/usd/usdLux/diskLight.h>
#include <pxr/usd/usdLux/distantLight.h>
#include <pxr/usd/usdLux/domeLight.h>
#include <pxr/usd/usdLux/rectLight.h>
#include <pxr/usd/usdLux/sphereLight.h>
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
            YBI_ASSERT(false); \
        } \
    }

static void CollectUSDCameras(const pxr::UsdPrim &root, std::vector<pxr::UsdGeomCamera> *out)
{
    YBI_ASSERT(out);
    out->clear();

    const pxr::Usd_PrimFlagsPredicate filterPredicate =
        pxr::UsdPrimIsActive && pxr::UsdPrimIsLoaded && !pxr::UsdPrimIsAbstract;
    std::vector<pxr::UsdPrim> stack;
    stack.push_back(root);
    while (!stack.empty())
    {
        const pxr::UsdPrim prim = stack.back();
        stack.pop_back();

        if (prim.IsA<pxr::UsdGeomCamera>())
        {
            out->push_back(pxr::UsdGeomCamera(prim));
        }

        for (const pxr::UsdPrim &child : prim.GetFilteredChildren(filterPredicate))
        {
            stack.push_back(child);
        }
    }
}

static std::string NormalizePurposeName(const std::string &input)
{
    std::string out = input;
    for (char &c : out)
    {
        c = static_cast<char>(std::tolower(static_cast<unsigned char>(c)));
    }
    return out;
}

static std::string PurposeTokenToName(const pxr::TfToken &purpose)
{
    if (purpose == pxr::UsdGeomTokens->default_)
    {
        return "default";
    }
    if (purpose == pxr::UsdGeomTokens->render)
    {
        return "render";
    }
    if (purpose == pxr::UsdGeomTokens->proxy)
    {
        return "proxy";
    }
    if (purpose == pxr::UsdGeomTokens->guide)
    {
        return "guide";
    }
    return "";
}

static bool IsPurposeAllowed(const pxr::UsdPrim &prim, const std::vector<std::string> &allowedPurposes)
{
    pxr::UsdGeomImageable imageable(prim);
    if (!imageable)
    {
        return true;
    }
    const std::string name = PurposeTokenToName(imageable.ComputePurpose());
    for (const std::string &allowed : allowedPurposes)
    {
        if (name == allowed)
        {
            return true;
        }
    }
    return false;
}

static bool IsVisibilityAllowed(const pxr::UsdPrim &prim)
{
    pxr::UsdGeomImageable imageable(prim);
    if (!imageable)
    {
        return true;
    }
    return imageable.ComputeVisibility() != pxr::UsdGeomTokens->invisible;
}

static bool
CreateRuntimeScenesFromBuildSceneDAG(const USDBuildSceneDAG &dag,
                                     ScenePool *scenePool,
                                     std::string *error)
{
    YBI_ASSERT(scenePool);
    if (dag.rootSceneIndex >= dag.scenes.size())
    {
        if (error)
        {
            *error = "build scene dag root index out of range";
        }
        return false;
    }

    scenePool->scenes.clear();
    scenePool->scenes.reserve(dag.scenes.size());
    for (size_t i = 0; i < dag.scenes.size(); i++)
    {
        scenePool->scenes.push_back(std::make_unique<Scene>());
    }

    for (uint32_t sceneIndex = 0; sceneIndex < dag.scenes.size(); sceneIndex++)
    {
        const USDBuildScene &buildScene = dag.scenes[sceneIndex];
        Scene *outScene = scenePool->scenes[sceneIndex].get();
        outScene->meshes.clear();
        outScene->curves.clear();
        outScene->instances.clear();
        outScene->childScenes.clear();
        outScene->meshes.reserve(buildScene.meshes.size());
        outScene->curves.reserve(buildScene.curves.size());
        outScene->instances.reserve(buildScene.instances.size());
        outScene->childScenes.reserve(buildScene.instances.size());

        for (const USDBuildSceneInstance &instance : buildScene.instances)
        {
            if (instance.childSceneIndex >= dag.scenes.size())
            {
                if (error)
                {
                    *error = "child scene index out of range while converting build scenes";
                }
                return false;
            }
            outScene->childScenes.push_back(scenePool->scenes[instance.childSceneIndex].get());
            outScene->instances.emplace_back(
                instance.parentFromLocal,
                static_cast<uint32_t>(outScene->childScenes.size() - 1));
        }
    }

    scenePool->rootSceneIndex = dag.rootSceneIndex;
    return true;
}

static pxr::UsdShadeMaterial GetPrimMaterial(const pxr::UsdPrim &prim)
{
    pxr::UsdShadeMaterialBindingAPI bindingApi(prim);
    pxr::UsdShadeMaterial material = bindingApi.ComputeBoundMaterial(pxr::UsdShadeTokens->full);
    if (!material)
    {
        material = bindingApi.ComputeBoundMaterial(pxr::UsdShadeTokens->preview);
    }
    if (!material)
    {
        material = bindingApi.ComputeBoundMaterial(pxr::UsdShadeTokens->allPurpose);
    }

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
static bool TryGetSingleConnectedShader(const pxr::UsdShadeInput &input,
                                        pxr::UsdShadeShader &shaderOut)
{
    auto sources = input.GetConnectedSources();
    if (sources.size() != 1)
    {
        // TODO: handle multiple or zero sources for a shading input.
        printf("Texture gather: unsupported source count (%zu) for input %s\n",
               sources.size(),
               input.GetBaseName().GetText());
        return false;
    }

    pxr::UsdPrim sourcePrim = sources[0].source.GetPrim();
    if (!sourcePrim.IsA<pxr::UsdShadeShader>())
    {
        // TODO: handle non-shader source nodes (e.g. nodegraphs/other source types).
        printf("Texture gather: unsupported source prim type %s for input %s\n",
               sourcePrim.GetTypeName().GetText(),
               input.GetBaseName().GetText());
        return false;
    }
    shaderOut = pxr::UsdShadeShader(sourcePrim);
    return true;
}

static std::string OutputNameToSwizzle(const pxr::TfToken &sourceName)
{
    const std::string name = sourceName.GetString();
    if (name == "r")
    {
        return "R";
    }
    if (name == "g")
    {
        return "G";
    }
    if (name == "b")
    {
        return "B";
    }
    if (name == "a")
    {
        return "A";
    }
    if (name == "rgb")
    {
        return "RGB";
    }
    if (name == "rgba")
    {
        return "RGBA";
    }
    return "";
}

static TextureWrapMode TextureWrapModeFromToken(const pxr::TfToken &token)
{
    if (token == pxr::TfToken("repeat"))
    {
        return TEXTURE_WRAP_MODE_REPEAT;
    }
    if (token == pxr::TfToken("clamp"))
    {
        return TEXTURE_WRAP_MODE_CLAMP;
    }
    if (token == pxr::TfToken("mirror"))
    {
        return TEXTURE_WRAP_MODE_MIRROR;
    }
    if (token == pxr::TfToken("black"))
    {
        return TEXTURE_WRAP_MODE_BLACK;
    }
    if (token == pxr::TfToken("useMetadata"))
    {
        return TEXTURE_WRAP_MODE_USE_METADATA;
    }
    return TEXTURE_WRAP_MODE_UNKNOWN;
}

static const char *TextureWrapModeToString(TextureWrapMode mode)
{
    switch (mode)
    {
        case TEXTURE_WRAP_MODE_REPEAT:
            return "repeat";
        case TEXTURE_WRAP_MODE_CLAMP:
            return "clamp";
        case TEXTURE_WRAP_MODE_MIRROR:
            return "mirror";
        case TEXTURE_WRAP_MODE_BLACK:
            return "black";
        case TEXTURE_WRAP_MODE_USE_METADATA:
            return "useMetadata";
        default:
            return "unknown";
    }
}

static bool TryGetImageTexturePath(const pxr::UsdShadeInput &input,
                                   std::string &outPath,
                                   std::string *outSwizzle,
                                   TextureWrapMode *outWrapS,
                                   TextureWrapMode *outWrapT);

static ybi::Vec3 ToFloat3(const pxr::GfVec3f &v)
{
    return ybi::Vec3(v[0], v[1], v[2]);
}

static ybi::Vec3 ToFloat3(const pxr::GfVec3d &v)
{
    return ybi::Vec3(static_cast<float>(v[0]),
                            static_cast<float>(v[1]),
                            static_cast<float>(v[2]));
}

static float ClampFinite(float value, float lo, float hi, float fallback)
{
    if (!std::isfinite(value))
    {
        return fallback;
    }
    if (value < lo)
    {
        return lo;
    }
    if (value > hi)
    {
        return hi;
    }
    return value;
}

static ybi::Vec3 ClampFiniteColor(const ybi::Vec3 &value, float lo, float hi)
{
    return ybi::Vec3(ClampFinite(value.x, lo, hi, lo),
                            ClampFinite(value.y, lo, hi, lo),
                            ClampFinite(value.z, lo, hi, lo));
}

static ybi::Vec3 ReadInputColor3f(const pxr::UsdShadeInput &input, const ybi::Vec3 &fallback)
{
    if (!input)
    {
        return fallback;
    }

    pxr::GfVec3f vec3fValue;
    if (input.Get(&vec3fValue))
    {
        return ToFloat3(vec3fValue);
    }

    pxr::GfVec3d vec3dValue;
    if (input.Get(&vec3dValue))
    {
        return ToFloat3(vec3dValue);
    }

    return fallback;
}

static float ReadInputFloat(const pxr::UsdShadeInput &input, float fallback)
{
    if (!input)
    {
        return fallback;
    }

    float floatValue = fallback;
    if (input.Get(&floatValue))
    {
        return floatValue;
    }

    double doubleValue = static_cast<double>(fallback);
    if (input.Get(&doubleValue))
    {
        return static_cast<float>(doubleValue);
    }

    return fallback;
}

static void FinalizePackedMaterial(PackedMaterial *packed)
{
    YBI_ASSERT(packed);
    packed->baseColor = ClampFiniteColor(packed->baseColor, 0.0f, 1.0f);
    packed->emissiveColor = ClampFiniteColor(packed->emissiveColor, 0.0f, 65504.0f);
    packed->roughness = ClampFinite(packed->roughness, 0.0f, 1.0f, 1.0f);
    packed->metallic = ClampFinite(packed->metallic, 0.0f, 1.0f, 0.0f);
    packed->ior = ClampFinite(packed->ior, 1.0f, 4.0f, 1.5f);
    packed->opacity = ClampFinite(packed->opacity, 0.0f, 1.0f, 1.0f);

    packed->flags = 0u;
    if (packed->emissiveColor.x > 0.0f || packed->emissiveColor.y > 0.0f ||
        packed->emissiveColor.z > 0.0f)
    {
        packed->flags |= MATERIAL_FLAG_HAS_EMISSION;
    }
}

static void ReadPreviewSurfaceMaterial(const pxr::UsdShadeShader &shader, MaterialInfo *outInfo)
{
    YBI_ASSERT(outInfo);
    outInfo->packed = {};
    outInfo->packed.baseColor = ReadInputColor3f(
        shader.GetInput(pxr::TfToken("diffuseColor")),
        ybi::Vec3(0.18f, 0.18f, 0.18f));
    outInfo->packed.emissiveColor = ReadInputColor3f(
        shader.GetInput(pxr::TfToken("emissiveColor")),
        ybi::Vec3(0.0f, 0.0f, 0.0f));
    outInfo->packed.roughness =
        ReadInputFloat(shader.GetInput(pxr::TfToken("roughness")), 0.5f);
    outInfo->packed.metallic =
        ReadInputFloat(shader.GetInput(pxr::TfToken("metallic")), 0.0f);
    outInfo->packed.ior = ReadInputFloat(shader.GetInput(pxr::TfToken("ior")), 1.5f);
    outInfo->packed.opacity =
        ReadInputFloat(shader.GetInput(pxr::TfToken("opacity")), 1.0f);

    for (const pxr::UsdShadeInput &input : shader.GetInputs())
    {
        std::string texturePath;
        std::string swizzle;
        TextureWrapMode wrapS = TEXTURE_WRAP_MODE_UNKNOWN;
        TextureWrapMode wrapT = TEXTURE_WRAP_MODE_UNKNOWN;
        if (!TryGetImageTexturePath(input, texturePath, &swizzle, &wrapS, &wrapT))
        {
            continue;
        }

        MaterialTextureInput textureInput = {};
        textureInput.inputName = input.GetBaseName().GetString();
        textureInput.texturePath = texturePath;
        textureInput.swizzle = swizzle;
        textureInput.wrapS = wrapS;
        textureInput.wrapT = wrapT;
        outInfo->textureInputs.push_back(std::move(textureInput));
    }

    FinalizePackedMaterial(&outInfo->packed);
}

template <typename LightSchemaT>
static void ReadCommonLightParams(const LightSchemaT &light, PackedLight *outLight)
{
    YBI_ASSERT(outLight);

    float intensity = 1.0f;
    light.GetIntensityAttr().Get(&intensity);
    float exposure = 0.0f;
    light.GetExposureAttr().Get(&exposure);

    pxr::GfVec3f color(1.0f, 1.0f, 1.0f);
    light.GetColorAttr().Get(&color);

    bool normalize = false;
    light.GetNormalizeAttr().Get(&normalize);

    outLight->color = ClampFiniteColor(ToFloat3(color), 0.0f, 65504.0f);
    outLight->emissionScale =
        std::max(0.0f, intensity) * std::pow(2.0f, std::max(-24.0f, std::min(exposure, 24.0f)));
    const uint32_t preserveFlags = outLight->flags & LIGHT_FLAG_ONE_SIDED;
    outLight->flags = preserveFlags | (normalize ? LIGHT_FLAG_NORMALIZED : 0u);
}

static ybi::Vec3 TransformPoint(const pxr::GfMatrix4d &m, const pxr::GfVec3d &p)
{
    return ToFloat3(m.Transform(p));
}

static ybi::Vec3 TransformDirection(const pxr::GfMatrix4d &m, const pxr::GfVec3d &v)
{
    return ybi::Normalize(ToFloat3(m.TransformDir(v)));
}

static ybi::Vec3 TransformVector(const pxr::GfMatrix4d &m, const pxr::GfVec3d &v)
{
    return ToFloat3(m.TransformDir(v));
}

static ybi::Vec3 SafeNormalizeOrDefault(const ybi::Vec3 &value, const ybi::Vec3 &fallback)
{
    const float len = ybi::Length(value);
    if (len <= 1.0e-8f)
    {
        return fallback;
    }
    return value / len;
}

static void FinalizeLightBasis(PackedLight *light)
{
    YBI_ASSERT(light);

    light->direction =
        SafeNormalizeOrDefault(light->direction, ybi::Vec3(0.0f, 0.0f, -1.0f));
    light->tangent =
        SafeNormalizeOrDefault(light->tangent, ybi::Vec3(1.0f, 0.0f, 0.0f));
    light->bitangent =
        SafeNormalizeOrDefault(light->bitangent, ybi::Cross(light->direction, light->tangent));

    if (std::fabs(ybi::Dot(light->direction, light->tangent)) > 0.999f)
    {
        const ybi::Vec3 up = std::fabs(light->direction.z) < 0.999f
                                   ? ybi::Vec3(0.0f, 0.0f, 1.0f)
                                   : ybi::Vec3(0.0f, 1.0f, 0.0f);
        light->tangent = ybi::Normalize(ybi::Cross(up, light->direction));
    }
    light->bitangent = ybi::Normalize(ybi::Cross(light->direction, light->tangent));
    light->tangent = ybi::Normalize(ybi::Cross(light->bitangent, light->direction));
}

static float ComputeAreaEmissionScale(float intensityScale, float area, bool normalized)
{
    if (normalized)
    {
        return intensityScale / std::max(area, 1.0e-4f);
    }
    return intensityScale;
}

static float ComputeLightLuminance(const ybi::Vec3 &rgb)
{
    return std::max(0.0f, rgb.x * 0.2126f + rgb.y * 0.7152f + rgb.z * 0.0722f);
}

static float ComputeLightSelectionWeight(const PackedLight &light)
{
    const float luminance = ComputeLightLuminance(light.color) * std::max(light.emissionScale, 0.0f);
    float weight = luminance;
    switch (static_cast<LightType>(light.type))
    {
        case LightType::Rect:
        case LightType::Disk:
        case LightType::Sphere:
        case LightType::Cylinder:
            weight *= std::max(light.areaScale, 1.0f);
            break;
        case LightType::Distant:
            weight *= 1024.0f;
            break;
        default:
            break;
    }
    return std::max(weight, 1.0e-4f);
}

static void ReadShadowLinkExcludes(const pxr::UsdPrim &prim, LightInfo *outInfo)
{
    YBI_ASSERT(outInfo);
    outInfo->shadowExcludePaths.clear();

    const pxr::UsdStagePtr stage = prim.GetStage();
    const pxr::UsdRelationship excludesRel =
        prim.GetRelationship(pxr::TfToken("collection:shadowLink:excludes"));
    if (!excludesRel)
    {
        return;
    }

    pxr::SdfPathVector targets;
    excludesRel.GetTargets(&targets);
    outInfo->shadowExcludePaths.reserve(targets.size());
    auto appendPath = [&](const std::string &path) {
        if (path.empty())
        {
            return;
        }
        if (std::find(outInfo->shadowExcludePaths.begin(),
                      outInfo->shadowExcludePaths.end(),
                      path) == outInfo->shadowExcludePaths.end())
        {
            outInfo->shadowExcludePaths.push_back(path);
        }
    };
    for (const pxr::SdfPath &target : targets)
    {
        if (!target.IsPrimPath())
        {
            continue;
        }
        appendPath(target.GetString());
        if (!stage)
        {
            continue;
        }

        const pxr::UsdPrim targetPrim = stage->GetPrimAtPath(target);
        if (!targetPrim)
        {
            continue;
        }
        if (targetPrim.IsInstance())
        {
            const pxr::UsdPrim prototype = targetPrim.GetPrototype();
            if (prototype)
            {
                appendPath(prototype.GetPath().GetString());
            }
        }
        else if (targetPrim.IsInstanceProxy())
        {
            const pxr::UsdPrim primInPrototype = targetPrim.GetPrimInPrototype();
            if (primInPrototype)
            {
                appendPath(primInPrototype.GetPath().GetString());
            }
        }
    }
}

static void CollectUsdLights(const pxr::UsdStageRefPtr &stage,
                             const pxr::UsdTimeCode &timeCode,
                             std::vector<LightInfo> *outLights)
{
    YBI_ASSERT(outLights);
    outLights->clear();

    const pxr::Usd_PrimFlagsPredicate filterPredicate =
        pxr::UsdPrimIsActive && pxr::UsdPrimIsLoaded && !pxr::UsdPrimIsAbstract;
    pxr::UsdGeomXformCache xformCache(timeCode);
    constexpr float kPi = 3.14159265358979323846f;

    for (const pxr::UsdPrim &prim : pxr::UsdPrimRange(stage->GetPseudoRoot(), filterPredicate))
    {
        if (!prim || !IsVisibilityAllowed(prim))
        {
            continue;
        }

        LightInfo info = {};
        const pxr::GfMatrix4d localToWorld = xformCache.GetLocalToWorldTransform(prim);
        info.lightPath = prim.GetPath().GetString();

        if (prim.IsA<pxr::UsdLuxDomeLight>())
        {
            pxr::UsdLuxDomeLight light(prim);
            info.packed.type = static_cast<uint32_t>(LightType::Dome);
            info.packed.position = TransformPoint(localToWorld, pxr::GfVec3d(0.0, 0.0, 0.0));
            info.packed.direction = TransformDirection(localToWorld, pxr::GfVec3d(0.0, 0.0, 1.0));
            ReadCommonLightParams(light, &info.packed);

            pxr::SdfAssetPath textureAsset;
            const pxr::UsdAttribute textureAttr = light.GetTextureFileAttr();
            if (textureAttr && textureAttr.Get(&textureAsset))
            {
                info.texturePath = ResolveUsdAssetPath(textureAttr, textureAsset);
            }
        }
        else if (prim.IsA<pxr::UsdLuxDistantLight>())
        {
            pxr::UsdLuxDistantLight light(prim);
            info.packed.type = static_cast<uint32_t>(LightType::Distant);
            info.packed.position = TransformPoint(localToWorld, pxr::GfVec3d(0.0, 0.0, 0.0));
            info.packed.direction = TransformDirection(localToWorld, pxr::GfVec3d(0.0, 0.0, -1.0));
            ReadCommonLightParams(light, &info.packed);
        }
        else if (prim.IsA<pxr::UsdLuxRectLight>())
        {
            pxr::UsdLuxRectLight light(prim);
            info.packed.type = static_cast<uint32_t>(LightType::Rect);
            info.packed.flags |= LIGHT_FLAG_ONE_SIDED;
            info.packed.position = TransformPoint(localToWorld, pxr::GfVec3d(0.0, 0.0, 0.0));
            const ybi::Vec3 axisX = TransformVector(localToWorld, pxr::GfVec3d(1.0, 0.0, 0.0));
            const ybi::Vec3 axisY = TransformVector(localToWorld, pxr::GfVec3d(0.0, 1.0, 0.0));
            const ybi::Vec3 axisZ = TransformVector(localToWorld, pxr::GfVec3d(0.0, 0.0, -1.0));
            info.packed.direction = SafeNormalizeOrDefault(axisZ, ybi::Vec3(0.0f, 0.0f, -1.0f));
            info.packed.tangent = SafeNormalizeOrDefault(axisX, ybi::Vec3(1.0f, 0.0f, 0.0f));
            info.packed.bitangent =
                SafeNormalizeOrDefault(axisY, ybi::Vec3(0.0f, 1.0f, 0.0f));
            ReadCommonLightParams(light, &info.packed);
            float localWidth = 0.0f;
            float localHeight = 0.0f;
            light.GetWidthAttr().Get(&localWidth);
            light.GetHeightAttr().Get(&localHeight);
            info.packed.width = std::max(localWidth, 0.0f) * ybi::Length(axisX);
            info.packed.height = std::max(localHeight, 0.0f) * ybi::Length(axisY);
            info.packed.areaScale = std::max(info.packed.width * info.packed.height, 1.0e-4f);
            const bool normalized = (info.packed.flags & LIGHT_FLAG_NORMALIZED) != 0u;
            info.packed.emissionScale =
                ComputeAreaEmissionScale(info.packed.emissionScale, info.packed.areaScale, normalized);
        }
        else if (prim.IsA<pxr::UsdLuxDiskLight>())
        {
            pxr::UsdLuxDiskLight light(prim);
            info.packed.type = static_cast<uint32_t>(LightType::Disk);
            info.packed.flags |= LIGHT_FLAG_ONE_SIDED;
            info.packed.position = TransformPoint(localToWorld, pxr::GfVec3d(0.0, 0.0, 0.0));
            const ybi::Vec3 axisX = TransformVector(localToWorld, pxr::GfVec3d(1.0, 0.0, 0.0));
            const ybi::Vec3 axisY = TransformVector(localToWorld, pxr::GfVec3d(0.0, 1.0, 0.0));
            const ybi::Vec3 axisZ = TransformVector(localToWorld, pxr::GfVec3d(0.0, 0.0, -1.0));
            info.packed.direction = SafeNormalizeOrDefault(axisZ, ybi::Vec3(0.0f, 0.0f, -1.0f));
            info.packed.tangent = SafeNormalizeOrDefault(axisX, ybi::Vec3(1.0f, 0.0f, 0.0f));
            info.packed.bitangent =
                SafeNormalizeOrDefault(axisY, ybi::Vec3(0.0f, 1.0f, 0.0f));
            ReadCommonLightParams(light, &info.packed);
            float localRadius = 0.0f;
            light.GetRadiusAttr().Get(&localRadius);
            const float worldScaleX = ybi::Length(axisX);
            const float worldScaleY = ybi::Length(axisY);
            info.packed.radius = std::max(localRadius, 0.0f) * 0.5f * (worldScaleX + worldScaleY);
            info.packed.areaScale =
                std::max(kPi * std::max(localRadius, 0.0f) * worldScaleX * std::max(localRadius, 0.0f) *
                             worldScaleY,
                         1.0e-4f);
            const bool normalized = (info.packed.flags & LIGHT_FLAG_NORMALIZED) != 0u;
            info.packed.emissionScale =
                ComputeAreaEmissionScale(info.packed.emissionScale, info.packed.areaScale, normalized);
        }
        else if (prim.IsA<pxr::UsdLuxSphereLight>())
        {
            pxr::UsdLuxSphereLight light(prim);
            info.packed.type = static_cast<uint32_t>(LightType::Sphere);
            info.packed.position = TransformPoint(localToWorld, pxr::GfVec3d(0.0, 0.0, 0.0));
            info.packed.direction = TransformDirection(localToWorld, pxr::GfVec3d(0.0, 0.0, -1.0));
            info.packed.tangent = TransformDirection(localToWorld, pxr::GfVec3d(1.0, 0.0, 0.0));
            info.packed.bitangent = TransformDirection(localToWorld, pxr::GfVec3d(0.0, 1.0, 0.0));
            ReadCommonLightParams(light, &info.packed);
            float localRadius = 0.0f;
            light.GetRadiusAttr().Get(&localRadius);
            const float worldScaleX = ybi::Length(TransformVector(localToWorld, pxr::GfVec3d(1.0, 0.0, 0.0)));
            const float worldScaleY = ybi::Length(TransformVector(localToWorld, pxr::GfVec3d(0.0, 1.0, 0.0)));
            const float worldScaleZ = ybi::Length(TransformVector(localToWorld, pxr::GfVec3d(0.0, 0.0, 1.0)));
            const float worldScale = (worldScaleX + worldScaleY + worldScaleZ) / 3.0f;
            info.packed.radius = std::max(localRadius, 0.0f) * worldScale;
            info.packed.areaScale =
                std::max(4.0f * kPi * info.packed.radius * info.packed.radius, 1.0e-4f);
            const bool normalized = (info.packed.flags & LIGHT_FLAG_NORMALIZED) != 0u;
            info.packed.emissionScale =
                ComputeAreaEmissionScale(info.packed.emissionScale, info.packed.areaScale, normalized);
        }
        else if (prim.IsA<pxr::UsdLuxCylinderLight>())
        {
            pxr::UsdLuxCylinderLight light(prim);
            info.packed.type = static_cast<uint32_t>(LightType::Cylinder);
            info.packed.position = TransformPoint(localToWorld, pxr::GfVec3d(0.0, 0.0, 0.0));
            const ybi::Vec3 axisX = TransformVector(localToWorld, pxr::GfVec3d(1.0, 0.0, 0.0));
            const ybi::Vec3 axisY = TransformVector(localToWorld, pxr::GfVec3d(0.0, 1.0, 0.0));
            const ybi::Vec3 axisZ = TransformVector(localToWorld, pxr::GfVec3d(0.0, 0.0, 1.0));
            info.packed.direction = SafeNormalizeOrDefault(axisZ, ybi::Vec3(0.0f, 0.0f, 1.0f));
            info.packed.tangent = SafeNormalizeOrDefault(axisX, ybi::Vec3(1.0f, 0.0f, 0.0f));
            info.packed.bitangent =
                SafeNormalizeOrDefault(axisY, ybi::Vec3(0.0f, 1.0f, 0.0f));
            ReadCommonLightParams(light, &info.packed);
            float localRadius = 0.0f;
            float localLength = 0.0f;
            light.GetRadiusAttr().Get(&localRadius);
            light.GetLengthAttr().Get(&localLength);
            const float worldRadiusScale =
                0.5f * (ybi::Length(axisX) + ybi::Length(axisY));
            info.packed.radius = std::max(localRadius, 0.0f) * worldRadiusScale;
            info.packed.length = std::max(localLength, 0.0f) * ybi::Length(axisZ);
            info.packed.areaScale =
                std::max(2.0f * kPi * info.packed.radius * info.packed.length, 1.0e-4f);
            const bool normalized = (info.packed.flags & LIGHT_FLAG_NORMALIZED) != 0u;
            info.packed.emissionScale =
                ComputeAreaEmissionScale(info.packed.emissionScale, info.packed.areaScale, normalized);
        }
        else
        {
            continue;
        }

        ReadShadowLinkExcludes(prim, &info);
        FinalizeLightBasis(&info.packed);
        info.packed.width = std::max(info.packed.width, 0.0f);
        info.packed.height = std::max(info.packed.height, 0.0f);
        info.packed.radius = std::max(info.packed.radius, 0.0f);
        info.packed.length = std::max(info.packed.length, 0.0f);
        info.packed.areaScale = std::max(info.packed.areaScale, 0.0f);
        info.packed.selectionWeight = ComputeLightSelectionWeight(info.packed);
        if (info.packed.emissionScale <= 0.0f ||
            (info.packed.color.x <= 0.0f && info.packed.color.y <= 0.0f &&
             info.packed.color.z <= 0.0f))
        {
            continue;
        }
        outLights->push_back(std::move(info));
    }
}

static bool TryGetImageTexturePath(const pxr::UsdShadeInput &input,
                                   std::string &outPath,
                                   std::string *outSwizzle,
                                   TextureWrapMode *outWrapS,
                                   TextureWrapMode *outWrapT)
{
    outPath.clear();
    if (outSwizzle)
    {
        outSwizzle->clear();
    }
    if (outWrapS)
    {
        *outWrapS = TEXTURE_WRAP_MODE_UNKNOWN;
    }
    if (outWrapT)
    {
        *outWrapT = TEXTURE_WRAP_MODE_UNKNOWN;
    }

    pxr::UsdShadeShader sourceShader;
    if (!TryGetSingleConnectedShader(input, sourceShader))
    {
        return false;
    }

    pxr::TfToken shaderId;
    sourceShader.GetShaderId(&shaderId);
    if (shaderId != pxr::TfToken("UsdUVTexture"))
    {
        // TODO: support non-UsdUVTexture image producers (MaterialX image nodes, etc.).
        printf("Texture gather: unsupported shader node type %s at %s\n",
               shaderId.GetText(),
               sourceShader.GetPath().GetText());
        return false;
    }

    pxr::UsdShadeInput fileInput = sourceShader.GetInput(pxr::TfToken("file"));
    if (!fileInput)
    {
        // TODO: support alternate file-bearing texture nodes/ports.
        printf("Texture gather: UsdUVTexture missing 'file' input at %s\n",
               sourceShader.GetPath().GetText());
        return false;
    }

    pxr::SdfAssetPath assetPath;
    if (!fileInput.Get(&assetPath))
    {
        // TODO: support connected/indirected file inputs.
        printf("Texture gather: could not read 'file' value at %s\n",
               sourceShader.GetPath().GetText());
        return false;
    }

    outPath = ResolveUsdAssetPath(fileInput.GetAttr(), assetPath);
    if (outPath.empty())
    {
        printf("Texture gather: empty texture file path at %s\n", sourceShader.GetPath().GetText());
        return false;
    }

    if (outSwizzle)
    {
        auto sources = input.GetConnectedSources();
        if (sources.size() == 1)
        {
            outSwizzle->assign(OutputNameToSwizzle(sources[0].sourceName));
        }
    }
    if (outWrapS)
    {
        pxr::TfToken wrapS;
        pxr::UsdShadeInput wrapSInput = sourceShader.GetInput(pxr::TfToken("wrapS"));
        if (wrapSInput && wrapSInput.Get(&wrapS))
        {
            *outWrapS = TextureWrapModeFromToken(wrapS);
        }
    }
    if (outWrapT)
    {
        pxr::TfToken wrapT;
        pxr::UsdShadeInput wrapTInput = sourceShader.GetInput(pxr::TfToken("wrapT"));
        if (wrapTInput && wrapTInput.Get(&wrapT))
        {
            *outWrapT = TextureWrapModeFromToken(wrapT);
        }
    }

    return true;
}

static void ReadNtcMaterialInfo(const pxr::UsdShadeMaterial &material,
                                MaterialInfo *outInfo)
{
    YBI_ASSERT(outInfo);
    outInfo->ntcDiffuseFile.clear();
    outInfo->ntcDiffuseTextureName.clear();

    const pxr::UsdPrim materialPrim = material.GetPrim();
    if (!materialPrim)
    {
        return;
    }

    const pxr::UsdAttribute ntcFileAttr =
        materialPrim.GetAttribute(pxr::TfToken("ybi:ntc:diffuseFile"));
    if (ntcFileAttr)
    {
        pxr::SdfAssetPath assetPath;
        if (ntcFileAttr.Get(&assetPath))
        {
            outInfo->ntcDiffuseFile = ResolveUsdAssetPath(ntcFileAttr, assetPath);
        }
        else
        {
            // Back-compat if authored as plain string.
            ntcFileAttr.Get(&outInfo->ntcDiffuseFile);
        }
    }

    const pxr::UsdAttribute ntcTextureNameAttr =
        materialPrim.GetAttribute(pxr::TfToken("ybi:ntc:diffuseTextureName"));
    if (ntcTextureNameAttr)
    {
        ntcTextureNameAttr.Get(&outInfo->ntcDiffuseTextureName);
    }
    if (!outInfo->ntcDiffuseFile.empty() && outInfo->ntcDiffuseTextureName.empty())
    {
        outInfo->ntcDiffuseTextureName = "diffuseColor";
    }
}

static void AppendUniqueTexturePath(std::vector<std::string> &paths, const std::string &path)
{
    if (std::find(paths.begin(), paths.end(), path) == paths.end())
    {
        paths.push_back(path);
    }
}

static int GetMaterialIndex(const std::unordered_map<std::string, int> &materialMap,
                            const pxr::UsdShadeMaterial &material)
{
    if (!material)
    {
        return -1;
    }
    const std::string materialPath = material.GetPath().GetString();
    auto it = materialMap.find(materialPath);
    if (it == materialMap.end())
    {
        return -1;
    }
    return it->second;
}

static bool BuildTriangulatedMeshSTAttribute(const pxr::UsdGeomMesh &mesh,
                                             const pxr::VtIntArray &faceCounts,
                                             const pxr::VtIntArray &faceIndices,
                                             int numTriangles,
                                             Attribute *outAttribute)
{
    YBI_ASSERT(outAttribute);

    pxr::UsdGeomPrimvarsAPI primvars(mesh);
    pxr::UsdGeomPrimvar stPrimvar = primvars.GetPrimvar(pxr::TfToken("st"));
    if (!stPrimvar)
    {
        return false;
    }

    pxr::VtVec2fArray stValues;
    if (!stPrimvar.Get(&stValues, 0.0))
    {
        return false;
    }
    if (stValues.empty())
    {
        return false;
    }

    pxr::VtIntArray stIndices;
    const bool hasExplicitIndices = stPrimvar.IsIndexed() && stPrimvar.GetIndices(&stIndices, 0.0);
    const pxr::TfToken interpolation = stPrimvar.GetInterpolation();

    Array<Vec2> texcoords(stValues);
    Array<int> triTexcoordIndices(static_cast<size_t>(numTriangles) * 3u);

    int faceIndexOffset = 0;
    int triIndexOffset = 0;
    for (int faceCount : faceCounts)
    {
        int faceTc[4] = {-1, -1, -1, -1};
        if (faceCount < 3 || faceCount > 4)
        {
            return false;
        }

        for (int corner = 0; corner < faceCount; corner++)
        {
            const int cornerIndex = faceIndexOffset + corner;
            int tcIndex = -1;
            if (hasExplicitIndices)
            {
                if (cornerIndex < 0 || cornerIndex >= static_cast<int>(stIndices.size()))
                {
                    return false;
                }
                tcIndex = stIndices[cornerIndex];
            }
            else if (interpolation == pxr::UsdGeomTokens->faceVarying)
            {
                tcIndex = cornerIndex;
            }
            else
            {
                if (cornerIndex < 0 || cornerIndex >= static_cast<int>(faceIndices.size()))
                {
                    return false;
                }
                tcIndex = faceIndices[cornerIndex];
            }

            if (tcIndex < 0 || tcIndex >= static_cast<int>(texcoords.size()))
            {
                return false;
            }
            faceTc[corner] = tcIndex;
        }

        if (faceCount == 3)
        {
            triTexcoordIndices[triIndexOffset++] = faceTc[0];
            triTexcoordIndices[triIndexOffset++] = faceTc[1];
            triTexcoordIndices[triIndexOffset++] = faceTc[2];
        }
        else
        {
            triTexcoordIndices[triIndexOffset++] = faceTc[0];
            triTexcoordIndices[triIndexOffset++] = faceTc[1];
            triTexcoordIndices[triIndexOffset++] = faceTc[2];
            triTexcoordIndices[triIndexOffset++] = faceTc[0];
            triTexcoordIndices[triIndexOffset++] = faceTc[2];
            triTexcoordIndices[triIndexOffset++] = faceTc[3];
        }
        faceIndexOffset += faceCount;
    }

    Array<uint8_t> dataBytes(sizeof(Vec2) * texcoords.size());
    if (texcoords.size() > 0)
    {
        memcpy(dataBytes.data(), texcoords.data(), sizeof(Vec2) * texcoords.size());
    }
    *outAttribute = Attribute(std::move(dataBytes),
                              std::move(triTexcoordIndices),
                              AttributeType::Float2,
                              PrimvarInterpolation::FaceVarying,
                              "st");
    return true;
}

static bool BuildSubdivisionMeshSTAttribute(const pxr::UsdGeomMesh &mesh,
                                            const pxr::VtIntArray &faceIndices,
                                            Attribute *outAttribute)
{
    YBI_ASSERT(outAttribute);

    pxr::UsdGeomPrimvarsAPI primvars(mesh);
    pxr::UsdGeomPrimvar stPrimvar = primvars.GetPrimvar(pxr::TfToken("st"));
    if (!stPrimvar)
    {
        return false;
    }

    pxr::VtVec2fArray stValues;
    if (!stPrimvar.Get(&stValues, 0.0))
    {
        return false;
    }
    if (stValues.empty())
    {
        return false;
    }

    pxr::VtIntArray stIndices;
    const bool hasExplicitIndices = stPrimvar.IsIndexed() && stPrimvar.GetIndices(&stIndices, 0.0);
    const pxr::TfToken interpolation = stPrimvar.GetInterpolation();

    Array<Vec2> texcoords(stValues);
    Array<int> cornerTexcoordIndices(faceIndices.size());
    for (size_t corner = 0; corner < faceIndices.size(); ++corner)
    {
        int tcIndex = -1;
        if (hasExplicitIndices)
        {
            if (corner >= stIndices.size())
            {
                return false;
            }
            tcIndex = stIndices[corner];
        }
        else if (interpolation == pxr::UsdGeomTokens->faceVarying)
        {
            tcIndex = int(corner);
        }
        else if (interpolation == pxr::UsdGeomTokens->vertex ||
                 interpolation == pxr::UsdGeomTokens->varying)
        {
            tcIndex = faceIndices[corner];
        }
        else
        {
            return false;
        }

        if (tcIndex < 0 || tcIndex >= int(texcoords.size()))
        {
            return false;
        }
        cornerTexcoordIndices[corner] = tcIndex;
    }

    Array<uint8_t> dataBytes(sizeof(Vec2) * texcoords.size());
    if (texcoords.size() > 0)
    {
        memcpy(dataBytes.data(), texcoords.data(), sizeof(Vec2) * texcoords.size());
    }
    *outAttribute = Attribute(std::move(dataBytes),
                              std::move(cornerTexcoordIndices),
                              AttributeType::Float2,
                              PrimvarInterpolation::FaceVarying,
                              "st");
    return true;
}

static PrimvarInterpolation ConvertPrimvarInterpolation(pxr::TfToken &token)
{
    if (token == "constant")
    {
        return PrimvarInterpolation::Constant;
    }
    if (token == "uniform")
    {
        return PrimvarInterpolation::Uniform;
    }
    if (token == "vertex")
    {
        return PrimvarInterpolation::Vertex;
    }
    if (token == "varying")
    {
        return PrimvarInterpolation::Varying;
    }
    if (token == "faceVarying")
    {
        return PrimvarInterpolation::FaceVarying;
    }
    return PrimvarInterpolation::Unknown;
}

static AttributeType ConvertPrimvarTypeName(pxr::SdfValueTypeName type)
{
    if (type == pxr::SdfValueTypeNames->FloatArray)
    {
        return AttributeType::Float;
    }
    if (type == pxr::SdfValueTypeNames->Double)
    {
        return AttributeType::Float;
    }
    if (type == pxr::SdfValueTypeNames->UCharArray)
    {
        return AttributeType::Int8;
    }
    if (type == pxr::SdfValueTypeNames->IntArray || type == pxr::SdfValueTypeNames->Int)
    {
        return AttributeType::Int;
    }
    if (type == pxr::SdfValueTypeNames->Float2Array)
    {
        return AttributeType::Float2;
    }
    if (type == pxr::SdfValueTypeNames->TexCoord2dArray)
    {
        return AttributeType::Float2;
    }
    if (type == pxr::SdfValueTypeNames->TexCoord2fArray)
    {
        return AttributeType::Float2;
    }
    if (type == pxr::SdfValueTypeNames->TexCoord2hArray)
    {
        return AttributeType::Float2;
    }
    if (type == pxr::SdfValueTypeNames->TexCoord3dArray)
    {
        return AttributeType::Float3;
    }
    if (type == pxr::SdfValueTypeNames->TexCoord3fArray)
    {
        return AttributeType::Float3;
    }
    if (type == pxr::SdfValueTypeNames->TexCoord3hArray)
    {
        return AttributeType::Float3;
    }
    if (type == pxr::SdfValueTypeNames->Float3Array)
    {
        return AttributeType::Float3;
    }
    if (type == pxr::SdfValueTypeNames->Point3fArray)
    {
        return AttributeType::Float3;
    }
    if (type == pxr::SdfValueTypeNames->Point3dArray)
    {
        return AttributeType::Float3;
    }
    if (type == pxr::SdfValueTypeNames->Point3hArray)
    {
        return AttributeType::Float3;
    }
    if (type == pxr::SdfValueTypeNames->Normal3fArray)
    {
        return AttributeType::Float3;
    }
    if (type == pxr::SdfValueTypeNames->Normal3dArray)
    {
        return AttributeType::Float3;
    }
    if (type == pxr::SdfValueTypeNames->Normal3hArray)
    {
        return AttributeType::Float3;
    }
    if (type == pxr::SdfValueTypeNames->Vector3fArray)
    {
        return AttributeType::Float3;
    }
    if (type == pxr::SdfValueTypeNames->Vector3hArray)
    {
        return AttributeType::Float3;
    }
    if (type == pxr::SdfValueTypeNames->Vector3dArray)
    {
        return AttributeType::Float3;
    }
    if (type == pxr::SdfValueTypeNames->Color3fArray)
    {
        return AttributeType::Float3;
    }
    if (type == pxr::SdfValueTypeNames->Color3hArray)
    {
        return AttributeType::Float3;
    }
    if (type == pxr::SdfValueTypeNames->Color3dArray)
    {
        return AttributeType::Float3;
    }
    if (type == pxr::SdfValueTypeNames->Color4fArray)
    {
        return AttributeType::Float4;
    }
    if (type == pxr::SdfValueTypeNames->Color4hArray)
    {
        return AttributeType::Float4;
    }
    if (type == pxr::SdfValueTypeNames->Color4dArray)
    {
        return AttributeType::Float4;
    }
    if (type == pxr::SdfValueTypeNames->BoolArray || type == pxr::SdfValueTypeNames->Bool)
    {
        return AttributeType::Bool;
    }
    if (type == pxr::SdfValueTypeNames->QuatfArray)
    {
        return AttributeType::Quaternion;
    }
    if (type == pxr::SdfValueTypeNames->QuatdArray)
    {
        return AttributeType::Quaternion;
    }
    if (type == pxr::SdfValueTypeNames->QuathArray)
    {
        return AttributeType::Quaternion;
    }
    if (type == pxr::SdfValueTypeNames->Token)
    {
        return AttributeType::Unknown;
    }
    YBI_LOGFATAL("Unsupported primvar type: %s", type.GetAsToken().GetText());
    return AttributeType::Unknown;
}

static MemoryView<uint8_t> ConvertPrimvarValues(Scene *scene, const pxr::VtValue &values)
{
    if (values.IsHolding<pxr::VtArray<pxr::GfVec2f>>())
    {
        pxr::VtArray<pxr::GfVec2f> array = values.Get<pxr::VtArray<pxr::GfVec2f>>();
        // MemoryView<Vec2> view = scene->arena.PushArray<Vec2>(array.size());
        // memcpy(view.data(), array.data(), array.size() * sizeof(Vec2));
        // return view.CastToBytes();
    }
    else
    {
        YBI_LOGFATAL("Unsupported primvar values type: %s", values.GetTypeName().c_str());
    }

    return MemoryView<uint8_t>();
}

static void ProcessPrimvars(pxr::UsdPrim prim, pxr::UsdTimeCode timeCode, Scene *scene)
{
    pxr::UsdGeomPrimvarsAPI primvarsAPI(prim);
    for (const pxr::UsdGeomPrimvar &primVar : primvarsAPI.GetPrimvars())
    {
        pxr::VtValue values;
        if (primVar.HasValue())
        {
            bool flattened = primVar.IsIndexed();
            USD_ASSERT(primVar.ComputeFlattened(&values, timeCode));
            pxr::TfToken interpolationToken = primVar.GetInterpolation();
            PrimvarInterpolation interpolation = ConvertPrimvarInterpolation(interpolationToken);

            pxr::SdfValueTypeName typeName = primVar.GetTypeName();
            AttributeType attrType = ConvertPrimvarTypeName(typeName);

            pxr::TfToken name = pxr::UsdGeomPrimvar::StripPrimvarsName(primVar.GetName());

            if (name == "st")
            {
                MemoryView<uint8_t> data = ConvertPrimvarValues(scene, values);
                scene->attributes.EmplaceBack(data, attrType, interpolation);
            }
            // printf("name: %s, typeName: %s, interp: %s\n",
            //        name.GetText(),
            //        typeName.GetAsToken().GetText(),
            //        interpolationToken.GetText());
        }
    }
}

static void
ProcessCatmullClarkMesh(pxr::UsdGeomMesh &mesh,
                        Scene *scene,
                        int materialIndex,
                        const Float3x4 &parentFromLocal,
                        pxr::UsdTimeCode timeCode = 0.0)
{
    pxr::VtVec3fArray positions;
    pxr::VtIntArray faceIndices;
    pxr::VtIntArray faceCounts;
    pxr::VtVec3fArray normals;

    pxr::TfToken interpolationBoundary;
    pxr::TfToken fvLinearInterpolation;
    pxr::TfToken triangleSubdivisionRule;

    pxr::VtIntArray cornerIndices;
    pxr::VtFloatArray cornerSharpnesses;
    pxr::VtIntArray creaseIndices;
    pxr::VtIntArray creaseLengths;
    pxr::VtFloatArray creaseSharpnesses;
    pxr::VtIntArray holeIndices;

    USD_ASSERT(mesh.GetPointsAttr().Get(&positions, timeCode));

    USD_ASSERT(mesh.GetFaceVertexIndicesAttr().Get(&faceIndices, timeCode));
    USD_ASSERT(mesh.GetFaceVertexCountsAttr().Get(&faceCounts, timeCode));

    if (mesh.GetNormalsAttr().Get(&normals, timeCode))
    {
        pxr::TfToken normalsInterpToken = mesh.GetNormalsInterpolation();
        for (auto n : normals)
        {
            // printf("help: %f %f %f\n", n[0], n[1], n[2]);
        }
    }

    USD_ASSERT(mesh.GetTriangleSubdivisionRuleAttr().Get(&triangleSubdivisionRule, timeCode));
    USD_ASSERT(mesh.GetInterpolateBoundaryAttr().Get(&interpolationBoundary, timeCode));
    USD_ASSERT(mesh.GetFaceVaryingLinearInterpolationAttr().Get(&fvLinearInterpolation, timeCode));

    mesh.GetCornerIndicesAttr().Get(&cornerIndices, timeCode);
    mesh.GetCornerSharpnessesAttr().Get(&cornerSharpnesses, timeCode);
    mesh.GetCreaseIndicesAttr().Get(&creaseIndices, timeCode);
    mesh.GetCreaseLengthsAttr().Get(&creaseLengths, timeCode);
    mesh.GetCreaseSharpnessesAttr().Get(&creaseSharpnesses, timeCode);
    mesh.GetHoleIndicesAttr().Get(&holeIndices, timeCode);

    pxr::UsdPrim prim = mesh.GetPrim();

    Array<Vec3> positionsArray(positions);
    Array<int> faceIndicesArray(faceIndices);
    Array<int> faceCountsArray(faceCounts);

    BoundaryInterpolation interpolation = BOUNDARY_INTERPOLATION_EDGE_AND_CORNER;
    if (interpolationBoundary == "none")
    {
        interpolation = BOUNDARY_INTERPOLATION_NONE;
    }
    else if (interpolationBoundary == "edgeOnly")
    {
        interpolation = BOUNDARY_INTERPOLATION_EDGE;
    }
    else if (interpolationBoundary == "edgeAndCorner")
    {
        interpolation = BOUNDARY_INTERPOLATION_EDGE_AND_CORNER;
    }

    FVarLinearInterpolation fvarLinear = FVAR_LINEAR_CORNERS_ONLY;
    if (fvLinearInterpolation == "none")
    {
        fvarLinear = FVAR_LINEAR_NONE;
    }
    else if (fvLinearInterpolation == "cornersOnly")
    {
        fvarLinear = FVAR_LINEAR_CORNERS_ONLY;
    }
    else if (fvLinearInterpolation == "cornersPlus1")
    {
        fvarLinear = FVAR_LINEAR_CORNERS_PLUS1;
    }
    else if (fvLinearInterpolation == "cornersPlus2")
    {
        fvarLinear = FVAR_LINEAR_CORNERS_PLUS2;
    }
    else if (fvLinearInterpolation == "boundaries")
    {
        fvarLinear = FVAR_LINEAR_BOUNDARIES;
    }
    else if (fvLinearInterpolation == "all")
    {
        fvarLinear = FVAR_LINEAR_ALL;
    }

    TriangleSubdivisionRule triangleSubdivision = TRIANGLE_SUBDIVISION_CATMULL_CLARK;
    if (triangleSubdivisionRule == "smooth")
    {
        triangleSubdivision = TRIANGLE_SUBDIVISION_SMOOTH;
    }

    Array<int> cornerIndicesArray(cornerIndices);
    Array<float> cornerSharpnessesArray(cornerSharpnesses);
    Array<int> creaseIndicesArray(creaseIndices);
    Array<int> creaseLengthsArray(creaseLengths);
    Array<float> creaseSharpnessesArray(creaseSharpnesses);

    scene->subdivisionMeshes.emplace_back(std::move(positionsArray),
                                          std::move(faceIndicesArray),
                                          std::move(faceCountsArray),
                                          std::move(cornerIndicesArray),
                                          std::move(cornerSharpnessesArray),
                                          std::move(creaseIndicesArray),
                                          std::move(creaseLengthsArray),
                                          std::move(creaseSharpnessesArray),
                                          std::move(holeIndices),
                                          interpolation,
                                          fvarLinear,
                                          triangleSubdivision);
    SubdivisionMesh &subdivMesh = scene->subdivisionMeshes.back();
    subdivMesh.primPath = prim.GetPath().GetString();
    subdivMesh.materialIndex = materialIndex;
    subdivMesh.parentFromLocal = parentFromLocal;
    Attribute stAttr;
    if (BuildSubdivisionMeshSTAttribute(mesh, faceIndices, &stAttr))
    {
        subdivMesh.attributes.push_back(std::move(stAttr));
    }
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

    Array<Vec3> positions(points);
    Array<float> curveWidths(widths);
    Array<int> curveOffsets(numCurves);
    int curveOffsetIndex = 0;

    for (int curveVertexCount : curveVertexCounts)
    {
        curveOffsets[curveOffsetIndex++] = totalNumVertices;
        totalNumVertices += curveVertexCount;
    }

    YBI_ASSERT(totalNumVertices == points.size());
    scene->curves.emplace_back(
        std::move(positions), std::move(curveWidths), std::move(curveOffsets));

    int curveFlags = 0;

    if (typeToken == "cubic")
    {
        curveFlags |= CurveFlags::CURVE_FLAGS_CUBIC;
    }
    else
    {
        YBI_ASSERT(typeToken == "linear");
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

void LoadUSDScene(ScenePool *scenePool, const std::string &filePath)
{
    const USDLoadOptions options = {};
    LoadUSDScene(scenePool, filePath, options);
}

void LoadUSDScene(ScenePool *scenePool, const std::string &filePath, const USDLoadOptions &options)
{
    YBI_ASSERT(scenePool);
    scenePool->camera = Camera();
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

    pxr::UsdPrim root = stage->GetPseudoRoot();

    std::vector<std::string> allowedPurposes;
    allowedPurposes.reserve(options.purposes.size());
    for (const std::string &purpose : options.purposes)
    {
        const std::string normalized = NormalizePurposeName(purpose);
        if (normalized == "default" || normalized == "render" || normalized == "proxy" ||
            normalized == "guide")
        {
            if (std::find(allowedPurposes.begin(), allowedPurposes.end(), normalized) ==
                allowedPurposes.end())
            {
                allowedPurposes.push_back(normalized);
            }
        }
    }
    if (allowedPurposes.empty())
    {
        allowedPurposes = {"default", "render"};
    }

    USDBuildSceneDAG buildSceneDAG = {};
    std::string buildSceneError;
    if (!BuildInstanceDAGFromUSD(stage, &buildSceneDAG, &buildSceneError))
    {
        printf("build instance dag failed: %s\n", buildSceneError.c_str());
        return;
    }
    std::vector<pxr::UsdGeomCamera> cameras;
    CollectUSDCameras(root, &cameras);

    std::string createScenesError;
    if (!CreateRuntimeScenesFromBuildSceneDAG(buildSceneDAG, scenePool, &createScenesError))
    {
        printf("build scene creation failed: %s\n", createScenesError.c_str());
        return;
    }

    if (cameras.size())
    {
        pxr::UsdGeomCamera &camera = cameras[0];
        const pxr::UsdTimeCode timeCode(0.0);
        pxr::GfCamera gfCam = camera.GetCamera(timeCode);
        pxr::GfFrustum frustum = gfCam.GetFrustum();
        pxr::GfMatrix4d viewM = frustum.ComputeViewMatrix().GetTranspose();
        pxr::GfMatrix4d projM = frustum.ComputeProjectionMatrix().GetTranspose();
        const pxr::UsdPrim cameraPrim = camera.GetPrim();
        pxr::UsdGeomXformCache xformCache(timeCode);
        const pxr::GfMatrix4d cameraLocalToWorld = xformCache.GetLocalToWorldTransform(cameraPrim);
        const pxr::GfVec3d cameraP = cameraLocalToWorld.Transform(pxr::GfVec3d(0.0, 0.0, 0.0));
        const pxr::GfVec3d cameraFwd = cameraLocalToWorld.TransformDir(pxr::GfVec3d(0.0, 0.0, -1.0));
        Vec3 worldPos = Vec3(float(cameraP[0]), float(cameraP[1]), float(cameraP[2]));
        Vec3 forward = Vec3(float(cameraFwd[0]), float(cameraFwd[1]), float(cameraFwd[2]));
        const float forwardLen = Length(forward);
        if (forwardLen > 1e-8f)
        {
            forward = forward * (1.0f / forwardLen);
        }
        else
        {
            forward = Vec3(0.0f, 0.0f, -1.0f);
        }
        const float verticalFovDegrees =
            float(gfCam.GetFieldOfView(pxr::GfCamera::FOVVertical));
        const pxr::GfRange1d nearFar = frustum.GetNearFar();
        const float nearPlane = float(nearFar.GetMin());

        Camera &uc = scenePool->camera;
        for (int i = 0; i < 4; i++)
        {
            for (int j = 0; j < 4; j++)
            {
                uc.cameraFromWorld.m[i][j] = static_cast<float>(viewM[i][j]);
                uc.clipFromCamera.m[i][j] = static_cast<float>(projM[i][j]);

                printf("%f ", projM[i][j]);
            }
        }

        if (settings)
        {
            pxr::GfVec2i resolution;
            if (settings.GetResolutionAttr().Get(&resolution))
            {
                uc.viewportWidth = resolution[0];
                uc.viewportHeight = resolution[1];
                printf("resolution: %i %i\n", resolution[0], resolution[1]);
            }
        }
        uc.worldPosition = worldPos;
        uc.forward = forward;
        uc.verticalFovDegrees =
            (std::isfinite(verticalFovDegrees) && verticalFovDegrees > 0.0f) ? verticalFovDegrees
                                                                              : 45.0f;
        uc.nearPlane = (std::isfinite(nearPlane) && nearPlane > 0.0f) ? nearPlane : 1.0f;
        uc.hasValidCamera = true;
        uc.path = cameraPrim.GetPath().GetString();

        double shutterOpen, shutterClose;
        camera.GetShutterOpenAttr().Get(&shutterOpen, 0.0);
        camera.GetShutterCloseAttr().Get(&shutterClose, 0.0);
        printf("open: %f close: %f\n", shutterOpen, shutterClose);
    }

    std::unordered_map<std::string, int> materialMap;
    std::vector<pxr::UsdShadeMaterial> materials;

    for (const USDBuildScene &buildScene : buildSceneDAG.scenes)
    {
        for (const USDBuildSceneCurve &curveRef : buildScene.curves)
        {
            const pxr::UsdPrim prim = stage->GetPrimAtPath(pxr::SdfPath(curveRef.path));
            if (!prim || !prim.IsA<pxr::UsdGeomBasisCurves>())
            {
                printf("curve prim missing or invalid: %s\n", curveRef.path.c_str());
                return;
            }
            if (!IsVisibilityAllowed(prim))
            {
                continue;
            }
            if (!IsPurposeAllowed(prim, allowedPurposes))
            {
                continue;
            }
            pxr::UsdShadeMaterial material = GetPrimMaterial(prim);
            AddMaterialToMap(materialMap, materials, material);
        }

        for (const USDBuildSceneMesh &meshRef : buildScene.meshes)
        {
            const pxr::UsdPrim prim = stage->GetPrimAtPath(pxr::SdfPath(meshRef.path));
            if (!prim || !prim.IsA<pxr::UsdGeomMesh>())
            {
                printf("mesh prim missing or invalid: %s\n", meshRef.path.c_str());
                return;
            }
            if (!IsVisibilityAllowed(prim))
            {
                continue;
            }
            if (!IsPurposeAllowed(prim, allowedPurposes))
            {
                continue;
            }

            pxr::UsdGeomMesh mesh(prim);
            const std::vector<pxr::UsdGeomSubset> subsets =
                pxr::UsdGeomSubset::GetAllGeomSubsets(mesh);
            if (subsets.size())
            {
                printf("subsets not handled yet\n");
            }
            else
            {
                pxr::UsdShadeMaterial material = GetPrimMaterial(mesh.GetPrim());
                AddMaterialToMap(materialMap, materials, material);
            }
        }
    }

    printf("num materials: %zi\n", materials.size());

    std::vector<MaterialInfo> materialTextures;
    materialTextures.reserve(materials.size());

    for (pxr::UsdShadeMaterial &material : materials)
    {
        MaterialInfo info = {};
        info.materialPath = material.GetPath().GetString();
        ReadNtcMaterialInfo(material, &info);
        info.packed = {};

        if (pxr::UsdShadeShader shader = material.ComputeSurfaceSource())
        {
            pxr::TfToken token;
            if (shader.GetShaderId(&token) && token == pxr::TfToken("UsdPreviewSurface"))
            {
                ReadPreviewSurfaceMaterial(shader, &info);
            }
            else
            {
                // TODO: support non-UsdPreviewSurface material networks.
                printf("Texture gather: material %s uses unsupported surface shader %s\n",
                       material.GetPath().GetText(),
                       token.GetText());
                FinalizePackedMaterial(&info.packed);
            }
        }
        else
        {
            printf("Texture gather: material %s has no surface source\n", material.GetPath().GetText());
            FinalizePackedMaterial(&info.packed);
        }
        materialTextures.push_back(std::move(info));
    }

    for (const MaterialInfo &info : materialTextures)
    {
        printf("material %s image textures: %zu\n",
               info.materialPath.c_str(),
               info.textureInputs.size());
        for (const MaterialTextureInput &textureInput : info.textureInputs)
        {
            printf("  %s (%s) wrapS=%s wrapT=%s\n",
                   textureInput.texturePath.c_str(),
                   textureInput.inputName.c_str(),
                   TextureWrapModeToString(textureInput.wrapS),
                   TextureWrapModeToString(textureInput.wrapT));
        }
        if (!info.ntcDiffuseFile.empty())
        {
            printf("  ntc diffuse: %s (%s)\n",
                   info.ntcDiffuseFile.c_str(),
                   info.ntcDiffuseTextureName.c_str());
        }
    }
    scenePool->materials = std::move(materialTextures);

    const pxr::UsdTimeCode lightTimeCode(0.0);
    CollectUsdLights(stage, lightTimeCode, &scenePool->lights);
    std::printf("num lights: %zu\n", scenePool->lights.size());
    for (const LightInfo &light : scenePool->lights)
    {
        std::printf("light %s type=%u color=(%.3f %.3f %.3f) scale=%.3f\n",
                    light.lightPath.c_str(),
                    light.packed.type,
                    light.packed.color.x,
                    light.packed.color.y,
                    light.packed.color.z,
                    light.packed.emissionScale);
        if (!light.shadowExcludePaths.empty())
        {
            std::printf("  shadowLink excludes: %zu\n", light.shadowExcludePaths.size());
        }
    }

    int total = 0;
    for (size_t sceneIndex = 0; sceneIndex < buildSceneDAG.scenes.size(); sceneIndex++)
    {
        const USDBuildScene &buildScene = buildSceneDAG.scenes[sceneIndex];
        Scene *outScene = scenePool->scenes[sceneIndex].get();

        for (const USDBuildSceneCurve &curveRef : buildScene.curves)
        {
            const pxr::UsdPrim prim = stage->GetPrimAtPath(pxr::SdfPath(curveRef.path));
            if (!prim || !prim.IsA<pxr::UsdGeomBasisCurves>())
            {
                printf("curve prim missing or invalid: %s\n", curveRef.path.c_str());
                return;
            }
            if (!IsPurposeAllowed(prim, allowedPurposes))
            {
                continue;
            }

            pxr::UsdGeomBasisCurves curve(prim);
            ProcessUSDBasisCurve(curve, outScene);
            outScene->curves.back().parentFromLocal = curveRef.parentFromLocal;
        }

        for (const USDBuildSceneMesh &meshRef : buildScene.meshes)
        {
            const pxr::UsdPrim prim = stage->GetPrimAtPath(pxr::SdfPath(meshRef.path));
            if (!prim || !prim.IsA<pxr::UsdGeomMesh>())
            {
                printf("mesh prim missing or invalid: %s\n", meshRef.path.c_str());
                return;
            }
            if (!IsPurposeAllowed(prim, allowedPurposes))
            {
                continue;
            }

            pxr::UsdGeomMesh mesh(prim);
            const pxr::UsdShadeMaterial material = GetPrimMaterial(mesh.GetPrim());
            const int materialIndex = GetMaterialIndex(materialMap, material);
            pxr::VtVec3fArray positions;
            pxr::VtIntArray faceIndices;
            pxr::VtIntArray faceCounts;
            pxr::TfToken scheme;

            mesh.GetSubdivisionSchemeAttr().Get(&scheme);
            if (scheme == pxr::UsdGeomTokens->catmullClark)
            {
                ProcessCatmullClarkMesh(mesh, outScene, materialIndex, meshRef.parentFromLocal);
                total++;
                continue;
            }
            else if (scheme == pxr::UsdGeomTokens->none)
            {
            }
            else
            {
                printf("%s\n", scheme.GetText());
            }

            USD_ASSERT(mesh.GetPointsAttr().Get(&positions, 0.0));
            USD_ASSERT(mesh.GetFaceVertexIndicesAttr().Get(&faceIndices, 0.0));
            USD_ASSERT(mesh.GetFaceVertexCountsAttr().Get(&faceCounts, 0.0));

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
                        YBI_ASSERT(0);
                    }
                }
                else
                {
                    numTriangles++;
                }
            }

            Array<Vec3> finalPositions(positions);
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
                        YBI_ASSERT(index < positions.size());
                        finalIndices[finalOffset++] = index;
                    }
                }
                else if (faceCount == 4)
                {
                    int tempIndices[4];
                    for (int i = 0; i < 4; i++)
                    {
                        int index = faceIndices[inputOffset++];
                        YBI_ASSERT(index < positions.size());
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
                    YBI_ASSERT(0);
                }
            }

            Attribute stAttr;
            const bool hasSt =
                BuildTriangulatedMeshSTAttribute(mesh, faceCounts, faceIndices, numTriangles, &stAttr);

            outScene->meshes.emplace_back(
                std::move(finalPositions), std::move(finalIndices), meshRef.parentFromLocal);
            Mesh &outMesh = outScene->meshes.back();
            outMesh.primPath = meshRef.path;
            outMesh.materialIndex = materialIndex;
            if (hasSt)
            {
                outMesh.attributes.push_back(std::move(stAttr));
            }
        }
    }

    printf("num cat clarks: %i\n", total);
}

YBI_NAMESPACE_END
