#include "application.h"
#include "core/arena.h"
#include "core/assert.h"
#include "core/math.h"
#include "core/platform.h"
#include "device/device.h"
#include "device/kernel_params.h"
#include "display_device/display_device.h"
#include "input/camera_controller.h"
#include "io/usd.h"
#include "render/progressive_state.h"
#include "render/virtual_texture_manager.h"
#include "scene/lut.h"
#include "scene/render_system.h"
#include "scene/scene.h"
#include <cerrno>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <exception>
#include <limits>
#include <memory>
#include <optional>
#include <string>
#include <utility>
#include <vector>

#include <GLFW/glfw3.h>
#ifdef WITH_VULKAN
#include <imgui.h>
#endif

using namespace ybi;

namespace
{
constexpr int RenderWidth = 1920;
constexpr int RenderHeight = 1080;

struct SceneSession
{
    Device *device;
    DisplayDevice *displayDevice;
    std::string filename;
    DeviceTexture dielectricEnergy;
    DeviceTexture dielectricAverage;
    std::unique_ptr<Arena> arena = std::make_unique<Arena>();
    DataContainer container;
    std::unique_ptr<DeviceArena> containerDeviceArena;
    std::unique_ptr<VirtualTextureManager> textures;
    std::unique_ptr<RenderSystem> renderSystem;
    Scene *scene = nullptr;
    size_t selectedCamera = 0;

    SceneSession(Device *device, DisplayDevice *displayDevice, std::string filename)
        : device(device), displayDevice(displayDevice), filename(std::move(filename))
    {
    }

    ~SceneSession()
    {
        device->DestroyTexture(dielectricAverage);
        device->DestroyTexture(dielectricEnergy);
    }

    bool Load(std::string *error)
    {
        std::unique_ptr<SceneLoader> loader = SceneLoader::CreateSceneLoader(filename);
        if (!loader)
        {
            *error = "Choose a .usd, .usda, .usdc, or .usdz file.";
            return false;
        }

        SceneImportContext importContext;
        scene = loader->LoadScene(arena.get(), &container, &importContext);
        if (!scene)
        {
            *error = "Unable to open a renderable USD stage.";
            return false;
        }

        const bool hasOpenPbr = !container.openPbrMaterials.empty();
        const bool hasMaterialX = !container.materialXMaterials.empty();
        if (hasOpenPbr && hasMaterialX)
        {
            *error = "Scenes containing both OpenPBR and MaterialX materials are not supported.";
            return false;
        }
        const MaterialBackend materialBackend =
            hasMaterialX ? MaterialBackend::MaterialX : MaterialBackend::OpenPbr;
        if (!device->CreatePipeline(materialBackend))
        {
            *error = "Unable to create the material rendering pipeline.";
            return false;
        }
        if (hasOpenPbr)
        {
            dielectricEnergy = LoadOpenPbrOpaqueDielectricEnergyComplementTexture(*device);
            dielectricAverage = LoadOpenPbrOpaqueDielectricAverageEnergyComplementTexture(*device);
            if (!dielectricEnergy.textureObject || !dielectricAverage.textureObject)
            {
                *error = "Unable to load the OpenPBR lookup tables.";
                return false;
            }
        }

        containerDeviceArena = device->CreateArena();
        container.Init(device, containerDeviceArena.get());
        textures = std::make_unique<VirtualTextureManager>(*arena, device, displayDevice);
        for (Material &material : container.materials)
        {
            MaterialTextureHandle handle = textures->AddTexture(material.diffuseUdimPaths);
            material.udimBase = handle.udimBase;
            material.udimTableOffset = handle.udimTableOffset;
            material.valid = material.udimBase != 0u;
        }
        renderSystem = std::make_unique<RenderSystem>(
            device, displayDevice, &container, scene, textures.get());
        renderSystem->Update(&importContext);
        textures->Finalize();

        if (container.GetBVHs(GeometryType::Scene).size() != 1)
        {
            *error = "The USD stage has no renderable scene.";
            return false;
        }
        AffineTransform renderFromCamera;
        Matrix4x4 cameraFromRaster;
        if (!container.camera.cameraFromWorld.TryInverse(&renderFromCamera) ||
            !container.camera.rasterFromCamera.TryInverse(&cameraFromRaster))
        {
            *error = "The USD camera is invalid.";
            return false;
        }
        return true;
    }
};

} // namespace

int main(int argc, char **argv)
{
    bool pathTrace = false;
    bool interactive = false;
    uint32_t samplesPerPixel = 0;
    std::string filename;
    std::string outputPath = "render.png";
    std::string cameraPath;
    std::optional<Vec2u> pixel;
    bool positionalOnly = false;
    const auto printUsage = []() {
        std::printf(
            "Usage: yeoubi [scene] [--path-trace | --interactive] [--spp 1..4096] [--output file.png] [--camera file.json] [--pixel X Y]\n"
            "A scene path selects offline rendering; --interactive opens the window instead.\n"
            "With no scene, open a USD file from File > Open USD.\n"
            "--spp sets offline samples (default 16). --output sets its PNG name (default render.png).\n"
            "Interactive mode ignores --output.\n"
            "--camera overrides USD cameras in either mode. Export one from the interactive File menu.\n"
            "--pixel renders one offline pixel and writes a 1x1 PNG.\n"
            "Esc: toggle camera controls; WASD/QE: move; mouse: look; R: reset accumulation.\n");
    };
    const auto parsePixelCoordinate = [](const char *text, uint32_t *coordinate) {
        if (!text[0] || text[0] == '-' || text[0] == '+')
            return false;
        char *end = nullptr;
        errno = 0;
        const unsigned long long value = std::strtoull(text, &end, 10);
        if (errno || end == text || *end || value > std::numeric_limits<uint32_t>::max())
            return false;
        *coordinate = static_cast<uint32_t>(value);
        return true;
    };
    for (int i = 1; i < argc; ++i)
    {
        if (!positionalOnly && std::strcmp(argv[i], "--help") == 0)
        {
            printUsage();
            return 0;
        }
        else if (!positionalOnly && std::strcmp(argv[i], "--") == 0)
            positionalOnly = true;
        else if (!positionalOnly && std::strcmp(argv[i], "--path-trace") == 0)
            pathTrace = true;
        else if (!positionalOnly && std::strcmp(argv[i], "--interactive") == 0)
        {
            interactive = true;
            pathTrace = true;
        }
        else if (!positionalOnly && std::strcmp(argv[i], "--spp") == 0 && i + 1 < argc)
        {
            char *end = nullptr;
            errno = 0;
            const long value = std::strtol(argv[++i], &end, 10);
            if (errno || end == argv[i] || *end || value < 1 || value > 4096)
            {
                std::fprintf(stderr, "--spp must be an integer from 1 to 4096.\n");
                return 1;
            }
            samplesPerPixel = static_cast<uint32_t>(value);
        }
        else if (!positionalOnly && std::strcmp(argv[i], "--output") == 0)
        {
            if (i + 1 >= argc || argv[i + 1][0] == '-')
            {
                std::fprintf(stderr, "--output requires a PNG filename.\n");
                return 1;
            }
            outputPath = argv[++i];
        }
        else if (!positionalOnly && std::strcmp(argv[i], "--camera") == 0)
        {
            if (i + 1 >= argc || argv[i + 1][0] == '-')
            {
                std::fprintf(stderr, "--camera requires a JSON filename.\n");
                return 1;
            }
            cameraPath = argv[++i];
        }
        else if (!positionalOnly && std::strcmp(argv[i], "--pixel") == 0)
        {
            uint32_t x = 0;
            uint32_t y = 0;
            if (i + 2 >= argc || !parsePixelCoordinate(argv[i + 1], &x) ||
                !parsePixelCoordinate(argv[i + 2], &y))
            {
                std::fprintf(stderr, "--pixel requires non-negative X and Y coordinates.\n");
                return 1;
            }
            pixel = Vec2u(x, y);
            i += 2;
        }
        else if (!positionalOnly && argv[i][0] == '-')
        {
            std::fprintf(stderr, "Unknown or incomplete option: %s\n", argv[i]);
            return 1;
        }
        else if (filename.empty())
            filename = argv[i];
        else
        {
            std::fprintf(stderr, "Expected one scene filename; extra argument: %s\n", argv[i]);
            return 1;
        }
    }
    if (samplesPerPixel == 0)
        samplesPerPixel = interactive ? 1u : 16u;
    if (pixel && (interactive || filename.empty()))
    {
        std::fprintf(stderr, "--pixel is available only for offline rendering with a scene.\n");
        return 1;
    }

    std::optional<CameraJson> cameraOverride;
    if (!cameraPath.empty())
    {
        CameraJson camera;
        std::string error;
        if (!LoadCameraJson(cameraPath, &camera, &error))
        {
            std::fprintf(stderr, "%s\n", error.c_str());
            return 1;
        }
        cameraOverride = camera;
    }

    if (interactive || filename.empty())
    {
        glfwInit();
        glfwWindowHint(GLFW_CLIENT_API, GLFW_NO_API);
        {
            InteractiveApplication app(DeviceKind::GPU, filename, cameraOverride);
            app.Run();
        }
        glfwTerminate();
    }
    else
    {
        OfflineApplication app(DeviceKind::GPU, filename, samplesPerPixel, outputPath,
                               cameraOverride, pixel);
        app.Run();
    }
    return 0;
}
