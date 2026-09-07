#include "core/arena.h"
#include "core/assert.h"
#include "device/device.h"
#include "device/kernel_params.h"
#include "display_device/display_device.h"
#include "input/camera_controller.h"
#include "io/usd.h"
#include "render/progressive_state.h"
#include "render/virtual_texture_manager.h"
#include "scene/render_system.h"
#include "scene/scene.h"
#include <cerrno>
#include <cstdlib>
#include <cstring>
#include <memory>
#include <string>
#include <vector>

#include <GLFW/glfw3.h>

#define STB_IMAGE_WRITE_IMPLEMENTATION
#include <stb_image_write.h>

using namespace ybi;

int main(int argc, char **argv)
{
    bool pathTrace = false;
    bool interactive = false;
    uint32_t samplesPerPixel = 0;
    std::string filename;
    bool positionalOnly = false;
    const auto printUsage = []() {
        std::printf(
            "Usage: yeoubi <scene> [--path-trace | --interactive] [--spp 1..4096]\n"
            "The scene filename and options may appear in any order.\n"
            "Default: base-color preview. --path-trace: 16 samples per frame.\n"
            "--interactive: progressive diffuse tracing, 1 sample per frame.\n"
            "--spp overrides samples per frame. Lighting is a constant environment.\n"
            "Esc: toggle camera controls; WASD/QE: move; mouse: look; R: reset accumulation.\n");
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
    if (filename.empty())
    {
        std::fprintf(stderr, "A scene filename is required.\n");
        printUsage();
        return 1;
    }
    if (samplesPerPixel == 0)
        samplesPerPixel = interactive ? 1u : 16u;

    Arena arena;
    DataContainer container;

    glfwInit();
    glfwWindowHint(GLFW_CLIENT_API, GLFW_NO_API);
    GLFWwindow *window = glfwCreateWindow(1920, 1080, "Yeoubi", nullptr, nullptr);
    YBI_ASSERT(window);
    glfwSetInputMode(window, GLFW_CURSOR, GLFW_CURSOR_NORMAL);

    std::unique_ptr<DisplayDevice> displayDevice = DisplayDevice::CreateDevice(window);
    YBI_ASSERT(displayDevice);

    std::unique_ptr<SceneLoader> sceneLoader = SceneLoader::CreateSceneLoader(filename);
    YBI_ASSERT(sceneLoader);
    SceneImportContext importContext;
    Scene *scene = sceneLoader->LoadScene(&arena, &container, &importContext);

    std::unique_ptr<Device> device = Device::CreateDevice(DeviceKind::GPU);

    // Virtual textures
    VirtualTextureManager manager(arena, device.get(), displayDevice.get());
    for (Material &mat : container.materials)
    {
        MaterialTextureHandle handle = manager.AddTexture(mat.diffuseUdimPaths);
        mat.udimBase = handle.udimBase;
        mat.udimTableOffset = handle.udimTableOffset;

        // TODO: this disables materials with no diffuse map. support actual materials
        mat.valid = mat.udimBase != 0u;
    }
    RenderSystem system(device.get(), displayDevice.get(), &container, scene, &manager);
    system.Update(&importContext);
    manager.Finalize();

    YBI_ASSERT(scene->container->GetBVHs(GeometryType::Scene).size() == 1);

    AffineTransform renderFromCamera;
    Matrix4x4 cameraFromRaster;

    const Camera &camera = scene->container->camera;
    bool result = camera.cameraFromWorld.TryInverse(&renderFromCamera);
    YBI_ASSERT(result);
    result = camera.rasterFromCamera.TryInverse(&cameraFromRaster);
    YBI_ASSERT(result);
    CameraController cameraController(window, renderFromCamera);

    const size_t imageSize =
        static_cast<size_t>(camera.width) * static_cast<size_t>(camera.height) * 4;
    const uint32_t numFramesInFlight = displayDevice->GetNumFramesInFlight();
    std::vector<std::unique_ptr<DisplayBuffer>> displayImageHandles;
    std::vector<DeviceExternalBuffer> deviceImageHandles;
    std::vector<std::unique_ptr<DisplaySemaphore>> displayFrameFinishedSemaphores;
    std::vector<DeviceExternalSemaphore> deviceFrameFinishedSemaphores;
    displayImageHandles.reserve(numFramesInFlight);
    deviceImageHandles.reserve(numFramesInFlight);
    displayFrameFinishedSemaphores.reserve(numFramesInFlight);
    deviceFrameFinishedSemaphores.reserve(numFramesInFlight);

    for (uint32_t frameIndex = 0; frameIndex < numFramesInFlight; ++frameIndex)
    {
        displayImageHandles.push_back(displayDevice->AllocateExternalBuffer(imageSize));
        deviceImageHandles.push_back(device->ImportExternalBuffer(*displayImageHandles.back()));
        displayFrameFinishedSemaphores.push_back(displayDevice->CreateDisplaySemaphore());
        deviceFrameFinishedSemaphores.push_back(
            device->ImportExternalSemaphore(*displayFrameFinishedSemaphores.back()));
    }

    KernelParams params{};
    params.pathTrace = pathTrace;
    params.samplesPerPixel = samplesPerPixel;
    params.bvhHandle = scene->container->GetBVHs(GeometryType::Scene)[0].handle.handle;
    params.cameraPosition = camera.position;
    params.fovRadians = camera.verticalFovRadians;
    params.cameraFromRaster = cameraFromRaster;
    params.renderFromCamera = renderFromCamera;
    params.maxDepth = 3;
    params.lensRadius = 0.f; // TODO: initialize properly
    params.imageWidth = camera.width;
    params.imageHeight = camera.height;
    params.scene = system.deviceScene;

    DeviceBuffer cudaParams = device->Alloc(sizeof(KernelParams));
    DeviceBuffer accumulationBuffer;
    if (interactive)
    {
        accumulationBuffer = device->Alloc(static_cast<size_t>(camera.width) *
                                           static_cast<size_t>(camera.height) * sizeof(Vec3f));
        params.accumulation = reinterpret_cast<Vec3f *>(accumulationBuffer.Data());
    }
    ProgressiveState progressiveState;

    uint32_t lastFrameIndex = 0;
    double lastTime = glfwGetTime();
    bool wasEscapePressed = false;
    bool wasResetPressed = false;
    bool resetPending = false;

    while (!glfwWindowShouldClose(window))
    {
        const double currentTime = glfwGetTime();
        const float deltaTime = static_cast<float>(currentTime - lastTime);
        lastTime = currentTime;

        glfwPollEvents();
        cameraController.UpdateCamera(window, deltaTime);
        params.renderFromCamera = cameraController.GetRenderFromCamera();
        Vec3f cameraPosition = params.renderFromCamera.GetTranslation();
        params.cameraPosition = cameraPosition;

        const bool escapePressed = glfwGetKey(window, GLFW_KEY_ESCAPE) == GLFW_PRESS;
        if (escapePressed && !wasEscapePressed)
        {
            int currentMode = glfwGetInputMode(window, GLFW_CURSOR);
            int cursorMode =
                currentMode == GLFW_CURSOR_DISABLED ? GLFW_CURSOR_NORMAL : GLFW_CURSOR_DISABLED;

            glfwSetInputMode(window, GLFW_CURSOR, cursorMode);
            cameraController.SetMovementEnabled(cursorMode == GLFW_CURSOR_DISABLED);
        }
        wasEscapePressed = escapePressed;
        const bool resetPressed = glfwGetKey(window, GLFW_KEY_R) == GLFW_PRESS;
        resetPending |= resetPressed && !wasResetPressed;
        wasResetPressed = resetPressed;

        if (!displayDevice->BeginFrame())
        {
            continue;
        }

        const uint32_t frameIndex = displayDevice->GetCurrentFrameIndex();
        lastFrameIndex = frameIndex;

        DeviceExternalBuffer &deviceImageHandle = deviceImageHandles[frameIndex];
        params.imageHandle = reinterpret_cast<uint64_t>(deviceImageHandle.buffer.Data());
        manager.UpdateParams(&params);
        manager.Update(&params);
        if (interactive)
        {
            progressiveState.BeginFrame(
                params.renderFromCamera,
                progressiveState.textureRevision, // manager.GetResidencyRevision(),
                samplesPerPixel,
                resetPending);
            params.accumulatedSamples = progressiveState.samples;
            params.sampleOffset = progressiveState.samples;
        }
        resetPending = false;
        device->CopyToDevice(cudaParams.Data(), &params, sizeof(KernelParams));
        device->DispatchKernel(camera.width, camera.height, cudaParams);
        if (interactive)
        {
            progressiveState.samples += samplesPerPixel;
            char title[128];
            std::snprintf(title,
                          sizeof(title),
                          "Yeoubi | Diffuse path tracing | %u spp | R: reset",
                          progressiveState.samples);
            glfwSetWindowTitle(window, title);
        }
        else if (params.pathTrace)
            params.sampleOffset += params.samplesPerPixel;
        manager.CopyFeedback();
        device->SignalExternalSemaphore(deviceFrameFinishedSemaphores[frameIndex]);

        displayDevice->EndFrame(*displayImageHandles[frameIndex],
                                *displayFrameFinishedSemaphores[frameIndex],
                                camera.width,
                                camera.height);
    }

    // Finish presentation before reading back or releasing shared CUDA/Vulkan resources.
    displayDevice->WaitIdle();

    const size_t outputSize =
        static_cast<size_t>(camera.width) * static_cast<size_t>(camera.height) * 4;
    void *outBytes = malloc(outputSize);
    device->CopyToHost(outBytes, deviceImageHandles[lastFrameIndex].buffer.Data(), outputSize);

    stbi_write_png("test.png", camera.width, camera.height, 4, outBytes, camera.width * 4);
    free(outBytes);

    for (DeviceExternalSemaphore &deviceFrameFinishedSemaphore : deviceFrameFinishedSemaphores)
    {
        device->FreeExternalSemaphore(deviceFrameFinishedSemaphore);
    }
    deviceFrameFinishedSemaphores.clear();
    displayFrameFinishedSemaphores.clear();

    for (DeviceExternalBuffer &deviceImageHandle : deviceImageHandles)
    {
        device->FreeExternalBuffer(deviceImageHandle);
    }
    deviceImageHandles.clear();
    displayImageHandles.clear();
    displayDevice.reset();
    if (accumulationBuffer)
        device->Free(accumulationBuffer);
    device->Free(cudaParams);

    glfwDestroyWindow(window);
    glfwTerminate();

    return 0;
}
