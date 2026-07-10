#include "core/arena.h"
#include "core/assert.h"
#include "device/device.h"
#include "device/kernel_params.h"
#include "display_device/display_device.h"
#include "input/camera_controller.h"
#include "io/usd.h"
#include "render/virtual_texture_manager.h"
#include "scene/scene.h"
#include <cstdlib>
#include <memory>
#include <string>
#include <vector>

#include <GLFW/glfw3.h>

#define STB_IMAGE_WRITE_IMPLEMENTATION
#include <stb_image_write.h>

using namespace ybi;

int main(int argc, char **argv)
{
    YBI_ASSERT(argc == 2);

    // char *filename = argv[0];

    std::string filename(argv[1]);
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
    Scene *scene = sceneLoader->LoadScene(&arena, &container);

    std::unique_ptr<Device> device = Device::CreateDevice(DeviceKind::GPU);

    // Virtual textures
    VirtualTextureManager manager(arena, device.get(), displayDevice.get());
    for (Material &mat : container.materials)
    {
        MaterialTextureHandle handle = manager.AddTexture(mat.diffuseUdimPaths);
        mat.udimBase = handle.udimBase;
        mat.udimTableOffset = handle.udimTableOffset;
    }
    manager.Finalize();

    RenderSystem system(device.get(), &container, scene);
    system.Update();

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

    KernelParams params;
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

    uint32_t lastFrameIndex = 0;
    double lastTime = glfwGetTime();
    bool wasEscapePressed = false;

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
        device->CopyToDevice(cudaParams.Data(), &params, sizeof(KernelParams));
        device->DispatchKernel(camera.width, camera.height, cudaParams);
        manager.CopyFeedback();
        device->SignalExternalSemaphore(deviceFrameFinishedSemaphores[frameIndex]);

        displayDevice->EndFrame(*displayImageHandles[frameIndex],
                                *displayFrameFinishedSemaphores[frameIndex],
                                camera.width,
                                camera.height);
    }

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
    device->Free(cudaParams);

    glfwDestroyWindow(window);
    glfwTerminate();

    return 0;
}
