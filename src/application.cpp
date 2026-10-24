#include "application.h"
#include "core/assert.h"
#include "core/math.h"
#include "core/platform.h"
#include "core/time.h"
#include "device/device.h"
#include "display_device/display_device.h"
#include "input/camera_controller.h"
#include "io/usd.h"
#include "kernel/launch/kernel_params.h"
#include "material/lut.h"
#include "render/render_system.h"
#include "render/virtual_texture/manager.h"
#include "render/wavefront_integrator.h"
#include "scene/scene.h"
#include <GLFW/glfw3.h>
#include <cstdio>
#include <cstdlib>
#include <exception>
#define STB_IMAGE_WRITE_IMPLEMENTATION
#include <stb_image_write.h>

#ifdef WITH_VULKAN
#include <imgui.h>
#endif

namespace ybi
{

GlfwWindow::GlfwWindow(int width, int height, const char *title)
    : window(glfwCreateWindow(width, height, title, nullptr, nullptr))
{
}

GlfwWindow::~GlfwWindow()
{
    if (window)
        glfwDestroyWindow(window);
}

SceneSession::SceneSession(Device *device, DisplayDevice *displayDevice, std::string filename)
    : device(device), displayDevice(displayDevice), filename(std::move(filename)),
      container(std::make_unique<DataContainer>())
{
}

SceneSession::~SceneSession()
{
    device->DestroyTexture(dielectricAverage);
    device->DestroyTexture(dielectricEnergy);
}

bool SceneSession::Load(std::string *error, bool validateUsdCamera)
{
    std::unique_ptr<SceneLoader> loader = SceneLoader::CreateSceneLoader(filename);
    if (!loader)
    {
        *error = "Choose a .usd, .usda, .usdc, or .usdz file.";
        return false;
    }

    SceneImportContext importContext;
    scene = loader->LoadScene(arena.get(), container.get(), &importContext);
    if (!scene)
    {
        *error = "Unable to open a renderable USD stage.";
        return false;
    }

    for (EnvmapLight &light : container->envmapLights)
    {
        if (!LoadEnvmap(light, error))
            return false;
    }

    const bool hasOpenPbr = !container->openPbrMaterials.empty();
    const bool hasMaterialX = !container->materialXMaterials.empty();
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
    container->Init(device, containerDeviceArena.get());
    textures = std::make_unique<VirtualTextureManager>(*arena, device, displayDevice);

    // TODO: this needs to be removed now
    for (Material &material : container->materials)
    {
        MaterialTextureHandle handle = textures->AddTexture(material.diffuseUdimPaths);
        material.udimBase = handle.udimBase;
        material.udimTableOffset = handle.udimTableOffset;
        material.valid = material.udimBase != 0u;
    }
    renderSystem = std::make_unique<RenderSystem>(
        device, displayDevice, container.get(), scene, textures.get());
    renderSystem->Update(&importContext);
    textures->Finalize();

    if (container->GetBVHs(GeometryType::Scene).size() != 1)
    {
        *error = "The USD stage has no renderable scene.";
        return false;
    }
    AffineTransform renderFromCamera;
    Matrix4x4 cameraFromRaster;
    if (validateUsdCamera && (!container->camera.cameraFromWorld.TryInverse(&renderFromCamera) ||
                              !container->camera.rasterFromCamera.TryInverse(&cameraFromRaster)))
    {
        *error = "The USD camera is invalid.";
        return false;
    }
    return true;
}

Application::Application(DeviceKind kind) : device(Device::CreateDevice(kind, false)) {}

// TODO: stop hardcoding width height
InteractiveApplication::InteractiveApplication(DeviceKind kind,
                                               const std::string &initialScenePath,
                                               std::optional<CameraJson> cameraOverride,
                                               size_t windowWidth,
                                               size_t windowHeight)
    : Application(kind), window(int(cameraOverride ? cameraOverride->width : windowWidth),
                                int(cameraOverride ? cameraOverride->height : windowHeight),
                                "yeoubi"),
      windowWidth(cameraOverride ? cameraOverride->width : windowWidth),
      windowHeight(cameraOverride ? cameraOverride->height : windowHeight),
      displayDevice(DisplayDevice::CreateDevice(window.Get())),
      cameraOverride(std::move(cameraOverride))
{
    params.samplesPerPixel = 1;
    params.maxDepth = 10;
    params.imageWidth = static_cast<int>(this->windowWidth);
    params.imageHeight = static_cast<int>(this->windowHeight);
    integrator = std::make_unique<WavefrontIntegrator>(device.get());
    kernelOutput = device->Alloc<Vec3f>(this->windowWidth * this->windowHeight);
    accumulation = device->Alloc<Vec3f>(this->windowWidth * this->windowHeight);
    params.kernelOutput = kernelOutput.Data();
    params.accumulation = accumulation.Data();
    DeviceKernelDesc kernelDesc;
    kernelDesc.filename = YBI_ACCUMULATE_KERNEL_PATH;
    kernelDesc.entryName = "AccumulateSamples";
    accumulateKernel = device->CreateDeviceKernel(kernelDesc);

    // TODO: this needs to be based on the format of the swapchain images
    const size_t imageSize = this->windowWidth * this->windowHeight * 4;
    const uint32_t framesInFlight = displayDevice->GetNumFramesInFlight();

    for (uint32_t i = 0; i < framesInFlight; ++i)
    {
        displayImages.push_back(displayDevice->AllocateExternalBuffer(imageSize));
        deviceImages.push_back(device->ImportExternalBuffer(*displayImages.back()));
        displaySemaphores.push_back(displayDevice->CreateDisplaySemaphore());
        deviceSemaphores.push_back(device->ImportExternalSemaphore(*displaySemaphores.back()));
    }
    if (!initialScenePath.empty())
        LoadScene(initialScenePath);
}

InteractiveApplication::~InteractiveApplication()
{
    displayDevice->WaitIdle();
    device->WaitIdle();
    device->Free(accumulation);
    device->Free(kernelOutput);
}

void InteractiveApplication::LoadScene(const std::string &path)
{
    displayDevice->WaitIdle();
    device->WaitIdle();
    session.reset();
    device->DestroyPipeline();

    if (cameraController)
    {
        cameraController->SetMovementEnabled(false);
        glfwSetInputMode(window.Get(), GLFW_CURSOR, GLFW_CURSOR_NORMAL);
    }
    auto next = std::make_unique<SceneSession>(device.get(), displayDevice.get(), path);
    std::string error;
    try
    {
        if (!next->Load(&error, !cameraOverride))
        {
            device->WaitIdle();
            next.reset();
            device->DestroyPipeline();
            loadError = error;
            std::fprintf(stderr, "USD load failed: %s\n", loadError.c_str());
            return;
        }
        session = std::move(next);
        CreateCamera();

        loadError.clear();
        std::fprintf(stdout, "Loaded USD scene: %s\n", path.c_str());
        std::fflush(stdout);
    }
    catch (const std::exception &exception)
    {
        device->WaitIdle();
        session.reset();
        next.reset();
        device->DestroyPipeline();
        loadError = "Unable to load the selected USD scene. See stderr for details.";
        std::fprintf(stderr, "USD load failed: %s\n", exception.what());
    }
}

void InteractiveApplication::CreateCamera()
{
    const Camera &camera = session->container->camera;
    AffineTransform renderFromCamera;
    Matrix4x4 cameraFromRaster;
    if (cameraOverride)
    {
        renderFromCamera = cameraOverride->renderFromCamera;
        cameraFromRaster = cameraOverride->cameraFromRaster;
    }
    else
    {
        const bool valid = camera.cameraFromWorld.TryInverse(&renderFromCamera) &&
                           camera.rasterFromCamera.TryInverse(&cameraFromRaster);
        YBI_ASSERT(valid);
    }

    cameraController = std::make_unique<CameraController>(
        window.Get(), renderFromCamera, cameraOverride ? cameraOverride->worldUp : camera.worldUp);

    glfwSetInputMode(window.Get(), GLFW_CURSOR, GLFW_CURSOR_NORMAL);
    params.renderFromCamera = renderFromCamera;
    params.cameraPosition = renderFromCamera.GetTranslation();
    params.cameraFromRaster = cameraFromRaster;
    params.fovRadians =
        cameraOverride ? cameraOverride->verticalFovRadians : camera.verticalFovRadians;
    params.lensRadius = cameraOverride ? cameraOverride->lensRadius : 0.f;
    params.bvhHandle = session->container->GetBVHs(GeometryType::Scene)[0].handle.handle;
    params.scene = session->renderSystem->deviceScene;
    params.opaqueDielectricEnergyComplementTexture = session->dielectricEnergy.textureObject;
    params.opaqueDielectricAverageEnergyComplementTexture =
        session->dielectricAverage.textureObject;
    params.sampleOffset = 0;
    params.accumulatedSamples = 0;
    progressiveState = ProgressiveState{};
}

void InteractiveApplication::SwitchCamera(size_t cameraIndex)
{
    // TODO: this feels weird
    session->selectedCamera = cameraIndex;
    session->container->camera = session->container->cameras[cameraIndex].camera;

    const Camera &camera = session->container->camera;
    AffineTransform renderFromCamera;
    Matrix4x4 cameraFromRaster;
    const bool valid = camera.cameraFromWorld.TryInverse(&renderFromCamera) &&
                       camera.rasterFromCamera.TryInverse(&cameraFromRaster);
    YBI_ASSERT(valid);
#if 0
    if (!cameraController)
        cameraController =
            std::make_unique<CameraController>(window.Get(), renderFromCamera, camera.worldUp);
    else
#endif
    cameraController->ResetCamera(renderFromCamera, camera.worldUp);
    cameraController->SetMovementEnabled(false);

#if 1
    glfwSetInputMode(window.Get(), GLFW_CURSOR, GLFW_CURSOR_NORMAL);
    params.renderFromCamera = renderFromCamera;
    params.cameraPosition = camera.position;
    params.cameraFromRaster = cameraFromRaster;
    params.fovRadians = camera.verticalFovRadians;
    params.lensRadius = 0.f;
    params.bvhHandle = session->container->GetBVHs(GeometryType::Scene)[0].handle.handle;
    params.scene = session->renderSystem->deviceScene;
    params.opaqueDielectricEnergyComplementTexture = session->dielectricEnergy.textureObject;
    params.opaqueDielectricAverageEnergyComplementTexture =
        session->dielectricAverage.textureObject;
    params.sampleOffset = 0;
    params.accumulatedSamples = 0;
    progressiveState = ProgressiveState{};
#endif
}

void InteractiveApplication::Run()
{
    double lastTime = glfwGetTime();
    bool wasEscapePressed = false;
    bool wasResetPressed = false;
    bool resetPending = false;
    while (!glfwWindowShouldClose(window.Get()))
    {
        const double currentTime = glfwGetTime();
        const float deltaTime = static_cast<float>(currentTime - lastTime);
        lastTime = currentTime;
        glfwPollEvents();

        if (session)
        {
            cameraController->UpdateCamera(window.Get(), deltaTime);
            params.renderFromCamera = cameraController->GetRenderFromCamera();
            params.cameraPosition = params.renderFromCamera.GetTranslation();
            const bool escapePressed = glfwGetKey(window.Get(), GLFW_KEY_ESCAPE) == GLFW_PRESS;
            if (escapePressed && !wasEscapePressed)
            {
                const bool capture =
                    glfwGetInputMode(window.Get(), GLFW_CURSOR) != GLFW_CURSOR_DISABLED;
                glfwSetInputMode(window.Get(),
                                 GLFW_CURSOR,
                                 capture ? GLFW_CURSOR_DISABLED : GLFW_CURSOR_NORMAL);
                cameraController->SetMovementEnabled(capture);
            }
            wasEscapePressed = escapePressed;
            const bool resetPressed = glfwGetKey(window.Get(), GLFW_KEY_R) == GLFW_PRESS;
            resetPending |= resetPressed && !wasResetPressed;
            wasResetPressed = resetPressed;
        }

        if (!displayDevice->BeginFrame())
            continue;

        const uint32_t frameIndex = displayDevice->GetCurrentFrameIndex();
        if (session)
        {
            const Timer frameTimer = StartTimer();
            const Timer vtTimer = StartTimer();
            params.imageHandle =
                reinterpret_cast<uint64_t>(deviceImages[frameIndex].buffer.Data());
            session->textures->UpdateParams(&params);
            device->StartProfiler();
            session->textures->Update(&params);
            device->StopProfiler();
            const double vtTime = StopTimer(vtTimer);
            printf("Vt time: %.3f ms\n", vtTime);

            progressiveState.BeginFrame(params.renderFromCamera,
                                        progressiveState.textureRevision,
                                        params.samplesPerPixel,
                                        resetPending);
            params.accumulatedSamples = progressiveState.samples;
            params.sampleOffset = progressiveState.samples;

            resetPending = false;
            device->StartProfiler();
            integrator->RenderSamples(params);
            device->StopProfiler();

            void *arguments[] = {&params};
            DeviceKernelLaunchDesc launch;
            launch.gridX = CeilDiv(windowWidth, 8);
            launch.gridY = CeilDiv(windowHeight, 8);
            launch.blockX = 8;
            launch.blockY = 8;
            launch.arguments = arguments;
            device->ExecuteKernel(*accumulateKernel, launch);
            progressiveState.samples += params.samplesPerPixel;
            session->textures->CopyFeedback();
            device->SignalExternalSemaphore(deviceSemaphores[frameIndex]);
            hasRenderedFrame = true;

            const double ms = StopTimer(frameTimer);
            printf("Frame time: %.3f ms\n", ms);
        }

        std::string pendingPath;
        size_t pendingCamera = ~size_t(0);
#ifdef WITH_VULKAN
        if (ImGui::BeginMainMenuBar())
        {
            if (ImGui::BeginMenu("File"))
            {
                if (ImGui::MenuItem("Open USD..."))
                {
#ifdef _WIN32
                    pendingPath = OpenFileDialog(
                        window.Get(), {"*.usd", "*.usda", "*.usdc", "*.usdz"}, &loadError);
#endif
                }
                if (ImGui::MenuItem("Export Camera JSON...", nullptr, false, session != nullptr))
                {
#ifdef _WIN32
                    const std::string path = SaveFileDialog(window.Get(), {"*.json"}, &loadError);
                    if (!path.empty())
                    {
                        CameraJson camera{};
                        camera.renderFromCamera = params.renderFromCamera;
                        camera.cameraFromRaster = params.cameraFromRaster;
                        camera.worldUp = cameraController->worldUp;
                        camera.verticalFovRadians = params.fovRadians;
                        camera.lensRadius = params.lensRadius;
                        camera.width = params.imageWidth;
                        camera.height = params.imageHeight;
                        if (SaveCameraJson(path, camera, &loadError))
                        {
                            loadError.clear();
                            std::fprintf(stdout, "Wrote camera JSON: %s\n", path.c_str());
                        }
#endif
                    }
                }
                ImGui::EndMenu();
            }
            ImGui::EndMainMenuBar();
        }

        ImGui::SetNextWindowPos(ImVec2(16.f, 40.f), ImGuiCond_FirstUseEver);
        ImGui::SetNextWindowBgAlpha(0.85f);
        ImGui::Begin(session ? "Camera" : "Scene", nullptr, ImGuiWindowFlags_AlwaysAutoResize);
        if (session)
        {
            if (cameraOverride)
                ImGui::TextUnformatted("Camera JSON override");
            else if (!session->container->cameras.empty())
            {
                const auto &cameras = session->container->cameras;
                const char *currentName = cameras[session->selectedCamera].name.c_str();
                if (ImGui::BeginCombo("USD camera", currentName))
                {
                    for (size_t i = 0; i < cameras.size(); ++i)
                    {
                        const bool selected = i == session->selectedCamera;
                        if (ImGui::Selectable(cameras[i].name.c_str(), selected))
                            pendingCamera = i;
                        if (selected)
                            ImGui::SetItemDefaultFocus();
                    }
                    ImGui::EndCombo();
                }
            }
            const Vec3f forward = params.renderFromCamera.TransformVector(Vec3f(0.f, 0.f, -1.f));
            const Vec3f up = params.renderFromCamera.TransformVector(Vec3f(0.f, 1.f, 0.f));
            const Vec3f position = params.cameraPosition;
            ImGui::Text("Position:  %.3f, %.3f, %.3f", position.x, position.y, position.z);
            ImGui::Text("Forward:   %.3f, %.3f, %.3f", forward.x, forward.y, forward.z);
            ImGui::Text("Up:        %.3f, %.3f, %.3f", up.x, up.y, up.z);
            ImGui::Text("Vertical FOV: %.2f deg",
                        params.fovRadians * (180.f / 3.14159265358979323846f));
            ImGui::Text("Resolution: %u x %u", params.imageWidth, params.imageHeight);
            ImGui::Text("Lens radius: %.3f", params.lensRadius);
            ImGui::Text("Move speed: %.2f", cameraController->moveSpeed);
            ImGui::Separator();
            ImGui::Text("SPP per frame: %u", params.samplesPerPixel);
            ImGui::Text("Accumulated SPP: %u", progressiveState.samples);
        }
        else
            ImGui::TextUnformatted("Open a USD file with File > Open USD...");
        if (!loadError.empty())
        {
            ImGui::Separator();
            ImGui::TextWrapped("Error: %s", loadError.c_str());
        }
        ImGui::End();
#endif

        if (session)
            displayDevice->EndFrame(displayImages[frameIndex].get(),
                                    displaySemaphores[frameIndex].get(),
                                    static_cast<int>(windowWidth),
                                    static_cast<int>(windowHeight));
        else
            displayDevice->EndFrame(nullptr, nullptr, 0, 0);

        if (pendingCamera != ~size_t(0) && session)
        {
            SwitchCamera(pendingCamera);
        }
        if (!pendingPath.empty())
            LoadScene(pendingPath);
    }
}

OfflineApplication::OfflineApplication(DeviceKind kind,
                                       const std::string &filePath,
                                       uint32_t sampleCount,
                                       const std::string &outputPath,
                                       std::optional<CameraJson> cameraOverride,
                                       std::optional<Vec2u> pixel)
    : Application(kind), filePath(filePath), outputPath(outputPath), sampleCount(sampleCount),
      singlePixel(pixel.has_value())
{
    session = std::make_unique<SceneSession>(device.get(), nullptr, filePath);
    std::string error;
    try
    {
        if (!session->Load(&error, !cameraOverride))
        {
            std::fprintf(stderr, "USD load failed: %s\n", error.c_str());
            std::exit(EXIT_FAILURE);
        }
    }
    catch (const std::exception &exception)
    {
        std::fprintf(stderr, "USD load failed: %s\n", exception.what());
        std::exit(EXIT_FAILURE);
    }

    if (cameraOverride)
    {
        params.renderFromCamera = cameraOverride->renderFromCamera;
        params.cameraFromRaster = cameraOverride->cameraFromRaster;
        params.cameraPosition = params.renderFromCamera.GetTranslation();
        params.fovRadians = cameraOverride->verticalFovRadians;
        params.lensRadius = cameraOverride->lensRadius;
        params.imageWidth = cameraOverride->width;
        params.imageHeight = cameraOverride->height;
    }
    else
    {
        const Camera &camera = session->container->camera;
        camera.cameraFromWorld.TryInverse(&params.renderFromCamera);
        camera.rasterFromCamera.TryInverse(&params.cameraFromRaster);
        params.cameraPosition = camera.position;
        params.fovRadians = camera.verticalFovRadians;
        params.imageWidth = camera.width;
        params.imageHeight = camera.height;
    }
    params.bvhHandle = session->container->GetBVHs(GeometryType::Scene)[0].handle.handle;
    params.scene = session->renderSystem->deviceScene;
    params.opaqueDielectricEnergyComplementTexture = session->dielectricEnergy.textureObject;
    params.opaqueDielectricAverageEnergyComplementTexture =
        session->dielectricAverage.textureObject;
    params.maxDepth = 10;
    params.samplesPerPixel = 1;
    if (pixel)
    {
        params.pixelOffsetX = pixel->x;
        params.pixelOffsetY = pixel->y;
    }

    const size_t pixelCount =
        singlePixel ? 1u : static_cast<size_t>(params.imageWidth) * params.imageHeight;
    integrator = std::make_unique<WavefrontIntegrator>(device.get());
    kernelOutput = device->Alloc<Vec3f>(pixelCount);
    accumulation = device->Alloc<Vec3f>(pixelCount);
    outputImage = device->Alloc<uint8_t>(pixelCount * 4);
    params.kernelOutput = kernelOutput.Data();
    params.accumulation = accumulation.Data();
    params.imageHandle = outputImage.Address();

    DeviceKernelDesc kernelDesc;
    kernelDesc.filename = YBI_ACCUMULATE_KERNEL_PATH;
    kernelDesc.entryName = "AccumulateSamples";
    accumulateKernel = device->CreateDeviceKernel(kernelDesc);
}

OfflineApplication::~OfflineApplication()
{
    device->WaitIdle();
    device->Free(outputImage);
    device->Free(accumulation);
    device->Free(kernelOutput);
}

void OfflineApplication::Run()
{
    const int outputWidth = singlePixel ? 1 : params.imageWidth;
    const int outputHeight = singlePixel ? 1 : params.imageHeight;
    for (uint32_t sampleIndex = 0; sampleIndex < sampleCount; sampleIndex++)
    {
        session->textures->UpdateParams(&params);
        device->StartProfiler();
        session->textures->UpdateSync(&params);
        device->StopProfiler();

        params.sampleOffset = sampleIndex;
        params.accumulatedSamples = sampleIndex;

        device->StartProfiler();
        const bool rendered = integrator->RenderSamples(params, Vec2u(outputWidth, outputHeight));
        device->StopProfiler();
        if (!rendered)
        {
            std::fprintf(stderr, "Offline rendering failed.\n");
            std::exit(EXIT_FAILURE);
        }

        KernelParams accumulationParams = params;
        accumulationParams.imageWidth = outputWidth;
        accumulationParams.imageHeight = outputHeight;
        void *arguments[] = {&accumulationParams};
        DeviceKernelLaunchDesc launch;
        launch.gridX = CeilDiv(outputWidth, 8);
        launch.gridY = CeilDiv(outputHeight, 8);
        launch.blockX = 8;
        launch.blockY = 8;
        launch.arguments = arguments;
        device->ExecuteKernel(*accumulateKernel, launch);
        session->textures->CopyFeedback();

        device->WaitIdle();
        printf("\n");
    }

    std::vector<unsigned char> pixels(outputImage.SizeInBytes());
    device->CopyToHost(pixels.data(), outputImage.Data(), pixels.size());
    if (!stbi_write_png(
            outputPath.c_str(), outputWidth, outputHeight, 4, pixels.data(), outputWidth * 4))
        std::fprintf(stderr, "Failed to write PNG: %s\n", outputPath.c_str());
    else
        std::fprintf(stdout, "Wrote PNG: %s\n", outputPath.c_str());
}

} // namespace ybi
