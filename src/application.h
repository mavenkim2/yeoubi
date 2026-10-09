#include "device/device_memory.h"
#include "device/device_texture.h"
#include "kernel/launch/kernel_params.h"
#include "core/vec2u.h"
#include "io/camera_json.h"
#include "render/progressive_state.h"
#include "render/path_guiding/lut.h"
#include <memory>
#include <optional>
#include <string>
#include <vector>

struct GLFWwindow;

namespace ybi
{

class Arena;
class CameraController;
class DeviceArena;
class DeviceKernel;
class VirtualTextureManager;
class WavefrontIntegrator;

enum class DeviceKind;

struct DataContainer;
struct Device;
struct DeviceExternalBuffer;
struct DeviceExternalSemaphore;
struct DisplayBuffer;
struct DisplayDevice;
struct DisplaySemaphore;
struct RenderSystem;
struct Scene;

class GlfwWindow
{
    GLFWwindow *window = nullptr;

public:
    GlfwWindow(int width, int height, const char *title);
    ~GlfwWindow();

    GlfwWindow(const GlfwWindow &) = delete;
    GlfwWindow &operator=(const GlfwWindow &) = delete;

    GLFWwindow *Get() const
    {
        return window;
    }
};

struct SceneSession
{
    Device *device;
    DisplayDevice *displayDevice;
    std::string filename;
    DeviceTexture dielectricEnergy;
    DeviceTexture dielectricAverage;
    VMFSplitTable vmfSplitTable;
    std::unique_ptr<Arena> arena = std::make_unique<Arena>();
    std::unique_ptr<DataContainer> container;
    std::unique_ptr<DeviceArena> containerDeviceArena;
    std::unique_ptr<VirtualTextureManager> textures;
    std::unique_ptr<RenderSystem> renderSystem;
    Scene *scene = nullptr;
    size_t selectedCamera = 0;

    SceneSession(Device *device, DisplayDevice *displayDevice, std::string filename);

    ~SceneSession();

    bool Load(std::string *error, bool validateUsdCamera = true);
    void PopulateKernelParams(KernelParams &params);
};

class Application
{
protected:
    std::unique_ptr<Device> device;

public:
    Application(DeviceKind kind);
    virtual void Run() = 0;
};

class InteractiveApplication : Application
{
private:
    // GUI display
    GlfwWindow window;
    size_t windowWidth;
    size_t windowHeight;

    std::unique_ptr<DisplayDevice> displayDevice;
    std::vector<std::unique_ptr<DisplayBuffer>> displayImages;
    std::vector<DeviceExternalBuffer> deviceImages;
    std::vector<std::unique_ptr<DisplaySemaphore>> displaySemaphores;
    std::vector<DeviceExternalSemaphore> deviceSemaphores;
    std::unique_ptr<CameraController> cameraController;
    std::unique_ptr<SceneSession> session;
    std::optional<CameraJson> cameraOverride;
    KernelParams params{};

    std::string loadError;

    ProgressiveState progressiveState;
    std::unique_ptr<WavefrontIntegrator> integrator;
    DeviceBuffer kernelOutput;
    DeviceBuffer accumulation;
    std::unique_ptr<DeviceKernel> accumulateKernel;
    bool hasRenderedFrame = false;

public:
    InteractiveApplication(DeviceKind kind, const std::string &initialScenePath = {},
                           std::optional<CameraJson> cameraOverride = {},
                           size_t windowWidth = 1920, size_t windowHeight = 1080);
    virtual ~InteractiveApplication();
    void Run() override;

private:
    void LoadScene(const std::string &path);
    void CreateCamera();
    void SwitchCamera(size_t cameraIndex);
};

class OfflineApplication : Application
{
private:
    std::string filePath;
    std::string outputPath;
    uint32_t sampleCount = 0;
    bool singlePixel = false;
    std::unique_ptr<SceneSession> session;
    KernelParams params{};
    DeviceBuffer cudaParams;
    DeviceBuffer kernelOutput;
    DeviceBuffer accumulation;
    DeviceBuffer outputImage;
    std::unique_ptr<DeviceKernel> accumulateKernel;

public:
    OfflineApplication(DeviceKind kind, const std::string &filePath, uint32_t sampleCount,
                       const std::string &outputPath,
                       std::optional<CameraJson> cameraOverride = {},
                       std::optional<Vec2u> pixel = {});
    virtual ~OfflineApplication();
    void Run() override;
};

} // namespace ybi
