#include <GLFW/glfw3.h>
#include "display_device/display_device.h"
#include "display_device/vulkan/vulkan.h"

namespace ybi
{

std::unique_ptr<DisplayDevice> DisplayDevice::CreateDevice(GLFWwindow *window)
{
#ifdef WITH_VULKAN
    return std::make_unique<Vulkan>(window);
#else
#error Only vulkan supported
#endif
}

}
