#pragma once

#include "core/assert.h"
#include "display_device/display_device.h"
#include <cstddef>
#include <cstdint>
#include <vector>
#if defined(_WIN32) && !defined(VK_USE_PLATFORM_WIN32_KHR)
#define VK_USE_PLATFORM_WIN32_KHR
#endif
#include <vulkan/vulkan.h>

struct GLFWwindow;

namespace ybi
{

extern PFN_vkCreateDebugUtilsMessengerEXT vkCreateDebugUtilsMessengerEXT;
extern PFN_vkDestroyDebugUtilsMessengerEXT vkDestroyDebugUtilsMessengerEXT;
#ifdef _WIN32
extern PFN_vkGetMemoryWin32HandleKHR vkGetMemoryWin32HandleKHR;
extern PFN_vkGetSemaphoreWin32HandleKHR vkGetSemaphoreWin32HandleKHR;
#endif

#define VK_CHECK(check) \
    do \
    { \
        VkResult result_ = (check); \
        YBI_ASSERT(result_ == VK_SUCCESS); \
    } while (0)

enum class ValidationMode
{
    Disabled,
    Enabled,
    Verbose,
};

enum class GPUDevicePreference
{
    Discrete,
    Integrated,
};

struct Vulkan : DisplayDevice
{
    GLFWwindow *window = nullptr;
    VkInstance instance = VK_NULL_HANDLE;
    VkDebugUtilsMessengerEXT debugMessenger = VK_NULL_HANDLE;
    VkPhysicalDevice physicalDevice = VK_NULL_HANDLE;
    VkPhysicalDeviceProperties physicalDeviceProperties = {};
    VkDevice device = VK_NULL_HANDLE;
    VkQueue graphicsQueue = VK_NULL_HANDLE;
    uint32_t graphicsQueueFamily = VK_QUEUE_FAMILY_IGNORED;
    VkSurfaceKHR windowSurface = VK_NULL_HANDLE;
    VkSwapchainKHR swapchain = VK_NULL_HANDLE;
    VkFormat swapchainImageFormat = VK_FORMAT_UNDEFINED;
    VkExtent2D swapchainExtent = {};
    std::vector<VkImage> swapchainImages;
    std::vector<VkImageLayout> swapchainImageLayouts;
    VkCommandPool commandPool = VK_NULL_HANDLE;
    std::vector<VkCommandBuffer> commandBuffers;
    std::vector<VkSemaphore> imageAvailableSemaphores;
    std::vector<VkSemaphore> renderFinishedSemaphores;
    std::vector<VkFence> inFlightFences;
    uint32_t currentSwapchainImageIndex = 0;
    bool debugUtilsEnabled = false;

    Vulkan(GLFWwindow *window,
           ValidationMode validationMode = ValidationMode::Enabled,
           GPUDevicePreference preference = GPUDevicePreference::Discrete);
    ~Vulkan();

    Vulkan(const Vulkan &) = delete;
    Vulkan &operator=(const Vulkan &) = delete;

    std::unique_ptr<DisplayBuffer> AllocateExternalBuffer(size_t size) override;
    std::unique_ptr<DisplaySemaphore> CreateDisplaySemaphore() override;
    bool BeginFrame() override;
    bool EndFrame(DisplayBuffer &framebuffer,
                  DisplaySemaphore &cudaFinishedSemaphore,
                  int width,
                  int height) override;
};

} // namespace ybi
