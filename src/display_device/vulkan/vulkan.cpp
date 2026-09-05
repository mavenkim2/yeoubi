#define GLFW_INCLUDE_VULKAN
#include "display_device/vulkan/vulkan.h"
#include "core/math.h"
#include <GLFW/glfw3.h>
#include <algorithm>
#include <cstdio>
#include <cstring>
#include <vector>
#include <vulkan/vulkan_core.h>

namespace ybi
{
PFN_vkCreateDebugUtilsMessengerEXT vkCreateDebugUtilsMessengerEXT;
PFN_vkDestroyDebugUtilsMessengerEXT vkDestroyDebugUtilsMessengerEXT;
#ifdef _WIN32
PFN_vkGetMemoryWin32HandleKHR vkGetMemoryWin32HandleKHR;
PFN_vkGetSemaphoreWin32HandleKHR vkGetSemaphoreWin32HandleKHR;
#endif

struct VulkanExternalBuffer : DisplayBuffer
{
    VkDevice device = VK_NULL_HANDLE;
    VkBuffer buffer = VK_NULL_HANDLE;
    VkDeviceMemory vkMemory = VK_NULL_HANDLE;

    ~VulkanExternalBuffer()
    {
#if defined(_WIN32)
        if (memory.handle)
        {
            CloseHandle(static_cast<HANDLE>(memory.handle));
            memory.handle = nullptr;
        }
#endif

        if (buffer)
        {
            vkDestroyBuffer(device, buffer, nullptr);
            buffer = VK_NULL_HANDLE;
        }

        if (vkMemory)
        {
            vkFreeMemory(device, vkMemory, nullptr);
            vkMemory = VK_NULL_HANDLE;
        }
    }
};

struct VulkanSemaphore : DisplaySemaphore
{
    VkDevice device = VK_NULL_HANDLE;
    VkSemaphore semaphore = VK_NULL_HANDLE;

    ~VulkanSemaphore()
    {
#if defined(_WIN32)
        if (handle.handle)
        {
            CloseHandle(static_cast<HANDLE>(handle.handle));
            handle.handle = nullptr;
        }
#endif

        if (semaphore)
        {
            vkDestroySemaphore(device, semaphore, nullptr);
            semaphore = VK_NULL_HANDLE;
        }
    }
};

static bool HasLayer(const std::vector<VkLayerProperties> &layers, const char *name)
{
    return std::any_of(layers.begin(), layers.end(), [name](const VkLayerProperties &layer) {
        return std::strcmp(layer.layerName, name) == 0;
    });
}

static bool HasExtension(const std::vector<VkExtensionProperties> &extensions, const char *name)
{
    return std::any_of(
        extensions.begin(), extensions.end(), [name](const VkExtensionProperties &extension) {
            return std::strcmp(extension.extensionName, name) == 0;
        });
}

static VKAPI_ATTR VkBool32 VKAPI_CALL
DebugUtilsMessengerCallback(VkDebugUtilsMessageSeverityFlagBitsEXT messageSeverity,
                            VkDebugUtilsMessageTypeFlagsEXT messageTypes,
                            const VkDebugUtilsMessengerCallbackDataEXT *callbackData,
                            void *userData)
{
    (void)messageSeverity;
    (void)messageTypes;
    (void)userData;

    fprintf(stderr, "[Vulkan] %s\n", callbackData->pMessage);
    return VK_FALSE;
}

static bool DeviceSupportsExtensions(VkPhysicalDevice device,
                                     const std::vector<const char *> &required)
{
    uint32_t extensionCount = 0;
    VK_CHECK(vkEnumerateDeviceExtensionProperties(device, nullptr, &extensionCount, nullptr));

    std::vector<VkExtensionProperties> extensions(extensionCount);
    VK_CHECK(
        vkEnumerateDeviceExtensionProperties(device, nullptr, &extensionCount, extensions.data()));

    for (const char *extension : required)
    {
        if (!HasExtension(extensions, extension))
        {
            return false;
        }
    }
    return true;
}

static uint32_t FindGraphicsQueueFamily(VkPhysicalDevice device, VkSurfaceKHR surface)
{
    uint32_t queueFamilyCount = 0;
    vkGetPhysicalDeviceQueueFamilyProperties(device, &queueFamilyCount, nullptr);

    std::vector<VkQueueFamilyProperties> queueFamilies(queueFamilyCount);
    vkGetPhysicalDeviceQueueFamilyProperties(device, &queueFamilyCount, queueFamilies.data());

    for (uint32_t index = 0; index < queueFamilyCount; ++index)
    {
        const VkQueueFamilyProperties &queueFamily = queueFamilies[index];
        if (queueFamily.queueCount == 0)
        {
            continue;
        }

        if ((queueFamily.queueFlags & VK_QUEUE_GRAPHICS_BIT) == 0)
        {
            continue;
        }

        VkBool32 surfaceSupported = VK_FALSE;
        VK_CHECK(vkGetPhysicalDeviceSurfaceSupportKHR(device, index, surface, &surfaceSupported));
        if (!surfaceSupported)
        {
            continue;
        }

#if defined(_WIN32)
        if (!vkGetPhysicalDeviceWin32PresentationSupportKHR(device, index))
        {
            continue;
        }
#endif

        return index;
    }

    return VK_QUEUE_FAMILY_IGNORED;
}

static int DevicePreferenceScore(VkPhysicalDeviceType type, GPUDevicePreference preference)
{
    if (preference == GPUDevicePreference::Discrete)
    {
        if (type == VK_PHYSICAL_DEVICE_TYPE_DISCRETE_GPU)
        {
            return 2;
        }
        if (type == VK_PHYSICAL_DEVICE_TYPE_INTEGRATED_GPU)
        {
            return 1;
        }
    }
    else
    {
        if (type == VK_PHYSICAL_DEVICE_TYPE_INTEGRATED_GPU)
        {
            return 2;
        }
        if (type == VK_PHYSICAL_DEVICE_TYPE_DISCRETE_GPU)
        {
            return 1;
        }
    }

    return type == VK_PHYSICAL_DEVICE_TYPE_CPU ? -1 : 0;
}

static uint32_t FindMemoryType(VkPhysicalDevice physicalDevice,
                               uint32_t memoryTypeBits,
                               VkMemoryPropertyFlags requiredFlags)
{
    VkPhysicalDeviceMemoryProperties memoryProperties = {};
    vkGetPhysicalDeviceMemoryProperties(physicalDevice, &memoryProperties);

    for (uint32_t index = 0; index < memoryProperties.memoryTypeCount; ++index)
    {
        const bool typeMatches = (memoryTypeBits & (1u << index)) != 0;
        const bool flagsMatch =
            (memoryProperties.memoryTypes[index].propertyFlags & requiredFlags) == requiredFlags;
        if (typeMatches && flagsMatch)
        {
            return index;
        }
    }

    return UINT32_MAX;
}

static VkSurfaceFormatKHR
ChooseSwapchainSurfaceFormat(const std::vector<VkSurfaceFormatKHR> &formats)
{
    YBI_ASSERT(!formats.empty());

    for (const VkSurfaceFormatKHR &format : formats)
    {
        if (format.format == VK_FORMAT_R8G8B8A8_UNORM &&
            format.colorSpace == VK_COLOR_SPACE_SRGB_NONLINEAR_KHR)
        {
            return format;
        }
    }

    return formats[0];
}

static VkPresentModeKHR
ChooseSwapchainPresentMode(const std::vector<VkPresentModeKHR> &presentModes)
{
    for (VkPresentModeKHR presentMode : presentModes)
    {
        if (presentMode == VK_PRESENT_MODE_MAILBOX_KHR)
        {
            return presentMode;
        }
    }

    return VK_PRESENT_MODE_FIFO_KHR;
}

static VkExtent2D ChooseSwapchainExtent(GLFWwindow *window,
                                        const VkSurfaceCapabilitiesKHR &capabilities)
{
    if (capabilities.currentExtent.width != UINT32_MAX)
    {
        return capabilities.currentExtent;
    }

    int framebufferWidth = 0;
    int framebufferHeight = 0;
    glfwGetFramebufferSize(window, &framebufferWidth, &framebufferHeight);

    VkExtent2D extent = {};
    extent.width = Clamp(static_cast<uint32_t>(framebufferWidth),
                         capabilities.minImageExtent.width,
                         capabilities.maxImageExtent.width);
    extent.height = Clamp(static_cast<uint32_t>(framebufferHeight),
                          capabilities.minImageExtent.height,
                          capabilities.maxImageExtent.height);
    return extent;
}

Vulkan::Vulkan(GLFWwindow *window, ValidationMode validationMode, GPUDevicePreference preference)
    : window(window)
{
    YBI_ASSERT(window);

    const uint32_t apiVersion = VK_API_VERSION_1_2;

    VkApplicationInfo appInfo = {};
    appInfo.sType = VK_STRUCTURE_TYPE_APPLICATION_INFO;
    appInfo.pApplicationName = "Yeoubi";
    appInfo.applicationVersion = VK_MAKE_API_VERSION(0, 1, 0, 0);
    appInfo.pEngineName = "Yeoubi";
    appInfo.engineVersion = VK_MAKE_API_VERSION(0, 0, 0, 1);
    appInfo.apiVersion = apiVersion;

    uint32_t layerCount = 0;
    VK_CHECK(vkEnumerateInstanceLayerProperties(&layerCount, nullptr));

    std::vector<VkLayerProperties> availableLayers(layerCount);
    VK_CHECK(vkEnumerateInstanceLayerProperties(&layerCount, availableLayers.data()));

    uint32_t instanceExtensionCount = 0;
    VK_CHECK(vkEnumerateInstanceExtensionProperties(nullptr, &instanceExtensionCount, nullptr));

    std::vector<VkExtensionProperties> availableInstanceExtensions(instanceExtensionCount);
    VK_CHECK(vkEnumerateInstanceExtensionProperties(
        nullptr, &instanceExtensionCount, availableInstanceExtensions.data()));

    std::vector<const char *> instanceLayers;
    if (validationMode != ValidationMode::Disabled &&
        HasLayer(availableLayers, "VK_LAYER_KHRONOS_validation"))
    {
        instanceLayers.push_back("VK_LAYER_KHRONOS_validation");
    }

    // GLFW
    YBI_ASSERT(glfwVulkanSupported());
    uint32_t glfwExtensionCount = 0;
    const char **glfwExtensions = glfwGetRequiredInstanceExtensions(&glfwExtensionCount);
    YBI_ASSERT(glfwExtensions);

    std::vector<const char *> instanceExtensions;
    for (uint32_t i = 0; i < glfwExtensionCount; ++i)
    {
        instanceExtensions.push_back(glfwExtensions[i]);
    }

    debugUtilsEnabled =
        validationMode != ValidationMode::Disabled &&
        HasExtension(availableInstanceExtensions, VK_EXT_DEBUG_UTILS_EXTENSION_NAME);
    if (debugUtilsEnabled)
    {
        instanceExtensions.push_back(VK_EXT_DEBUG_UTILS_EXTENSION_NAME);
    }

    VkDebugUtilsMessengerCreateInfoEXT debugCreateInfo = {};
    debugCreateInfo.sType = VK_STRUCTURE_TYPE_DEBUG_UTILS_MESSENGER_CREATE_INFO_EXT;
    debugCreateInfo.messageSeverity = VK_DEBUG_UTILS_MESSAGE_SEVERITY_ERROR_BIT_EXT |
                                      VK_DEBUG_UTILS_MESSAGE_SEVERITY_WARNING_BIT_EXT;
    if (validationMode == ValidationMode::Verbose)
    {
        debugCreateInfo.messageSeverity |= VK_DEBUG_UTILS_MESSAGE_SEVERITY_INFO_BIT_EXT |
                                           VK_DEBUG_UTILS_MESSAGE_SEVERITY_VERBOSE_BIT_EXT;
    }
    debugCreateInfo.messageType = VK_DEBUG_UTILS_MESSAGE_TYPE_GENERAL_BIT_EXT |
                                  VK_DEBUG_UTILS_MESSAGE_TYPE_VALIDATION_BIT_EXT |
                                  VK_DEBUG_UTILS_MESSAGE_TYPE_PERFORMANCE_BIT_EXT;
    debugCreateInfo.pfnUserCallback = DebugUtilsMessengerCallback;

    VkInstanceCreateInfo instanceCreateInfo = {};
    instanceCreateInfo.sType = VK_STRUCTURE_TYPE_INSTANCE_CREATE_INFO;
    instanceCreateInfo.pApplicationInfo = &appInfo;
    instanceCreateInfo.enabledLayerCount = static_cast<uint32_t>(instanceLayers.size());
    instanceCreateInfo.ppEnabledLayerNames = instanceLayers.data();
    instanceCreateInfo.enabledExtensionCount = static_cast<uint32_t>(instanceExtensions.size());
    instanceCreateInfo.ppEnabledExtensionNames = instanceExtensions.data();
    instanceCreateInfo.pNext = debugUtilsEnabled ? &debugCreateInfo : nullptr;

    VK_CHECK(vkCreateInstance(&instanceCreateInfo, nullptr, &instance));

    // Create window surface
    VK_CHECK(glfwCreateWindowSurface(instance, window, nullptr, &windowSurface));

    // Load instance functions
    vkCreateDebugUtilsMessengerEXT = reinterpret_cast<PFN_vkCreateDebugUtilsMessengerEXT>(
        vkGetInstanceProcAddr(instance, "vkCreateDebugUtilsMessengerEXT"));
    vkDestroyDebugUtilsMessengerEXT = reinterpret_cast<PFN_vkDestroyDebugUtilsMessengerEXT>(
        vkGetInstanceProcAddr(instance, "vkDestroyDebugUtilsMessengerEXT"));

    if (debugUtilsEnabled)
    {
        YBI_ASSERT(vkCreateDebugUtilsMessengerEXT);
        YBI_ASSERT(vkDestroyDebugUtilsMessengerEXT);
        VK_CHECK(
            vkCreateDebugUtilsMessengerEXT(instance, &debugCreateInfo, nullptr, &debugMessenger));
    }

    uint32_t deviceCount = 0;
    VK_CHECK(vkEnumeratePhysicalDevices(instance, &deviceCount, nullptr));
    YBI_ASSERT(deviceCount > 0);

    std::vector<VkPhysicalDevice> devices(deviceCount);
    VK_CHECK(vkEnumeratePhysicalDevices(instance, &deviceCount, devices.data()));

    const std::vector<const char *> deviceExtensions = {
        VK_KHR_SWAPCHAIN_EXTENSION_NAME,
        VK_KHR_EXTERNAL_MEMORY_EXTENSION_NAME,
        VK_KHR_EXTERNAL_SEMAPHORE_EXTENSION_NAME,
#if defined(_WIN32)
        VK_KHR_EXTERNAL_MEMORY_WIN32_EXTENSION_NAME,
        VK_KHR_EXTERNAL_SEMAPHORE_WIN32_EXTENSION_NAME,
#endif
    };

    int bestScore = -1;
    for (VkPhysicalDevice candidate : devices)
    {
        VkPhysicalDeviceProperties properties = {};
        vkGetPhysicalDeviceProperties(candidate, &properties);

        const uint32_t queueFamily = FindGraphicsQueueFamily(candidate, windowSurface);
        if (queueFamily == VK_QUEUE_FAMILY_IGNORED)
        {
            continue;
        }

        if (!DeviceSupportsExtensions(candidate, deviceExtensions))
        {
            continue;
        }

        const int score = DevicePreferenceScore(properties.deviceType, preference);
        if (score > bestScore)
        {
            bestScore = score;
            physicalDevice = candidate;
            physicalDeviceProperties = properties;
            graphicsQueueFamily = queueFamily;
        }
    }

    YBI_ASSERT(physicalDevice != VK_NULL_HANDLE);
    YBI_ASSERT(graphicsQueueFamily != VK_QUEUE_FAMILY_IGNORED);

    const float queuePriority = 1.0f;
    VkDeviceQueueCreateInfo queueCreateInfo = {};
    queueCreateInfo.sType = VK_STRUCTURE_TYPE_DEVICE_QUEUE_CREATE_INFO;
    queueCreateInfo.queueFamilyIndex = graphicsQueueFamily;
    queueCreateInfo.queueCount = 1;
    queueCreateInfo.pQueuePriorities = &queuePriority;

    VkPhysicalDeviceFeatures deviceFeatures = {};

    VkDeviceCreateInfo deviceCreateInfo = {};
    deviceCreateInfo.sType = VK_STRUCTURE_TYPE_DEVICE_CREATE_INFO;
    deviceCreateInfo.queueCreateInfoCount = 1;
    deviceCreateInfo.pQueueCreateInfos = &queueCreateInfo;
    deviceCreateInfo.pEnabledFeatures = &deviceFeatures;
    deviceCreateInfo.enabledExtensionCount = static_cast<uint32_t>(deviceExtensions.size());
    deviceCreateInfo.ppEnabledExtensionNames = deviceExtensions.data();

    VK_CHECK(vkCreateDevice(physicalDevice, &deviceCreateInfo, nullptr, &device));
    vkGetDeviceQueue(device, graphicsQueueFamily, 0, &graphicsQueue);

    // Load device functions
#if defined(_WIN32)
    vkGetMemoryWin32HandleKHR = reinterpret_cast<PFN_vkGetMemoryWin32HandleKHR>(
        vkGetDeviceProcAddr(device, "vkGetMemoryWin32HandleKHR"));
    YBI_ASSERT(vkGetMemoryWin32HandleKHR);
    vkGetSemaphoreWin32HandleKHR = reinterpret_cast<PFN_vkGetSemaphoreWin32HandleKHR>(
        vkGetDeviceProcAddr(device, "vkGetSemaphoreWin32HandleKHR"));
    YBI_ASSERT(vkGetSemaphoreWin32HandleKHR);
#endif

    VkSurfaceCapabilitiesKHR surfaceCapabilities = {};
    VK_CHECK(vkGetPhysicalDeviceSurfaceCapabilitiesKHR(
        physicalDevice, windowSurface, &surfaceCapabilities));

    uint32_t surfaceFormatCount = 0;
    VK_CHECK(vkGetPhysicalDeviceSurfaceFormatsKHR(
        physicalDevice, windowSurface, &surfaceFormatCount, nullptr));
    YBI_ASSERT(surfaceFormatCount > 0);
    std::vector<VkSurfaceFormatKHR> surfaceFormats(surfaceFormatCount);
    VK_CHECK(vkGetPhysicalDeviceSurfaceFormatsKHR(
        physicalDevice, windowSurface, &surfaceFormatCount, surfaceFormats.data()));

    uint32_t presentModeCount = 0;
    VK_CHECK(vkGetPhysicalDeviceSurfacePresentModesKHR(
        physicalDevice, windowSurface, &presentModeCount, nullptr));
    YBI_ASSERT(presentModeCount > 0);
    std::vector<VkPresentModeKHR> presentModes(presentModeCount);
    VK_CHECK(vkGetPhysicalDeviceSurfacePresentModesKHR(
        physicalDevice, windowSurface, &presentModeCount, presentModes.data()));

    const VkSurfaceFormatKHR surfaceFormat = ChooseSwapchainSurfaceFormat(surfaceFormats);
    const VkPresentModeKHR presentMode = ChooseSwapchainPresentMode(presentModes);
    swapchainExtent = ChooseSwapchainExtent(window, surfaceCapabilities);
    swapchainImageFormat = surfaceFormat.format;

    uint32_t imageCount = surfaceCapabilities.minImageCount + 1;
    if (surfaceCapabilities.maxImageCount > 0 && imageCount > surfaceCapabilities.maxImageCount)
    {
        imageCount = surfaceCapabilities.maxImageCount;
    }

    VkSwapchainCreateInfoKHR swapchainInfo = {};
    swapchainInfo.sType = VK_STRUCTURE_TYPE_SWAPCHAIN_CREATE_INFO_KHR;
    swapchainInfo.surface = windowSurface;
    swapchainInfo.minImageCount = imageCount;
    swapchainInfo.imageFormat = swapchainImageFormat;
    swapchainInfo.imageColorSpace = surfaceFormat.colorSpace;
    swapchainInfo.imageExtent = swapchainExtent;
    swapchainInfo.imageArrayLayers = 1;
    swapchainInfo.imageUsage =
        VK_IMAGE_USAGE_TRANSFER_DST_BIT | VK_IMAGE_USAGE_COLOR_ATTACHMENT_BIT;
    swapchainInfo.imageSharingMode = VK_SHARING_MODE_EXCLUSIVE;
    swapchainInfo.preTransform = surfaceCapabilities.currentTransform;
    swapchainInfo.compositeAlpha = VK_COMPOSITE_ALPHA_OPAQUE_BIT_KHR;
    swapchainInfo.presentMode = presentMode;
    swapchainInfo.clipped = VK_TRUE;
    swapchainInfo.oldSwapchain = VK_NULL_HANDLE;

    VK_CHECK(vkCreateSwapchainKHR(device, &swapchainInfo, nullptr, &swapchain));

    uint32_t swapchainImageCount = 0;
    VK_CHECK(vkGetSwapchainImagesKHR(device, swapchain, &swapchainImageCount, nullptr));
    swapchainImages.resize(swapchainImageCount);
    VK_CHECK(
        vkGetSwapchainImagesKHR(device, swapchain, &swapchainImageCount, swapchainImages.data()));
    swapchainImageLayouts.resize(swapchainImageCount, VK_IMAGE_LAYOUT_UNDEFINED);

    VkCommandPoolCreateInfo commandPoolInfo = {};
    commandPoolInfo.sType = VK_STRUCTURE_TYPE_COMMAND_POOL_CREATE_INFO;
    commandPoolInfo.flags = VK_COMMAND_POOL_CREATE_RESET_COMMAND_BUFFER_BIT;
    commandPoolInfo.queueFamilyIndex = graphicsQueueFamily;
    VK_CHECK(vkCreateCommandPool(device, &commandPoolInfo, nullptr, &commandPool));

    commandBuffers.resize(numFramesInFlight);
    VkCommandBufferAllocateInfo commandBufferInfo = {};
    commandBufferInfo.sType = VK_STRUCTURE_TYPE_COMMAND_BUFFER_ALLOCATE_INFO;
    commandBufferInfo.commandPool = commandPool;
    commandBufferInfo.level = VK_COMMAND_BUFFER_LEVEL_PRIMARY;
    commandBufferInfo.commandBufferCount = static_cast<uint32_t>(commandBuffers.size());
    VK_CHECK(vkAllocateCommandBuffers(device, &commandBufferInfo, commandBuffers.data()));

    imageAvailableSemaphores.resize(numFramesInFlight);
    renderFinishedSemaphores.resize(swapchainImages.size());
    inFlightFences.resize(numFramesInFlight);

    VkSemaphoreCreateInfo semaphoreInfo = {};
    semaphoreInfo.sType = VK_STRUCTURE_TYPE_SEMAPHORE_CREATE_INFO;

    VkFenceCreateInfo fenceInfo = {};
    fenceInfo.sType = VK_STRUCTURE_TYPE_FENCE_CREATE_INFO;
    fenceInfo.flags = VK_FENCE_CREATE_SIGNALED_BIT;

    for (uint32_t index = 0; index < numFramesInFlight; ++index)
    {
        VK_CHECK(
            vkCreateSemaphore(device, &semaphoreInfo, nullptr, &imageAvailableSemaphores[index]));
        VK_CHECK(vkCreateFence(device, &fenceInfo, nullptr, &inFlightFences[index]));
    }

    for (VkSemaphore &semaphore : renderFinishedSemaphores)
    {
        VK_CHECK(vkCreateSemaphore(device, &semaphoreInfo, nullptr, &semaphore));
    }
}

void Vulkan::WaitIdle()
{
    if (device)
        VK_CHECK(vkDeviceWaitIdle(device));
}

Vulkan::~Vulkan()
{
    if (device)
    {
        WaitIdle();
        for (VkFence fence : inFlightFences)
        {
            if (fence)
            {
                vkDestroyFence(device, fence, nullptr);
            }
        }
        inFlightFences.clear();

        for (VkSemaphore semaphore : imageAvailableSemaphores)
        {
            if (semaphore)
            {
                vkDestroySemaphore(device, semaphore, nullptr);
            }
        }
        imageAvailableSemaphores.clear();

        for (VkSemaphore semaphore : renderFinishedSemaphores)
        {
            if (semaphore)
            {
                vkDestroySemaphore(device, semaphore, nullptr);
            }
        }
        renderFinishedSemaphores.clear();

        if (commandPool)
        {
            vkDestroyCommandPool(device, commandPool, nullptr);
            commandPool = VK_NULL_HANDLE;
            commandBuffers.clear();
        }

        if (swapchain)
        {
            vkDestroySwapchainKHR(device, swapchain, nullptr);
            swapchain = VK_NULL_HANDLE;
            swapchainImages.clear();
            swapchainImageLayouts.clear();
        }
        vkDestroyDevice(device, nullptr);
        device = VK_NULL_HANDLE;
    }

    if (windowSurface)
    {
        vkDestroySurfaceKHR(instance, windowSurface, nullptr);
        windowSurface = VK_NULL_HANDLE;
    }

    if (debugMessenger)
    {
        YBI_ASSERT(vkDestroyDebugUtilsMessengerEXT);
        vkDestroyDebugUtilsMessengerEXT(instance, debugMessenger, nullptr);
        debugMessenger = VK_NULL_HANDLE;
    }

    if (instance)
    {
        vkDestroyInstance(instance, nullptr);
        instance = VK_NULL_HANDLE;
    }
}

std::unique_ptr<DisplayBuffer> Vulkan::AllocateExternalBuffer(size_t size)
{
    YBI_ASSERT(size > 0);

    const VkBufferUsageFlags usage =
        VK_BUFFER_USAGE_TRANSFER_SRC_BIT | VK_BUFFER_USAGE_STORAGE_BUFFER_BIT;

    auto result = std::make_unique<VulkanExternalBuffer>();
    result->device = device;
    result->size = size;

    VkExternalMemoryBufferCreateInfo externalBufferInfo = {};
    externalBufferInfo.sType = VK_STRUCTURE_TYPE_EXTERNAL_MEMORY_BUFFER_CREATE_INFO;
#if defined(_WIN32)
    externalBufferInfo.handleTypes = VK_EXTERNAL_MEMORY_HANDLE_TYPE_OPAQUE_WIN32_BIT;
#else
#error External Vulkan memory export is only implemented for Win32 right now.
#endif

    VkBufferCreateInfo bufferInfo = {};
    bufferInfo.sType = VK_STRUCTURE_TYPE_BUFFER_CREATE_INFO;
    bufferInfo.pNext = &externalBufferInfo;
    bufferInfo.size = static_cast<VkDeviceSize>(size);
    bufferInfo.usage = usage;
    bufferInfo.sharingMode = VK_SHARING_MODE_EXCLUSIVE;

    VK_CHECK(vkCreateBuffer(device, &bufferInfo, nullptr, &result->buffer));

    VkMemoryDedicatedRequirements dedicatedRequirements = {};
    dedicatedRequirements.sType = VK_STRUCTURE_TYPE_MEMORY_DEDICATED_REQUIREMENTS;

    VkMemoryRequirements2 memoryRequirements = {};
    memoryRequirements.sType = VK_STRUCTURE_TYPE_MEMORY_REQUIREMENTS_2;
    memoryRequirements.pNext = &dedicatedRequirements;

    VkBufferMemoryRequirementsInfo2 memoryRequirementsInfo = {};
    memoryRequirementsInfo.sType = VK_STRUCTURE_TYPE_BUFFER_MEMORY_REQUIREMENTS_INFO_2;
    memoryRequirementsInfo.buffer = result->buffer;
    vkGetBufferMemoryRequirements2(device, &memoryRequirementsInfo, &memoryRequirements);

    const uint32_t memoryTypeIndex =
        FindMemoryType(physicalDevice,
                       memoryRequirements.memoryRequirements.memoryTypeBits,
                       VK_MEMORY_PROPERTY_DEVICE_LOCAL_BIT);
    YBI_ASSERT(memoryTypeIndex != UINT32_MAX);

    VkMemoryDedicatedAllocateInfo dedicatedAllocateInfo = {};
    dedicatedAllocateInfo.sType = VK_STRUCTURE_TYPE_MEMORY_DEDICATED_ALLOCATE_INFO;
    dedicatedAllocateInfo.buffer = result->buffer;

    VkExportMemoryAllocateInfo exportInfo = {};
    exportInfo.sType = VK_STRUCTURE_TYPE_EXPORT_MEMORY_ALLOCATE_INFO;
    exportInfo.pNext = &dedicatedAllocateInfo;
#if defined(_WIN32)
    exportInfo.handleTypes = VK_EXTERNAL_MEMORY_HANDLE_TYPE_OPAQUE_WIN32_BIT;
#endif

    VkMemoryAllocateInfo allocateInfo = {};
    allocateInfo.sType = VK_STRUCTURE_TYPE_MEMORY_ALLOCATE_INFO;
    allocateInfo.pNext = &exportInfo;
    allocateInfo.allocationSize = memoryRequirements.memoryRequirements.size;
    allocateInfo.memoryTypeIndex = memoryTypeIndex;

    VK_CHECK(vkAllocateMemory(device, &allocateInfo, nullptr, &result->vkMemory));
    VK_CHECK(vkBindBufferMemory(device, result->buffer, result->vkMemory, 0));

#if defined(_WIN32)
    HANDLE handle = nullptr;
    VkMemoryGetWin32HandleInfoKHR handleInfo = {};
    handleInfo.sType = VK_STRUCTURE_TYPE_MEMORY_GET_WIN32_HANDLE_INFO_KHR;
    handleInfo.memory = result->vkMemory;
    handleInfo.handleType = VK_EXTERNAL_MEMORY_HANDLE_TYPE_OPAQUE_WIN32_BIT;
    YBI_ASSERT(vkGetMemoryWin32HandleKHR);
    VK_CHECK(vkGetMemoryWin32HandleKHR(device, &handleInfo, &handle));

    result->memory.type = ExternalMemoryHandleType::OpaqueWin32;
    result->memory.handle = handle;
    result->memory.size = static_cast<size_t>(memoryRequirements.memoryRequirements.size);
#endif

    return result;
}

std::unique_ptr<DisplaySemaphore> Vulkan::CreateDisplaySemaphore()
{
    auto result = std::make_unique<VulkanSemaphore>();
    result->device = device;

    VkExportSemaphoreCreateInfo exportInfo = {};
    exportInfo.sType = VK_STRUCTURE_TYPE_EXPORT_SEMAPHORE_CREATE_INFO;
#if defined(_WIN32)
    exportInfo.handleTypes = VK_EXTERNAL_SEMAPHORE_HANDLE_TYPE_OPAQUE_WIN32_BIT;
#else
#error External Vulkan semaphore export is only implemented for Win32 right now.
#endif

    VkSemaphoreCreateInfo semaphoreInfo = {};
    semaphoreInfo.sType = VK_STRUCTURE_TYPE_SEMAPHORE_CREATE_INFO;
    semaphoreInfo.pNext = &exportInfo;
    VK_CHECK(vkCreateSemaphore(device, &semaphoreInfo, nullptr, &result->semaphore));

#if defined(_WIN32)
    HANDLE handle = nullptr;
    VkSemaphoreGetWin32HandleInfoKHR handleInfo = {};
    handleInfo.sType = VK_STRUCTURE_TYPE_SEMAPHORE_GET_WIN32_HANDLE_INFO_KHR;
    handleInfo.semaphore = result->semaphore;
    handleInfo.handleType = VK_EXTERNAL_SEMAPHORE_HANDLE_TYPE_OPAQUE_WIN32_BIT;
    YBI_ASSERT(vkGetSemaphoreWin32HandleKHR);
    VK_CHECK(vkGetSemaphoreWin32HandleKHR(device, &handleInfo, &handle));

    result->handle.type = ExternalMemoryHandleType::OpaqueWin32;
    result->handle.handle = handle;
#endif

    return result;
}

bool Vulkan::BeginFrame()
{
    YBI_ASSERT(!commandBuffers.empty());
    YBI_ASSERT(!imageAvailableSemaphores.empty());
    YBI_ASSERT(!inFlightFences.empty());

    VkFence fence = inFlightFences[currentFrame];
    VK_CHECK(vkWaitForFences(device, 1, &fence, VK_TRUE, UINT64_MAX));

    VkSemaphore imageAvailableSemaphore = imageAvailableSemaphores[currentFrame];
    const VkResult acquireResult = vkAcquireNextImageKHR(device,
                                                         swapchain,
                                                         UINT64_MAX,
                                                         imageAvailableSemaphore,
                                                         VK_NULL_HANDLE,
                                                         &currentSwapchainImageIndex);
    if (acquireResult == VK_ERROR_OUT_OF_DATE_KHR)
    {
        return false;
    }
    YBI_ASSERT(acquireResult == VK_SUCCESS || acquireResult == VK_SUBOPTIMAL_KHR);

    VK_CHECK(vkResetFences(device, 1, &fence));

    VkCommandBuffer commandBuffer = commandBuffers[currentFrame];
    VK_CHECK(vkResetCommandBuffer(commandBuffer, 0));

    VkCommandBufferBeginInfo beginInfo = {};
    beginInfo.sType = VK_STRUCTURE_TYPE_COMMAND_BUFFER_BEGIN_INFO;
    beginInfo.flags = VK_COMMAND_BUFFER_USAGE_ONE_TIME_SUBMIT_BIT;
    VK_CHECK(vkBeginCommandBuffer(commandBuffer, &beginInfo));

    return true;
}

bool Vulkan::EndFrame(DisplayBuffer &framebuffer,
                      DisplaySemaphore &cudaFinishedSemaphore,
                      int width,
                      int height)
{
    YBI_ASSERT(width > 0);
    YBI_ASSERT(height > 0);
    YBI_ASSERT(static_cast<uint32_t>(width) == swapchainExtent.width);
    YBI_ASSERT(static_cast<uint32_t>(height) == swapchainExtent.height);

    VulkanExternalBuffer &vulkanFramebuffer = static_cast<VulkanExternalBuffer &>(framebuffer);
    VulkanSemaphore &vulkanCudaFinishedSemaphore =
        static_cast<VulkanSemaphore &>(cudaFinishedSemaphore);

    VkCommandBuffer commandBuffer = commandBuffers[currentFrame];
    VkImage swapchainImage = swapchainImages[currentSwapchainImageIndex];

    VkImageMemoryBarrier transferBarrier = {};
    transferBarrier.sType = VK_STRUCTURE_TYPE_IMAGE_MEMORY_BARRIER;
    transferBarrier.srcAccessMask = 0;
    transferBarrier.dstAccessMask = VK_ACCESS_TRANSFER_WRITE_BIT;
    transferBarrier.oldLayout = swapchainImageLayouts[currentSwapchainImageIndex];
    transferBarrier.newLayout = VK_IMAGE_LAYOUT_TRANSFER_DST_OPTIMAL;
    transferBarrier.srcQueueFamilyIndex = VK_QUEUE_FAMILY_IGNORED;
    transferBarrier.dstQueueFamilyIndex = VK_QUEUE_FAMILY_IGNORED;
    transferBarrier.image = swapchainImage;
    transferBarrier.subresourceRange.aspectMask = VK_IMAGE_ASPECT_COLOR_BIT;
    transferBarrier.subresourceRange.baseMipLevel = 0;
    transferBarrier.subresourceRange.levelCount = 1;
    transferBarrier.subresourceRange.baseArrayLayer = 0;
    transferBarrier.subresourceRange.layerCount = 1;

    vkCmdPipelineBarrier(commandBuffer,
                         VK_PIPELINE_STAGE_TOP_OF_PIPE_BIT,
                         VK_PIPELINE_STAGE_TRANSFER_BIT,
                         0,
                         0,
                         nullptr,
                         0,
                         nullptr,
                         1,
                         &transferBarrier);

    VkBufferImageCopy copyRegion = {};
    copyRegion.bufferOffset = 0;
    copyRegion.bufferRowLength = static_cast<uint32_t>(width);
    copyRegion.bufferImageHeight = static_cast<uint32_t>(height);
    copyRegion.imageSubresource.aspectMask = VK_IMAGE_ASPECT_COLOR_BIT;
    copyRegion.imageSubresource.mipLevel = 0;
    copyRegion.imageSubresource.baseArrayLayer = 0;
    copyRegion.imageSubresource.layerCount = 1;
    copyRegion.imageOffset = {0, 0, 0};
    copyRegion.imageExtent = {static_cast<uint32_t>(width), static_cast<uint32_t>(height), 1};

    vkCmdCopyBufferToImage(commandBuffer,
                           vulkanFramebuffer.buffer,
                           swapchainImage,
                           VK_IMAGE_LAYOUT_TRANSFER_DST_OPTIMAL,
                           1,
                           &copyRegion);

    VkImageMemoryBarrier presentBarrier = {};
    presentBarrier.sType = VK_STRUCTURE_TYPE_IMAGE_MEMORY_BARRIER;
    presentBarrier.srcAccessMask = VK_ACCESS_TRANSFER_WRITE_BIT;
    presentBarrier.dstAccessMask = 0;
    presentBarrier.oldLayout = VK_IMAGE_LAYOUT_TRANSFER_DST_OPTIMAL;
    presentBarrier.newLayout = VK_IMAGE_LAYOUT_PRESENT_SRC_KHR;
    presentBarrier.srcQueueFamilyIndex = VK_QUEUE_FAMILY_IGNORED;
    presentBarrier.dstQueueFamilyIndex = VK_QUEUE_FAMILY_IGNORED;
    presentBarrier.image = swapchainImage;
    presentBarrier.subresourceRange.aspectMask = VK_IMAGE_ASPECT_COLOR_BIT;
    presentBarrier.subresourceRange.baseMipLevel = 0;
    presentBarrier.subresourceRange.levelCount = 1;
    presentBarrier.subresourceRange.baseArrayLayer = 0;
    presentBarrier.subresourceRange.layerCount = 1;

    vkCmdPipelineBarrier(commandBuffer,
                         VK_PIPELINE_STAGE_TRANSFER_BIT,
                         VK_PIPELINE_STAGE_BOTTOM_OF_PIPE_BIT,
                         0,
                         0,
                         nullptr,
                         0,
                         nullptr,
                         1,
                         &presentBarrier);

    VK_CHECK(vkEndCommandBuffer(commandBuffer));

    VkSemaphore waitSemaphores[] = {
        imageAvailableSemaphores[currentFrame],
        vulkanCudaFinishedSemaphore.semaphore,
    };
    VkPipelineStageFlags waitStages[] = {
        VK_PIPELINE_STAGE_TRANSFER_BIT,
        VK_PIPELINE_STAGE_TRANSFER_BIT,
    };
    VkSemaphore signalSemaphore = renderFinishedSemaphores[currentSwapchainImageIndex];

    VkSubmitInfo submitInfo = {};
    submitInfo.sType = VK_STRUCTURE_TYPE_SUBMIT_INFO;
    submitInfo.waitSemaphoreCount = 2;
    submitInfo.pWaitSemaphores = waitSemaphores;
    submitInfo.pWaitDstStageMask = waitStages;
    submitInfo.commandBufferCount = 1;
    submitInfo.pCommandBuffers = &commandBuffer;
    submitInfo.signalSemaphoreCount = 1;
    submitInfo.pSignalSemaphores = &signalSemaphore;

    VkFence fence = inFlightFences[currentFrame];
    VK_CHECK(vkQueueSubmit(graphicsQueue, 1, &submitInfo, fence));

    VkPresentInfoKHR presentInfo = {};
    presentInfo.sType = VK_STRUCTURE_TYPE_PRESENT_INFO_KHR;
    presentInfo.waitSemaphoreCount = 1;
    presentInfo.pWaitSemaphores = &signalSemaphore;
    presentInfo.swapchainCount = 1;
    presentInfo.pSwapchains = &swapchain;
    presentInfo.pImageIndices = &currentSwapchainImageIndex;

    const VkResult presentResult = vkQueuePresentKHR(graphicsQueue, &presentInfo);
    if (presentResult == VK_ERROR_OUT_OF_DATE_KHR)
    {
        return false;
    }
    YBI_ERROR(presentResult == VK_SUCCESS || presentResult == VK_SUBOPTIMAL_KHR,
              "Error: %i\n",
              presentResult);

    swapchainImageLayouts[currentSwapchainImageIndex] = VK_IMAGE_LAYOUT_PRESENT_SRC_KHR;
    currentFrame = (currentFrame + 1) % numFramesInFlight;
    frameCount++;
    return true;
}

} // namespace ybi
