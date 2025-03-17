#include "Renderer.hpp"

#include <iostream>

#include <vulkan/vulkan_win32.h>


Renderer::Renderer(HWND hwnd, HINSTANCE hInstance)
{
    CreateInstance();

    //// Create surface
    //VkWin32SurfaceCreateInfoKHR win32SurfaceCreateInfo{};
    //win32SurfaceCreateInfo.sType = VK_STRUCTURE_TYPE_WIN32_SURFACE_CREATE_INFO_KHR;
    //win32SurfaceCreateInfo.hinstance = hInstance;
    //win32SurfaceCreateInfo.hwnd = hwnd;

    //VkSurfaceKHR surface{};
    //VkResult result = vkCreateWin32SurfaceKHR(mInstance, &win32SurfaceCreateInfo, nullptr, &surface);
    //if (result != VK_SUCCESS) {
    //    throw std::runtime_error("Failed to create Win32 surface: " + std::to_string(result));
    //}
    //mSurface = surface;

    //// Create device
    //vk::PhysicalDevice physicalDevice = mInstance.enumeratePhysicalDevices().front();
    //float queuePriority = 1.0f;
    //vk::DeviceQueueCreateInfo queueCreateInfo{ {}, 0, 1, &queuePriority };
    //vk::DeviceCreateInfo deviceCreateInfo{ {}, queueCreateInfo };
    //mDevice = physicalDevice.createDevice(deviceCreateInfo);
    //mQueue = mDevice.getQueue(0, 0);

    //// Create render pass
    //vk::RenderPassCreateInfo renderPassInfo{};
    //mRenderPass = mDevice.createRenderPass(renderPassInfo);

    //// Create graphics pipeline
    //vk::PipelineLayoutCreateInfo pipelineLayoutInfo{};
    //mPipelineLayout = mDevice.createPipelineLayout(pipelineLayoutInfo);

    //vk::GraphicsPipelineCreateInfo pipelineInfo{};
    //pipelineInfo.layout = mPipelineLayout;
    //mPipeline = mDevice.createGraphicsPipeline(nullptr, pipelineInfo).value;

    //// Create command buffer
    //vk::CommandPoolCreateInfo poolInfo{};
    //mCommandPool = mDevice.createCommandPool(poolInfo);

    //vk::CommandBufferAllocateInfo allocInfo{ mCommandPool, vk::CommandBufferLevel::ePrimary, 1 };
    //mCommandBuffer = mDevice.allocateCommandBuffers(allocInfo).front();

    //// Create clear color
    //mClearColor = vk::ClearValue({0.0f, 0.0f, 0.0f, 1.0f});
}


Renderer::~Renderer()
{
    mDevice.destroyPipeline(mPipeline);
    mDevice.destroyPipelineLayout(mPipelineLayout);
    mDevice.destroyRenderPass(mRenderPass);
    mDevice.destroyCommandPool(mCommandPool);
    mDevice.destroy();
    mInstance.destroySurfaceKHR(mSurface);
    mInstance.destroy();
}


void Renderer::Update(uint32_t width, uint32_t height)
{
    mWidth = width;
    mHeight = height;
}


void Renderer::Render()
{
    //vk::CommandBufferBeginInfo beginInfo{};
    //mCommandBuffer.begin(beginInfo);

    //vk::RenderPassBeginInfo renderPassInfo{};
    //renderPassInfo.renderPass = mRenderPass;
    //renderPassInfo.framebuffer = mFramebuffer;
    //renderPassInfo.renderArea.extent.width = mWidth;
    //renderPassInfo.renderArea.extent.height = mHeight;

    //renderPassInfo.clearValueCount = 1;
    //renderPassInfo.pClearValues = &mClearColor;

    //mCommandBuffer.beginRenderPass(renderPassInfo, vk::SubpassContents::eInline);
    //mCommandBuffer.bindPipeline(vk::PipelineBindPoint::eGraphics, mPipeline);
    //mCommandBuffer.draw(3, 1, 0, 0);
    //mCommandBuffer.endRenderPass();
    //mCommandBuffer.end();

    //vk::SubmitInfo submitInfo{};
    //submitInfo.commandBufferCount = 1;
    //submitInfo.pCommandBuffers = &mCommandBuffer;
    //mQueue.submit(submitInfo);
    //mQueue.waitIdle();
}


void Renderer::CreateInstance()
{
    // Set application info
    vk::ApplicationInfo appInfo {
        "Hello Triangle", VK_MAKE_VERSION(1,0,0), "No Engine", VK_MAKE_VERSION(1,0,0), VK_API_VERSION_1_0
    };

    // Enable validation layers in debug mode
#ifdef NDEBUG
    const bool enableValidationLayers = false;
#else
    const bool enableValidationLayers = true;
#endif

    // Determine validation layers
    auto layers = (enableValidationLayers) ? ValidationLayers() : std::vector<std::string>{};
    std::vector<const char*> enabledLayerNames{ layers.size() };
    std::transform(layers.begin(), layers.end(), enabledLayerNames.begin(), [](const std::string& str) {
        return str.c_str();
    });

    // Determine extensions
    auto extensions = InstanceExtensions();
    std::vector<const char*> enabledExtensionNames{ extensions.size() };
    std::transform(extensions.begin(), extensions.end(), enabledExtensionNames.begin(), [](const std::string& str) {
        return str.c_str();
    });

    // Set instance info
    vk::InstanceCreateInfo instanceCreateInfo {
        {}, &appInfo, enabledLayerNames, enabledExtensionNames
    };

    // Create instance
    mInstance = vk::createInstance(instanceCreateInfo);
    if (!mInstance) {
        throw std::runtime_error("Failed to create Vulkan instance.");
    }
}


std::vector<std::string> Renderer::ValidationLayers()
{
    const std::vector<std::string> requiredLayers = {
        "VK_LAYER_KHRONOS_validation"
    };

    auto availableLayers = vk::enumerateInstanceLayerProperties();

    for (std::string layerName : requiredLayers) {
        bool layerFound = false;

        for (const auto& layerProperties : availableLayers) {
            if (layerName == layerProperties.layerName) {
                layerFound = true;
                break;
            }
        }

        if (!layerFound) {
            throw std::runtime_error("Unable to find required validation layer " + std::string(layerName));
        }
    }

    return requiredLayers;
}


std::vector<std::string> Renderer::InstanceExtensions()
{
    std::vector<std::string> extensionNames;

    // Query available extensions
    auto availableExtensions = vk::enumerateInstanceExtensionProperties();

    // Check if required extensions are available
    for (const auto& extension : availableExtensions) {
        if (std::strcmp(extension.extensionName, VK_KHR_SURFACE_EXTENSION_NAME) == 0 ||
            std::strcmp(extension.extensionName, VK_KHR_WIN32_SURFACE_EXTENSION_NAME) == 0) {
            extensionNames.push_back(extension.extensionName);
        }
    }

    return extensionNames;
}
