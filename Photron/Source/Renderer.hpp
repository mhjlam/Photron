#pragma once

#include <Windows.h>
#include <vulkan/vulkan.hpp>


class Renderer
{
private:
    vk::Instance mInstance;
    vk::SurfaceKHR mSurface;
    vk::Device mDevice;
    vk::Queue mQueue;
    vk::SwapchainKHR mSwapchain;
    vk::RenderPass mRenderPass;
    vk::PipelineLayout mPipelineLayout;
    vk::Pipeline mPipeline;
    vk::CommandPool mCommandPool;
    vk::CommandBuffer mCommandBuffer;
    vk::Framebuffer mFramebuffer;

    vk::ClearValue mClearColor;

    uint32_t mWidth;
    uint32_t mHeight;

public:
    Renderer(HWND hwnd, HINSTANCE hInstance);
    ~Renderer();

    //void Initialize();

    void Update(uint32_t width, uint32_t height);
    void Render();

private:
    void CreateInstance();

    std::vector<std::string> ValidationLayers();
    std::vector<std::string> InstanceExtensions();
};
