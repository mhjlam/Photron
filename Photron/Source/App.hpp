#pragma once

#include <string>
#include <memory>

#include <Windows.h>
#include <vulkan/vulkan.hpp>

class Renderer;

class App
{
public:
    static LRESULT CALLBACK WindowProc(HWND hwnd, UINT uMsg, WPARAM wParam, LPARAM lParam);

public:
    App(HINSTANCE hInstance);
    ~App();

    void Run() const;

private:
    HWND mHwnd;
    HINSTANCE mHinstance;

    std::shared_ptr<Renderer> mRenderer;

private:
    void SetupConsole(std::wstring title);
};
