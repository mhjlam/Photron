#pragma once

#include <iostream>
#include <filesystem>

#include "App.hpp"
#include "Renderer.hpp"


constexpr uint32_t WINDOW_WIDTH = 800;
constexpr uint32_t WINDOW_HEIGHT = 600;

constexpr LPCWSTR WINDOW_NAME = L"Photran";
constexpr LPCWSTR WINDOW_CLASS = L"WINDOW_CLASS";


LRESULT CALLBACK App::WindowProc(HWND hwnd, UINT uMsg, WPARAM wParam, LPARAM lParam)
{
    if (uMsg == WM_CLOSE) {
        ::PostQuitMessage(0);
    }
    return ::DefWindowProc(hwnd, uMsg, wParam, lParam);
}


App::App(HINSTANCE hInstance) : mHwnd{}, mHinstance{hInstance}
{
    SetupConsole(WINDOW_NAME);

    // Create window class
    WNDCLASSEX wcex{};
    ::ZeroMemory(&wcex, sizeof(WNDCLASSEX));
    wcex.lpszClassName = WINDOW_CLASS;
    wcex.cbSize = sizeof(WNDCLASSEX);
    wcex.style = CS_HREDRAW | CS_VREDRAW;
    wcex.lpfnWndProc = WindowProc;
    wcex.hInstance = mHinstance;
    wcex.hIcon = ::LoadIcon(nullptr, IDI_APPLICATION);
    wcex.hIconSm = nullptr;
    wcex.hCursor = ::LoadCursor(nullptr, IDC_ARROW);
    wcex.hbrBackground = nullptr;
    wcex.lpszMenuName = nullptr;
    wcex.cbClsExtra = 0;
    wcex.cbWndExtra = 0;

    if (!::RegisterClassEx(&wcex)) {
        throw std::runtime_error("Window registration failed!");
    }

    // Determine the position of the window
    RECT taskbarRect{};
    LONG taskbarW = 0L;
    LONG taskbarH = 0L;
    
    // Find taskbar location and dimensions
    if (::GetWindowRect(::FindWindow(L"Shell_traywnd", nullptr), &taskbarRect)) {
        taskbarW = taskbarRect.right - taskbarRect.left;
        taskbarH = taskbarRect.bottom - taskbarRect.top;

        // Taskbar on top or bottom
        if (taskbarW > taskbarH) {
            taskbarW = 0;

            // Taskbar on top
            if (taskbarRect.top == 0) {
                taskbarH = -taskbarH;
            }
        }
        // Taskbar on left or right
        else {
            taskbarH = 0;

            // Taskbar on left
            if (taskbarRect.left == 0) {
                taskbarW = -taskbarW;
            }
        }
    }

    uint32_t winX = (::GetSystemMetrics(SM_CXSCREEN) - WINDOW_WIDTH  - taskbarW) / 2;
    uint32_t winY = (::GetSystemMetrics(SM_CYSCREEN) - WINDOW_HEIGHT - taskbarH) / 2;

    // Create window
    mHwnd = CreateWindowEx(WS_EX_CLIENTEDGE, WINDOW_CLASS, WINDOW_NAME, WS_OVERLAPPEDWINDOW,
                               winX, winY, WINDOW_WIDTH, WINDOW_HEIGHT, 0, 0, hInstance, 0);

    if (!mHwnd) {
        throw std::runtime_error("Window creation failed.");
    }

    // Find executable directory
    CHAR exeFilePath[MAX_PATH] = { 0 };
    ::GetModuleFileNameA(NULL, &exeFilePath[0], MAX_PATH);
    std::string executablePath = std::filesystem::path(std::string(exeFilePath)).remove_filename().string();

    // Find resource directory
    std::string resourcePath = std::string(executablePath) + "Resources\\";
    if (::GetFileAttributesA(resourcePath.c_str()) == INVALID_FILE_ATTRIBUTES) {
        throw std::runtime_error("Unable to locate resource directory " + resourcePath);
    }

    mRenderer = std::make_shared<Renderer>(mHwnd, mHinstance);
}

App::~App()
{
    mRenderer.reset();

    ::DestroyWindow(mHwnd);
    ::UnregisterClass(WINDOW_CLASS, mHinstance);
}

void App::Run() const
{
    ::ShowWindow(mHwnd, SW_SHOW);

    MSG msg = {};
    while (msg.message != WM_QUIT) {
        if (::PeekMessage(&msg, nullptr, 0, 0, PM_REMOVE)) {
            ::TranslateMessage(&msg);
            ::DispatchMessage(&msg);
        }

        mRenderer->Update(WINDOW_WIDTH, WINDOW_HEIGHT);
        mRenderer->Render();
    }
}


void App::SetupConsole(std::wstring title)
{
    // Allocate console
    ::AllocConsole();
    ::AttachConsole(::GetCurrentProcessId());
    FILE* stream;
    std::ignore = freopen_s(&stream, "CONIN$", "r", stdin);
    std::ignore = freopen_s(&stream, "CONOUT$", "w+", stdout);
    std::ignore = freopen_s(&stream, "CONOUT$", "w+", stderr);

    // Enable flags so we can color the output
    HANDLE consoleHandle = GetStdHandle(STD_OUTPUT_HANDLE);
    DWORD dwMode = 0;
    ::GetConsoleMode(consoleHandle, &dwMode);
    dwMode |= ENABLE_VIRTUAL_TERMINAL_PROCESSING;
    ::SetConsoleMode(consoleHandle, dwMode);
    ::SetConsoleTitle(title.c_str());
}
