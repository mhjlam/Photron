#include <iostream>
#include <stdexcept>

#include <windows.h>
#include <vulkan/vulkan.h>

#include "App.hpp"


int WINAPI WinMain(_In_ HINSTANCE hInstance, _In_opt_ HINSTANCE hPrevInstance, _In_ LPSTR lpCmdLine, _In_ int nShowCmd)
{
    try {
        App app{hInstance};
        app.Run();
    }
    catch (const std::exception& e) {
        std::cerr << "Error: " << e.what() << std::endl;
        ::MessageBoxA(nullptr, e.what(), "ERROR", MB_ICONEXCLAMATION | MB_OK);
        return EXIT_FAILURE;
    }
    return EXIT_SUCCESS;
}
