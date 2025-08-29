# Photron: Monte Carlo Photon Transport in Multi-Layered Tissues

Photron is a program that simulates Monte Carlo photon transport in multi-layered materials (such as skin tissue), and was initially developed as part of my experimentation project into subsurface scattering at Utrecht University.

The original [MCML](https://omlc.org/software/mc/mcml/index.html) program was created in 1994 by Lihong Wang and Steven Jacques while they worked together at the University of Texas M. D. Anderson Cancer Center, in the Laser Biology Research Laboratory. The program allows multiple planar layers of tissue, each with different optical properties (μa, μs, g, n) and thickness. The same authors released an improved and time-resolved [MCML 2.0](https://github.com/lhvwang/MCML) version in 1996.

The C++ reimplementation presented here is based on this 2.0 version. It simulates shooting an infinitely small photon on a multi-layered surface and renders the result in 3D using an OpenGL context.

## Build Requirements

* C++20
* Windows SDK 10.0
* [Visual Studio 2022](https://visualstudio.microsoft.com/vs/)

Open the provided Visual Studio solution to build the solution. Required libraries (except the Vulkan SDK) are included.

## Libraries

* [VulkanSDK](https://vulkan.lunarg.com) (1.4.304.1)
* [Dear ImGui](https://github.com/ocornut/imgui) (1.91.8)
* [JSON for Modern C++](https://github.com/nlohmann/json) (3.11.3)

## Configuration

## Quick Guide

## License

This source code is released under the [CC BY-SA 4.0](https://creativecommons.org/licenses/by-sa/4.0) license.

## References

* Wang, L. V.; Jacques, S. L.; Zheng, L.-Q.; "CONV — convolution for responses to a finite diameter photon beam incident on multi-layered tissues," Computer Methods and Programs in Biomedicine 54(3) 141-150 (1997)
* Wang, L. V.; Jacques, S. L.; Zheng, L. Q.; "MCML — Monte Carlo modeling of light transport in multilayered tissues," Computer Methods and Programs in Biomedicine 47(2) 131-146 (1995)
* Jacques, S. L.; Wang, L. V.; "Monte Carlo modeling of light transport in tissues," Optical Thermal Response of Laser Irradiated Tissue 73–100 (1995)
* Wang, L. V.; Jacques, S. L.; "Optimized radial and angular positions in Monte Carlo modeling," Medical Physics 21(7) 1081-1083 (1994)
