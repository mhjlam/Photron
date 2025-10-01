<!-- markdownlint-disable MD033 -->

# Photron

<p align="center">Monte Carlo Photon Transport Simulation and Real-Time 3D Visualization</p>

<p align="center">
   <img src="media/photron-demo.gif" alt="Photron Demo" width="600" style="border-radius: 8px;"/>
   <br>
   <em>Real-time Monte Carlo photon transport simulation with interactive 3D visualization.</em>
</p>

## Overview

**Photron** is a high-performance Monte Carlo photon transport simulator with real-time 3D visualization, designed for subsurface light transport analysis in multi-layered biological tissues and materials. Built on modern C++20 and OpenGL 4.5, it combines the proven MCML (Monte Carlo Multi-Layered) algorithm with advanced rendering techniques to provide both accurate simulation and intuitive visualization of photon behavior in complex geometries.

The application features a dual-architecture design separating the simulation engine from the rendering system, enabling both interactive GUI mode for research and visualization, and headless batch processing for high-throughput computational workflows.

## Features

<table align="center">
   <tr>
      <td align="center"><img src="media/voxel-rendering.png" alt="Voxel Rendering" width="300"/><br><strong>Energy Deposition Visualization</strong></td>
      <td align="center"><img src="media/photon-paths.png" alt="Photon Paths" width="300"/><br><strong>Interactive Path Tracing</strong></td>
      <td align="center"><img src="media/material-layers.png" alt="Material Layers" width="300"/><br><strong>Multi-Layer Geometry</strong></td>
   </tr>
</table>

### Simulation

- **MCML 3.0 Algorithm**: Advanced Monte Carlo photon transport based on proven MCML methodology with modern optimizations.
- **Multi-Layered Materials**: Support for complex geometries with varying optical properties (absorption, scattering, anisotropy, refractive index).
- **Voxelized Discretization**: High-resolution 3D voxel grids for accurate energy deposition tracking and spatial analysis.
- **Interactive Photon Launching**: Launching of additional photons after initial simulation and monitoring energy conservation and performance metrics.

### Visualization

- **Real-Time Visualization**: Modern OpenGL 4.5 pipeline with shader-based rendering and interactive camera.
- **Energy Mapping**: Color-coded visualization of absorbed energy, photon density, and scattering events.
- **Photon Paths**: Interactive display of photon trajectories, scattering points, and material interfaces.
- **Render Modes**: Switchable rendering modes including absorption, emittance, volumes, and photon paths.
- **Performance Monitoring**: Real-time performance metrics and simulation progress tracking.

### Algorithms

- **3D DDA Traversal**: Digital Differential Analyzer for robust, precision-safe voxel traversal without floating-point errors.
- **BVH Acceleration**: Bounding Volume Hierarchy with Surface Area Heuristic (SAH) for logarithmic-complexity ray-triangle intersection.
- **Instanced Rendering**: GPU-optimized instanced rendering for thousands of voxels and even more photon path segments.
- **Spatial Optimization**: Octree-like spatial data structures and caching systems for performance-critical operations.

## Scene Configuration

Photron uses **TOML** configuration files to define simulation parameters like geometry, layers, and material properties. These config files support a range of scenarios from simple educational demonstrations to more complex research simulations.

### Configuration Files

A number of configuration files are already defined. These showcase several of the capabilities and scenarios that Photron can be used for.

#### Basic Geometry

- **`default.toml`** - Simple rectangular box for basic use-case testing.
- **`pyramid.toml`** - Inverted pyramidal structure for directional scattering studies.
- **`sphere.toml`** - Spherical geometry using a crude icosahedral approximation.

#### Material Properties

- **`absorption.toml`** - High-absorption materials for energy deposition analysis.
- **`scatter-*.toml`** - Various scattering configurations (forward, backward, high-scatter).
- **`anisotropy.toml`** - Anisotropic scattering parameter studies.
- **`interface.toml`** - Refractive index interface investigations.

#### Multiple Layers

- **`nested.toml`** - Nested geometries with multiple material boundaries.
- **`mixed.toml`** - Heterogeneous materials with varying properties.
- **`stack-*.toml`** - Multi-layer mediums, where the layers are stacked in different directions.

### Configuration Structure

The Photron TOML configuration file is kept intentionally simple and has a general and (light) source section, and the material and geometry definitions of as many layers as are desired.

#### Optical Properties

- **Refractive Index (eta)**: `1.0-3.0`, defines material interface behavior.
- **Absorption Coefficient (mua)**: `0.0+`, controls energy absorption rate.
- **Scattering Coefficient (mus)**: `0.0+`, determines scattering frequency.
- **Anisotropy Factor (ani)**: `-1.0 to +1.0`, controls scattering directionality.

#### Example file

```toml
[general]
photons = 100                # Number of Monte Carlo photons
voxel_size = 0.01            # Voxel edge length in cm
log = false                  # Enable/disable detailed logging
deterministic = false        # Random or fixed seed

[source]
position = [0.0, 0.2, 0.0]   # Light source position (cm)
direction = [0, -1, -0.5]    # Emission direction (normalized)

[[layer]]
eta = 1.37                   # Refractive index
mua = 1.0                    # Absorption coefficient (1/cm)
mus = 10.0                   # Scattering coefficient (1/cm) 
ani = 0.1                    # Anisotropy factor (-1 to +1)

# 3D mesh geometry definition
vertices = [
    [-1.0,  0.0,  1.0],      # Vertex coordinates in cm
    [ 1.0,  0.0,  1.0],
    # ... additional vertices
]

faces = [
    [0, 1, 2],               # Triangle face indices
    [0, 2, 3],
    # ... additional faces
]
```

## User Guide

Photron has both an interactive and headless mode. Starting the program without any arguments will show an empty screen with a config loader/selector in the middle of the window:

```bash
./Photron.exe
```

For a short explanation of the controls of the GUI, hover over the (?) tooltip on the top-right of the screen.

Supplying an argument to a valid config .toml file will run its simulation and then show the graphical user interface window with the results:

```bash
./Photron.exe config/default.toml
```

Another way to run Photron is in headless mode, which will only show output in the console and won't display a GUI. Detailed results of the simulation are automatically written to output files:

```bash
./Photron.exe config/test.toml --headless
```

Run `./Photron.exe --help` for more information.

## Technical Implementation

### Algorithms

#### Monte Carlo Photon Transport (MCML)

- **Photon Launching**: Configurable source geometries and emission patterns.
- **Step Size Calculation**: Statistically accurate transport distance sampling.
- **Scattering Events**: Henyey-Greenberg phase function with anisotropy support.
- **Interface Handling**: Fresnel reflection/transmission with angle-dependent probabilities.
- **Energy Deposition**: Accurate tracking in discrete voxel grid with energy conservation validation.

#### 3D Digital Differential Analyzer (DDA)

- **Precision-Safe Traversal**: Integer-based algorithm avoiding floating-point precision issues.
- **Voxel-Perfect Sampling**: Guarantees no missed voxels along photon paths.
- **Performance Optimized**: `O(n)` complexity where `n` is traversed voxel count.
- **Multi-Medium Support**: Handles transitions between materials with different properties.

#### Bounding Volume Hierarchy (BVH)

- **SAH-Optimized Construction**: Surface Area Heuristic for optimal tree partitioning.
- **Logarithmic Ray Intersection**: `O(log n)` complexity for ray-triangle intersection.
- **Memory Efficient**: Compact node representation with triangle index lists.
- **Dynamic Rebuilding**: Automatic reconstruction when geometry changes.

#### GPU Acceleration

- **Shader Pipeline**: Modern GLSL shaders with programmable vertex/fragment processing.
- **Instanced Rendering**: Efficient batch rendering of millions of voxels and path segments.
- **Buffer Optimization**: Vertex buffer object (VBO) caching and reuse strategies.

### Applications

- **Material Science**: Subsurface scattering analysis in translucent materials and composites.
- **Computer Graphics**: Accurate subsurface scattering for realistic material rendering.
- **Optical Engineering**: Light guide design and optical component optimization.
- **Physics Simulation**: Interactive demonstration of photon transport phenomena.
- **Parameter Studies**: Interactive exploration of material optical properties and their effects.

## Building and Dependencies

### Requirements

- **C++ Compiler**: C++20 support (MSVC 2022, GCC 11+, Clang 13+).
- **CMake**: 3.20 or higher for modern CMake features.
- **vcpkg**: Package manager for dependency resolution.
- **OpenGL**: 4.5 or later and a supported GPU.

### Libraries

- [GLFW 3.4](https://www.glfw.org/)
- [GLEW 2.2](http://glew.sourceforge.net/)
- [GLM 1.0](https://github.com/g-truc/glm)
- [ImGui 1.91](https://github.com/ocornut/imgui)
- [cxxopts 3.3](https://github.com/jarro2783/cxxopts)
- [toml++ 3.4](https://github.com/marzer/tomlplusplus)

### Build Instructions

#### Windows (Visual Studio)

```powershell
# Clone repository
git clone https://github.com/mhjlam/photron.git
cd photron

# Configure with vcpkg
cmake -S . -B build -DCMAKE_TOOLCHAIN_FILE="$ENV:VCPKG_ROOT/scripts/buildsystems/vcpkg.cmake"

# Build debug version  
cmake --build build --config Debug

# Build release version
cmake --build build --config Release
```

#### Linux (vcpkg)

```bash
# Install dependencies via vcpkg
vcpkg install glfw3 glew glm imgui cxxopts tomlplusplus

# Configure and build
cmake -S . -B build -DCMAKE_TOOLCHAIN_FILE=$VCPKG_ROOT/scripts/buildsystems/vcpkg.cmake
cmake --build build --config Release

# Run from bin directory (required for shader loading)
cd bin && ./Photron config/default.toml
```

## File Structure

```txt
Photron/
├── bin/                # Executable output and runtime files
│   ├── config/         # Configuration files copied for runtime
│   ├── out/            # Output directory for simulation results
│   └── shaders/        # GLSL shader files copied for runtime
├── build/              # CMake build outputs and intermediate files
├── config/             # Source configuration files (TOML format)
├── shaders/            # Source GLSL shader programs
├── src/                # C++ source code (coordination layer)
│   ├── common/         # Shared utilities (error handling, file I/O)
│   ├── math/           # Mathematical algorithms (DDA, BVH, geometry)
│   ├── renderer/       # OpenGL rendering and photon visualizer
│   └── simulator/      # Monte Carlo simulation engine
├── CMakeLists.txt      # CMake build configuration
└── vcpkg.json          # Dependency manifest for vcpkg
```

## Acknowledgments

- The MCML algorithm is based on the Monte Carlo Multi-Layered algorithm by Wang & Jacques (1996).
- DDA Implementation inspired by "A Fast Voxel Traversal Algorithm" by Amanatides & Woo (1987).
- Bounding Volume Hierarchy and Surface Area Heuristic optimizations from modern ray tracing literature.

## License

This software is licensed under the [CC BY-NC-SA](https://creativecommons.org/licenses/by-nc-sa/4.0/) license.
