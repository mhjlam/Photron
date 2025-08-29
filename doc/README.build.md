# Photron Build System

This directory contains both **Makefile** and **CMake** build systems for the Photron Monte Carlo Photon Transport Renderer.

## Project Structure

The project is organized as follows:
```
Photron/
├── build/           # Build files (you are here)
│   ├── Makefile     # Unified cross-platform Makefile
│   └── CMakeLists.txt # Modern CMake build configuration
├── src/             # Source code
│   ├── main.cpp     # Main application entry point
│   ├── structs/     # Data structures and types
│   ├── renderer/    # OpenGL rendering code
│   ├── simulator/   # Monte Carlo simulation engine
│   └── utilities/   # Utility functions and helpers
├── config/          # Configuration files (*.in)
├── bin/             # Compiled executables (generated)
└── doc/             # Documentation
```

## Dependencies

This project requires:
- **OpenGL**: Graphics rendering  
- **GLFW**: Modern windowing and input handling (replaces GLUT)
- **GLEW**: OpenGL extension loading
- **C++17 compiler**: GCC, Clang, or MSVC

### Installing Dependencies

**Recommended: Using vcpkg (for CMake)**

The project includes a `vcpkg.json` file that automatically manages dependencies when using CMake with vcpkg:

1. **Install vcpkg** (if not already installed):
   ```bash
   git clone https://github.com/Microsoft/vcpkg.git
   cd vcpkg
   ./bootstrap-vcpkg.sh  # On Linux/macOS
   ./bootstrap-vcpkg.bat # On Windows
   ```

2. **Set environment variable**:
   ```bash
   export VCPKG_ROOT=/path/to/vcpkg  # Linux/macOS
   set VCPKG_ROOT=C:\path\to\vcpkg   # Windows
   ```

3. **Build with CMake** (dependencies will be installed automatically):
   ```bash
   cmake . -B output
   cmake --build output
   ```

**Alternative: Manual Installation**

**Windows:**
- Install MinGW-w64 or Visual Studio
- Download and install GLFW and GLEW libraries
- Update library paths in Makefile if needed

**Ubuntu/Debian:**
```bash
sudo apt install build-essential libgl1-mesa-dev libglu1-mesa-dev libglfw3-dev libglew-dev
```

**macOS:**
```bash
brew install glfw glew
```

## Build Options

### Option 1: Using Make (Recommended)

The Makefile provides a simple, unified build system that works across platforms.

**Basic commands:**
```bash
# Build release version (default)
make

# Build debug version
make debug

# Build with verbose output
make VERBOSE=1

# Clean build artifacts
make clean

# Deep clean (remove all generated files)
make distclean

# Show help
make help

# Build and run
make run
```

**Configuration options:**
- `BUILD_TYPE=debug|release` - Build configuration (default: release)
- `CXX=g++|clang++|cl` - Compiler selection (default: g++)
- `VERBOSE=0|1` - Verbose compilation output (default: 0)

**Examples:**
```bash
# Debug build with Clang
make debug CXX=clang++

# Release build with verbose output
make BUILD_TYPE=release VERBOSE=1

# Using MSVC on Windows
make CXX=cl
```

### Option 2: Using CMake

CMake provides more advanced features and better IDE integration. With vcpkg integration, dependencies are automatically managed.

**Basic commands (with vcpkg):**
```bash
# Make sure VCPKG_ROOT is set in your environment
# Configure and build (from build directory)
cmake . -B output
cmake --build output

# Or specify build type
cmake -DCMAKE_BUILD_TYPE=Release . -B output
cmake --build output

# For debug builds
cmake -DCMAKE_BUILD_TYPE=Debug . -B output
cmake --build output
```

**Advanced CMake usage:**
```bash
# Generate Visual Studio project files (Windows)
cmake -G "Visual Studio 17 2022" . -B output

# Generate Xcode project files (macOS)  
cmake -G Xcode . -B output

# Build with specific compiler
cmake -DCMAKE_CXX_COMPILER=clang++ . -B output

# Install to system directories
cmake --build output --target install
```

**Without vcpkg (manual dependencies):**
If you prefer to manage dependencies manually, remove or ignore the `vcpkg.json` file and install libraries using your system package manager.

## Platform-Specific Notes

### Windows
- Both MinGW/GCC and MSVC are supported
- Library paths may need adjustment in Makefile
- CMake will try to auto-detect OpenGL libraries

### Linux/Unix
- Standard package manager installations should work
- Use `sudo make install` to install system-wide

### macOS
- OpenGL frameworks are used instead of libraries
- Install dependencies via Homebrew

## Troubleshooting

**Library not found errors:**
- Update library paths in the Makefile `LIBDIRS` and `CXXFLAGS` sections
- For CMake, set environment variables like `GLUT_ROOT_PATH`

**Compilation errors:**
- Ensure you have C++17 compatible compiler
- Check that all dependencies are installed
- Try with verbose output: `make VERBOSE=1`

**Runtime errors:**
- Make sure config files are in `bin/config/` directory
- Check that OpenGL drivers are properly installed

## Output

Both build systems will:
1. Create the `../bin/` directory
2. Compile the `photron` executable (or `photron.exe` on Windows)
3. Copy configuration files to `../bin/config/`

## Development

For development, the Makefile includes additional targets:
- `make deps` - Generate dependency files
- `make clean` - Remove build artifacts only
- `make distclean` - Remove all generated files

The CMake system automatically handles dependencies and provides better IDE integration for development.

## License

See the main project documentation for license information.
