# Photron v1.0.0 - Monte Carlo Photon Transport Renderer
# Cross-platform Makefile with comprehensive build targets

# Project metadata
PROJECT = photron
VERSION = 1.0.0
SRCDIR = ../src
BINDIR = ../bin
$(BUILDDIR)/renderer/window.o: $(SRCDIR)/renderer/window.cpp $(HEADERS)
	$(ECHO) "Compiling $<..."
ifeq ($(CXX),cl)
	$(Q)$(CXX) $(CXXFLAGS) /c $(SRCDIR)/renderer/window.cpp /Fo$@
else
	$(Q)$(CXX) $(CXXFLAGS) $(INCLUDES) -c $(SRCDIR)/renderer/window.cpp -o $@
endif

# Copy configuration and shader files to bin directory
.PHONY: copy-configs output
CONFIGDIR = ../config
SHADERSDIR = ../shaders

# Build configuration
BUILD_TYPE ?= release
VERBOSE ?= 0

# Platform detection and cross-platform commands
ifeq ($(OS),Windows_NT)
	SHELL = cmd.exe
	EXEEXT = .exe
	MKDIR_P = if not exist "$(1)" mkdir "$(1)"
	RM = del /Q
	RMDIR = rmdir /S /Q
	DEVNULL = NUL
	PATHSEP = \\
	CP = copy
	XCOPY = xcopy /E /I /Y
	
	# For Windows, use g++ by default (user can override with make CXX=clang++)
	CXX ?= g++
	ifeq ($(CXX),cl)
		CXXFLAGS = /EHsc /std:c++20 /W3
		LIBS = opengl32.lib glu32.lib glew32.lib glfw3dll.lib gdi32.lib user32.lib
		COMPILER_FOUND = cl
		ifeq ($(BUILD_TYPE),debug)
			CXXFLAGS += /Od /Zi /DEBUG
		else
			CXXFLAGS += /O2 /DNDEBUG
		endif
	else
		# For g++, clang++, and other GCC-compatible compilers
		CXXFLAGS = -std=c++20 -Wall -Wextra -pedantic -Wconversion -Wsign-conversion
		LIBS = -lopengl32 -lglu32 -lglew32 -lglfw3dll -lgdi32 -luser32
		COMPILER_FOUND = $(CXX)
		ifeq ($(BUILD_TYPE),debug)
			CXXFLAGS += -g -O0 -DDEBUG
		else
			CXXFLAGS += -O2 -DNDEBUG
		endif
	endif
	
	# vcpkg paths for Windows
	VCPKG_DIR = ../vcpkg_installed/x64-windows
	INCLUDES = -I$(SRCDIR) -I$(SRCDIR)/structs -I$(SRCDIR)/renderer -I$(SRCDIR)/simulator -I$(SRCDIR)/utilities -I$(SRCDIR)/mcml -I$(VCPKG_DIR)/include
	LIBDIRS = -L$(VCPKG_DIR)/lib
else
	# Unix-like environment (Linux, macOS, WSL when not Windows_NT)
	EXEEXT =
	MKDIR_P = mkdir -p $(1)
	RM = rm -f
	RMDIR = rm -rf
	DEVNULL = /dev/null
	PATHSEP = /
	CP = cp
	XCOPY = cp -r
	# Unix-like systems - use c++ by default
	CXX ?= c++
	CXXFLAGS = -std=c++20 -Wall -Wextra -pedantic -Wconversion -Wsign-conversion
	LIBS = -lGL -lGLU -lGLEW -lglfw -lm
	COMPILER_FOUND = $(CXX)
	ifeq ($(BUILD_TYPE),debug)
		CXXFLAGS += -g -O0 -DDEBUG
	else
		CXXFLAGS += -O2 -DNDEBUG
	endif
	
	# Standard system paths for Unix-like systems
	INCLUDES = -I$(SRCDIR) -I$(SRCDIR)/structs -I$(SRCDIR)/renderer -I$(SRCDIR)/simulator -I$(SRCDIR)/utilities -I$(SRCDIR)/mcml
	LIBDIRS =
endif

# Verbose output control
ifeq ($(VERBOSE),1)
	Q =
	ECHO = @echo
else
	Q = @
	ECHO = @echo
endif

TARGET = $(BINDIR)/$(PROJECT)$(EXEEXT)
TEST_TARGET = $(BINDIR)/$(PROJECT)_test$(EXEEXT)

# Source files (equivalent to CMake GLOB_RECURSE, matching exact CMake configuration)
PHOTRON_SOURCES = $(SRCDIR)/main.cpp \
                  $(SRCDIR)/renderer/renderer.cpp \
                  $(SRCDIR)/renderer/window.cpp \
                  $(SRCDIR)/simulator/initializer.cpp \
                  $(SRCDIR)/simulator/simulator.cpp \
                  $(SRCDIR)/utilities/experimenter.cpp \
                  $(SRCDIR)/utilities/utilities.cpp

# MCML sources are no longer needed - only random.hpp is used via header include
MCML_SOURCES = 

# Combine all sources
ALL_SOURCES = $(PHOTRON_SOURCES) $(MCML_SOURCES)
HEADERS = $(wildcard $(SRCDIR)/*.hpp) $(wildcard $(SRCDIR)/*/*.hpp)

# Object files (put in build/output directory as requested)
OBJECTS = $(patsubst $(SRCDIR)/%.cpp, $(BUILDDIR)/%.o, $(ALL_SOURCES))

# Default target
.PHONY: all
all: info $(TARGET) copy-configs
	$(ECHO) "Build completed successfully with $(COMPILER_FOUND) compiler in $(BUILD_TYPE) mode"

# Information about build configuration
.PHONY: info
info:
	$(ECHO) "Photron $(VERSION) Build Information:"
	$(ECHO) "   Platform:     $(OS)"
	$(ECHO) "   Compiler:     $(COMPILER_FOUND)"
	$(ECHO) "   Build Type:   $(BUILD_TYPE)"
	$(ECHO) "   C++ Standard: C++20"
	$(ECHO) "   Source Dir:   $(SRCDIR)"
	$(ECHO) "   Binary Dir:   $(BINDIR)"
	$(ECHO) "   Build Dir:    $(BUILDDIR)"
	$(ECHO) "   Flags:        $(CXXFLAGS)"
	$(ECHO)

# Create necessary directories
$(BINDIR):
	$(Q)$(call MKDIR_P,$(BINDIR))

$(BUILDDIR):
	$(Q)$(call MKDIR_P,$(BUILDDIR))
	$(Q)$(call MKDIR_P,$(BUILDDIR)/renderer)
	$(Q)$(call MKDIR_P,$(BUILDDIR)/simulator)
	$(Q)$(call MKDIR_P,$(BUILDDIR)/utilities)

# Build Photron executable
.PHONY: photron
photron: $(TARGET) copy-configs

$(TARGET): $(BINDIR) $(BUILDDIR) $(OBJECTS)
	$(ECHO) "Linking $(PROJECT) v$(VERSION)..."
ifeq ($(CXX),cl)
	$(Q)$(CXX) $(CXXFLAGS) $(OBJECTS) /Fe$@ $(LIBS)
else
	$(Q)$(CXX) $(CXXFLAGS) $(LDFLAGS) -o $@ $(OBJECTS) $(LIBDIRS) $(LIBS)
endif

# Object file compilation rules (object files go in output directory)
$(BUILDDIR)/main.o: $(SRCDIR)/main.cpp $(HEADERS)
	$(ECHO) "Compiling $<..."
ifeq ($(CXX),cl)
	$(Q)$(CXX) $(CXXFLAGS) /c $(SRCDIR)/main.cpp /Fo$@
else
	$(Q)$(CXX) $(CXXFLAGS) $(INCLUDES) -c $(SRCDIR)/main.cpp -o $@
endif

$(BUILDDIR)/renderer/renderer.o: $(SRCDIR)/renderer/renderer.cpp $(HEADERS)
	$(ECHO) "Compiling $<..."
ifeq ($(CXX),cl)
	$(Q)$(CXX) $(CXXFLAGS) /c $(SRCDIR)/renderer/renderer.cpp /Fo$@
else
	$(Q)$(CXX) $(CXXFLAGS) $(INCLUDES) -c $(SRCDIR)/renderer/renderer.cpp -o $@
endif

$(BUILDDIR)/simulator/initializer.o: $(SRCDIR)/simulator/initializer.cpp $(HEADERS)
	$(ECHO) "Compiling $<..."
ifeq ($(CXX),cl)
	$(Q)$(CXX) $(CXXFLAGS) /c $(SRCDIR)/simulator/initializer.cpp /Fo$@
else
	$(Q)$(CXX) $(CXXFLAGS) $(INCLUDES) -c $(SRCDIR)/simulator/initializer.cpp -o $@
endif

$(BUILDDIR)/simulator/mcml_integration.o: $(SRCDIR)/simulator/mcml_integration.cpp $(HEADERS)
	$(ECHO) "Compiling $<..."
ifeq ($(CXX),cl)
	$(Q)$(CXX) $(CXXFLAGS) /c $(SRCDIR)/simulator/mcml_integration.cpp /Fo$@
else
	$(Q)$(CXX) $(CXXFLAGS) $(INCLUDES) -c $(SRCDIR)/simulator/mcml_integration.cpp -o $@
endif

$(BUILDDIR)/simulator/simulator.o: $(SRCDIR)/simulator/simulator.cpp $(HEADERS)
	$(ECHO) "Compiling $<..."
ifeq ($(CXX),cl)
	$(Q)$(CXX) $(CXXFLAGS) /c $(SRCDIR)/simulator/simulator.cpp /Fo$@
else
	$(Q)$(CXX) $(CXXFLAGS) $(INCLUDES) -c $(SRCDIR)/simulator/simulator.cpp -o $@
endif

$(BUILDDIR)/utilities/experimenter.o: $(SRCDIR)/utilities/experimenter.cpp $(HEADERS)
	$(ECHO) "Compiling $<..."
ifeq ($(CXX),cl)
	$(Q)$(CXX) $(CXXFLAGS) /c $(SRCDIR)/utilities/experimenter.cpp /Fo$@
else
	$(Q)$(CXX) $(CXXFLAGS) $(INCLUDES) -c $(SRCDIR)/utilities/experimenter.cpp -o $@
endif

$(BUILDDIR)/utilities/utilities.o: $(SRCDIR)/utilities/utilities.cpp $(HEADERS)
	$(ECHO) "Compiling $<..."
ifeq ($(CXX),cl)
	$(Q)$(CXX) $(CXXFLAGS) /c $(SRCDIR)/utilities/utilities.cpp /Fo$@
else
	$(Q)$(CXX) $(CXXFLAGS) $(INCLUDES) -c $(SRCDIR)/utilities/utilities.cpp -o $@
endif

# Copy configuration and shader files to bin directory
.PHONY: copy-configs
copy-configs: $(BINDIR)
	$(ECHO) "Copying configuration files to target directory..."
ifeq ($(OS),Windows_NT)
	$(Q)$(call MKDIR_P,$(BINDIR)$(PATHSEP)config)
	$(Q)$(XCOPY) "$(CONFIGDIR)$(PATHSEP)*" "$(BINDIR)$(PATHSEP)config$(PATHSEP)" >$(DEVNULL) 2>&1 || echo Config files not found
	$(ECHO) "Copying shader files to target directory..."
	$(Q)$(call MKDIR_P,$(BINDIR)$(PATHSEP)shaders)
	$(Q)$(XCOPY) "$(SHADERSDIR)$(PATHSEP)*" "$(BINDIR)$(PATHSEP)shaders$(PATHSEP)" >$(DEVNULL) 2>&1 || echo Shader files not found
else
	$(Q)$(call MKDIR_P,$(BINDIR)/config)
	$(Q)$(XCOPY) "$(CONFIGDIR)/"* "$(BINDIR)/config/" 2>$(DEVNULL) || echo "Config files not found"
	$(ECHO) "Copying shader files to target directory..."
	$(Q)$(call MKDIR_P,$(BINDIR)/shaders)
	$(Q)$(XCOPY) "$(SHADERSDIR)/"* "$(BINDIR)/shaders/" 2>$(DEVNULL) || echo "Shader files not found"
endif

# Testing targets
.PHONY: test
test: $(TARGET)
	$(ECHO) "Running Photron tests..."
ifeq ($(OS),Windows_NT)
	$(ECHO) "Running: $(PROJECT)$(EXEEXT) config$(PATHSEP)config1.in"
	$(Q)cd $(BINDIR) && $(PROJECT)$(EXEEXT) config$(PATHSEP)config1.in
else
	$(ECHO) "Running: $(PROJECT)$(EXEEXT) config/config1.in"
	$(Q)cd $(BINDIR) && ./$(PROJECT)$(EXEEXT) config/config1.in
endif

# Debug build
.PHONY: debug
debug:
	$(ECHO) "Building in debug mode..."
	$(Q)$(MAKE) BUILD_TYPE=debug all

# Release build (default)
.PHONY: release  
release:
	$(ECHO) "Building in release mode..."
	$(Q)$(MAKE) BUILD_TYPE=release all

# Static analysis with clang-tidy (if available)
.PHONY: analyze
analyze:
	$(ECHO) "Running static analysis..."
ifeq ($(OS),Windows_NT)
	$(Q)where clang-tidy >$(DEVNULL) 2>&1 && ( \
		clang-tidy $(ALL_SOURCES) -- $(CXXFLAGS) $(INCLUDES) \
	) || echo "clang-tidy not found - install it for static analysis"
else
	$(Q)command -v clang-tidy >$(DEVNULL) 2>&1 && ( \
		clang-tidy $(ALL_SOURCES) -- $(CXXFLAGS) $(INCLUDES) \
	) || echo "clang-tidy not found - install it for static analysis"
endif

# Format code with clang-format (if available)
.PHONY: format
format:
	$(ECHO) "Formatting source code..."
ifeq ($(OS),Windows_NT)
	$(Q)where clang-format >$(DEVNULL) 2>&1 && ( \
		clang-format -i $(ALL_SOURCES) $(HEADERS) \
	) || echo "clang-format not found - install it for code formatting"
else
	$(Q)command -v clang-format >$(DEVNULL) 2>&1 && ( \
		clang-format -i $(ALL_SOURCES) $(HEADERS) \
	) || echo "clang-format not found - install it for code formatting"
endif

# Install target
.PHONY: install
install: $(TARGET)
	$(ECHO) "Installing Photron v$(VERSION)..."
ifdef PREFIX
	$(call MKDIR_P,$(PREFIX)/bin)
ifeq ($(OS),Windows_NT)
	$(CP) "$(TARGET)" "$(PREFIX)$(PATHSEP)bin$(PATHSEP)"
else
	$(CP) $(TARGET) $(PREFIX)/bin/
endif
	$(ECHO) "Installed to $(PREFIX)/bin/"
else
	$(ECHO) "PREFIX not set. Use: make install PREFIX=/usr/local"
endif

# Clean targets  
.PHONY: clean
clean:
	$(ECHO) "Cleaning build artifacts..."
ifeq ($(OS),Windows_NT)
	$(Q)if exist "$(BUILDDIR)" $(RMDIR) "$(BUILDDIR)" >$(DEVNULL) 2>&1 || echo.
	$(Q)if exist "$(TARGET)" $(RM) "$(TARGET)" >$(DEVNULL) 2>&1
	$(Q)if exist "$(BINDIR)$(PATHSEP)config" $(RMDIR) "$(BINDIR)$(PATHSEP)config" >$(DEVNULL) 2>&1 || echo.
	$(Q)if exist "$(BINDIR)$(PATHSEP)shaders" $(RMDIR) "$(BINDIR)$(PATHSEP)shaders" >$(DEVNULL) 2>&1 || echo.
else
	$(Q)$(RMDIR) $(BUILDDIR) 2>$(DEVNULL) || true
	$(Q)$(RM) $(TARGET) 2>$(DEVNULL) || true
	$(Q)$(RMDIR) "$(BINDIR)/config" "$(BINDIR)/shaders" 2>$(DEVNULL) || true
endif

.PHONY: clean-all
clean-all: clean
	$(ECHO) "Cleaning all generated files..."

.PHONY: distclean
distclean: clean-all
	$(ECHO) "Distribution cleanup completed."

# Print discovered source files (equivalent to CMake debug messages)
.PHONY: sources
sources:
	$(ECHO) "Discovered Photron source files:"
ifeq ($(OS),Windows_NT)
	$(Q)for %%f in ($(PHOTRON_SOURCES)) do $(ECHO) "  %%f"
	$(ECHO) ""
	$(ECHO) "Discovered MCML source files:"
	$(Q)for %%f in ($(MCML_SOURCES)) do $(ECHO) "  %%f"
else
	$(Q)for src in $(PHOTRON_SOURCES); do echo "  $$src"; done
	$(ECHO) ""
	$(ECHO) "Discovered MCML source files:"
	$(Q)for src in $(MCML_SOURCES); do echo "  $$src"; done
endif

# Help target
.PHONY: help
help:
	$(ECHO) "Photron $(VERSION) Makefile Help"
	$(ECHO) ""
	$(ECHO) "Build Targets:"
	$(ECHO) "  all          - Build Photron executable (default)"
	$(ECHO) "  photron      - Build Photron executable" 
	$(ECHO) "  debug        - Build with debug flags"
	$(ECHO) "  release      - Build with optimization flags"
	$(ECHO) ""
	$(ECHO) "Testing Targets:"
	$(ECHO) "  test         - Run test with sample config"
	$(ECHO) ""
	$(ECHO) "Quality Targets:"
	$(ECHO) "  analyze      - Run static analysis (clang-tidy)"
	$(ECHO) "  format       - Format source code (clang-format)"
	$(ECHO) ""
	$(ECHO) "Maintenance Targets:"
	$(ECHO) "  clean        - Remove build artifacts"
	$(ECHO) "  clean-all    - Remove all generated files"
	$(ECHO) "  distclean    - Full cleanup for distribution"
	$(ECHO) "  install      - Install Photron (set PREFIX)"
	$(ECHO) ""
	$(ECHO) "Utility Targets:"
	$(ECHO) "  info         - Show build configuration"
	$(ECHO) "  sources      - Show discovered source files"
	$(ECHO) "  help         - Show this help message"
	$(ECHO) ""
	$(ECHO) "Variables:"
	$(ECHO) "  CXX          - C++ compiler (default: g++/c++)"
	$(ECHO) "  BUILD_TYPE   - debug or release (default: release)"
	$(ECHO) "  VERBOSE      - Show detailed build output (0 or 1)"  
	$(ECHO) "  PREFIX       - Installation prefix"

# Prevent deletion of intermediate files
.SECONDARY: $(OBJECTS)
