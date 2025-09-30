/**
 * @file shader_utils.cpp
 * @brief OpenGL shader compilation and loading utilities
 *
 * Provides essential utilities for OpenGL shader management:
 * - GLSL source code compilation with error reporting
 * - Shader file loading from filesystem
 * - Compilation error handling and diagnostics
 * - Cross-platform shader resource management
 *
 * These utilities support the rendering pipeline for Monte Carlo photon
 * visualization with robust error handling for development workflows.
 */

#include "shader_utils.hpp"

#include <array>
#include <fstream>
#include <iostream>
#include <sstream>

namespace ShaderUtils
{

/**
 * @brief Compile GLSL shader from source code
 * @param source GLSL source code as string
 * @param shader_type OpenGL shader type (vertex, fragment, etc.)
 * @return OpenGL shader object ID, or 0 on failure
 *
 * Compiles GLSL source code into OpenGL shader object with comprehensive
 * error reporting. Essential for dynamic shader loading and development.
 */
GLuint compile_shader(const std::string& source, GLenum shader_type) {
	// Create shader object for specified type
	GLuint shader = glCreateShader(shader_type);
	const char* source_cstr = source.c_str();

	// Load source code and compile shader
	glShaderSource(shader, 1, &source_cstr, nullptr);
	glCompileShader(shader);

	// Check compilation status and handle errors
	GLint success;
	glGetShaderiv(shader, GL_COMPILE_STATUS, &success);
	if (!success) {
		// Extract and report compilation error details
		std::array<char, 512> info_log {};
		glGetShaderInfoLog(shader, static_cast<GLsizei>(info_log.size()), nullptr, info_log.data());
		std::cerr << "Shader compilation failed: " << info_log.data() << std::endl;

		// Clean up failed shader object
		glDeleteShader(shader);
		return 0;
	}

	return shader;
}

/**
 * @brief Load shader source code from file
 * @param file_path Path to GLSL shader file
 * @return Shader source code as string, empty on failure
 *
 * Loads GLSL source code from filesystem with error handling.
 * Supports the development workflow for shader iteration and testing.
 */
std::string load_shader_source(const std::string& file_path) {
	// Open shader file for reading
	std::ifstream file(file_path);
	if (!file.is_open()) {
		std::cerr << "Failed to open shader file: " << file_path << std::endl;
		return "";
	}

	// Read entire file content into string
	std::stringstream buffer;
	buffer << file.rdbuf();
	return buffer.str();
}

} // namespace ShaderUtils
