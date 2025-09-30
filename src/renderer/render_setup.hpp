/**
 * @file render_setup.hpp
 * @brief Centralized rendering pipeline setup and configuration
 *
 * Provides utilities for initializing and configuring the OpenGL rendering
 * pipeline, including shader programs, vertex arrays, and rendering state.
 */

#pragma once

#include <functional>
#include <iostream>
#include <string>
#include <vector>

#include <GL/glew.h>

#include "shader_utils.hpp"

namespace RenderSetup
{
/**
 * @struct ShaderConfig
 * @brief Configuration structure for basic shader setup
 *
 * Contains paths to shader files and error identification for the rendering pipeline.
 */
struct ShaderConfig
{
	std::string vertex_shader_path;   ///< Path to vertex shader source file
	std::string fragment_shader_path; ///< Path to fragment shader source file
	std::string error_name;           ///< Identifier for error messages and debugging
};

/**
 * @struct UniformConfig
 * @brief Configuration structure for uniform caching
 *
 * Maps uniform names to location storage pointers for efficient uniform access.
 */
struct UniformConfig
{
	std::string name;    ///< Name of the uniform variable in the shader
	GLint* location_ptr; ///< Pointer to store the cached uniform location
};

/**
 * @struct VertexAttributeConfig
 * @brief Configuration structure for vertex attribute setup
 *
 * Defines vertex attribute parameters for OpenGL vertex array configuration.
 */
struct VertexAttributeConfig
{
	GLuint location;      ///< Vertex attribute location (layout location in shader)
	GLint size;           ///< Number of components per vertex attribute (1-4)
	GLenum type;          ///< Data type of each component (GL_FLOAT, GL_INT, etc.)
	GLboolean normalized; ///< Whether fixed-point data should be normalized
	GLsizei stride;       ///< Byte offset between consecutive vertex attributes
	const void* pointer;  ///< Offset of first component in the vertex buffer
	GLuint divisor {0};   ///< Instance divisor: 0 = per vertex, 1+ = per instance
};

/**
 * @typedef GeometrySetupFunction
 * @brief Function type for base geometry setup in instanced rendering
 *
 * Callback function that configures the base geometry VBO for instanced rendering.
 *
 * @param vbo The vertex buffer object ID to configure with base geometry data
 * @return bool True if setup succeeded, false on error
 */
using GeometrySetupFunction = std::function<bool(GLuint vbo)>;

/**
 * @brief Template function for basic rendering setup (lines, points, triangles)
 *
 * Handles complete OpenGL rendering pipeline setup including shader loading,
 * program creation, uniform caching, VAO/VBO setup, and vertex attributes.
 *
 * @tparam VertexType The vertex data structure type for template specialization
 * @param shader_config Configuration containing shader paths and error name
 * @param uniforms Vector of uniform configurations for caching locations
 * @param vertex_attributes Vector of vertex attribute configurations
 * @param[out] shader_program Reference to store the created shader program ID
 * @param[out] vao Reference to store the created vertex array object ID
 * @param[out] vbo Reference to store the created vertex buffer object ID
 * @return bool True if setup succeeded, false on any error
 */
template<typename VertexType>
inline bool setup_basic_rendering(const ShaderConfig& shader_config,
								  const std::vector<UniformConfig>& uniforms,
								  const std::vector<VertexAttributeConfig>& vertex_attributes,
								  GLuint& shader_program,
								  GLuint& vao,
								  GLuint& vbo) {
	// Load shaders
	std::string vertex_source = ShaderUtils::load_shader_source(shader_config.vertex_shader_path);
	std::string fragment_source = ShaderUtils::load_shader_source(shader_config.fragment_shader_path);

	if (vertex_source.empty() || fragment_source.empty()) {
		std::cerr << "Failed to load " << shader_config.error_name << " shaders" << std::endl;
		return false;
	}

	// Create shader program
	GLuint vertex_shader = ShaderUtils::compile_shader(vertex_source, GL_VERTEX_SHADER);
	GLuint fragment_shader = ShaderUtils::compile_shader(fragment_source, GL_FRAGMENT_SHADER);

	if (vertex_shader == 0 || fragment_shader == 0) {
		if (vertex_shader) {
			glDeleteShader(vertex_shader);
		}
		if (fragment_shader) {
			glDeleteShader(fragment_shader);
		}
		return false;
	}

	shader_program = glCreateProgram();
	glAttachShader(shader_program, vertex_shader);
	glAttachShader(shader_program, fragment_shader);
	glLinkProgram(shader_program);

	GLint success;
	glGetProgramiv(shader_program, GL_LINK_STATUS, &success);
	if (!success) {
		char info_log[512];
		glGetProgramInfoLog(shader_program, 512, nullptr, info_log);
		std::cerr << "Program linking failed for " << shader_config.error_name << ": " << info_log << std::endl;
		glDeleteProgram(shader_program);
		shader_program = 0;
		glDeleteShader(vertex_shader);
		glDeleteShader(fragment_shader);
		return false;
	}

	glDeleteShader(vertex_shader);
	glDeleteShader(fragment_shader);

	// Cache uniform locations
	for (const auto& uniform : uniforms) {
		if (uniform.location_ptr) {
			*uniform.location_ptr = glGetUniformLocation(shader_program, uniform.name.c_str());
		}
	}

	// Create VAO and VBO
	glGenVertexArrays(1, &vao);
	glGenBuffers(1, &vbo);

	glBindVertexArray(vao);
	glBindBuffer(GL_ARRAY_BUFFER, vbo);

	// Setup vertex attributes
	for (const auto& attr : vertex_attributes) {
		glVertexAttribPointer(attr.location, attr.size, attr.type, attr.normalized, attr.stride, attr.pointer);
		glEnableVertexAttribArray(attr.location);
		if (attr.divisor > 0) {
			glVertexAttribDivisor(attr.location, attr.divisor);
		}
	}

	glBindVertexArray(0);
	return true;
}

/**
 * @brief Template function for instanced rendering setup (point instanced, line instanced)
 *
 * Handles complete instanced rendering pipeline setup including shader loading,
 * base geometry configuration, and per-instance attribute setup.
 *
 * @tparam InstanceType The instance data structure type for template specialization
 * @param shader_config Configuration containing shader paths and error name
 * @param uniforms Vector of uniform configurations for caching locations
 * @param geometry_setup Function callback to configure base geometry VBO
 * @param instance_attributes Vector of per-instance attribute configurations
 * @param[out] shader_program Reference to store the created shader program ID
 * @param[out] vao Reference to store the created vertex array object ID
 * @param[out] base_vbo Reference to store the base geometry VBO ID
 * @param[out] instance_vbo Reference to store the instance data VBO ID
 * @return bool True if setup succeeded, false on any error
 */
template<typename InstanceType>
inline bool setup_instanced_rendering(const ShaderConfig& shader_config,
									  const std::vector<UniformConfig>& uniforms,
									  const GeometrySetupFunction& geometry_setup,
									  const std::vector<VertexAttributeConfig>& instance_attributes,
									  GLuint& shader_program,
									  GLuint& vao,
									  GLuint& base_vbo,
									  GLuint& instance_vbo) {
	// Load shaders
	std::string vertex_source = ShaderUtils::load_shader_source(shader_config.vertex_shader_path);
	std::string fragment_source = ShaderUtils::load_shader_source(shader_config.fragment_shader_path);

	if (vertex_source.empty() || fragment_source.empty()) {
		std::cerr << "Failed to load " << shader_config.error_name << " shaders" << std::endl;
		return false;
	}

	// Create shader program
	GLuint vertex_shader = ShaderUtils::compile_shader(vertex_source, GL_VERTEX_SHADER);
	GLuint fragment_shader = ShaderUtils::compile_shader(fragment_source, GL_FRAGMENT_SHADER);

	if (vertex_shader == 0 || fragment_shader == 0) {
		if (vertex_shader) {
			glDeleteShader(vertex_shader);
		}
		if (fragment_shader) {
			glDeleteShader(fragment_shader);
		}
		return false;
	}

	shader_program = glCreateProgram();
	glAttachShader(shader_program, vertex_shader);
	glAttachShader(shader_program, fragment_shader);
	glLinkProgram(shader_program);

	GLint success;
	glGetProgramiv(shader_program, GL_LINK_STATUS, &success);
	if (!success) {
		char info_log[512];
		glGetProgramInfoLog(shader_program, 512, nullptr, info_log);
		std::cerr << "Program linking failed for " << shader_config.error_name << ": " << info_log << std::endl;
		glDeleteProgram(shader_program);
		shader_program = 0;
		glDeleteShader(vertex_shader);
		glDeleteShader(fragment_shader);
		return false;
	}

	glDeleteShader(vertex_shader);
	glDeleteShader(fragment_shader);

	// Cache uniform locations
	for (const auto& uniform : uniforms) {
		if (uniform.location_ptr) {
			*uniform.location_ptr = glGetUniformLocation(shader_program, uniform.name.c_str());
		}
	}

	// Create VAO and VBOs
	glGenVertexArrays(1, &vao);
	glGenBuffers(1, &base_vbo);
	glGenBuffers(1, &instance_vbo);

	glBindVertexArray(vao);

	// Setup base geometry using the provided function
	if (!geometry_setup(base_vbo)) {
		std::cerr << "Failed to setup base geometry for " << shader_config.error_name << std::endl;
		glBindVertexArray(0);
		return false;
	}

	// Setup instance data buffer
	glBindBuffer(GL_ARRAY_BUFFER, instance_vbo);

	// Setup instance attributes
	for (const auto& attr : instance_attributes) {
		glVertexAttribPointer(attr.location, attr.size, attr.type, attr.normalized, attr.stride, attr.pointer);
		glEnableVertexAttribArray(attr.location);
		if (attr.divisor > 0) {
			glVertexAttribDivisor(attr.location, attr.divisor);
		}
	}

	glBindVertexArray(0);
	return true;
}

} // namespace RenderSetup
