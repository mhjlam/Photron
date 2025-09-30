/**
 * @file shader.hpp
 * @brief OpenGL shader program wrapper with automatic resource management
 *
 * Provides a high-level interface for GLSL shader compilation, linking,
 * and uniform variable management. Features performance optimizations
 * including uniform location caching and automatic error handling.
 */

#pragma once

#include <string>
#include <unordered_map>

#include <GL/glew.h>
#include <glm/glm.hpp>

/**
 * @class Shader
 * @brief RAII wrapper for OpenGL shader programs with uniform management
 *
 * The Shader class provides a complete abstraction over OpenGL shader programs,
 * handling compilation, linking, and uniform variable management. Key features:
 *
 * **Shader Loading:**
 * - Load from files with automatic source code reading
 * - Load from source strings for runtime generation
 * - Automatic error handling with detailed compilation messages
 *
 * **Performance Optimizations:**
 * - Uniform location caching to avoid repeated glGetUniformLocation calls
 * - Efficient uniform setter methods for all GLM types
 * - Minimal OpenGL state changes with use/unuse pattern
 *
 * **Resource Management:**
 * - RAII-based automatic shader program cleanup
 * - Exception-safe resource allocation
 * - Proper OpenGL context management
 *
 * **Usage Pattern:**
 * ```cpp
 * Shader basic_shader;
 * if (basic_shader.load_from_file("shaders/basic.vert", "shaders/basic.frag")) {
 *     basic_shader.use();
 *     basic_shader.set_uniform("mvp_matrix", mvp);
 *     basic_shader.set_uniform("color", glm::vec4(1.0f, 0.0f, 0.0f, 1.0f));
 *     // ... render geometry
 *     basic_shader.unuse();
 * }
 * ```
 */
class Shader
{
public:
	/**
	 * @brief Default constructor creates empty shader program
	 */
	Shader() = default;

	/**
	 * @brief Destructor automatically cleans up OpenGL shader program
	 */
	~Shader();

	// Shader loading methods

	/**
	 * @brief Load and compile shader program from vertex and fragment shader files
	 *
	 * Reads shader source code from files, compiles both shaders, and links
	 * them into a complete program. Handles all compilation and linking errors.
	 *
	 * @param vertex_path Path to vertex shader file (.vert)
	 * @param fragment_path Path to fragment shader file (.frag)
	 * @return true if compilation and linking succeeded, false otherwise
	 *
	 * @note Paths are relative to the executable's working directory.
	 *       Error messages are written to stderr on failure.
	 */
	bool load_from_file(const std::string& vertex_path, const std::string& fragment_path);

	/**
	 * @brief Load and compile shader program from source code strings
	 *
	 * Compiles shader program directly from provided GLSL source code.
	 * Useful for runtime-generated shaders or embedded shader code.
	 *
	 * @param vertex_source Complete GLSL vertex shader source code
	 * @param fragment_source Complete GLSL fragment shader source code
	 * @return true if compilation and linking succeeded, false otherwise
	 */
	bool load_from_source(const std::string& vertex_source, const std::string& fragment_source);

	// OpenGL state management

	/**
	 * @brief Activate shader program for rendering
	 *
	 * Makes this shader program active for subsequent draw calls.
	 * Must be called before setting uniforms or issuing draw commands.
	 */
	void use() const;

	/**
	 * @brief Deactivate shader program (restore default)
	 *
	 * Unbinds shader program, restoring OpenGL default state.
	 * Good practice for clean state management.
	 */
	void unuse() const;

	// Uniform variable setters

	/**
	 * @brief Set boolean uniform variable
	 * @param name Uniform variable name in shader
	 * @param value Boolean value to set
	 */
	void set_uniform(const std::string& name, bool value);

	/**
	 * @brief Set integer uniform variable
	 * @param name Uniform variable name in shader
	 * @param value Integer value to set
	 */
	void set_uniform(const std::string& name, int value);

	/**
	 * @brief Set floating-point uniform variable
	 * @param name Uniform variable name in shader
	 * @param value Float value to set
	 */
	void set_uniform(const std::string& name, float value);

	/**
	 * @brief Set 2D vector uniform variable
	 * @param name Uniform variable name in shader
	 * @param value 2D vector value to set
	 */
	void set_uniform(const std::string& name, const glm::vec2& value);

	/**
	 * @brief Set 3D vector uniform variable
	 * @param name Uniform variable name in shader
	 * @param value 3D vector value to set
	 */
	void set_uniform(const std::string& name, const glm::vec3& value);

	/**
	 * @brief Set 4D vector uniform variable
	 * @param name Uniform variable name in shader
	 * @param value 4D vector (RGBA color or homogeneous coordinate)
	 */
	void set_uniform(const std::string& name, const glm::vec4& value);

	/**
	 * @brief Set 3x3 matrix uniform variable
	 * @param name Uniform variable name in shader
	 * @param value 3x3 transformation matrix
	 */
	void set_uniform(const std::string& name, const glm::mat3& value);

	/**
	 * @brief Set 4x4 matrix uniform variable
	 * @param name Uniform variable name in shader
	 * @param value 4x4 transformation matrix (MVP, model, view, projection)
	 */
	void set_uniform(const std::string& name, const glm::mat4& value);

	// Query methods

	/**
	 * @brief Get OpenGL program ID for advanced operations
	 * @return GLuint OpenGL program object ID (0 if invalid)
	 */
	GLuint id() const { return program_id_; }

	/**
	 * @brief Check if shader program is valid and ready for use
	 * @return true if program compiled and linked successfully
	 */
	bool is_valid() const { return program_id_ != 0; }

private:
	/// OpenGL shader program object ID (0 if not created)
	GLuint program_id_ {0};

	/// Cache uniform locations to avoid repeated glGetUniformLocation calls
	mutable std::unordered_map<std::string, GLint> uniform_location_cache_;

	// Helper methods for shader compilation

	/**
	 * @brief Read complete file contents into string
	 * @param path File path to read
	 * @return std::string Complete file contents
	 */
	std::string read_file(const std::string& path);

	/**
	 * @brief Compile individual shader from source code
	 * @param source GLSL source code
	 * @param type OpenGL shader type (GL_VERTEX_SHADER, GL_FRAGMENT_SHADER)
	 * @return GLuint Compiled shader object ID (0 on failure)
	 */
	GLuint compile_shader(const std::string& source, GLenum type);

	/**
	 * @brief Link vertex and fragment shaders into complete program
	 * @param vertex_shader Compiled vertex shader object
	 * @param fragment_shader Compiled fragment shader object
	 * @return true if linking succeeded, false otherwise
	 */
	bool link_program(GLuint vertex_shader, GLuint fragment_shader);

	/**
	 * @brief Get cached uniform location (with performance optimization)
	 * @param name Uniform variable name
	 * @return GLint Uniform location (-1 if not found)
	 */
	GLint get_uniform_location(const std::string& name);
};
