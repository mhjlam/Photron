/**
 * @file vertex_array.hpp
 * @brief OpenGL Vertex Array Object (VAO) management and configuration
 *
 * Provides RAII-based VAO management for organizing vertex attribute data
 * and streamlining OpenGL rendering state management.
 */

#pragma once

#include <GL/glew.h>

/**
 * @class VertexArray
 * @brief RAII-based OpenGL Vertex Array Object (VAO) wrapper
 *
 * Manages OpenGL vertex array objects with automatic resource management
 * and provides convenient methods for vertex attribute configuration.
 * Uses lazy initialization to create the VAO only when needed.
 */
class VertexArray
{
public:
	/**
	 * @brief Default constructor (lazy initialization)
	 */
	VertexArray() = default;

	/**
	 * @brief Destructor - automatically cleans up OpenGL VAO
	 */
	~VertexArray();

	/**
	 * @brief Bind this VAO for subsequent OpenGL operations
	 */
	void bind() const;

	/**
	 * @brief Unbind the current VAO (bind VAO 0)
	 */
	void unbind() const;

	/**
	 * @brief Enable a vertex attribute array
	 * @param index Vertex attribute index to enable
	 */
	void enable_attribute(GLuint index) const;

	/**
	 * @brief Disable a vertex attribute array
	 * @param index Vertex attribute index to disable
	 */
	void disable_attribute(GLuint index) const;

	/**
	 * @brief Configure a vertex attribute pointer
	 * @param index Vertex attribute index
	 * @param size Number of components per vertex (1-4)
	 * @param type Data type of each component
	 * @param normalized Whether fixed-point data should be normalized
	 * @param stride Byte offset between consecutive vertex attributes
	 * @param offset Offset of first component in the buffer
	 */
	void set_attribute_pointer(
		GLuint index, GLint size, GLenum type, GLboolean normalized, GLsizei stride, const void* offset) const;

	/**
	 * @brief Configure a floating-point vertex attribute (convenience method)
	 * @param index Vertex attribute index
	 * @param size Number of float components per vertex (1-4)
	 * @param stride Byte offset between consecutive vertex attributes
	 * @param offset Offset of first component in the buffer
	 * @param normalized Whether to normalize the data (default: false)
	 */
	void set_float_attribute(
		GLuint index, GLint size, GLsizei stride, const void* offset, GLboolean normalized = GL_FALSE) const;

	/**
	 * @brief Get the OpenGL VAO ID
	 * @return GLuint VAO identifier (0 if not yet created)
	 */
	GLuint id() const { return vao_id_; }

	/**
	 * @brief Check if the VAO has been created and is valid
	 * @return bool True if VAO is valid, false otherwise
	 */
	bool is_valid() const { return vao_id_ != 0; }

private:
	/**
	 * @brief Lazy initialization - create VAO if not yet created
	 */
	void ensure_created() const;

	mutable GLuint vao_id_ {0}; ///< OpenGL VAO identifier (0 = not created)
};
