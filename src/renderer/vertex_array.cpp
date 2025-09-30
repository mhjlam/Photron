/**
 * @file vertex_array.cpp
 * @brief OpenGL Vertex Array Object (VAO) management and attribute configuration
 *
 * Implements RAII-based VAO management for OpenGL vertex data organization:
 * - Automatic VAO creation and cleanup
 * - Vertex attribute pointer configuration
 * - Attribute enable/disable state management
 * - Type-safe vertex layout specification
 *
 * VAOs efficiently encapsulate vertex attribute state for rendering
 * Monte Carlo photon paths, voxel geometry, and UI elements.
 */

#include "vertex_array.hpp"

/**
 * @brief Destructor ensuring proper OpenGL VAO cleanup
 *
 * Automatically deletes OpenGL vertex array object when VertexArray
 * instance goes out of scope, preventing resource leaks.
 */
VertexArray::~VertexArray() {
	// Clean up OpenGL vertex array if it was created
	if (vao_id_ != 0) {
		glDeleteVertexArrays(1, &vao_id_);
	}
}

/**
 * @brief Bind vertex array to OpenGL context
 *
 * Makes this VAO the active vertex array, enabling all configured
 * vertex attribute pointers and states for rendering operations.
 */
void VertexArray::bind() const {
	// Ensure VAO exists before binding
	ensure_created();
	if (vao_id_ != 0) {
		glBindVertexArray(vao_id_);
	}
}

/**
 * @brief Unbind vertex array from OpenGL context
 *
 * Clears active VAO binding to prevent accidental state modifications
 * and ensure clean OpenGL state for subsequent operations.
 */
void VertexArray::unbind() const {
	glBindVertexArray(0);
}

/**
 * @brief Enable vertex attribute at specified index
 * @param index Vertex attribute index to enable
 *
 * Activates vertex attribute for use in shader programs.
 * Must be called before configuring attribute pointer.
 */
void VertexArray::enable_attribute(GLuint index) const {
	glEnableVertexAttribArray(index);
}

/**
 * @brief Disable vertex attribute at specified index
 * @param index Vertex attribute index to disable
 *
 * Deactivates vertex attribute to prevent unintended data access
 * in shader programs.
 */
void VertexArray::disable_attribute(GLuint index) const {
	glDisableVertexAttribArray(index);
}

/**
 * @brief Configure vertex attribute pointer
 * @param index Vertex attribute index
 * @param size Number of components per attribute (1-4)
 * @param type OpenGL data type (GL_FLOAT, GL_INT, etc.)
 * @param normalized Whether fixed-point data should be normalized
 * @param stride Byte stride between consecutive attributes
 * @param offset Byte offset to first component
 *
 * Low-level attribute configuration for maximum flexibility.
 */
void VertexArray::set_attribute_pointer(GLuint index,
										GLint size,
										GLenum type,
										GLboolean normalized,
										GLsizei stride,
										const void* offset) const {
	glVertexAttribPointer(index, size, type, normalized, stride, offset);
}

/**
 * @brief Configure floating-point vertex attribute with automatic enable
 * @param index Vertex attribute index
 * @param size Number of float components per attribute (1-4)
 * @param stride Byte stride between consecutive attributes
 * @param offset Byte offset to first component
 * @param normalized Whether data should be normalized to [0,1] or [-1,1]
 *
 * Convenient helper for common floating-point attribute configuration.
 * Automatically enables the attribute for immediate use.
 */
void VertexArray::set_float_attribute(GLuint index,
									  GLint size,
									  GLsizei stride,
									  const void* offset,
									  GLboolean normalized) const {
	// Enable attribute and configure as float array
	enable_attribute(index);
	set_attribute_pointer(index, size, GL_FLOAT, normalized, stride, offset);
}

/**
 * @brief Ensure OpenGL vertex array object exists
 *
 * Lazy initialization helper that creates OpenGL VAO only when needed,
 * optimizing resource usage and supporting deferred initialization.
 */
void VertexArray::ensure_created() const {
	// Create VAO if it doesn't exist yet
	if (vao_id_ == 0) {
		glGenVertexArrays(1, &const_cast<VertexArray*>(this)->vao_id_);
	}
}
