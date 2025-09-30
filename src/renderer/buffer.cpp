/**
 * @file buffer.cpp
 * @brief OpenGL buffer object management and data operations
 *
 * Implements RAII-based OpenGL buffer management for rendering operations:
 * - Automatic buffer creation and cleanup
 * - Type-safe buffer binding and operations
 * - Efficient data upload and update mechanisms
 * - Memory management with proper resource cleanup
 *
 * Essential for managing vertex data, index data, and uniform buffers
 * in the Monte Carlo photon visualization pipeline.
 */

#include "buffer.hpp"

/**
 * @brief Destructor ensuring proper OpenGL resource cleanup
 *
 * Automatically deletes OpenGL buffer object when Buffer instance
 * goes out of scope, preventing resource leaks.
 */
Buffer::~Buffer() {
	// Clean up OpenGL buffer if it was created
	if (buffer_id_ != 0) {
		glDeleteBuffers(1, &buffer_id_);
	}
}

/**
 * @brief Bind buffer to OpenGL context
 *
 * Makes this buffer the active buffer for its type, enabling
 * subsequent OpenGL operations to target this buffer.
 */
void Buffer::bind() const {
	// Ensure buffer exists before binding
	ensure_created();
	if (buffer_id_ != 0) {
		glBindBuffer(static_cast<GLenum>(type_), buffer_id_);
	}
}

/**
 * @brief Unbind buffer from OpenGL context
 *
 * Clears the active buffer binding for this buffer type,
 * preventing accidental operations on this buffer.
 */
void Buffer::unbind() const {
	glBindBuffer(static_cast<GLenum>(type_), 0);
}

/**
 * @brief Upload data to GPU buffer memory
 * @param data Pointer to source data
 * @param size_bytes Size of data in bytes
 * @param usage Buffer usage hint for OpenGL optimization
 *
 * Allocates GPU memory and uploads data with usage hints for
 * optimal performance based on access patterns.
 */
void Buffer::upload_data(const void* data, size_t size_bytes, BufferUsage usage) {
	// Ensure buffer exists and bind for data operations
	ensure_created();
	bind();

	// Upload data with appropriate usage hint
	glBufferData(static_cast<GLenum>(type_), static_cast<GLsizeiptr>(size_bytes), data, static_cast<GLenum>(usage));
	size_bytes_ = size_bytes;

	// Unbind after operation to prevent accidental modifications
	unbind();
}

/**
 * @brief Update subset of buffer data
 * @param data Pointer to new data
 * @param size_bytes Size of data to update
 * @param offset Byte offset into buffer
 *
 * Efficiently updates portion of existing buffer without reallocation.
 * Useful for dynamic data like particle positions or animation.
 */
void Buffer::update_data(const void* data, size_t size_bytes, size_t offset) {
	// Bind buffer for sub-data update
	bind();
	glBufferSubData(
		static_cast<GLenum>(type_), static_cast<GLintptr>(offset), static_cast<GLsizeiptr>(size_bytes), data);

	// Unbind after update operation
	unbind();
}

/**
 * @brief Ensure OpenGL buffer object exists
 *
 * Lazy initialization helper that creates OpenGL buffer object
 * only when needed, optimizing resource usage.
 */
void Buffer::ensure_created() const {
	// Create buffer if it doesn't exist yet
	if (buffer_id_ == 0) {
		glGenBuffers(1, &const_cast<Buffer*>(this)->buffer_id_);
	}
}
