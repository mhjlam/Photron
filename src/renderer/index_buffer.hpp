/**
 * @file index_buffer.hpp
 * @brief OpenGL index buffer specialization for efficient mesh rendering
 *
 * Provides a specialized buffer class for managing vertex indices in OpenGL
 * rendering operations, optimized for triangle mesh data.
 */

#pragma once

#include "buffer.hpp"

/**
 * @class IndexBuffer
 * @brief Specialized buffer for OpenGL index data management
 * 
 * Inherits from Buffer and provides index-specific functionality for uploading
 * and managing vertex indices in OpenGL rendering operations. Optimized for
 * triangle mesh rendering with 32-bit unsigned integer indices.
 */
class IndexBuffer : public Buffer
{
public:
	/**
	 * @brief Construct an index buffer
	 * 
	 * Creates an index buffer specialized for GL_ELEMENT_ARRAY_BUFFER usage
	 * for efficient indexed triangle rendering.
	 */
	IndexBuffer() : Buffer(BufferType::Index) {}
	
	/**
	 * @brief Default destructor
	 * 
	 * Resource cleanup is handled by the base Buffer class.
	 */
	~IndexBuffer() = default;

	/**
	 * @brief Upload index data to the GPU buffer
	 * 
	 * Uploads a vector of 32-bit unsigned integer indices to the OpenGL buffer.
	 * Indices are used for efficient triangle mesh rendering by referencing
	 * vertices in a vertex buffer.
	 * 
	 * @param indices Vector of vertex indices (uint32_t) to upload
	 * @param usage Buffer usage hint for OpenGL optimization (default: Static)
	 */
	void upload_indices(const std::vector<uint32_t>& indices, BufferUsage usage = BufferUsage::Static) {
		upload_data(indices, usage);
	}

	/**
	 * @brief Update existing index data in the buffer
	 * 
	 * Updates a portion of the existing buffer data without reallocating.
	 * Useful for dynamic mesh updates during rendering operations.
	 * 
	 * @param indices Vector of new index data
	 * @param offset Index offset (not byte offset) where update begins (default: 0)
	 */
	void update_indices(const std::vector<uint32_t>& indices, size_t offset = 0) {
		update_data(indices.data(), indices.size() * sizeof(uint32_t), offset * sizeof(uint32_t));
	}

	/**
	 * @brief Get the number of indices currently stored in the buffer
	 * 
	 * Calculates index count by dividing total buffer size by uint32_t size.
	 * For triangle meshes, divide by 3 to get triangle count.
	 * 
	 * @return size_t Number of indices in the buffer
	 */
	size_t index_count() const { return size() / sizeof(uint32_t); }
};
