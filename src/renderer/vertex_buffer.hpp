/**
 * @file vertex_buffer.hpp
 * @brief OpenGL vertex buffer specialization for vertex attribute data
 *
 * Provides a specialized buffer class for managing vertex attribute data
 * in OpenGL rendering operations, with type-safe vertex format handling.
 */

#pragma once

#include "renderer/buffer.hpp"
#include "renderer/vertex.hpp"

/**
 * @class VertexBuffer
 * @brief Specialized buffer for OpenGL vertex attribute data management
 *
 * Inherits from Buffer and provides vertex-specific functionality for uploading
 * and managing RenderVertex data in OpenGL rendering operations. Automatically
 * handles vertex data formatting and size calculations.
 */
class VertexBuffer : public Buffer
{
public:
	/**
	 * @brief Construct a vertex buffer
	 *
	 * Creates a vertex buffer specialized for GL_ARRAY_BUFFER usage.
	 */
	VertexBuffer() : Buffer(BufferType::Vertex) {}

	/**
	 * @brief Default destructor
	 *
	 * Resource cleanup is handled by the base Buffer class.
	 */
	~VertexBuffer() = default;

	/**
	 * @brief Upload vertex data to the GPU buffer
	 *
	 * Uploads a vector of RenderVertex structures to the OpenGL buffer.
	 * Automatically calculates the correct buffer size and data layout.
	 *
	 * @param vertices Vector of RenderVertex data to upload
	 * @param usage Buffer usage hint for OpenGL optimization (default: Static)
	 */
	void upload_vertices(const std::vector<RenderVertex>& vertices, BufferUsage usage = BufferUsage::Static) {
		upload_data(vertices, usage);
	}

	/**
	 * @brief Update existing vertex data in the buffer
	 *
	 * Updates a portion of the existing buffer data without reallocating.
	 * Useful for dynamic vertex data updates during rendering.
	 *
	 * @param vertices Vector of new vertex data
	 * @param offset Vertex offset (not byte offset) where update begins (default: 0)
	 */
	void update_vertices(const std::vector<RenderVertex>& vertices, size_t offset = 0) {
		update_data(vertices.data(), vertices.size() * sizeof(RenderVertex), offset * sizeof(RenderVertex));
	}

	/**
	 * @brief Get the number of vertices currently stored in the buffer
	 *
	 * Calculates vertex count by dividing total buffer size by RenderVertex size.
	 *
	 * @return size_t Number of vertices in the buffer
	 */
	size_t vertex_count() const { return size() / sizeof(RenderVertex); }
};
