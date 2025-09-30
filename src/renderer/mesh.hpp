/**
 * @file mesh.hpp
 * @brief OpenGL mesh abstraction for efficient 3D rendering
 *
 * Provides a high-level interface for managing vertex data, index buffers,
 * and OpenGL rendering state. Supports multiple primitive types and optimized
 * GPU memory management for real-time 3D visualization.
 */

#pragma once

// Standard library includes
#include <memory>
#include <vector>

#include <GL/glew.h>
#include <glm/glm.hpp>

#include "renderer/index_buffer.hpp"
#include "renderer/vertex.hpp"
#include "renderer/vertex_array.hpp"
#include "renderer/vertex_buffer.hpp"

// Forward declare OpenGL types and constants
using GLuint = unsigned int;

/**
 * @enum PrimitiveType
 * @brief OpenGL primitive rendering modes
 *
 * Defines supported primitive types for mesh rendering with direct
 * mapping to OpenGL primitive constants for efficient rendering.
 */
enum class PrimitiveType
{
	Points = GL_POINTS,               ///< Individual point sprites
	Lines = GL_LINES,                 ///< Disconnected line segments
	LineStrip = GL_LINE_STRIP,        ///< Connected line sequence
	Triangles = GL_TRIANGLES,         ///< Individual triangles (most common)
	TriangleStrip = GL_TRIANGLE_STRIP ///< Connected triangle strip (optimized)
};

/**
 * @class Mesh
 * @brief High-level OpenGL mesh abstraction with automatic resource management
 *
 * The Mesh class provides a complete abstraction over OpenGL vertex data management,
 * including vertex buffers, index buffers, and vertex array objects. Key features:
 *
 * **Vertex Data Management:**
 * - Support for multiple vertex data formats (position-only, position+color)
 * - Automatic vertex attribute configuration
 * - Efficient GPU memory updates with dirty tracking
 *
 * **Index Buffer Support:**
 * - Optional indexed rendering for memory efficiency
 * - Automatic index buffer management
 * - Support for 32-bit indices for large meshes
 *
 * **Rendering Optimization:**
 * - Lazy GPU uploads to minimize bandwidth usage
 * - Dirty state tracking to avoid unnecessary uploads
 * - Multiple primitive type support for different visualization needs
 *
 * **Resource Management:**
 * - RAII-based OpenGL resource cleanup
 * - Smart pointer management for automatic memory handling
 * - Exception-safe resource allocation
 *
 * **Usage Examples:**
 * ```cpp
 * Mesh voxel_mesh;
 * voxel_mesh.set_vertices(positions, colors);
 * voxel_mesh.set_primitive_type(PrimitiveType::Triangles);
 * voxel_mesh.upload_to_gpu();
 * voxel_mesh.render();
 * ```
 */
class Mesh
{
public:
	/**
	 * @brief Construct new Mesh with default settings
	 *
	 * Initializes mesh with triangle primitive type and creates
	 * OpenGL resource wrappers (VAO, VBO, IBO).
	 */
	Mesh();

	/**
	 * @brief Destructor handles automatic OpenGL resource cleanup
	 */
	~Mesh() = default;

	// Vertex data configuration methods

	/**
	 * @brief Set vertex data using complete RenderVertex structures
	 *
	 * Most flexible method supporting position, color, normal, and texture data.
	 * Marks GPU data as dirty for next render call.
	 *
	 * @param vertices Vector of complete vertex data structures
	 */
	void set_vertices(const std::vector<RenderVertex>& vertices);

	/**
	 * @brief Set vertex data using position-only data
	 *
	 * Convenience method for simple geometry. Colors default to white,
	 * normals to (0,0,1), and texture coordinates to (0,0).
	 *
	 * @param positions Vector of 3D vertex positions
	 */
	void set_vertices(const std::vector<glm::vec3>& positions);

	/**
	 * @brief Set vertex data with positions and per-vertex colors
	 *
	 * Common method for colored geometry visualization. Normals default
	 * to (0,0,1) and texture coordinates to (0,0).
	 *
	 * @param positions Vector of 3D vertex positions
	 * @param colors Vector of RGBA colors (must match positions size)
	 */
	void set_vertices(const std::vector<glm::vec3>& positions, const std::vector<glm::vec4>& colors);

	/**
	 * @brief Set index buffer for indexed rendering
	 *
	 * Enables indexed rendering for memory efficiency. Indices reference
	 * vertices in the vertex buffer. Marks GPU data as dirty.
	 *
	 * @param indices Vector of vertex indices (32-bit unsigned integers)
	 */
	void set_indices(const std::vector<uint32_t>& indices);

	// GPU memory management

	/**
	 * @brief Upload all mesh data to GPU buffers
	 *
	 * Uploads vertex data and index data (if present) to GPU memory.
	 * Call after setting vertex/index data and before rendering.
	 * Automatically handles buffer creation and binding.
	 */
	void upload_to_gpu();

	/**
	 * @brief Update vertex data on GPU without recreating buffers
	 *
	 * Efficient method for updating vertex data when buffer size hasn't changed.
	 * Faster than full upload_to_gpu() for dynamic meshes.
	 */
	void update_vertices();

	/**
	 * @brief Update index data on GPU without recreating buffers
	 *
	 * Efficient method for updating index data when buffer size hasn't changed.
	 */
	void update_indices();

	// Rendering control

	/**
	 * @brief Bind mesh for rendering (activate VAO)
	 *
	 * Must be called before issuing draw commands. Sets up all vertex
	 * attribute pointers and binds appropriate buffers.
	 */
	void bind() const;

	/**
	 * @brief Unbind mesh after rendering (restore default VAO)
	 */
	void unbind() const;

	/**
	 * @brief Render the mesh using current OpenGL state
	 *
	 * Issues appropriate draw call based on whether indices are present.
	 * Mesh must be bound and shaders must be active before calling.
	 */
	void render() const;

	// Configuration and query methods

	/**
	 * @brief Set primitive rendering type
	 * @param type Primitive type for rendering (points, lines, triangles, etc.)
	 */
	void set_primitive_type(PrimitiveType type) { primitive_type_ = type; }

	/**
	 * @brief Get current primitive rendering type
	 * @return PrimitiveType Current primitive type
	 */
	PrimitiveType get_primitive_type() const { return primitive_type_; }

	/**
	 * @brief Get number of vertices in mesh
	 * @return size_t Vertex count
	 */
	size_t vertex_count() const { return vertices_.size(); }

	/**
	 * @brief Get number of indices in mesh
	 * @return size_t Index count (0 if no indices)
	 */
	size_t index_count() const { return indices_.size(); }

	/**
	 * @brief Check if mesh uses indexed rendering
	 * @return true if mesh has index buffer, false for array rendering
	 */
	bool has_indices() const { return !indices_.empty(); }

private:
	// CPU-side data storage
	std::vector<RenderVertex> vertices_; ///< Vertex data array
	std::vector<uint32_t> indices_;      ///< Index data array (optional)

	// OpenGL resource wrappers
	std::unique_ptr<VertexArray> vao_;            ///< Vertex Array Object - manages vertex state
	std::unique_ptr<VertexBuffer> vertex_buffer_; ///< Vertex Buffer Object - stores vertex data on GPU
	std::unique_ptr<IndexBuffer> index_buffer_;   ///< Index Buffer Object - stores indices on GPU

	// Rendering configuration
	PrimitiveType primitive_type_; ///< Current primitive rendering type
	bool gpu_data_dirty_;          ///< Flag indicating if GPU data needs upload
};
