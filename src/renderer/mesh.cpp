/**
 * @file mesh.cpp
 * @brief 3D mesh representation and GPU upload for OpenGL rendering
 *
 * Implements the Mesh class for managing 3D geometry in rendering pipeline:
 * - Vertex and index data management with dirty flag optimization
 * - Multiple vertex format support (positions, colors, normals)
 * - Automatic GPU buffer synchronization and VAO configuration
 * - Primitive type support (triangles, lines, points)
 *
 * Essential for rendering Monte Carlo photon paths, voxel geometry,
 * and material layer visualization in the 3D scene.
 */

#include "mesh.hpp"

#include <GL/glew.h>

/**
 * @brief Construct mesh with default triangle primitive type
 *
 * Initializes mesh with OpenGL buffer objects and sets up for
 * triangle rendering. Marks GPU data as dirty for initial upload.
 */
Mesh::Mesh() : primitive_type_(PrimitiveType::Triangles), gpu_data_dirty_(true) {
	// Initialize OpenGL buffer objects for rendering
	vao_ = std::make_unique<VertexArray>();
	vertex_buffer_ = std::make_unique<VertexBuffer>();
	index_buffer_ = std::make_unique<IndexBuffer>();
}

/**
 * @brief Set complete vertex data with RenderVertex structures
 * @param vertices Vector of complete vertex data including position, color, normal
 *
 * Updates mesh with full vertex data and marks GPU buffers for re-upload.
 * Most efficient method when complete vertex data is available.
 */
void Mesh::set_vertices(const std::vector<RenderVertex>& vertices) {
	vertices_ = vertices;
	gpu_data_dirty_ = true;
}

/**
 * @brief Set vertices from position data only
 * @param positions Vector of 3D vertex positions
 *
 * Constructs vertex data using only positions with default colors and normals.
 * Suitable for simple geometry where only position matters.
 */
void Mesh::set_vertices(const std::vector<glm::vec3>& positions) {
	// Clear and reserve space for efficient construction
	vertices_.clear();
	vertices_.reserve(positions.size());

	// Convert positions to complete vertex structures
	for (const auto& position : positions) {
		vertices_.emplace_back(position);
	}

	gpu_data_dirty_ = true;
}

/**
 * @brief Set vertices from positions and colors
 * @param positions Vector of 3D vertex positions
 * @param colors Vector of RGBA vertex colors
 *
 * Constructs vertex data with positions and colors, using default normals.
 * Colors are matched by index with white fallback for missing colors.
 */
void Mesh::set_vertices(const std::vector<glm::vec3>& positions, const std::vector<glm::vec4>& colors) {
	// Prepare vertex array with known size
	vertices_.clear();
	vertices_.reserve(positions.size());

	// Combine positions and colors into vertex structures
	for (size_t i = 0; i < positions.size(); ++i) {
		glm::vec4 color = (i < colors.size()) ? colors[i] : glm::vec4(1.0f);
		vertices_.emplace_back(positions[i], color);
	}

	gpu_data_dirty_ = true;
}

/**
 * @brief Set index data for indexed rendering
 * @param indices Vector of vertex indices for triangle/line connectivity
 *
 * Updates mesh with index data for indexed rendering and marks GPU
 * buffers for re-upload. Enables efficient vertex reuse in complex geometry.
 */
void Mesh::set_indices(const std::vector<uint32_t>& indices) {
	indices_ = indices;
	gpu_data_dirty_ = true;
}

/**
 * @brief Upload mesh data to GPU buffers and configure VAO
 *
 * Synchronizes CPU mesh data with GPU buffers using dirty flag optimization.
 * Configures vertex attribute pointers for position and color data.
 * Only performs upload when mesh data has been modified.
 */
void Mesh::upload_to_gpu() {
	// Skip upload if GPU data is current
	if (!gpu_data_dirty_) {
		return;
	}

	// Upload vertex data to GPU buffer
	vertex_buffer_->upload_vertices(vertices_);

	// Upload index data if present
	if (!indices_.empty()) {
		index_buffer_->upload_indices(indices_);
	}

	// Configure vertex array object with attribute pointers
	vao_->bind();

	// Bind vertex buffer for attribute configuration
	vertex_buffer_->bind();

	// Configure position attribute (shader location 0)
	vao_->set_float_attribute(
		0, 3, sizeof(RenderVertex), reinterpret_cast<const void*>(offsetof(RenderVertex, position)));

	// Configure color attribute (shader location 1)
	vao_->set_float_attribute(1, 4, sizeof(RenderVertex), reinterpret_cast<const void*>(offsetof(RenderVertex, color)));

	// Bind index buffer for indexed rendering if available
	if (!indices_.empty()) {
		index_buffer_->bind();
	}

	// Clean up OpenGL state
	vao_->unbind();
	vertex_buffer_->unbind();
	if (!indices_.empty()) {
		index_buffer_->unbind();
	}

	gpu_data_dirty_ = false;
}

void Mesh::update_vertices() {
	if (vertices_.empty()) {
		return;
	}
	vertex_buffer_->update_vertices(vertices_);
}

void Mesh::update_indices() {
	if (indices_.empty()) {
		return;
	}
	index_buffer_->update_indices(indices_);
}

void Mesh::bind() const {
	vao_->bind();
}

void Mesh::unbind() const {
	vao_->unbind();
}

void Mesh::render() const {
	if (vertices_.empty()) {
		return;
	}

	bind();

	if (has_indices()) {
		glDrawElements(
			static_cast<GLenum>(primitive_type_), static_cast<GLsizei>(indices_.size()), GL_UNSIGNED_INT, nullptr);
	}
	else {
		glDrawArrays(static_cast<GLenum>(primitive_type_), 0, static_cast<GLsizei>(vertices_.size()));
	}

	unbind();
}
