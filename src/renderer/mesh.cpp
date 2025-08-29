#include "mesh.hpp"
#include <GL/glew.h>

Mesh::Mesh() : primitive_type_(PrimitiveType::TRIANGLES), gpu_data_dirty_(true) {
    vao_ = std::make_unique<VertexArray>();
    vertex_buffer_ = std::make_unique<VertexBuffer>();
    index_buffer_ = std::make_unique<IndexBuffer>();
}

void Mesh::set_vertices(const std::vector<RenderVertex>& vertices) {
    vertices_ = vertices;
    gpu_data_dirty_ = true;
}

void Mesh::set_vertices(const std::vector<glm::vec3>& positions) {
    vertices_.clear();
    vertices_.reserve(positions.size());
    
    for (const auto& pos : positions) {
        vertices_.emplace_back(pos);
    }

    gpu_data_dirty_ = true;
}

void Mesh::set_vertices(const std::vector<glm::vec3>& positions, const std::vector<glm::vec4>& colors) {
    vertices_.clear();
    vertices_.reserve(positions.size());

    for (size_t i = 0; i < positions.size(); ++i) {
        glm::vec4 color = (i < colors.size()) ? colors[i] : glm::vec4(1.0f);
        vertices_.emplace_back(positions[i], color);
    }

    gpu_data_dirty_ = true;
}

void Mesh::set_indices(const std::vector<unsigned int>& indices) {
    indices_ = indices;
    gpu_data_dirty_ = true;
}

void Mesh::upload_to_gpu() {
    if (!gpu_data_dirty_) return;

    // Upload vertices
    vertex_buffer_->upload_vertices(vertices_);

    // Upload indices if we have them
    if (!indices_.empty()) {
        index_buffer_->upload_indices(indices_);
    }

    // Setup VAO
    vao_->bind();

    // Bind vertex buffer and setup attributes
    vertex_buffer_->bind();
    
    // Position attribute (location 0)
    vao_->set_float_attribute(0, 3, sizeof(RenderVertex), reinterpret_cast<const void*>(offsetof(RenderVertex, position)));
    
    // Normal attribute (location 1) - we'll use this for color for now
    vao_->set_float_attribute(1, 4, sizeof(RenderVertex), reinterpret_cast<const void*>(offsetof(RenderVertex, color)));

    // Bind index buffer if we have indices
    if (!indices_.empty()) {
        index_buffer_->bind();
    }

    vao_->unbind();
    vertex_buffer_->unbind();
    if (!indices_.empty()) {
        index_buffer_->unbind();
    }

    gpu_data_dirty_ = false;
}

void Mesh::update_vertices() {
    if (vertices_.empty()) return;
    vertex_buffer_->update_vertices(vertices_);
}

void Mesh::update_indices() {
    if (indices_.empty()) return;
    index_buffer_->update_indices(indices_);
}

void Mesh::bind() const {
    vao_->bind();
}

void Mesh::unbind() const {
    vao_->unbind();
}

void Mesh::render() const {
    if (vertices_.empty()) return;

    bind();

    if (has_indices()) {
        glDrawElements(static_cast<GLenum>(primitive_type_), static_cast<GLsizei>(indices_.size()), GL_UNSIGNED_INT, nullptr);
    } else {
        glDrawArrays(static_cast<GLenum>(primitive_type_), 0, static_cast<GLsizei>(vertices_.size()));
    }

    unbind();
}
