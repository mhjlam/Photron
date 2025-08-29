#pragma once

#include <vector>
#include <memory>
#include <glm/glm.hpp>
#include "vertex.hpp"
#include "vertex_array.hpp"
#include "vertex_buffer.hpp"
#include "index_buffer.hpp"

enum class PrimitiveType {
    POINTS = GL_POINTS,
    LINES = GL_LINES,
    LINE_STRIP = GL_LINE_STRIP,
    TRIANGLES = GL_TRIANGLES,
    TRIANGLE_STRIP = GL_TRIANGLE_STRIP
};

/**
 * @class Mesh
 * @brief Modern OpenGL mesh for rendering
 */
class Mesh {
public:
    Mesh();
    ~Mesh() = default;

    // Vertex data setup
    void set_vertices(const std::vector<RenderVertex>& vertices);
    void set_vertices(const std::vector<glm::vec3>& positions);
    void set_vertices(const std::vector<glm::vec3>& positions, const std::vector<glm::vec4>& colors);
    void set_indices(const std::vector<unsigned int>& indices);

    // GPU upload/update
    void upload_to_gpu();
    void update_vertices();
    void update_indices();

    // Rendering
    void bind() const;
    void unbind() const;
    void render() const;

    // Getters/Setters
    void set_primitive_type(PrimitiveType type) { primitive_type_ = type; }
    PrimitiveType get_primitive_type() const { return primitive_type_; }
    size_t vertex_count() const { return vertices_.size(); }
    size_t index_count() const { return indices_.size(); }
    bool has_indices() const { return !indices_.empty(); }

private:
    std::vector<RenderVertex> vertices_;
    std::vector<unsigned int> indices_;
    
    std::unique_ptr<VertexArray> vao_;
    std::unique_ptr<VertexBuffer> vertex_buffer_;
    std::unique_ptr<IndexBuffer> index_buffer_;
    
    PrimitiveType primitive_type_;
    bool gpu_data_dirty_;
};
