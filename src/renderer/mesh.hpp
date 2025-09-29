#pragma once

// Standard library includes
#include <memory>
#include <vector>

// Third-party includes for interface types
#include <GL/glew.h>  // For GL constants
#include <glm/glm.hpp>

// Project includes for member types
#include "renderer/index_buffer.hpp"
#include "renderer/vertex.hpp"  // For RenderVertex type
#include "renderer/vertex_array.hpp"
#include "renderer/vertex_buffer.hpp"

// Forward declare OpenGL types and constants
using GLuint = unsigned int;

enum class PrimitiveType
{
	Points = GL_POINTS,
	Lines = GL_LINES,
	LineStrip = GL_LINE_STRIP,
	Triangles = GL_TRIANGLES,
	TriangleStrip = GL_TRIANGLE_STRIP
};

class Mesh
{
public:
	Mesh();
	~Mesh() = default;

	// Vertex data setup
	void set_vertices(const std::vector<RenderVertex>& vertices);
	void set_vertices(const std::vector<glm::vec3>& positions);
	void set_vertices(const std::vector<glm::vec3>& positions, const std::vector<glm::vec4>& colors);
	void set_indices(const std::vector<uint32_t>& indices);

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
	std::vector<uint32_t> indices_;

	std::unique_ptr<VertexArray> vao_;
	std::unique_ptr<VertexBuffer> vertex_buffer_;
	std::unique_ptr<IndexBuffer> index_buffer_;

	PrimitiveType primitive_type_;
	bool gpu_data_dirty_;
};
