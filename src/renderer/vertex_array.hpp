#pragma once

#include <GL/glew.h>

class VertexArray
{
public:
	VertexArray() = default;
	~VertexArray();

	void bind() const;
	void unbind() const;

	void enable_attribute(GLuint index) const;
	void disable_attribute(GLuint index) const;
	void set_attribute_pointer(GLuint index, GLint size, GLenum type, GLboolean normalized, GLsizei stride,
							   const void* offset) const;
	void set_float_attribute(GLuint index, GLint size, GLsizei stride, const void* offset,
							 GLboolean normalized = GL_FALSE) const;

	GLuint id() const { return vao_id_; }
	bool is_valid() const { return vao_id_ != 0; }

private:
	void ensure_created() const;

	mutable GLuint vao_id_ {0};
};
