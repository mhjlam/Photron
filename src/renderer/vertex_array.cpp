#include "vertex_array.hpp"

VertexArray::~VertexArray() {
    if (vao_id_ != 0) {
        glDeleteVertexArrays(1, &vao_id_);
    }
}

void VertexArray::bind() const {
    ensure_created();
    if (vao_id_ != 0) {
        glBindVertexArray(vao_id_);
    }
}

void VertexArray::unbind() const {
    glBindVertexArray(0);
}

void VertexArray::enable_attribute(GLuint index) const {
    glEnableVertexAttribArray(index);
}

void VertexArray::disable_attribute(GLuint index) const {
    glDisableVertexAttribArray(index);
}

void VertexArray::set_attribute_pointer(GLuint index, GLint size, GLenum type, GLboolean normalized, GLsizei stride, const void* offset) const {
    glVertexAttribPointer(index, size, type, normalized, stride, offset);
}

void VertexArray::set_float_attribute(GLuint index, GLint size, GLsizei stride, const void* offset, GLboolean normalized) const {
    enable_attribute(index);
    set_attribute_pointer(index, size, GL_FLOAT, normalized, stride, offset);
}

void VertexArray::ensure_created() const {
    if (vao_id_ == 0) {
        glGenVertexArrays(1, &const_cast<VertexArray*>(this)->vao_id_);
    }
}
