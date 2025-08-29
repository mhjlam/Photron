#include "buffer.hpp"

Buffer::~Buffer() {
    if (buffer_id_ != 0) {
        glDeleteBuffers(1, &buffer_id_);
    }
}

void Buffer::bind() const {
    ensure_created();
    if (buffer_id_ != 0) {
        glBindBuffer(static_cast<GLenum>(type_), buffer_id_);
    }
}

void Buffer::unbind() const {
    glBindBuffer(static_cast<GLenum>(type_), 0);
}

void Buffer::upload_data(const void* data, size_t size_bytes, BufferUsage usage) {
    ensure_created();
    bind();
    glBufferData(static_cast<GLenum>(type_), static_cast<GLsizeiptr>(size_bytes), data, static_cast<GLenum>(usage));
    size_bytes_ = size_bytes;
    unbind();
}

void Buffer::update_data(const void* data, size_t size_bytes, size_t offset) {
    bind();
    glBufferSubData(static_cast<GLenum>(type_), static_cast<GLintptr>(offset), static_cast<GLsizeiptr>(size_bytes), data);
    unbind();
}

void Buffer::ensure_created() const {
    if (buffer_id_ == 0) {
        glGenBuffers(1, &const_cast<Buffer*>(this)->buffer_id_);
    }
}
