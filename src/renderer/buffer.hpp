#pragma once

#include <GL/glew.h>
#include <vector>

enum class BufferType {
    VERTEX_BUFFER = GL_ARRAY_BUFFER,
    INDEX_BUFFER = GL_ELEMENT_ARRAY_BUFFER
};

enum class BufferUsage {
    STATIC_DRAW = GL_STATIC_DRAW,
    DYNAMIC_DRAW = GL_DYNAMIC_DRAW,
    STREAM_DRAW = GL_STREAM_DRAW
};

/**
 * @class Buffer
 * @brief Base class for OpenGL buffer objects
 */
class Buffer {
public:
    explicit Buffer(BufferType type) : buffer_id_(0), type_(type) {}
    ~Buffer();

    void bind() const;
    void unbind() const;
    
    void upload_data(const void* data, size_t size_bytes, BufferUsage usage);
    void update_data(const void* data, size_t size_bytes, size_t offset = 0);

    template<typename T>
    void upload_data(const std::vector<T>& data, BufferUsage usage) {
        upload_data(data.data(), data.size() * sizeof(T), usage);
    }

    template<typename T>
    void upload_data(const T* data, size_t count, BufferUsage usage) {
        upload_data(static_cast<const void*>(data), count * sizeof(T), usage);
    }

    GLuint id() const { return buffer_id_; }
    bool is_valid() const { return buffer_id_ != 0; }
    size_t size() const { return size_bytes_; }

private:
    void ensure_created() const;

    mutable GLuint buffer_id_;
    BufferType type_;
    size_t size_bytes_ = 0;
};
