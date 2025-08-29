#pragma once

#include "buffer.hpp"
#include "vertex.hpp"

/**
 * @class VertexBuffer
 * @brief Vertex buffer object wrapper
 */
class VertexBuffer : public Buffer {
public:
    VertexBuffer() : Buffer(BufferType::VERTEX_BUFFER) {}
    ~VertexBuffer() = default;

    void upload_vertices(const std::vector<RenderVertex>& vertices, BufferUsage usage = BufferUsage::STATIC_DRAW) {
        upload_data(vertices, usage);
    }

    void update_vertices(const std::vector<RenderVertex>& vertices, size_t offset = 0) {
        update_data(vertices.data(), vertices.size() * sizeof(RenderVertex), offset * sizeof(RenderVertex));
    }

    size_t vertex_count() const {
        return size() / sizeof(RenderVertex);
    }
};
