#pragma once

#include "buffer.hpp"

/**
 * @class IndexBuffer
 * @brief Index buffer object wrapper
 */
class IndexBuffer : public Buffer {
public:
    IndexBuffer() : Buffer(BufferType::INDEX_BUFFER) {}
    ~IndexBuffer() = default;

    void upload_indices(const std::vector<unsigned int>& indices, BufferUsage usage = BufferUsage::STATIC_DRAW) {
        upload_data(indices, usage);
    }

    void update_indices(const std::vector<unsigned int>& indices, size_t offset = 0) {
        update_data(indices.data(), indices.size() * sizeof(unsigned int), offset * sizeof(unsigned int));
    }

    size_t index_count() const {
        return size() / sizeof(unsigned int);
    }
};
