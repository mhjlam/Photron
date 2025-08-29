#pragma once

#include "buffer.hpp"

class IndexBuffer : public Buffer
{
public:
	IndexBuffer() : Buffer(BufferType::Index) {}
	~IndexBuffer() = default;

	void upload_indices(const std::vector<uint32_t>& indices, BufferUsage usage = BufferUsage::Static) {
		upload_data(indices, usage);
	}

	void update_indices(const std::vector<uint32_t>& indices, size_t offset = 0) {
		update_data(indices.data(), indices.size() * sizeof(uint32_t), offset * sizeof(uint32_t));
	}

	size_t index_count() const { return size() / sizeof(uint32_t); }
};
