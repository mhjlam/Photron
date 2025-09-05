#pragma once

#include <vector>
#include <span>
#include <concepts>
#include <type_traits>

#include <GL/glew.h>

enum class BufferType
{
	Vertex = GL_ARRAY_BUFFER,
	Index = GL_ELEMENT_ARRAY_BUFFER
};

enum class BufferUsage
{
	Static = GL_STATIC_DRAW,
	Dynamic = GL_DYNAMIC_DRAW,
	Stream = GL_STREAM_DRAW
};

// C++20 concepts for type safety
template<typename T>
concept TriviallyCopiable = std::is_trivially_copyable_v<T>;

template<typename T>
concept ContiguousRange = std::ranges::contiguous_range<T> && std::ranges::sized_range<T>;

class Buffer
{
public:
	explicit Buffer(BufferType type) noexcept : type_(type) {}
	~Buffer();

	void bind() const;
	void unbind() const;

	void upload_data(const void* data, size_t size_bytes, BufferUsage usage);
	void update_data(const void* data, size_t size_bytes, size_t offset = 0);

	// C++20: Modern template with concepts for type safety
	template<TriviallyCopiable T>
	void upload_data(const std::vector<T>& data, BufferUsage usage) {
		upload_data(data.data(), data.size() * sizeof(T), usage);
	}

	// C++20: Span interface for zero-copy operations
	template<TriviallyCopiable T>
	void upload_data(std::span<const T> data, BufferUsage usage) {
		upload_data(data.data(), data.size() * sizeof(T), usage);
	}

	// C++20: Concepts-constrained raw pointer interface
	template<TriviallyCopiable T>
	void upload_data(const T* data, size_t count, BufferUsage usage) {
		upload_data(static_cast<const void*>(data), count * sizeof(T), usage);
	}

	// C++20: Generic contiguous range support
	template<ContiguousRange R>
		requires TriviallyCopiable<std::ranges::range_value_t<R>>
	void upload_data(const R& range, BufferUsage usage) {
		using ValueType = std::ranges::range_value_t<R>;
		upload_data(std::ranges::data(range), std::ranges::size(range) * sizeof(ValueType), usage);
	}

	GLuint id() const noexcept { return buffer_id_; }
	bool is_valid() const noexcept { return buffer_id_ != 0; }
	size_t size() const noexcept { return size_bytes_; }

	private:
	void ensure_created() const;

	mutable GLuint buffer_id_ {0};
	BufferType type_;
	size_t size_bytes_ {0};
};
