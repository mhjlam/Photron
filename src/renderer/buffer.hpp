/**
 * @file buffer.hpp
 * @brief RAII-based OpenGL buffer object management
 *
 * Provides modern C++ wrapper for OpenGL buffer objects with automatic
 * resource management, type safety, and efficient data operations.
 */

#pragma once

#include <concepts>
#include <span>
#include <type_traits>
#include <vector>

#include <GL/glew.h>

/**
 * @enum BufferType
 * @brief OpenGL buffer object type enumeration
 */
enum class BufferType
{
	Vertex = GL_ARRAY_BUFFER,       ///< Vertex buffer for vertex attribute data
	Index = GL_ELEMENT_ARRAY_BUFFER ///< Index buffer for indexed drawing
};

/**
 * @enum BufferUsage
 * @brief OpenGL buffer usage pattern hint enumeration
 */
enum class BufferUsage
{
	Static = GL_STATIC_DRAW,   ///< Data will be modified once and used many times
	Dynamic = GL_DYNAMIC_DRAW, ///< Data will be modified repeatedly and used many times
	Stream = GL_STREAM_DRAW    ///< Data will be modified once and used at most a few times
};

// C++20 concepts for type safety
template<typename T>
concept TriviallyCopiable = std::is_trivially_copyable_v<T>;

template<typename T>
concept ContiguousRange = std::ranges::contiguous_range<T> && std::ranges::sized_range<T>;

/**
 * @class Buffer
 * @brief RAII-based OpenGL buffer object wrapper with modern C++20 features
 *
 * Provides type-safe, efficient buffer management with automatic resource cleanup,
 * C++20 concepts for type safety, and multiple data upload interfaces.
 */
class Buffer
{
public:
	/**
	 * @brief Construct buffer with specified type
	 * @param type Buffer type (vertex or index buffer)
	 */
	explicit Buffer(BufferType type) noexcept : type_(type) {}

	/**
	 * @brief Destructor - automatically cleans up OpenGL buffer
	 */
	~Buffer();

	/**
	 * @brief Bind this buffer for subsequent OpenGL operations
	 */
	void bind() const;

	/**
	 * @brief Unbind the current buffer (bind 0)
	 */
	void unbind() const;

	/**
	 * @brief Upload raw data to the buffer
	 * @param data Pointer to data to upload
	 * @param size_bytes Size of data in bytes
	 * @param usage Buffer usage hint for optimization
	 */
	void upload_data(const void* data, size_t size_bytes, BufferUsage usage);

	/**
	 * @brief Update existing buffer data
	 * @param data Pointer to new data
	 * @param size_bytes Size of data in bytes
	 * @param offset Byte offset into buffer (default: 0)
	 */
	void update_data(const void* data, size_t size_bytes, size_t offset = 0);

	/**
	 * @brief Upload vector data to buffer (C++20 concepts)
	 * @tparam T Type that satisfies TriviallyCopiable concept
	 * @param data Vector of data to upload
	 * @param usage Buffer usage hint for optimization
	 */
	template<TriviallyCopiable T>
	void upload_data(const std::vector<T>& data, BufferUsage usage) {
		upload_data(data.data(), data.size() * sizeof(T), usage);
	}

	/**
	 * @brief Upload span data to buffer (C++20 zero-copy interface)
	 * @tparam T Type that satisfies TriviallyCopiable concept
	 * @param data Span of data to upload
	 * @param usage Buffer usage hint for optimization
	 */
	template<TriviallyCopiable T>
	void upload_data(std::span<const T> data, BufferUsage usage) {
		upload_data(data.data(), data.size() * sizeof(T), usage);
	}

	/**
	 * @brief Upload array data to buffer (C++20 concepts-constrained)
	 * @tparam T Type that satisfies TriviallyCopiable concept
	 * @param data Pointer to array data
	 * @param count Number of elements to upload
	 * @param usage Buffer usage hint for optimization
	 */
	template<TriviallyCopiable T>
	void upload_data(const T* data, size_t count, BufferUsage usage) {
		upload_data(static_cast<const void*>(data), count * sizeof(T), usage);
	}

	/**
	 * @brief Upload contiguous range data to buffer (C++20 generic ranges)
	 * @tparam R Range type that satisfies ContiguousRange concept
	 * @param range Contiguous range of data to upload
	 * @param usage Buffer usage hint for optimization
	 */
	template<ContiguousRange R> 
	requires TriviallyCopiable<std::ranges::range_value_t<R>>
	void upload_data(const R& range, BufferUsage usage) {
		using ValueType = std::ranges::range_value_t<R>;
		upload_data(std::ranges::data(range), std::ranges::size(range) * sizeof(ValueType), usage);
	}

	/**
	 * @brief Get the OpenGL buffer ID
	 * @return GLuint Buffer identifier (0 if not yet created)
	 */
	GLuint id() const noexcept { return buffer_id_; }

	/**
	 * @brief Check if the buffer has been created and is valid
	 * @return bool True if buffer is valid, false otherwise
	 */
	bool is_valid() const noexcept { return buffer_id_ != 0; }

	/**
	 * @brief Get the current buffer size in bytes
	 * @return size_t Buffer size in bytes
	 */
	size_t size() const noexcept { return size_bytes_; }

private:
	/**
	 * @brief Lazy initialization - create buffer if not yet created
	 */
	void ensure_created() const;

	mutable GLuint buffer_id_ {0}; ///< OpenGL buffer identifier (0 = not created)
	BufferType type_;              ///< Buffer type (vertex or index)
	size_t size_bytes_ {0};        ///< Current buffer size in bytes
};
