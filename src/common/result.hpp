/**
 * @file result.hpp
 * @brief Modern C++ Result type for structured error handling
 *
 * Provides a Rust-style Result<T, E> type that eliminates the need for
 * exceptions, mixed return types, or output parameters. Enables explicit
 * and type-safe error handling throughout the codebase.
 */

#pragma once

#include <string>
#include <utility>
#include <variant>

/**
 * @class Result
 * @brief Type-safe result container for operations that can succeed or fail
 *
 * The Result class provides a modern alternative to exception-based error handling
 * by explicitly representing success and failure states in the type system.
 * This approach offers several advantages:
 *
 * - **Explicit Error Handling**: Callers must explicitly handle both success and failure cases
 * - **No Exception Overhead**: Zero-cost abstraction with no exception handling overhead
 * - **Composable**: Easy to chain operations and propagate errors upward
 * - **Type Safe**: Compile-time guarantees about error handling completeness
 *
 * Based on Rust's Result<T, E> and C++23's proposed std::expected.
 *
 * Usage examples:
 * @code
 * Result<int, std::string> divide(int a, int b) {
 *     if (b == 0) return Result<int, std::string>::error("Division by zero");
 *     return Result<int, std::string>::ok(a / b);
 * }
 *
 * auto result = divide(10, 2);
 * if (result.is_ok()) {
 *     std::cout << "Result: " << result.value() << std::endl;
 * } else {
 *     std::cout << "Error: " << result.error() << std::endl;
 * }
 * @endcode
 *
 * @tparam T Type of the success value
 * @tparam E Type of the error value
 */
template<typename T, typename E>
class Result
{
public:
	// Private constructors - use static factory methods for clarity
private:
	template<typename U>
	explicit Result(U&& value) : data_(std::forward<U>(value)) {}

public:
	/**
	 * @brief Create a successful result containing a value
	 *
	 * Static factory method that constructs a Result in the success state.
	 * Preferred over direct construction for API clarity.
	 *
	 * @tparam U Type of the success value (deduced)
	 * @param value The success value to store
	 * @return Result<T,E> Result in success state
	 */
	template<typename U = T>
	static Result ok(U&& value) {
		return Result(std::forward<U>(value));
	}

	/**
	 * @brief Create an error result containing an error value
	 *
	 * Static factory method that constructs a Result in the error state.
	 * Preferred over direct construction for API clarity.
	 *
	 * @tparam F Type of the error value (deduced)
	 * @param err The error value to store
	 * @return Result<T,E> Result in error state
	 */
	template<typename F = E>
	static Result error(F&& err) {
		return Result(std::forward<F>(err));
	}

	/**
	 * @brief Check if the result represents success
	 * @return true if result contains a success value, false if error
	 */
	bool is_ok() const { return std::holds_alternative<T>(data_); }

	/**
	 * @brief Check if the result represents an error
	 * @return true if result contains an error value, false if success
	 */
	bool is_error() const { return std::holds_alternative<E>(data_); }

	/**
	 * @brief Access the success value (const lvalue reference)
	 *
	 * @return const T& Reference to the success value
	 * @throws std::bad_variant_access if result is in error state
	 */
	const T& value() const& { return std::get<T>(data_); }

	/**
	 * @brief Access the success value (mutable lvalue reference)
	 *
	 * @return T& Mutable reference to the success value
	 * @throws std::bad_variant_access if result is in error state
	 */
	T& value() & { return std::get<T>(data_); }

	/**
	 * @brief Access the success value (rvalue reference)
	 *
	 * @return T&& Rvalue reference to the success value for move semantics
	 * @throws std::bad_variant_access if result is in error state
	 */
	T&& value() && { return std::get<T>(std::move(data_)); }

	/**
	 * @brief Access the error value (const lvalue reference)
	 *
	 * @return const E& Reference to the error value
	 * @throws std::bad_variant_access if result is in success state
	 */
	const E& error() const& { return std::get<E>(data_); }

	/**
	 * @brief Access the error value (mutable lvalue reference)
	 *
	 * @return E& Mutable reference to the error value
	 * @throws std::bad_variant_access if result is in success state
	 */
	E& error() & { return std::get<E>(data_); }

	/**
	 * @brief Access the error value (rvalue reference)
	 *
	 * @return E&& Rvalue reference to the error value for move semantics
	 * @throws std::bad_variant_access if result is in success state
	 */
	E&& error() && { return std::get<E>(std::move(data_)); }

	/**
	 * @brief Get success value or return default if error
	 *
	 * Provides safe access to the success value with a fallback.
	 * Useful for cases where a reasonable default exists.
	 *
	 * @tparam U Type of the default value
	 * @param default_value Default value to return if result is error
	 * @return T Success value if available, otherwise the default value
	 */
	template<typename U>
	T value_or(U&& default_value) const& {
		return is_ok() ? value() : static_cast<T>(std::forward<U>(default_value));
	}

	template<typename U>
	T value_or(U&& default_value) && {
		return is_ok() ? std::move(*this).value() : static_cast<T>(std::forward<U>(default_value));
	}

	// Convenience operators
	explicit operator bool() const { return is_ok(); }

private:
	std::variant<T, E> data_;
};

// Specialization for void success type
template<typename E>
class Result<void, E>
{
public:
	// Constructors
	Result() : data_(std::monostate {}) {}

	template<typename F = E>
	Result(F&& error) : data_(std::forward<F>(error)) {}

	// Static factory methods
	static Result ok() { return Result(); }

	template<typename F = E>
	static Result error(F&& err) {
		return Result(std::forward<F>(err));
	}

	// Check if result is success or error
	bool is_ok() const { return std::holds_alternative<std::monostate>(data_); }
	bool is_error() const { return std::holds_alternative<E>(data_); }

	// Access error (throws if success)
	const E& error() const& { return std::get<E>(data_); }
	E& error() & { return std::get<E>(data_); }
	E&& error() && { return std::get<E>(std::move(data_)); }

	// Convenience operators
	explicit operator bool() const { return is_ok(); }

private:
	std::variant<std::monostate, E> data_;
};