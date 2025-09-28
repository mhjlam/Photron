#pragma once

#include <variant>
#include <string>
#include <utility>

/**
 * @brief Result type for structured error handling
 * 
 * Similar to Rust's Result<T, E> or C++23's std::expected.
 * Replaces mixed bool/exception/void return patterns with consistent error handling.
 */
template<typename T, typename E>
class Result {
public:
    // Private constructors - use static factory methods for clarity
private:
    template<typename U>
    explicit Result(U&& value) : data_(std::forward<U>(value)) {}
    
public:
    // Static factory methods for clarity
    template<typename U = T>
    static Result ok(U&& value) {
        return Result(std::forward<U>(value));
    }
    
    template<typename F = E>  
    static Result error(F&& err) {
        return Result(std::forward<F>(err));
    }
    
    // Check if result is success or error
    bool is_ok() const { return std::holds_alternative<T>(data_); }
    bool is_error() const { return std::holds_alternative<E>(data_); }
    
    // Access value (throws if error)
    const T& value() const & { return std::get<T>(data_); }
    T& value() & { return std::get<T>(data_); }
    T&& value() && { return std::get<T>(std::move(data_)); }
    
    // Access error (throws if success)
    const E& error() const & { return std::get<E>(data_); }
    E& error() & { return std::get<E>(data_); }
    E&& error() && { return std::get<E>(std::move(data_)); }
    
    // Safe access with defaults
    template<typename U>
    T value_or(U&& default_value) const & {
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
class Result<void, E> {
public:
    // Constructors
    Result() : data_(std::monostate{}) {}
    
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
    const E& error() const & { return std::get<E>(data_); }
    E& error() & { return std::get<E>(data_); }
    E&& error() && { return std::get<E>(std::move(data_)); }
    
    // Convenience operators
    explicit operator bool() const { return is_ok(); }
    
private:
    std::variant<std::monostate, E> data_;
};