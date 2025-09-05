#pragma once

#include "glm_types.hpp"
#include "numerical_constants.hpp"
#include <algorithm>
#include <cmath>
#include <concepts>
#include <type_traits>

/**
 * @file geometric_utils.hpp
 * @brief Modern geometric utility functions for robust and efficient 3D calculations
 * 
 * This file provides optimized geometric operations using GLM vectorization,
 * consistent numerical precision handling, and C++20 concepts for type safety.
 */

// C++20 concepts for type safety
template<typename T>
concept FloatingPoint = std::floating_point<T>;

template<typename T>
concept Arithmetic = std::is_arithmetic_v<T>;

namespace GeometricUtils {
    
    /**
     * @brief Safe floating-point equality comparison with epsilon tolerance
     * @note Now uses C++20 concepts for compile-time type checking
     */
    template<FloatingPoint T>
    constexpr bool approximately_equal(T a, T b, T epsilon = static_cast<T>(NumericalConstants::GEOMETRIC_EPSILON)) noexcept {
        return std::abs(a - b) <= epsilon;
    }
    
    /**
     * @brief Safe floating-point zero comparison
     * @note Enhanced with C++20 concepts and constexpr
     */
    template<FloatingPoint T>
    constexpr bool approximately_zero(T value, T epsilon = static_cast<T>(NumericalConstants::GEOMETRIC_EPSILON)) noexcept {
        return std::abs(value) <= epsilon;
    }
    
    /**
     * @brief Robust vector normalization with zero-length protection
     */
    inline glm::dvec3 safe_normalize(const glm::dvec3& v) {
        const double length_sq = glm::dot(v, v);
        if (length_sq < NumericalConstants::GEOMETRIC_EPSILON * NumericalConstants::GEOMETRIC_EPSILON) {
            return glm::dvec3(1.0, 0.0, 0.0); // Default direction for zero-length vectors
        }
        return v / std::sqrt(length_sq);
    }
    
    /**
     * @brief Fast squared distance calculation (avoids expensive sqrt)
     */
    inline double distance_squared(const glm::dvec3& a, const glm::dvec3& b) {
        const glm::dvec3 diff = b - a;
        return glm::dot(diff, diff);
    }
    
    /**
     * @brief Robust barycentric coordinate clamping for triangle intersections
     */
    inline bool valid_barycentric_coordinates(double u, double v) {
        return u >= 0.0 && v >= 0.0 && (u + v) <= 1.0;
    }
    
    /**
     * @brief Modern reflection vector calculation using GLM
     */
    inline glm::dvec3 reflect(const glm::dvec3& incident, const glm::dvec3& normal) {
        return incident - 2.0 * glm::dot(incident, normal) * normal;
    }
    
    /**
     * @brief Robust refraction vector calculation using Snell's law
     */
    inline glm::dvec3 refract(const glm::dvec3& incident, const glm::dvec3& normal, double eta_ratio) {
        const double cos_i = -glm::dot(incident, normal);
        const double sin_t_sq = eta_ratio * eta_ratio * (1.0 - cos_i * cos_i);
        
        if (sin_t_sq >= 1.0) {
            // Total internal reflection
            return glm::dvec3(0.0);
        }
        
        const double cos_t = std::sqrt(1.0 - sin_t_sq);
        return eta_ratio * incident + (eta_ratio * cos_i - cos_t) * normal;
    }
    
    /**
     * @brief Create orthonormal basis from a single direction vector
     */
    inline void create_orthonormal_basis(const glm::dvec3& w, glm::dvec3& u, glm::dvec3& v) {
        // Choose a vector not parallel to w
        if (std::abs(w.z) < 0.9) {
            u = glm::normalize(glm::cross(w, glm::dvec3(0, 0, 1)));
        } else {
            u = glm::normalize(glm::cross(w, glm::dvec3(1, 0, 0)));
        }
        v = glm::cross(w, u);
    }
    
    /**
     * @brief Clamp value to range with epsilon tolerance for boundaries
     */
    template<typename T>
    inline T clamp_with_epsilon(T value, T min_val, T max_val, T epsilon = static_cast<T>(NumericalConstants::GEOMETRIC_EPSILON)) {
        if (value < min_val + epsilon) return min_val;
        if (value > max_val - epsilon) return max_val;
        return value;
    }
    
    /**
     * @brief Calculate triangle area using cross product
     */
    inline double triangle_area(const glm::dvec3& v0, const glm::dvec3& v1, const glm::dvec3& v2) {
        const glm::dvec3 edge1 = v1 - v0;
        const glm::dvec3 edge2 = v2 - v0;
        return 0.5 * glm::length(glm::cross(edge1, edge2));
    }
    
    /**
     * @brief Calculate triangle normal using cross product
     */
    inline glm::dvec3 triangle_normal(const glm::dvec3& v0, const glm::dvec3& v1, const glm::dvec3& v2) {
        const glm::dvec3 edge1 = v1 - v0;
        const glm::dvec3 edge2 = v2 - v0;
        return safe_normalize(glm::cross(edge1, edge2));
    }
}
