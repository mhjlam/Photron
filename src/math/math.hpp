/**
 * @file math.hpp
 * @brief Modern geometric utility functions for robust and efficient 3D calculations
 *
 * This file provides optimized geometric operations using GLM vectorization,
 * consistent numerical precision handling, and C++20 concepts for type safety.
 */

#pragma once

#include <algorithm>
#include <cmath>
#include <concepts>
#include <numbers>
#include <type_traits>

#define GLM_ENABLE_EXPERIMENTAL
#include <glm/glm.hpp>
#include <glm/gtc/matrix_transform.hpp>
#include <glm/gtc/type_ptr.hpp>
#include <glm/gtx/norm.hpp>

// C++20 concepts for type safety
template<typename T>
concept FloatingPoint = std::floating_point<T>;

// C++20 concept for arithmetic types
template<typename T>
concept Arithmetic = std::is_arithmetic_v<T>;

// Conversion helpers between float and double precision
inline glm::vec3 to_float(const glm::dvec3& v) {
	return glm::vec3(static_cast<float>(v.x), static_cast<float>(v.y), static_cast<float>(v.z));
}

inline glm::vec4 to_float(const glm::dvec4& v) {
	return glm::vec4(
		static_cast<float>(v.x), static_cast<float>(v.y), static_cast<float>(v.z), static_cast<float>(v.w));
}

/**
 * @namespace MathConstants
 * @brief Centralized mathematical and numerical constants for consistent precision
 *
 * Provides all mathematical constants and epsilon values used throughout the codebase.
 * Using these centralized constants instead of local definitions ensures numerical
 * consistency across different modules and prevents precision-related bugs.
 *
 * **Mathematical Constants:**
 * - Standard trigonometric constants (π, 2π, π/2, etc.)
 * - Inverse constants for efficient division avoidance
 * - Degree/radian conversion factors
 *
 * **Precision Constants:**
 * - Geometric epsilon for spatial calculations
 * - Boundary epsilon for photon-surface interactions
 * - Simulation-specific thresholds and tolerances
 *
 * All constants are `constexpr` for compile-time evaluation and optimal performance.
 */
namespace MathConstants
{
/// Mathematical constant π (3.14159...)
constexpr double PI = std::numbers::pi;

/// Mathematical constant 2π for full rotations
constexpr double TWO_PI = 2.0 * PI;

/// Mathematical constant π/2 for right angles
constexpr double HALF_PI = PI / 2.0;

/// Inverse of π for efficient multiplication instead of division
constexpr double INV_PI = 1.0 / PI;

/// Inverse of 2π for efficient angle normalization
constexpr double INV_TWO_PI = 1.0 / TWO_PI;

/// Conversion factor from radians to degrees (180/π)
constexpr double RAD_TO_DEG = 180.0 / PI;

/// Conversion factor from degrees to radians (π/180)
constexpr double DEG_TO_RAD = PI / 180.0;

/// Numerical epsilon for geometric calculations and comparisons
constexpr double GEOMETRIC_EPSILON = 1e-12;

/// Larger epsilon for photon boundary nudging to prevent precision issues
constexpr double BOUNDARY_EPSILON = 2e-9;

/// MCML photon weight termination threshold for simulation efficiency
constexpr double MCML_WEIGHT_THRESHOLD = 1e-4;

/// Standard epsilon for photon position nudging to avoid surface precision issues
constexpr double PHOTON_NUDGE_EPSILON = 1e-6;

/// Epsilon for surface interaction calculations
constexpr double SURFACE_NUDGE_EPSILON = 1e-6;

/// Minimum photon step size threshold for simulation stability
constexpr double STEP_SIZE_THRESHOLD = 1e-10;

/// High precision threshold for normalization and distance calculations
constexpr double PRECISION_THRESHOLD = 1e-12;

/// Energy threshold for voxel rendering and emittance detection (float precision)
constexpr float ENERGY_THRESHOLD = 1e-8f;

/// UV coordinate matching epsilon for renderer texture operations
constexpr float UV_MATCH_EPSILON = 0.001f;

/// Vertex proximity threshold for geometric calculations
constexpr double VERTEX_THRESHOLD = 1e-6;

/// Ray intersection epsilon for avoiding self-intersection artifacts
constexpr double RAY_INTERSECTION_EPSILON = 1e-12;
} // namespace MathConstants

/**
 * @namespace GeometricUtils
 * @brief Robust geometric utility functions with C++20 concepts
 *
 * Provides type-safe geometric operations using modern C++20 concepts for
 * compile-time type checking. All functions are designed to handle floating-point
 * precision issues gracefully and consistently.
 */
namespace GeometricUtils
{
/**
 * @brief Safe floating-point equality comparison with configurable tolerance
 *
 * Compares two floating-point values for approximate equality using an epsilon
 * tolerance to handle floating-point precision limitations. Uses C++20 concepts
 * to ensure type safety at compile time.
 *
 * @tparam T Floating-point type (enforced by FloatingPoint concept)
 * @param a First value to compare
 * @param b Second value to compare
 * @param epsilon Tolerance for equality comparison
 * @return true if values are approximately equal within tolerance
 */
template<FloatingPoint T>
constexpr bool approximately_equal(T a, T b, T epsilon = static_cast<T>(MathConstants::GEOMETRIC_EPSILON)) noexcept {
	return std::abs(a - b) <= epsilon;
}

/**
 * @brief Safe floating-point zero comparison with configurable tolerance
 *
 * Tests if a floating-point value is approximately zero using epsilon tolerance.
 * Essential for robust geometric calculations where exact zero comparisons fail.
 *
 * @tparam T Floating-point type (enforced by FloatingPoint concept)
 * @param value Value to test against zero
 * @param epsilon Tolerance for zero comparison
 * @return true if value is approximately zero within tolerance
 */
template<FloatingPoint T>
constexpr bool approximately_zero(T value, T epsilon = static_cast<T>(MathConstants::GEOMETRIC_EPSILON)) noexcept {
	return std::abs(value) <= epsilon;
}

/**
 * @brief Robust vector normalization with zero-length protection
 */
inline glm::dvec3 safe_normalize(const glm::dvec3& v) {
	const double length_sq = glm::dot(v, v);
	if (length_sq < MathConstants::GEOMETRIC_EPSILON * MathConstants::GEOMETRIC_EPSILON) {
		return glm::dvec3(1.0, 0.0, 0.0); // Default direction for zero-length vectors
	}
	return v / std::sqrt(length_sq);
}

/**
 * @brief Fast squared distance calculation (avoids expensive sqrt)
 */
constexpr inline double distance_squared(const glm::dvec3& a, const glm::dvec3& b) {
	const glm::dvec3 diff = b - a;
	return glm::dot(diff, diff);
}

/**
 * @brief Robust barycentric coordinate clamping for triangle intersections
 */
constexpr inline bool valid_barycentric_coordinates(double u, double v) {
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
	}
	else {
		u = glm::normalize(glm::cross(w, glm::dvec3(1, 0, 0)));
	}
	v = glm::cross(w, u);
}

/**
 * @brief Clamp value to range with epsilon tolerance for boundaries
 */
template<typename T>
inline T clamp_with_epsilon(T value, T min, T max, T epsilon = static_cast<T>(MathConstants::GEOMETRIC_EPSILON)) {
	if (value < min + epsilon) {
		return min;
	}
	if (value > max - epsilon) {
		return max;
	}
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

} // namespace GeometricUtils
