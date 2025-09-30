/**
 * @file range.hpp
 * @brief Templated N-dimensional range/bounding box utilities
 *
 * Provides efficient N-dimensional bounding box operations with GLM integration
 * for 1D, 2D, 3D, and 4D coordinate spaces commonly used in geometric computations.
 */

#pragma once

#include "math.hpp"

/**
 * @struct Range
 * @brief N-dimensional templated range/bounding box with GLM vector integration
 *
 * A flexible template class for representing axis-aligned bounding boxes
 * in 1D through 4D coordinate spaces. Uses appropriate GLM vector types
 * for efficient vectorized operations and maintains consistency with
 * the rest of the rendering and simulation pipeline.
 *
 * @tparam N Dimensionality (1, 2, 3, or 4)
 * @tparam T Arithmetic type (typically float or double)
 *
 * Key features:
 * - Automatic GLM vector type selection based on dimensionality
 * - Efficient bounds checking with vectorized operations
 * - Consistent interface across all dimensions
 * - C++20 concepts for type safety
 *
 * Usage examples:
 * @code
 * Range3 bbox(glm::dvec3(0,0,0), glm::dvec3(10,10,10));
 * if (bbox.includes(glm::dvec3(5,5,5))) { ... }
 * auto center = bbox.center();
 * @endcode
 */
template<int N, Arithmetic T = double>
struct Range
{
	/// GLM vector type automatically selected based on dimensionality N
	using VecType = std::conditional_t<
		N == 1,
		T,
		std::conditional_t<N == 2, glm::vec<2, T>, std::conditional_t<N == 3, glm::vec<3, T>, glm::vec<4, T>>>>;

	VecType min_bounds {T(0)}; ///< Minimum bounds in all dimensions
	VecType max_bounds {T(0)}; ///< Maximum bounds in all dimensions

	/**
	 * @brief Default constructor initializes bounds to zero
	 */
	Range() = default;

	/**
	 * @brief Construct range from minimum and maximum points
	 *
	 * @param min_pt Minimum bounds in all dimensions
	 * @param max_pt Maximum bounds in all dimensions
	 */
	Range(const VecType& min_pt, const VecType& max_pt) noexcept : min_bounds(min_pt), max_bounds(max_pt) {}

	// Accessor methods

	/**
	 * @brief Get minimum bounds (const reference)
	 * @return const VecType& Minimum bounds vector
	 */
	const VecType& min_bound() const noexcept { return min_bounds; }

	/**
	 * @brief Get maximum bounds (const reference)
	 * @return const VecType& Maximum bounds vector
	 */
	const VecType& max_bound() const noexcept { return max_bounds; }

	/**
	 * @brief Get minimum bounds (mutable reference)
	 * @return VecType& Minimum bounds vector for modification
	 */
	VecType& min_bound() noexcept { return min_bounds; }

	/**
	 * @brief Get maximum bounds (mutable reference)
	 * @return VecType& Maximum bounds vector for modification
	 */
	VecType& max_bound() noexcept { return max_bounds; }

	/**
	 * @brief Calculate geometric center of the range
	 * @return VecType Center point coordinates
	 */
	VecType center() const noexcept { return (min_bounds + max_bounds) * T(0.5); }

	/**
	 * @brief Calculate size/extent in each dimension
	 * @return VecType Size vector (max - min in each dimension)
	 */
	VecType size() const noexcept { return max_bounds - min_bounds; }

	/**
	 * @brief Test if a point lies within the range (inclusive bounds)
	 *
	 * Efficiently tests point inclusion using appropriate comparison
	 * operations for the given dimensionality.
	 *
	 * @param point Point coordinates to test
	 * @return true if point is within bounds (inclusive), false otherwise
	 */
	bool includes(const VecType& point) const noexcept {
		if constexpr (N == 1) {
			return point >= min_bounds && point <= max_bounds;
		}
		else {
			for (int i = 0; i < N; ++i) {
				if (point[i] < min_bounds[i] || point[i] > max_bounds[i])
					return false;
			}
			return true;
		}
	}
};

// Commonly used type aliases

/**
 * @typedef Range1
 * @brief 1D range for scalar intervals (e.g., time ranges, parameter bounds)
 */
using Range1 = Range<1, double>;

/**
 * @typedef Range2
 * @brief 2D range for rectangular regions (e.g., texture coordinates, screen areas)
 */
using Range2 = Range<2, double>;

/**
 * @typedef Range3
 * @brief 3D range for axis-aligned bounding boxes (primary use case in 3D graphics)
 */
using Range3 = Range<3, double>;
using Range4 = Range<4, double>;

using Range1f = Range<1, float>;
using Range2f = Range<2, float>;
using Range3f = Range<3, float>;
using Range4f = Range<4, float>;
