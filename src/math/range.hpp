#pragma once

#include <cmath>
#include <concepts>

#include "glm_types.hpp"

// C++20 concept for arithmetic types
template<typename T>
concept Arithmetic = std::is_arithmetic_v<T>;

// Simple templated range class with C++20 concepts
template<int N, Arithmetic T = double>
struct Range
{
	// We'll use the appropriate GLM vector type based on dimension
	using VecType = std::conditional_t<
		N == 1, T,
		std::conditional_t<N == 2, glm::vec<2, T>, std::conditional_t<N == 3, glm::vec<3, T>, glm::vec<4, T>>>>;

	VecType min_bounds{T(0)};
	VecType max_bounds{T(0)};

	// Default constructor uses default member initialization
	Range() = default;

	// Vector constructor
	Range(const VecType& min_pt, const VecType& max_pt) noexcept : min_bounds(min_pt), max_bounds(max_pt) {}

	// Convenience methods
	const VecType& min_bound() const noexcept { return min_bounds; }
	const VecType& max_bound() const noexcept { return max_bounds; }
	VecType& min_bound() noexcept { return min_bounds; }
	VecType& max_bound() noexcept { return max_bounds; }

	VecType center() const noexcept { return (min_bounds + max_bounds) * T(0.5); }
	VecType size() const noexcept { return max_bounds - min_bounds; }

	// Vector-based includes
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

// Type aliases for common use cases
using Range1 = Range<1, double>;
using Range2 = Range<2, double>;
using Range3 = Range<3, double>;
using Range4 = Range<4, double>;

using Range1f = Range<1, float>;
using Range2f = Range<2, float>;
using Range3f = Range<3, float>;
using Range4f = Range<4, float>;
