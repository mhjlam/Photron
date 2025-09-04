#pragma once

#include <cmath>

#include "glm_types.hpp"

// Simple templated range class
template<int N, typename T = double>
struct Range
{
	// We'll use the appropriate GLM vector type based on dimension
	using VecType = std::conditional_t<
		N == 1, T,
		std::conditional_t<N == 2, glm::vec<2, T>, std::conditional_t<N == 3, glm::vec<3, T>, glm::vec<4, T>>>>;

	VecType min_bounds;
	VecType max_bounds;

	// Default constructor
	Range() : min_bounds(T(0)), max_bounds(T(0)) {}

	// Vector constructor
	Range(const VecType& min_pt, const VecType& max_pt) : min_bounds(min_pt), max_bounds(max_pt) {}

	// Convenience methods
	const VecType& min_bound() const { return min_bounds; }
	const VecType& max_bound() const { return max_bounds; }
	VecType& min_bound() { return min_bounds; }
	VecType& max_bound() { return max_bounds; }

	VecType center() const { return (min_bounds + max_bounds) * T(0.5); }
	VecType size() const { return max_bounds - min_bounds; }

	// Vector-based includes
	bool includes(const VecType& point) const {
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

// Range3 specialization for backwards compatibility
template<>
struct Range<3, double>
{
	glm::dvec3 min_bounds;
	glm::dvec3 max_bounds;

	Range() : min_bounds(0.0), max_bounds(0.0) {}

	Range(double x0, double x1, double y0, double y1, double z0, double z1) :
		min_bounds(x0, y0, z0), max_bounds(x1, y1, z1) {}

	Range(const glm::dvec3& min_pt, const glm::dvec3& max_pt) : min_bounds(min_pt), max_bounds(max_pt) {}

	bool includes(double x, double y, double z) const {
		return (x >= min_bounds.x && y >= min_bounds.y && z >= min_bounds.z && x <= max_bounds.x && y <= max_bounds.y
				&& z <= max_bounds.z);
	}

	bool includes(const glm::dvec3& point) const { return includes(point.x, point.y, point.z); }

	// GLM convenience methods
	const glm::dvec3& min_bound() const { return min_bounds; }
	const glm::dvec3& max_bound() const { return max_bounds; }
	glm::dvec3& min_bound() { return min_bounds; }
	glm::dvec3& max_bound() { return max_bounds; }
	glm::dvec3 center() const { return (min_bounds + max_bounds) * 0.5; }
	glm::dvec3 size() const { return max_bounds - min_bounds; }

	// Computed properties for convenience
	double width() const { return std::fabs(max_bounds.x - min_bounds.x); }
	double height() const { return std::fabs(max_bounds.y - min_bounds.y); }
	double depth() const { return std::fabs(max_bounds.z - min_bounds.z); }
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
