/**
 * @file cuboid.cpp
 * @brief Axis-aligned bounding box geometry and ray intersection algorithms
 *
 * Implements high-performance AABB (Axis-Aligned Bounding Box) operations for:
 * - Ray-box intersection using the slab method
 * - Boolean operations on bounding volumes
 * - Geometric queries and spatial relationships
 * - BVH-compatible intersection tests
 *
 * All algorithms are optimized for Monte Carlo simulations and 3D rendering
 * with consistent epsilon handling for numerical stability.
 */

#include "cuboid.hpp"

#include <algorithm>
#include <array>
#include <iostream>
#include <limits>
#include <ranges>

#include "math/ray.hpp"

/**
 * @brief Construct cuboid from individual coordinates
 * @param x1,y1,z1 Minimum corner coordinates
 * @param x2,y2,z2 Maximum corner coordinates
 *
 * Creates an axis-aligned bounding box from explicit coordinate values.
 * No validation is performed - caller must ensure min <= max for each axis.
 */
Cuboid::Cuboid(double x1, double y1, double z1, double x2, double y2, double z2) noexcept :
	min_point_(x1, y1, z1), max_point_(x2, y2, z2) {
}

/**
 * @brief Construct cuboid from corner points
 * @param min_pt Minimum corner point
 * @param max_pt Maximum corner point
 *
 * Creates AABB from GLM vector points for convenient geometric operations.
 */
Cuboid::Cuboid(const glm::dvec3& min_pt, const glm::dvec3& max_pt) noexcept : min_point_(min_pt), max_point_(max_pt) {
}

/**
 * @brief Update cuboid bounds
 * @param min_pt New minimum corner
 * @param max_pt New maximum corner
 *
 * Allows dynamic resizing of bounding boxes during spatial operations.
 */
void Cuboid::set_bounds(const glm::dvec3& min_pt, const glm::dvec3& max_pt) noexcept {
	min_point_ = min_pt;
	max_point_ = max_pt;
}

/**
 * @brief Calculate geometric center of cuboid
 * @return Center point as 3D vector
 *
 * Computes centroid for spatial partitioning and BVH construction.
 */
glm::dvec3 Cuboid::center() const noexcept {
	return (min_point_ + max_point_) * 0.5;
}

/**
 * @brief Calculate cuboid dimensions
 * @return Size vector with width, height, depth
 *
 * Provides extent information for spatial analysis and partitioning.
 */
glm::dvec3 Cuboid::size() const noexcept {
	return max_point_ - min_point_;
}

/**
 * @brief Calculate cuboid volume
 * @return Volume as scalar value
 *
 * Used for spatial density calculations and BVH cost metrics.
 */
double Cuboid::volume() const noexcept {
	const auto sz = max_point_ - min_point_;
	return sz.x * sz.y * sz.z;
}

/**
 * @brief Test if point lies within cuboid
 * @param point 3D point to test
 * @return True if point is inside or on boundary
 *
 * Performs inclusive containment test using standard AABB algorithm.
 */
bool Cuboid::contains(const glm::dvec3& point) const noexcept {
	return point.x >= min_point_.x && point.x <= max_point_.x && point.y >= min_point_.y && point.y <= max_point_.y
		   && point.z >= min_point_.z && point.z <= max_point_.z;
}

/**
 * @brief Compute intersection of two cuboids
 * @param other Second cuboid for intersection
 * @return New cuboid representing intersection volume
 *
 * Uses vectorized min/max operations for efficient boolean geometry.
 * Returns degenerate cuboid if no intersection exists.
 */
Cuboid Cuboid::intersection(const Cuboid& other) const noexcept {
	// Calculate intersection bounds using component-wise operations
	const glm::dvec3 new_min = glm::max(min_point_, other.min_point_);
	const glm::dvec3 new_max = glm::min(max_point_, other.max_point_);

	// Check for valid intersection (min <= max for all axes)
	if (new_min.x > new_max.x || new_min.y > new_max.y || new_min.z > new_max.z) {
		return Cuboid {0, 0, 0, 0, 0, 0}; // Return degenerate cuboid
	}

	return Cuboid {new_min, new_max};
}

/**
 * @brief Compute bounding box union of two cuboids
 * @param other Second cuboid for union
 * @return New cuboid encompassing both input cuboids
 *
 * Essential for BVH construction and spatial data structure updates.
 */
Cuboid Cuboid::union_bounds(const Cuboid& other) const noexcept {
	return Cuboid {glm::min(min_point_, other.min_point_), glm::max(max_point_, other.max_point_)};
}

/**
 * @brief Test if two cuboids intersect
 * @param other Second cuboid for intersection test
 * @return True if cuboids have non-zero overlap
 *
 * Fast separating axis test for collision detection and spatial queries.
 */
bool Cuboid::intersects(const Cuboid& other) const noexcept {
	return !(max_point_.x < other.min_point_.x || min_point_.x > other.max_point_.x || max_point_.y < other.min_point_.y
			 || min_point_.y > other.max_point_.y || max_point_.z < other.min_point_.z
			 || min_point_.z > other.max_point_.z);
}

/**
 * @brief Generate all 8 corner points of cuboid
 * @return Array containing all corner vertices
 *
 * Provides complete vertex set for mesh generation and geometric analysis.
 * Uses C++20 aggregate initialization for efficient construction.
 */
std::array<glm::dvec3, 8> Cuboid::corners() const noexcept {
	return {{{min_point_.x, min_point_.y, min_point_.z},
			 {max_point_.x, min_point_.y, min_point_.z},
			 {min_point_.x, max_point_.y, min_point_.z},
			 {max_point_.x, max_point_.y, min_point_.z},
			 {min_point_.x, min_point_.y, max_point_.z},
			 {max_point_.x, min_point_.y, max_point_.z},
			 {min_point_.x, max_point_.y, max_point_.z},
			 {max_point_.x, max_point_.y, max_point_.z}}};
}

/**
 * @brief Ray-AABB intersection with internal surface handling
 * @param ray Input ray for intersection test
 * @param intersection Output intersection point
 * @param normal Output surface normal at intersection
 * @return Distance to intersection, or max double if no intersection
 *
 * Implements the slab method for ray-box intersection with special handling
 * for rays originating inside the volume. Critical for Monte Carlo photon
 * transport where particles traverse voxel boundaries from inside.
 */
double Cuboid::intersect_ray_internal(const Ray& ray, glm::dvec3& intersection, glm::dvec3& normal) const {
	// Extract ray parameters for slab intersection tests
	const glm::dvec3& origin = ray.origin();
	const glm::dvec3& direction = ray.direction();

	// Validate ray direction to avoid degenerate cases
	if (glm::length(direction) < 1e-12) {
		return std::numeric_limits<double>::max(); // Invalid direction
	}

	// Validate cuboid bounds to ensure proper geometry
	if (min_point_.x > max_point_.x || min_point_.y > max_point_.y || min_point_.z > max_point_.z) {
		return std::numeric_limits<double>::max(); // Invalid bounds
	}

	// Initialize slab intersection parameters
	double t_near = -std::numeric_limits<double>::max();
	double t_far = std::numeric_limits<double>::max();

	// Test intersection with each axis-aligned slab
	for (int i = 0; i < 3; ++i) {
		double dir_component = (i == 0) ? direction.x : (i == 1) ? direction.y : direction.z;
		double orig_component = (i == 0) ? origin.x : (i == 1) ? origin.y : origin.z;
		double min_component = (i == 0) ? min_point_.x : (i == 1) ? min_point_.y : min_point_.z;
		double max_component = (i == 0) ? max_point_.x : (i == 1) ? max_point_.y : max_point_.z;

		// Handle parallel ray case (direction component near zero)
		if (std::abs(dir_component) < 1e-9) {
			// Ray parallel to slab - check if origin is within bounds
			if (orig_component < min_component - 1e-9 || orig_component > max_component + 1e-9) {
				return std::numeric_limits<double>::max(); // No intersection possible
			}
		}
		else {
			// Calculate parametric distances to slab faces
			double t1 = (min_component - orig_component) / dir_component;
			double t2 = (max_component - orig_component) / dir_component;

			// Ensure t1 represents near intersection, t2 represents far
			if (t1 > t2) {
				std::swap(t1, t2);
			}

			// Set overall intersection interval
			t_near = std::max(t_near, t1);
			t_far = std::min(t_far, t2);

			// Early exit if intervals don't overlap
			if (t_near > t_far) {
				return std::numeric_limits<double>::max(); // No intersection
			}
		}
	}

	// Select appropriate intersection point based on ray origin
	// Internal rays need exit point, external rays need entry point
	double t;

	// Determine ray origin location and select intersection accordingly
	if (std::abs(t_near) < 1e-9) {
		// Ray starts on or very near boundary - use exit point
		t = t_far;
	}
	else if (t_near <= 1e-9) {
		// Ray starts inside volume - use exit point for internal intersection
		t = t_far;
	}
	else {
		// Ray starts outside volume - use entry point
		t = t_near;
	}

	// Validate intersection is forward and within reasonable bounds
	if (t <= 1e-12 || t > 1e6) {
		return std::numeric_limits<double>::max();
	}

	// Compute intersection point using parametric ray equation
	intersection = origin + direction * t;

	// Ensure intersection point has valid floating-point values
	if (!std::isfinite(intersection.x) || !std::isfinite(intersection.y) || !std::isfinite(intersection.z)) {
		return std::numeric_limits<double>::max();
	}

	// Determine surface normal based on intersection location
	glm::dvec3 center = (min_point_ + max_point_) * 0.5;
	glm::dvec3 relative = intersection - center;
	glm::dvec3 size = max_point_ - min_point_;

	// Normalize relative position to [-1,1] range for each axis
	relative.x /= size.x * 0.5;
	relative.y /= size.y * 0.5;
	relative.z /= size.z * 0.5;

	// Select normal based on axis with largest relative component
	if (std::abs(relative.x) >= std::abs(relative.y) && std::abs(relative.x) >= std::abs(relative.z)) {
		normal = (relative.x > 0) ? glm::dvec3(1, 0, 0) : glm::dvec3(-1, 0, 0);
	}
	else if (std::abs(relative.y) >= std::abs(relative.z)) {
		normal = (relative.y > 0) ? glm::dvec3(0, 1, 0) : glm::dvec3(0, -1, 0);
	}
	else {
		normal = (relative.z > 0) ? glm::dvec3(0, 0, 1) : glm::dvec3(0, 0, -1);
	}

	return t;
}

/**
 * @brief BVH-compatible ray-AABB intersection test
 * @param ray Input ray for intersection testing
 * @param t_min Output parameter for near intersection distance
 * @param t_max Output parameter for far intersection distance
 * @return True if ray intersects the AABB
 *
 * Optimized slab method implementation for BVH traversal algorithms.
 * Returns intersection interval [t_min, t_max] for use in acceleration
 * structures and spatial queries.
 */
bool Cuboid::intersect_ray(const Ray& ray, double& t_min, double& t_max) const {
	const glm::dvec3& ray_origin = ray.origin();
	const glm::dvec3& ray_dir = ray.direction(); // Direction need not be normalized for slab test

	// Initialize intersection interval
	t_min = 0.0;
	t_max = std::numeric_limits<double>::max();

	// Perform slab intersection test for each coordinate axis
	for (int axis = 0; axis < 3; ++axis) {
		const double axis_min = min_point_[axis];
		const double axis_max = max_point_[axis];
		const double ray_orig_axis = ray_origin[axis];
		const double ray_dir_axis = ray_dir[axis];

		// Handle parallel ray case (direction component near zero)
		if (std::abs(ray_dir_axis) < MathConstants::GEOMETRIC_EPSILON) {
			// Check if ray origin lies within slab bounds
			if (ray_orig_axis < axis_min || ray_orig_axis > axis_max) {
				return false; // Ray misses slab entirely
			}
			// Ray is within slab bounds, continue to next axis
		}
		else {
			// Calculate parametric intersection distances
			const double inv_ray_dir = 1.0 / ray_dir_axis;
			double t1 = (axis_min - ray_orig_axis) * inv_ray_dir;
			double t2 = (axis_max - ray_orig_axis) * inv_ray_dir;

			// Ensure t1 represents near intersection, t2 represents far
			if (t1 > t2) {
				std::swap(t1, t2);
			}

			// Set intersection interval with current slab
			t_min = std::max(t_min, t1);
			t_max = std::min(t_max, t2);

			// Early termination if intervals become disjoint
			if (t_min > t_max) {
				return false;
			}
		}
	}

	// Ray intersects AABB if final t_max is non-negative
	return t_max >= 0.0;
}
