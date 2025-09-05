#include "cuboid.hpp"

#include <algorithm>
#include <array>
#include <iostream>
#include <limits>
#include <ranges>

#include "ray.hpp"

// C++20 constructors
Cuboid::Cuboid(double x1, double y1, double z1, double x2, double y2, double z2) noexcept :
	min_point_(x1, y1, z1), max_point_(x2, y2, z2) {
}

Cuboid::Cuboid(const glm::dvec3& min_pt, const glm::dvec3& max_pt) noexcept : min_point_(min_pt), max_point_(max_pt) {
}

// C++20 setter
void Cuboid::set_bounds(const glm::dvec3& min_pt, const glm::dvec3& max_pt) noexcept {
	min_point_ = min_pt;
	max_point_ = max_pt;
}

// Modern utility methods
glm::dvec3 Cuboid::center() const noexcept {
	return (min_point_ + max_point_) * 0.5;
}

glm::dvec3 Cuboid::size() const noexcept {
	return max_point_ - min_point_;
}

double Cuboid::volume() const noexcept {
	const auto sz = max_point_ - min_point_;
	return sz.x * sz.y * sz.z;
}

bool Cuboid::contains(const glm::dvec3& point) const noexcept {
	return point.x >= min_point_.x && point.x <= max_point_.x && point.y >= min_point_.y && point.y <= max_point_.y
		   && point.z >= min_point_.z && point.z <= max_point_.z;
}

// Advanced geometric operations with C++20
Cuboid Cuboid::intersection(const Cuboid& other) const noexcept {
	const glm::dvec3 new_min = glm::max(min_point_, other.min_point_);
	const glm::dvec3 new_max = glm::min(max_point_, other.max_point_);

	// If no intersection, return degenerate cuboid
	if (new_min.x > new_max.x || new_min.y > new_max.y || new_min.z > new_max.z) {
		return Cuboid {0, 0, 0, 0, 0, 0};
	}

	return Cuboid {new_min, new_max};
}

Cuboid Cuboid::union_bounds(const Cuboid& other) const noexcept {
	return Cuboid {glm::min(min_point_, other.min_point_), glm::max(max_point_, other.max_point_)};
}

bool Cuboid::intersects(const Cuboid& other) const noexcept {
	return !(max_point_.x < other.min_point_.x || min_point_.x > other.max_point_.x || max_point_.y < other.min_point_.y
			 || min_point_.y > other.max_point_.y || max_point_.z < other.min_point_.z
			 || min_point_.z > other.max_point_.z);
}

// C++20 array-based corner generation
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

double Cuboid::intersect_ray_internal(const Ray& ray, glm::dvec3& intersection, glm::dvec3& normal) const {
	// Classic ray-AABB intersection using the slab method
	const glm::dvec3& origin = ray.origin();
	const glm::dvec3& direction = ray.direction();

	// Check for degenerate cases
	if (glm::length(direction) < 1e-12) {
		return std::numeric_limits<double>::max(); // Invalid direction
	}

	// Check if the voxel bounds are valid
	if (min_point_.x > max_point_.x || min_point_.y > max_point_.y || min_point_.z > max_point_.z) {
		return std::numeric_limits<double>::max(); // Invalid bounds
	}

	// Calculate t values for each slab
	double t_near = -std::numeric_limits<double>::max();
	double t_far = std::numeric_limits<double>::max();

	// For each axis
	for (int i = 0; i < 3; ++i) {
		double dir_component = (i == 0) ? direction.x : (i == 1) ? direction.y : direction.z;
		double orig_component = (i == 0) ? origin.x : (i == 1) ? origin.y : origin.z;
		double min_component = (i == 0) ? min_point_.x : (i == 1) ? min_point_.y : min_point_.z;
		double max_component = (i == 0) ? max_point_.x : (i == 1) ? max_point_.y : max_point_.z;

		if (std::abs(dir_component) < 1e-9) {
			// Ray is parallel to this slab
			// Add a small tolerance for boundary conditions
			if (orig_component < min_component - 1e-9 || orig_component > max_component + 1e-9) {
				return std::numeric_limits<double>::max(); // No intersection
			}
		}
		else {
			// Calculate intersection distances
			double t1 = (min_component - orig_component) / dir_component;
			double t2 = (max_component - orig_component) / dir_component;

			if (t1 > t2)
				std::swap(t1, t2);

			t_near = std::max(t_near, t1);
			t_far = std::min(t_far, t2);

			if (t_near > t_far) {
				return std::numeric_limits<double>::max(); // No intersection
			}
		}
	}

	// Choose the appropriate t value for internal intersection
	// For ray starting inside, we want the exit point (t_far)
	// For ray starting outside, we want the entry point (t_near)
	double t;

	// Special handling for boundary cases - if we're very close to t_near = 0,
	// it means we're starting on or very near the boundary
	if (std::abs(t_near) < 1e-9) {
		// Ray starts on or very near boundary - always use exit point
		t = t_far;
	}
	else if (t_near <= 1e-9) {
		// Ray starts inside - use exit point
		t = t_far;
	}
	else {
		// Ray starts outside - use entry point
		t = t_near;
	}

	// Must be a forward intersection with reasonable distance
	if (t <= 1e-12 || t > 1e6) {
		return std::numeric_limits<double>::max();
	}

	// Calculate intersection point
	intersection = origin + direction * t;

	// Validate intersection point is within reasonable bounds
	if (!std::isfinite(intersection.x) || !std::isfinite(intersection.y) || !std::isfinite(intersection.z)) {
		return std::numeric_limits<double>::max();
	}

	// Calculate normal based on which face we hit
	glm::dvec3 center = (min_point_ + max_point_) * 0.5;
	glm::dvec3 relative = intersection - center;
	glm::dvec3 size = max_point_ - min_point_;

	// Normalize relative position
	relative.x /= size.x * 0.5;
	relative.y /= size.y * 0.5;
	relative.z /= size.z * 0.5;

	// Find the axis with the largest relative component
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

// BVH-style ray intersection method (from AABB)
bool Cuboid::intersect_ray(const Ray& ray, double& t_min, double& t_max) const {
	const glm::dvec3& ray_origin = ray.origin();
	const glm::dvec3& ray_dir = ray.direction(); // Use direction directly, no need to normalize for slab test

	t_min = 0.0;
	t_max = std::numeric_limits<double>::max();

	// Test intersection with each axis-aligned slab
	for (int axis = 0; axis < 3; ++axis) {
		const double axis_min = min_point_[axis];
		const double axis_max = max_point_[axis];
		const double ray_orig_axis = ray_origin[axis];
		const double ray_dir_axis = ray_dir[axis];

		if (std::abs(ray_dir_axis) < MathConstants::GEOMETRIC_EPSILON) {
			// Ray is parallel to the slab - check if ray origin is within slab bounds
			if (ray_orig_axis < axis_min || ray_orig_axis > axis_max) {
				return false;
			}
			// Ray is within the slab, continue to next axis
		}
		else {
			// Calculate intersection distances with slab faces
			const double inv_ray_dir = 1.0 / ray_dir_axis;
			double t1 = (axis_min - ray_orig_axis) * inv_ray_dir;
			double t2 = (axis_max - ray_orig_axis) * inv_ray_dir;

			// Ensure t1 <= t2
			if (t1 > t2) {
				std::swap(t1, t2);
			}

			// Update intersection interval
			t_min = std::max(t_min, t1);
			t_max = std::min(t_max, t2);

			// Early exit if intervals don't overlap
			if (t_min > t_max) {
				return false;
			}
		}
	}

	// Ray intersects AABB if t_max is non-negative
	return t_max >= 0.0;
}
