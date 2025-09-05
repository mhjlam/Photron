#pragma once

#include <limits>
#include <span>

#include "concepts.hpp"
#include "math.hpp"

// Forward declaration
class Ray;

/**
 * @brief Modern C++20 Axis-Aligned Bounding Box (AABB/Cuboid) implementation
 * @details Unified class combining AABB functionality for BVH operations with
 *          advanced geometric operations. Optimized with C++20 concepts, constexpr support,
 *          and template-based generic interfaces for point-like types.
 *          Replaces both AABB struct and original Cuboid class.
 */
class Cuboid
{
private:
	glm::dvec3 min_point_;
	glm::dvec3 max_point_;

public:
	// Modern C++20 constructors with concepts
	Cuboid() noexcept :
		min_point_(std::numeric_limits<double>::max()), max_point_(-std::numeric_limits<double>::max()) {}

	Cuboid(double x1, double y1, double z1, double x2, double y2, double z2) noexcept;

	template<Point3D P>
	Cuboid(const P& min_pt, const P& max_pt) noexcept :
		min_point_(min_pt.x, min_pt.y, min_pt.z), max_point_(max_pt.x, max_pt.y, max_pt.z) {}

	Cuboid(const glm::dvec3& min_pt, const glm::dvec3& max_pt) noexcept;

	// C++20 constexpr getters
	[[nodiscard]] const glm::dvec3& min_point() const noexcept { return min_point_; }
	[[nodiscard]] const glm::dvec3& max_point() const noexcept { return max_point_; }

	// Modern C++20 setters with concepts
	void set_bounds(const glm::dvec3& min_pt, const glm::dvec3& max_pt) noexcept;

	template<Point3D P>
	void set_bounds(const P& min_pt, const P& max_pt) noexcept {
		min_point_ = glm::dvec3(min_pt.x, min_pt.y, min_pt.z);
		max_point_ = glm::dvec3(max_pt.x, max_pt.y, max_pt.z);
	}

	// BVH-style expansion operations (from AABB)
	void expand(const Cuboid& other) noexcept {
		min_point_ = glm::min(min_point_, other.min_point_);
		max_point_ = glm::max(max_point_, other.max_point_);
	}

	template<Point3D P>
	void expand(const P& point) noexcept {
		min_point_ = glm::min(min_point_, static_cast<glm::dvec3>(point));
		max_point_ = glm::max(max_point_, static_cast<glm::dvec3>(point));
	}

	// C++20 utility methods with [[nodiscard]]
	[[nodiscard]] glm::dvec3 center() const noexcept;
	[[nodiscard]] glm::dvec3 size() const noexcept;
	[[nodiscard]] double volume() const noexcept;

	// Surface area calculation (for BVH SAH heuristic)
	[[nodiscard]] constexpr double surface_area() const noexcept {
		const auto extent = max_point_ - min_point_;
		return 2.0 * (extent.x * extent.y + extent.y * extent.z + extent.z * extent.x);
	}

	// Generic point containment with C++20 concepts
	[[nodiscard]] bool contains(const glm::dvec3& point) const noexcept;

	template<Point3D P>
	[[nodiscard]] bool contains(const P& point) const noexcept {
		return point.x >= min_point_.x && point.x <= max_point_.x && point.y >= min_point_.y && point.y <= max_point_.y
			   && point.z >= min_point_.z && point.z <= max_point_.z;
	}

	// Batch operations with C++20 ranges
	template<std::ranges::range Range>
	[[nodiscard]] bool contains_all(const Range& points) const noexcept
		requires Point3D<std::ranges::range_value_t<Range>>
	{
		return std::ranges::all_of(points, [this](const auto& point) { return contains(point); });
	}

	template<std::ranges::range Range>
	[[nodiscard]] bool contains_any(const Range& points) const noexcept
		requires Point3D<std::ranges::range_value_t<Range>>
	{
		return std::ranges::any_of(points, [this](const auto& point) { return contains(point); });
	}

	// Advanced geometric operations
	[[nodiscard]] Cuboid intersection(const Cuboid& other) const noexcept;
	[[nodiscard]] Cuboid union_bounds(const Cuboid& other) const noexcept;
	[[nodiscard]] bool intersects(const Cuboid& other) const noexcept;

	// C++20 span-based corner access
	[[nodiscard]] std::array<glm::dvec3, 8> corners() const noexcept;

	// Ray intersection methods (enhanced from both AABB and original Cuboid)
	[[nodiscard]] bool intersect_ray(const Ray& ray, double& t_min, double& t_max) const;
	[[nodiscard]] double intersect_ray_internal(const Ray& ray, glm::dvec3& intersection, glm::dvec3& normal) const;

	// C++20 comparison operators
	[[nodiscard]] bool operator==(const Cuboid& other) const noexcept = default;
	[[nodiscard]] auto operator<=>(const Cuboid& other) const noexcept = default;
};
