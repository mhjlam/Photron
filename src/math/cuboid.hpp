/**
 * @file cuboid.hpp
 * @brief Axis-aligned bounding box operations and ray intersection
 *
 * Provides efficient AABB (Axis-Aligned Bounding Box) operations including
 * ray intersection tests, volume calculations, and spatial queries.
 */

#pragma once

#include <limits>
#include <span>

#include "math/concepts.hpp"
#include "math/math.hpp"

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
	/**
	 * @brief Default constructor creating empty/invalid bounding box
	 *
	 * Initializes an empty cuboid with inverted bounds that can be expanded
	 * to contain geometry through subsequent expand() operations.
	 */
	Cuboid() noexcept :
		min_point_(std::numeric_limits<double>::max()), max_point_(-std::numeric_limits<double>::max()) {}

	/**
	 * @brief Construct cuboid from coordinate components
	 * @param x1 Minimum x coordinate
	 * @param y1 Minimum y coordinate
	 * @param z1 Minimum z coordinate
	 * @param x2 Maximum x coordinate
	 * @param y2 Maximum y coordinate
	 * @param z2 Maximum z coordinate
	 */
	Cuboid(double x1, double y1, double z1, double x2, double y2, double z2) noexcept;

	/**
	 * @brief Construct cuboid from generic 3D points
	 * @tparam P Point type satisfying Point3D concept
	 * @param min_pt Minimum corner point
	 * @param max_pt Maximum corner point
	 *
	 * Template constructor accepting any point-like type with x, y, z members.
	 */
	template<Point3D P>
	Cuboid(const P& min_pt, const P& max_pt) noexcept :
		min_point_(min_pt.x, min_pt.y, min_pt.z), max_point_(max_pt.x, max_pt.y, max_pt.z) {}

	/**
	 * @brief Construct cuboid from GLM vectors
	 * @param min_pt Minimum corner as GLM vector
	 * @param max_pt Maximum corner as GLM vector
	 */
	Cuboid(const glm::dvec3& min_pt, const glm::dvec3& max_pt) noexcept;

	// C++20 constexpr getters
	/**
	 * @brief Get minimum corner point
	 * @return Reference to minimum corner coordinates
	 */
	[[nodiscard]] const glm::dvec3& min_point() const noexcept { return min_point_; }

	/**
	 * @brief Get maximum corner point
	 * @return Reference to maximum corner coordinates
	 */
	[[nodiscard]] const glm::dvec3& max_point() const noexcept { return max_point_; }

	// Modern C++20 setters with concepts
	/**
	 * @brief Set cuboid bounds from GLM vectors
	 * @param min_pt New minimum corner coordinates
	 * @param max_pt New maximum corner coordinates
	 */
	void set_bounds(const glm::dvec3& min_pt, const glm::dvec3& max_pt) noexcept;

	/**
	 * @brief Set cuboid bounds from generic 3D points
	 * @tparam P Point type satisfying Point3D concept
	 * @param min_pt New minimum corner point
	 * @param max_pt New maximum corner point
	 */
	template<Point3D P>
	void set_bounds(const P& min_pt, const P& max_pt) noexcept {
		min_point_ = glm::dvec3(min_pt.x, min_pt.y, min_pt.z);
		max_point_ = glm::dvec3(max_pt.x, max_pt.y, max_pt.z);
	}

	// BVH-style expansion operations (from AABB)
	/**
	 * @brief Expand bounds to include another cuboid
	 * @param other Cuboid to include in expansion
	 *
	 * Modifies this cuboid to be the union of both bounding boxes,
	 * essential for BVH construction and spatial partitioning.
	 */
	void expand(const Cuboid& other) noexcept {
		min_point_ = glm::min(min_point_, other.min_point_);
		max_point_ = glm::max(max_point_, other.max_point_);
	}

	/**
	 * @brief Expand bounds to include a point
	 * @tparam P Point type satisfying Point3D concept
	 * @param point Point to include in expansion
	 *
	 * Grows the bounding box as needed to contain the specified point.
	 */
	template<Point3D P>
	void expand(const P& point) noexcept {
		min_point_ = glm::min(min_point_, static_cast<glm::dvec3>(point));
		max_point_ = glm::max(max_point_, static_cast<glm::dvec3>(point));
	}

	// C++20 utility methods with [[nodiscard]]
	/**
	 * @brief Calculate geometric center of cuboid
	 * @return Center point coordinates
	 */
	[[nodiscard]] glm::dvec3 center() const noexcept;

	/**
	 * @brief Calculate size dimensions of cuboid
	 * @return Extent vector (max - min for each axis)
	 */
	[[nodiscard]] glm::dvec3 size() const noexcept;

	/**
	 * @brief Calculate volume of cuboid
	 * @return Volume in cubic units
	 */
	[[nodiscard]] double volume() const noexcept;

	/**
	 * @brief Calculate surface area for SAH heuristic
	 * @return Total surface area of all six faces
	 *
	 * Used in BVH construction's Surface Area Heuristic (SAH)
	 * for optimal spatial partitioning decisions.
	 */
	[[nodiscard]] constexpr double surface_area() const noexcept {
		const auto extent = max_point_ - min_point_;
		return 2.0 * (extent.x * extent.y + extent.y * extent.z + extent.z * extent.x);
	}

	// Generic point containment with C++20 concepts
	/**
	 * @brief Test if point is inside cuboid bounds
	 * @param point Point to test for containment
	 * @return True if point is within bounds (inclusive)
	 */
	[[nodiscard]] bool contains(const glm::dvec3& point) const noexcept;

	/**
	 * @brief Test if generic point is inside cuboid bounds
	 * @tparam P Point type satisfying Point3D concept
	 * @param point Point to test for containment
	 * @return True if point is within bounds (inclusive)
	 */
	template<Point3D P>
	[[nodiscard]] bool contains(const P& point) const noexcept {
		return point.x >= min_point_.x && point.x <= max_point_.x && point.y >= min_point_.y && point.y <= max_point_.y
			   && point.z >= min_point_.z && point.z <= max_point_.z;
	}

	// Batch operations with C++20 ranges
	/**
	 * @brief Test if all points in range are contained
	 * @tparam Range C++20 range type containing Point3D elements
	 * @param points Range of points to test
	 * @return True if all points are within bounds
	 *
	 * Efficient batch containment test using C++20 ranges and algorithms.
	 */
	template<std::ranges::range Range>
	[[nodiscard]] bool contains_all(const Range& points) const noexcept
		requires Point3D<std::ranges::range_value_t<Range>>
	{
		return std::ranges::all_of(points, [this](const auto& point) { return contains(point); });
	}

	/**
	 * @brief Test if any point in range is contained
	 * @tparam Range C++20 range type containing Point3D elements
	 * @param points Range of points to test
	 * @return True if at least one point is within bounds
	 *
	 * Efficient batch containment test with early termination.
	 */
	template<std::ranges::range Range>
	[[nodiscard]] bool contains_any(const Range& points) const noexcept
		requires Point3D<std::ranges::range_value_t<Range>>
	{
		return std::ranges::any_of(points, [this](const auto& point) { return contains(point); });
	}

	// Advanced geometric operations
	/**
	 * @brief Compute intersection of two cuboids
	 * @param other Cuboid to intersect with
	 * @return New cuboid representing the overlapping volume
	 *
	 * Returns empty cuboid if no intersection exists.
	 */
	[[nodiscard]] Cuboid intersection(const Cuboid& other) const noexcept;

	/**
	 * @brief Compute union of two cuboids
	 * @param other Cuboid to union with
	 * @return New cuboid encompassing both volumes
	 *
	 * Creates minimal bounding box containing both cuboids.
	 */
	[[nodiscard]] Cuboid union_bounds(const Cuboid& other) const noexcept;

	/**
	 * @brief Test if two cuboids intersect
	 * @param other Cuboid to test intersection with
	 * @return True if cuboids overlap in any dimension
	 */
	[[nodiscard]] bool intersects(const Cuboid& other) const noexcept;

	/**
	 * @brief Get all 8 corner vertices of the cuboid
	 * @return Array of corner positions in consistent ordering
	 *
	 * Returns corners in standard ordering for visualization and
	 * geometric processing applications.
	 */
	[[nodiscard]] std::array<glm::dvec3, 8> corners() const noexcept;

	// Ray intersection methods (enhanced from both AABB and original Cuboid)
	/**
	 * @brief Test ray intersection with cuboid
	 * @param ray Ray to test for intersection
	 * @param t_min Output parameter for near intersection distance
	 * @param t_max Output parameter for far intersection distance
	 * @return True if ray intersects the cuboid
	 *
	 * Efficient AABB ray intersection using slab method.
	 */
	[[nodiscard]] bool intersect_ray(const Ray& ray, double& t_min, double& t_max) const;

	/**
	 * @brief Detailed ray intersection with surface information
	 * @param ray Ray to intersect
	 * @param intersection Output intersection point coordinates
	 * @param normal Output surface normal at intersection
	 * @return Distance along ray to intersection, or negative if no hit
	 *
	 * Provides complete intersection data including surface normal
	 * for lighting calculations and Monte Carlo photon transport.
	 */
	[[nodiscard]] double intersect_ray_internal(const Ray& ray, glm::dvec3& intersection, glm::dvec3& normal) const;

	// C++20 comparison operators
	/**
	 * @brief Equality comparison of cuboids
	 * @param other Cuboid to compare with
	 * @return True if both min and max points are equal
	 */
	[[nodiscard]] bool operator==(const Cuboid& other) const noexcept = default;

	/**
	 * @brief Three-way comparison of cuboids
	 * @param other Cuboid to compare with
	 * @return Comparison result for ordering operations
	 *
	 * Enables automatic generation of <, <=, >, >= operators.
	 */
	[[nodiscard]] auto operator<=>(const Cuboid& other) const noexcept = default;
};
