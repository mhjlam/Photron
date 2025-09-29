#pragma once

#include <cstdint>
#include <optional>
#include <span>
#include <vector>

// Include types needed for member variables and interface
#include "math/bvh.hpp"    // Needed for BVH member
#include "math/range.hpp"  // Needed for Range3 member
#include "math/triangle.hpp"  // Needed for Triangle vector
#include "math/concepts.hpp"  // Needed for template concepts

// Forward declarations for types only used in method signatures
class Ray;

class Layer
{
public:
	uint8_t id {0};
	uint8_t tissue_id {0};
	std::vector<Triangle> mesh;

	// Default constructor uses default member initialization
	Layer() = default;

	// Delete copy constructor and copy assignment (non-copyable due to BVH)
	Layer(const Layer&) = delete;
	Layer& operator=(const Layer&) = delete;

	// Enable move constructor and move assignment
	Layer(Layer&&) = default;
	Layer& operator=(Layer&&) = default;

	// Update the mesh geometry based on the mesh triangles
	void update_geometry() noexcept;

	// Add a triangle to the mesh
	void add_triangle(const Triangle& triangle);

	// Set all triangles at once
	template<TriangleContainer Container>
	void set_triangles(const Container& triangles)
		requires std::ranges::sized_range<Container>
	{
		mesh.assign(std::ranges::begin(triangles), std::ranges::end(triangles));
		calculate_bounds();
		bvh_.build(mesh);
	}

	// Get the triangles in the mesh (zero-copy access)
	[[nodiscard]] std::span<const Triangle> triangles() const noexcept { return std::span<const Triangle>(mesh); }

	// Get the bounding box of the mesh
	[[nodiscard]] const Range3& bounds() const noexcept { return bounds_; }

	// Fast check if a point might be inside (bounding box test)
	[[nodiscard]] bool might_contain(const glm::dvec3& point) const noexcept {
		return has_bounds_ && bounds_.includes(point);
	}

	// Check if a point is inside the mesh using BVH-accelerated ray casting
	[[nodiscard]] bool contains_point(const glm::dvec3& point) const;

	// Find the first intersection of a ray with the mesh using BVH
	// Returns the distance to intersection and the triangle hit
	[[nodiscard]] std::pair<double, Triangle> find_first_intersection(const Ray& ray) const;

	// Find all intersections of a ray with the mesh (uses BVH for initial filtering)
	[[nodiscard]] std::vector<std::pair<double, Triangle>> find_all_intersections(const Ray& ray) const;

	// Check if the mesh has any triangles
	[[nodiscard]] constexpr bool is_empty() const noexcept { return mesh.empty(); }

	// Check if the mesh bounds have been calculated
	[[nodiscard]] constexpr bool has_bounds() const noexcept { return has_bounds_; }

	// Validate and fix normal orientations to ensure they point outward
	void validate_and_fix_normals();

	bool operator==(const Layer& other) const noexcept { return other.id == id; }

private:
	// Geometry data
	Range3 bounds_;
	bool has_bounds_ {false};
	BVH bvh_; // Spatial acceleration structure

	// Calculate the bounding box of the mesh
	void calculate_bounds();

	// Extend the current bounds to include a triangle
	void extend_bounds(const Triangle& triangle);
};
