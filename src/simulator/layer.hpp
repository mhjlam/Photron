/**
 * @file layer.hpp
 * @brief Geometric layer representation for multi-layered material simulation
 *
 * Defines the Layer class which represents individual material layers within
 * a multi-layered medium. Each layer consists of triangular mesh geometry with
 * associated material properties and optimized spatial acceleration structures.
 */

#pragma once

#include <cstdint>
#include <optional>
#include <span>
#include <vector>

#include "math/bvh.hpp"
#include "math/concepts.hpp"
#include "math/range.hpp"
#include "math/triangle.hpp"

// Forward declarations
class Ray;

/**
 * @class Layer
 * @brief Geometric material layer with triangular mesh representation
 *
 * The Layer class represents a single material layer within a multi-layered
 * medium used for Monte Carlo photon transport simulation. Each layer consists of:
 *
 * **Geometric Representation:**
 * - **Triangular Mesh**: Defines layer boundaries and interfaces
 * - **Bounding Volume Hierarchy**: Accelerates ray-mesh intersection queries
 * - **Spatial Bounds**: Axis-aligned bounding box for culling operations
 *
 * **Material Association:**
 * - **Layer ID**: Unique identifier within the medium
 * - **Material ID**: Links to material properties (absorption, scattering, etc.)
 *
 * **Key Features:**
 * - **Move-only Semantics**: Prevents expensive copying of large mesh data
 * - **Automatic BVH Updates**: Spatial acceleration structures rebuilt on geometry changes
 * - **Template-based Mesh Loading**: Efficient bulk triangle insertion
 * - **Intersection Queries**: Fast ray-layer intersection with surface normals
 *
 * **Typical Use Cases:**
 * - Skin layers in dermatological simulations
 * - Tissue interfaces in medical imaging
 * - Material boundaries in industrial applications
 * - Optical component surfaces
 *
 * The layer system supports complex multi-layered geometries where photons
 * must navigate between materials with different optical properties.
 */
class Layer
{
public:
	uint8_t id {0};             ///< Unique layer identifier within the medium
	uint8_t tissue_id {0};      ///< Associated material/tissue type identifier
	std::vector<Triangle> mesh; ///< Triangular mesh defining layer geometry

	/**
	 * @brief Default constructor creates empty layer
	 *
	 * Initializes layer with default IDs and empty mesh.
	 * Call add_triangle() or set_triangles() to define geometry.
	 */
	Layer() = default;

	// Move-only semantics (BVH makes copying expensive)
	Layer(const Layer&) = delete;            ///< Copy constructor disabled
	Layer& operator=(const Layer&) = delete; ///< Copy assignment disabled
	Layer(Layer&&) = default;                ///< Move constructor
	Layer& operator=(Layer&&) = default;     ///< Move assignment

	/**
	 * @brief Update spatial acceleration structures after geometry changes
	 *
	 * Rebuilds bounding volume hierarchy and recalculates spatial bounds
	 * after mesh modifications. Should be called after adding/modifying triangles.
	 */
	void update_geometry() noexcept;

	/**
	 * @brief Add a single triangle to the layer mesh
	 *
	 * Appends triangle to mesh and updates spatial bounds.
	 * For bulk operations, prefer set_triangles() for better performance.
	 *
	 * @param triangle Triangle to add to the mesh
	 */
	void add_triangle(const Triangle& triangle);

	/**
	 * @brief Set complete triangle mesh from container
	 *
	 * Efficiently replaces entire mesh with triangles from any container
	 * that satisfies the TriangleContainer concept. Automatically updates
	 * bounds and rebuilds spatial acceleration structures.
	 *
	 * @tparam Container Type satisfying TriangleContainer concept
	 * @param triangles Container of Triangle objects to set as mesh
	 */
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

	// Validate and correct normal orientations to ensure they point outward
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
