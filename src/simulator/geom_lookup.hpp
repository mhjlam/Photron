/**
 * @file geom_lookup.hpp
 * @brief Spatial query and geometry management for Monte Carlo simulation
 *
 * The GeomLookup class handles all spatial operations including medium lookup,
 * voxel access, geometric containment tests, and coordinate transformations.
 * It provides efficient spatial queries using DDA acceleration structures.
 */

#pragma once

#include <memory>
#include <vector>

#include <glm/glm.hpp>

#include "math/range.hpp"

// Forward declarations
class Medium;
class Voxel;
class DDA;
class Photon;
class Cuboid;

/**
 * @class GeomLookup
 * @brief Manages spatial queries and geometry operations for photon transport simulation
 *
 * The GeomLookup provides efficient spatial query operations for the Monte Carlo
 * photon transport simulation. It handles medium lookup, voxel access, and geometric
 * containment tests using DDA (Digital Differential Analyzer) acceleration structures.
 *
 * Key responsibilities:
 * - Medium identification at 3D positions
 * - Voxel access and coordinate conversion
 * - Geometric containment tests
 * - Path tracking assistance for DDA traversal
 *
 * The class maintains references to simulation data but does not own it, ensuring
 * proper lifetime management and allowing for efficient reinitialization.
 */
class GeomLookup
{
public:
	/**
	 * @brief Construct GeomLookup with references to simulation data
	 *
	 * @param mediums Reference to vector of mediums (not owned)
	 * @param ddas Reference to vector of DDA instances (not owned)
	 */
	explicit GeomLookup(std::vector<Medium>& mediums, const std::vector<std::unique_ptr<DDA>>& ddas);

	/**
	 * @brief Find the medium containing a given 3D position
	 *
	 * Performs spatial query to determine which medium (if any) contains
	 * the specified point using geometry intersection tests.
	 *
	 * @param position 3D world coordinate to test
	 * @return Pointer to containing medium, or nullptr if outside all media
	 */
	Medium* find_medium_at(const glm::dvec3& position) const;

	/**
	 * @brief Find medium at position using DDA-optimized spatial queries
	 *
	 * Uses DDA acceleration structures for fast medium lookup at
	 * specified 3D coordinates with improved performance and robustness.
	 *
	 * @param position 3D world coordinate for medium query
	 * @return Pointer to containing medium, or nullptr if outside geometry
	 */
	Medium* find_medium_at_with_dda(const glm::dvec3& position) const;

	/**
	 * @brief Check if position is inside any simulation medium
	 *
	 * Fast test to determine if a point is within the simulation domain
	 * without identifying the specific medium.
	 *
	 * @param position 3D world coordinate to test
	 * @return true if inside any medium, false if in ambient space
	 */
	bool is_inside_any_medium(const glm::dvec3& position) const;

	/**
	 * @brief Check if 3D position is inside simulation geometry
	 *
	 * Tests whether the given position is within any medium boundary
	 * using geometric intersection tests.
	 *
	 * @param position 3D world position to test
	 * @return bool True if position is within any medium boundary
	 */
	bool is_point_inside_geometry(const glm::dvec3& position) const;

	/**
	 * @brief Find voxel at specified 3D position across all media
	 *
	 * Delegates to appropriate medium to find voxel containing the
	 * given position, handling multi-medium scenarios efficiently.
	 *
	 * @param position 3D world coordinate for voxel lookup
	 * @return Pointer to containing voxel, or nullptr if outside voxelized regions
	 */
	Voxel* voxel_at(const glm::dvec3& position) const;

	/**
	 * @brief Access voxel by discrete grid coordinates
	 *
	 * Provides direct access to voxel data using integer grid coordinates.
	 * Searches across all mediums for the specified voxel coordinates.
	 *
	 * @param x Grid X coordinate
	 * @param y Grid Y coordinate
	 * @param z Grid Z coordinate
	 * @return Voxel* Pointer to voxel at coordinates, nullptr if invalid
	 */
	Voxel* voxel_grid(uint32_t x, uint32_t y, uint32_t z) const;

	/**
	 * @brief Compute 3D bounding box corners for given voxel
	 *
	 * Calculates world-space corner coordinates of voxel for geometric
	 * operations and visualization. Delegates to owning medium.
	 *
	 * @param voxel Pointer to voxel for corner calculation
	 * @return Cuboid containing min and max corner coordinates
	 */
	Cuboid voxel_corners(Voxel* voxel) const;

	/**
	 * @brief Find last surface voxel using DDA traversal for accurate exit recording
	 *
	 * Traces backward from exit point to find the last voxel containing
	 * geometry for proper energy classification and recording.
	 *
	 * @param photon Photon that has exited the geometry
	 * @param exit_direction Direction of photon exit
	 * @return Pointer to last surface voxel, or nullptr if not found
	 */
	Voxel* find_last_surface_voxel_with_dda(const Photon& photon, const glm::dvec3& exit_direction);

private:
	std::vector<Medium>* mediums_;                  ///< Reference to mediums vector (not owned)
	const std::vector<std::unique_ptr<DDA>>* ddas_; ///< Reference to DDA instances (not owned)
};
