/**
 * @file dda.hpp
 * @brief 3D Digital Differential Analyzer for robust voxel traversal
 *
 * Implements the 3D DDA (Digital Differential Analyzer) algorithm for accurate
 * ray-voxel intersection and traversal. This is critical for Monte Carlo photon
 * transport simulation where photons must accurately traverse discrete voxel grids
 * without missing voxels or introducing floating-point precision errors.
 *
 * The DDA algorithm is significantly more robust than traditional geometric
 * ray-box intersection methods, especially for long rays or small voxels where
 * floating-point precision becomes problematic.
 */

#pragma once

#include <algorithm>
#include <cmath>
#include <vector>

#include <glm/glm.hpp>

#include "math/math.hpp"

/**
 * @class VoxelDDA3D
 * @brief 3D Digital Differential Analyzer for precise voxel grid traversal
 *
 * The VoxelDDA3D class implements the 3D extension of Bresenham's line algorithm
 * for discrete voxel grid traversal. This algorithm is essential for Monte Carlo
 * photon transport simulation where photons must traverse voxelized media with
 * perfect accuracy.
 *
 * **Key Features:**
 * - **Precision**: Uses integer arithmetic to avoid floating-point errors
 * - **Completeness**: Guarantees no voxels are missed along ray path
 * - **Efficiency**: O(n) complexity where n is number of traversed voxels
 * - **Robustness**: Handles edge cases and boundary conditions correctly
 *
 * **Algorithm Overview:**
 * 1. Convert ray to discrete grid coordinates
 * 2. Calculate step direction and delta distances
 * 3. Step through voxels one at a time using integer increments
 * 4. Track entry/exit points and normals for each voxel
 *
 * **Applications:**
 * - Photon path tracking in Monte Carlo simulation
 * - Energy deposition calculation in discrete voxel grids
 * - Ray casting for visualization and collision detection
 * - Volume rendering traversal
 *
 * Based on:
 * "A Fast Voxel Traversal Algorithm for Ray Tracing" by Amanatides & Woo (1987)
 */
class DDA
{
public:
	/**
	 * @struct StepResult
	 * @brief Results of a single voxel traversal step
	 *
	 * Contains all information about a voxel that was traversed during
	 * ray marching, including spatial coordinates, entry information,
	 * and distance metrics.
	 */
	struct StepResult
	{
		glm::ivec3 voxel_coords;   ///< Discrete voxel coordinates (ix, iy, iz)
		glm::dvec3 world_position; ///< World space position at voxel entry
		glm::dvec3 entry_normal;   ///< Surface normal of entry face (±1 in one axis, 0 in others)
		double distance_traveled;  ///< Cumulative distance from ray origin to this voxel
		bool valid;                ///< True if this represents a valid voxel intersection

		/**
		 * @brief Default constructor creates invalid step result
		 */
		StepResult() : voxel_coords(0), world_position(0.0), entry_normal(0.0), distance_traveled(0.0), valid(false) {}
	};

	/**
	 * @struct TraversalResult
	 * @brief Complete results of ray traversal through voxel grid
	 *
	 * Aggregates all voxels intersected by a ray, along with summary
	 * information about the traversal path and exit conditions.
	 * Essential for photon transport calculations and energy deposition.
	 */
	struct TraversalResult
	{
		std::vector<StepResult> voxels; ///< Ordered sequence of all traversed voxels
		double total_distance;          ///< Total distance traveled through the grid
		bool hit_boundary;              ///< True if ray exited grid boundaries (not absorbed)
		glm::ivec3 exit_voxel;          ///< Coordinates of last valid voxel before exit
		glm::dvec3 exit_position;       ///< World space position where ray exited grid
		glm::dvec3 exit_normal;         ///< Surface normal of grid face where ray exited
	};

private:
	// Grid configuration parameters
	glm::ivec3 grid_dimensions_; ///< Grid dimensions in voxels (nx, ny, nz)
	glm::dvec3 grid_origin_;     ///< World coordinates of grid corner (minimum bounds)
	double voxel_size_;          ///< Uniform edge length of each cubic voxel
	glm::dvec3 grid_bounds_max_; // Maximum bounds of the grid

	// Current state
	glm::dvec3 ray_origin_;    // Ray starting point
	glm::dvec3 ray_direction_; // Ray direction (normalized)

	// DDA state variables
	glm::ivec3 current_voxel_; // Current voxel coordinates
	glm::dvec3 delta_dist_;    // Distance ray travels for each unit step in each axis
	glm::dvec3 side_dist_;     // Distance from current position to next grid line
	glm::ivec3 step_;          // Direction to step (+1 or -1 for each axis)
	glm::ivec3 side_;          // Which axis had the shortest distance to next grid line
	double last_distance_;     // Track cumulative distance for per-voxel calculations

public:
	/**
	 * @brief Construct DDA with grid parameters
	 */
	DDA(const glm::ivec3& grid_dimensions, const glm::dvec3& grid_origin, double voxel_size);

	/**
	 * @brief Initialize DDA for a specific ray
	 */
	void initialize_ray(const glm::dvec3& origin, const glm::dvec3& direction);

	/**
	 * @brief Perform single DDA step
	 * @return StepResult for the next voxel, or invalid result if ray exits grid
	 */
	StepResult step();

	/**
	 * @brief Traverse entire ray path through voxel grid
	 * @param max_distance Maximum distance to traverse (optional)
	 * @return Complete traversal result
	 */
	TraversalResult traverse(double max_distance = std::numeric_limits<double>::max());

	/**
	 * @brief Find voxel containing a specific world position
	 * @param position World coordinates
	 * @return Voxel coordinates, or (-1,-1,-1) if outside grid
	 */
	glm::ivec3 world_to_voxel(const glm::dvec3& position) const;

	/**
	 * @brief Convert voxel coordinates to world position (center of voxel)
	 * @param voxel_coords Voxel coordinates
	 * @return World coordinates of voxel center
	 */
	glm::dvec3 voxel_to_world(const glm::ivec3& voxel_coords) const;

	/**
	 * @brief Check if voxel coordinates are within grid bounds
	 */
	bool is_valid_voxel(const glm::ivec3& voxel_coords) const;

	/**
	 * @brief Get the world bounds of a specific voxel
	 */
	std::pair<glm::dvec3, glm::dvec3> get_voxel_bounds(const glm::ivec3& voxel_coords) const;

private:
	/**
	 * @brief Calculate entry normal based on which face was crossed
	 */
	glm::dvec3 calculate_entry_normal(int face_axis, int step_direction) const;

	/**
	 * @brief Calculate world position from current DDA state
	 */
	glm::dvec3 calculate_world_position() const;
};
