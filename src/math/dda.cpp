/**
 * @file dda.cpp
 * @brief Implementation of 3D Digital Differential Analyzer for voxel traversal
 *
 * Implements the classic 3D DDA algorithm for efficient ray-voxel grid traversal.
 * Used extensively in Monte Carlo photon transport for determining which voxels
 * a photon path intersects and the distances traveled in each voxel.
 */

#include "dda.hpp"

#include <cmath>
#include <limits>

/**
 * @brief Construct DDA for 3D voxel grid traversal
 * @param grid_dimensions Number of voxels in each dimension (x, y, z)
 * @param grid_origin World coordinate origin of the voxel grid
 * @param voxel_size Size of each voxel in world units
 *
 * Initializes the DDA algorithm with grid parameters for efficient
 * ray-voxel traversal in Monte Carlo photon transport simulations.
 */
DDA::DDA(const glm::ivec3& grid_dimensions, const glm::dvec3& grid_origin, double voxel_size) :
	grid_dimensions_(grid_dimensions), grid_origin_(grid_origin), voxel_size_(voxel_size) {
	// Calculate maximum grid bounds for boundary detection
	grid_bounds_max_ = grid_origin_ + glm::dvec3(grid_dimensions_) * voxel_size_;
}

/**
 * @brief Initialize ray for voxel traversal
 * @param origin Ray starting point in world coordinates
 * @param direction Ray direction vector (will be normalized)
 *
 * Sets up DDA state for ray traversal including step directions,
 * delta distances, and initial voxel coordinates. Must be called
 * before using step() or traverse() methods.
 */
void DDA::initialize_ray(const glm::dvec3& origin, const glm::dvec3& direction) {
	ray_origin_ = origin;
	ray_direction_ = glm::normalize(direction);
	last_distance_ = 0.0; // Reset distance tracking for new ray

	// Transform ray origin to grid coordinate system
	glm::dvec3 ray_pos = (ray_origin_ - grid_origin_) / voxel_size_;

	// Determine starting voxel coordinates (floor to get integer indices)
	current_voxel_ = glm::ivec3(floor(ray_pos.x), floor(ray_pos.y), floor(ray_pos.z));

	// Calculate delta distances for DDA algorithm
	// Delta distance: how far along ray to cross one voxel boundary per axis
	for (int i = 0; i < 3; ++i) {
		if (std::abs(ray_direction_[i]) < MathConstants::GEOMETRIC_EPSILON) {
			// Ray is parallel to this axis - will never cross boundaries
			delta_dist_[i] = std::numeric_limits<double>::max();
		}
		else {
			// Distance to cross one voxel width along this axis
			delta_dist_[i] = std::abs(1.0 / ray_direction_[i]) * voxel_size_;
		}
	}

	// Calculate initial step directions and side distances
	for (int i = 0; i < 3; ++i) {
		if (ray_direction_[i] < 0) {
			// Moving in negative direction along this axis
			step_[i] = -1;
			side_dist_[i] = (ray_pos[i] - current_voxel_[i]) * delta_dist_[i];
		}
		else {
			// Moving in positive direction along this axis
			step_[i] = 1;
			side_dist_[i] = (current_voxel_[i] + 1.0 - ray_pos[i]) * delta_dist_[i];
		}
	}
}

/**
 * @brief Perform single DDA step to next voxel
 * @return StepResult containing voxel coordinates, world position, and distance
 *
 * Advances the DDA algorithm by one step to the next voxel boundary.
 * Returns information about the current voxel before stepping, including
 * entry normal and distance traveled within the voxel.
 */
DDA::StepResult DDA::step() {
	StepResult result;

	// Validate current voxel is within grid bounds
	if (!is_valid_voxel(current_voxel_)) {
		result.valid = false;
		return result;
	}

	// Record current voxel information before stepping
	result.voxel_coords = current_voxel_;
	result.world_position = calculate_world_position();
	result.valid = true;

	// Core DDA algorithm: find next grid boundary
	// Determine which axis has shortest distance to next grid line
	int min_axis = 0;
	double min_dist = side_dist_[0];

	for (int i = 1; i < 3; ++i) {
		if (side_dist_[i] < min_dist) {
			min_dist = side_dist_[i];
			min_axis = i;
		}
	}

	// Calculate voxel entry properties
	result.entry_normal = calculate_entry_normal(min_axis, -step_[min_axis]);

	// Calculate distance traveled within THIS voxel (not cumulative)
	double previous_distance = last_distance_;
	result.distance_traveled = min_dist - previous_distance;
	last_distance_ = min_dist;

	// Advance to next voxel using DDA step
	side_dist_[min_axis] += delta_dist_[min_axis]; // Set distance to next boundary
	current_voxel_[min_axis] += step_[min_axis];   // Move to adjacent voxel
	side_[min_axis] = 1;                           // Mark which boundary face was crossed

	return result;
}

/**
 * @brief Traverse voxel grid up to maximum distance
 * @param max_distance Maximum distance to travel along ray
 * @return TraversalResult containing all visited voxels and exit information
 *
 * Performs complete DDA traversal collecting all voxels intersected
 * by the ray within the specified maximum distance. Handles grid
 * boundary detection and provides exit position/normal information.
 */
DDA::TraversalResult DDA::traverse(double max_distance) {
	TraversalResult result;
	result.total_distance = 0.0;
	result.hit_boundary = false;

	while (result.total_distance < max_distance) {
		StepResult step_result = step();

		if (!step_result.valid) {
			// Ray has exited the grid
			result.hit_boundary = true;
			if (!result.voxels.empty()) {
				const auto& last_voxel = result.voxels.back();
				result.exit_voxel = last_voxel.voxel_coords;
				result.exit_position = last_voxel.world_position;

				// Calculate exit normal (opposite of last entry normal)
				result.exit_normal = -last_voxel.entry_normal;
			}
			break;
		}

		// Check if adding this voxel would exceed max_distance
		if (result.total_distance + step_result.distance_traveled > max_distance) {
			// Trim the distance to not exceed max_distance
			step_result.distance_traveled = max_distance - result.total_distance;
			result.total_distance = max_distance;
			result.voxels.push_back(step_result);
			break;
		}

		result.total_distance += step_result.distance_traveled;
		result.voxels.push_back(step_result);
	}

	return result;
}

/**
 * @brief Convert world coordinates to voxel indices
 * @param position World position to convert
 * @return Voxel coordinates, or (-1,-1,-1) if position is outside grid
 *
 * Transforms world coordinates to discrete voxel indices using
 * floor operations. Returns invalid coordinates for out-of-bounds positions.
 */
glm::ivec3 DDA::world_to_voxel(const glm::dvec3& position) const {
	glm::dvec3 grid_pos = (position - grid_origin_) / voxel_size_;

	glm::ivec3 voxel_coords(
		static_cast<int>(floor(grid_pos.x)), static_cast<int>(floor(grid_pos.y)), static_cast<int>(floor(grid_pos.z)));

	// Check bounds
	if (!is_valid_voxel(voxel_coords)) {
		return glm::ivec3(-1, -1, -1); // Invalid voxel
	}

	return voxel_coords;
}

/**
 * @brief Convert voxel indices to world coordinates
 * @param voxel_coords Voxel indices to convert
 * @return World position at center of specified voxel
 *
 * Transforms discrete voxel indices to world coordinates,
 * returning the center point of the specified voxel.
 */
glm::dvec3 DDA::voxel_to_world(const glm::ivec3& voxel_coords) const {
	// Return center of voxel
	return grid_origin_ + (glm::dvec3(voxel_coords) + 0.5) * voxel_size_;
}

/**
 * @brief Check if voxel coordinates are within grid bounds
 * @param voxel_coords Voxel indices to validate
 * @return True if coordinates are within valid grid range
 *
 * Validates that voxel coordinates are within the defined grid dimensions.
 * Used for boundary checking during traversal operations.
 */
bool DDA::is_valid_voxel(const glm::ivec3& voxel_coords) const {
	return (voxel_coords.x >= 0 && voxel_coords.x < grid_dimensions_.x && voxel_coords.y >= 0
			&& voxel_coords.y < grid_dimensions_.y && voxel_coords.z >= 0 && voxel_coords.z < grid_dimensions_.z);
}

/**
 * @brief Get axis-aligned bounding box of specified voxel
 * @param voxel_coords Voxel indices to get bounds for
 * @return Pair containing minimum and maximum bound coordinates
 *
 * Returns the world-space axis-aligned bounding box (AABB) for
 * the specified voxel, useful for intersection testing and visualization.
 */
std::pair<glm::dvec3, glm::dvec3> DDA::get_voxel_bounds(const glm::ivec3& voxel_coords) const {
	glm::dvec3 min_bound = grid_origin_ + glm::dvec3(voxel_coords) * voxel_size_;
	glm::dvec3 max_bound = min_bound + glm::dvec3(voxel_size_);

	return std::make_pair(min_bound, max_bound);
}

/**
 * @brief Calculate surface normal for voxel face entry
 * @param face_axis Axis perpendicular to the entered face (0=x, 1=y, 2=z)
 * @param step_direction Direction of entry (-1 or +1)
 * @return Unit normal vector pointing into the voxel
 *
 * Computes the surface normal vector for the voxel face that was
 * crossed during DDA traversal, used for lighting calculations.
 */
glm::dvec3 DDA::calculate_entry_normal(int face_axis, int step_direction) const {
	glm::dvec3 normal(0.0);
	normal[face_axis] = static_cast<double>(step_direction);
	return normal;
}

/**
 * @brief Calculate current world position along ray
 * @return World coordinates of current DDA position
 *
 * Computes the current world position along the ray based on
 * the DDA algorithm state, used for intersection point calculations.
 */
glm::dvec3 DDA::calculate_world_position() const {
	// Calculate current position along the ray
	double min_dist =
		std::min({side_dist_[0] - delta_dist_[0], side_dist_[1] - delta_dist_[1], side_dist_[2] - delta_dist_[2]});

	return ray_origin_ + ray_direction_ * min_dist;
}
