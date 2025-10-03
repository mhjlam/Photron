/**
 * @file geom_lookup.cpp
 * @brief Implementation of spatial query and geometry management
 *
 * Contains complete implementations of all spatial operations for Monte Carlo
 * photon transport simulation, including medium lookup, voxel access, and
 * geometric containment tests using DDA acceleration structures.
 */

#include "geom_lookup.hpp"

#include <algorithm>
#include <iostream>

#include "common/config.hpp"
#include "common/error_handler.hpp"
#include "math/cuboid.hpp"
#include "math/dda.hpp"
#include "math/range.hpp"
#include "simulator/medium.hpp"
#include "simulator/photon.hpp"
#include "simulator/voxel.hpp"

// Initialize geometry lookup with medium collection and DDA structures
GeomLookup::GeomLookup(std::vector<Medium>& mediums, const std::vector<std::unique_ptr<DDA>>& ddas) :
	mediums_(&mediums), ddas_(&ddas) {
}

// Find which medium contains the given point using simple containment test
Medium* GeomLookup::find_medium_at(const glm::dvec3& position) const {
	for (auto& medium : *mediums_) {
		if (medium.contains_point(position)) {
			return const_cast<Medium*>(&medium);
		}
	}
	return nullptr; // Point is outside all mediums
}

// Find medium using DDA acceleration for improved boundary precision
Medium* GeomLookup::find_medium_at_with_dda(const glm::dvec3& position) const {
	// Epsilon tolerance for floating-point boundary calculations
	static const double EPSILON = 1e-9;

	// Test each medium using its DDA structure
	for (size_t i = 0; i < mediums_->size() && i < ddas_->size(); ++i) {
		const auto& medium = (*mediums_)[i];
		DDA* dda = (*ddas_)[i].get();
		Medium* mutable_medium = const_cast<Medium*>(&medium);

		// Check if position maps to valid voxel in DDA grid
		glm::ivec3 voxel_coords = dda->world_to_voxel(position);
		if (dda->is_valid_voxel(voxel_coords)) {
			glm::dvec3 mutable_pos = position;
			Voxel* voxel = mutable_medium->voxel_at(mutable_pos);
			if (voxel && voxel->material) {
				return mutable_medium;
			}
		}

		// Stage 2: Epsilon nudging for boundary cases
		glm::dvec3 robust_position = position;
		bool dda_success = false;

		for (double eps_scale = 1.0; eps_scale <= 1000.0 && !dda_success; eps_scale *= 10.0) {
			// Test 8-directional nudging with increasing epsilon
			std::vector<glm::dvec3> nudge_directions = {
				glm::dvec3(-EPSILON * eps_scale, -EPSILON * eps_scale, -EPSILON * eps_scale),
				glm::dvec3(-EPSILON * eps_scale, -EPSILON * eps_scale, EPSILON * eps_scale),
				glm::dvec3(-EPSILON * eps_scale, EPSILON * eps_scale, -EPSILON * eps_scale),
				glm::dvec3(-EPSILON * eps_scale, EPSILON * eps_scale, EPSILON * eps_scale),
				glm::dvec3(EPSILON * eps_scale, -EPSILON * eps_scale, -EPSILON * eps_scale),
				glm::dvec3(EPSILON * eps_scale, -EPSILON * eps_scale, EPSILON * eps_scale),
				glm::dvec3(EPSILON * eps_scale, EPSILON * eps_scale, -EPSILON * eps_scale),
				glm::dvec3(EPSILON * eps_scale, EPSILON * eps_scale, EPSILON * eps_scale)};

			for (const auto& nudge : nudge_directions) {
				glm::dvec3 nudged_pos = position + nudge;
				glm::ivec3 nudged_coords = dda->world_to_voxel(nudged_pos);

				if (dda->is_valid_voxel(nudged_coords)) {
					glm::dvec3 test_pos = nudged_pos;
					Voxel* voxel = mutable_medium->voxel_at(test_pos);
					if (voxel && voxel->material) {
						return mutable_medium;
					}
				}
			}
		}

		// Stage 3: Direct voxel access (bypass DDA)
		glm::dvec3 bypass_pos = position;
		Voxel* direct_voxel = mutable_medium->voxel_at(bypass_pos);
		if (direct_voxel && direct_voxel->material) {
			// Found material voxel via direct access
			return mutable_medium;
		}

		// Stage 4: Geometric fallback with local search
		if (medium.contains_point(position)) {
			// Point is geometrically inside - search nearby for material voxels
			for (double search_radius = EPSILON; search_radius <= 0.01; search_radius *= 10.0) {
				for (int dx = -1; dx <= 1; dx++) {
					for (int dy = -1; dy <= 1; dy++) {
						for (int dz = -1; dz <= 1; dz++) {
							glm::dvec3 search_pos =
								position + glm::dvec3(dx * search_radius, dy * search_radius, dz * search_radius);
							glm::dvec3 test_pos = search_pos;
							Voxel* search_voxel = mutable_medium->voxel_at(test_pos);
							if (search_voxel && search_voxel->material) {
								return mutable_medium;
							}
						}
					}
				}
			}
		}
	}

	return nullptr; // All fallback stages failed - truly in ambient space
}

// Check if position is inside any medium
bool GeomLookup::is_inside_any_medium(const glm::dvec3& position) const {
	return find_medium_at(position) != nullptr;
}

// Test geometric containment across all mediums
bool GeomLookup::is_point_inside_geometry(const glm::dvec3& position) const {
	for (const auto& medium : *mediums_) {
		if (medium.contains_point(position)) {
			return true;
		}
	}
	return false;
}

// Get voxel at world position by finding containing medium
Voxel* GeomLookup::voxel_at(const glm::dvec3& position) const {
	Medium* medium = find_medium_at(position);
	if (medium) {
		glm::dvec3 pos = position; // Copy for non-const medium method
		return medium->voxel_at(pos);
	}
	return nullptr;
}

// Map grid coordinates to voxel via combined world bounds
Voxel* GeomLookup::voxel_grid(uint32_t x, uint32_t y, uint32_t z) const {
	// Convert grid indices to world position
	double voxel_size = Config::get().vox_size();

	// Compute unified bounds across all mediums
	Range3 bounds;
	bool first_medium = true;

	for (const auto& medium : *mediums_) {
		Range3 medium_bounds = medium.get_bounds();
		if (first_medium) {
			bounds = medium_bounds;
			first_medium = false;
		}
		else {
			// Expand bounds to include this medium
			bounds.min_bounds.x = std::min(bounds.min_bounds.x, medium_bounds.min_bounds.x);
			bounds.min_bounds.y = std::min(bounds.min_bounds.y, medium_bounds.min_bounds.y);
			bounds.min_bounds.z = std::min(bounds.min_bounds.z, medium_bounds.min_bounds.z);
			bounds.max_bounds.x = std::max(bounds.max_bounds.x, medium_bounds.max_bounds.x);
			bounds.max_bounds.y = std::max(bounds.max_bounds.y, medium_bounds.max_bounds.y);
			bounds.max_bounds.z = std::max(bounds.max_bounds.z, medium_bounds.max_bounds.z);
		}
	}

	glm::dvec3 position(bounds.min_bounds.x + (x + 0.5) * voxel_size,
						bounds.min_bounds.y + (y + 0.5) * voxel_size,
						bounds.min_bounds.z + (z + 0.5) * voxel_size);

	// Find which medium contains this position and get the voxel
	return voxel_at(position);
}

// Get world-space corners of voxel cuboid
Cuboid GeomLookup::voxel_corners(Voxel* voxel) const {
	if (!voxel) {
		return Cuboid(); // Invalid voxel
	}

	// Use first medium for corner calculation (TODO: multi-medium support)
	if (!mediums_->empty()) {
		return (*mediums_)[0].voxel_corners(voxel);
	}

	return Cuboid(); // No mediums available
}

// Find last surface voxel before photon exit using DDA backtracking
Voxel* GeomLookup::find_last_surface_voxel_with_dda(const Photon& photon, const glm::dvec3& exit_direction) {
	// Locate photon's current medium
	Medium* current_medium = find_medium_at(photon.position);
	if (!current_medium) {
		return nullptr;
	}

	// Find corresponding DDA index
	size_t medium_index = 0;
	for (size_t i = 0; i < mediums_->size(); ++i) {
		if (&(*mediums_)[i] == current_medium) {
			medium_index = i;
			break;
		}
	}

	if (medium_index >= ddas_->size()) {
		// DDA unavailable - use step-back method
		glm::dvec3 step_back_pos = photon.position - exit_direction * 1e-6;
		Medium* last_medium = find_medium_at_with_dda(step_back_pos);
		if (last_medium) {
			glm::dvec3 mutable_pos = step_back_pos;
			return last_medium->voxel_at(mutable_pos);
		}
		return nullptr;
	}

	DDA* dda = (*ddas_)[medium_index].get();

	// Trace backward from current position to find the path
	glm::dvec3 start_pos = photon.position - exit_direction * 0.01; // Start slightly inside medium
	glm::dvec3 end_pos = photon.position;
	glm::dvec3 direction = glm::normalize(end_pos - start_pos);
	double total_distance = glm::length(end_pos - start_pos);

	if (total_distance < 1e-12) {
		// Very short distance, use current voxel
		return photon.voxel;
	}

	// Initialize DDA for this ray
	dda->initialize_ray(start_pos, direction);

	// Traverse voxels using DDA
	DDA::TraversalResult result = dda->traverse(total_distance);

	// Find last surface voxel along DDA traversal path
	Voxel* last_surface_voxel = nullptr;

	for (size_t i = 0; i < result.voxels.size(); ++i) {
		const auto& step = result.voxels[i];
		// Access voxel at DDA step position
		glm::dvec3 mutable_pos = step.world_position;
		Voxel* voxel = current_medium->voxel_at(mutable_pos);

		if (voxel && voxel->material) {
			// Check if voxel is marked as external surface
			if (voxel->is_surface_voxel) {
				last_surface_voxel = voxel;
			}
		}
	}

	// Return found surface voxel or fallback to current photon voxel
	return last_surface_voxel ? last_surface_voxel : photon.voxel;
}
