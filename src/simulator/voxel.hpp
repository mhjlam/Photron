#pragma once

#include "math/glm_types.hpp"
#include "simulator/tissue.hpp"

struct Voxel
{
private:
	// Volume fraction thresholds for classification
	static constexpr double FULLY_INSIDE_THRESHOLD = 0.999;
	static constexpr double FULLY_OUTSIDE_THRESHOLD = 0.999;
	static constexpr double PARTIAL_THRESHOLD = 0.001;

public:
	// Use GLM internally but expose original interface
	glm::uvec3 coords; // (ix, iy, iz) coordinates
	double absorption;
	double emittance;
	Tissue* tissue;
	
	// Partial volume support for boundary physics
	double volume_fraction_inside;  // Fraction of voxel volume inside the geometry [0.0, 1.0]
	double volume_fraction_outside; // Fraction of voxel volume outside the geometry [0.0, 1.0]
	bool is_boundary_voxel;         // True if voxel intersects geometry boundary
	
	Voxel() : coords(0), absorption(0), emittance(0), tissue(nullptr), 
			  volume_fraction_inside(0.0), volume_fraction_outside(0.0), is_boundary_voxel(false) {}

	Voxel(uint32_t x, uint32_t y, uint32_t z) : coords(x, y, z), absorption(0), emittance(0), tissue(nullptr),
												 volume_fraction_inside(0.0), volume_fraction_outside(0.0), is_boundary_voxel(false) {}

	const uint32_t& ix() const { return coords.x; }
	const uint32_t& iy() const { return coords.y; }
	const uint32_t& iz() const { return coords.z; }

	// GLM convenience methods
	const glm::uvec3& indices() const { return coords; }
	void set_indices(const glm::uvec3& idx) { coords = idx; }
	
	// Volume fraction utilities
	bool is_fully_inside() const { return volume_fraction_inside >= FULLY_INSIDE_THRESHOLD; }
	bool is_fully_outside() const { return volume_fraction_outside >= FULLY_OUTSIDE_THRESHOLD; }
	bool is_partial() const { return is_boundary_voxel && volume_fraction_inside > PARTIAL_THRESHOLD && volume_fraction_outside > PARTIAL_THRESHOLD; }

	bool operator==(const Voxel& other) const { return (other.coords == coords && other.tissue == tissue); }
};
