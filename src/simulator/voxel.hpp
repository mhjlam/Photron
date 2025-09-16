#pragma once

#include <glm/glm.hpp>

#include "simulator/tissue.hpp"

struct Voxel
{
private:
	// Volume fraction thresholds for classification
	static constexpr double FULLY_INSIDE_THRESHOLD = 0.999;
	static constexpr double FULLY_OUTSIDE_THRESHOLD = 0.999;
	static constexpr double PARTIAL_THRESHOLD = 0.001;
	static constexpr double EMITTANCE_THRESHOLD = 0.01;  // Minimum 1% volume overlap for emittance recording

public:
	// Use GLM internally but expose original interface
	glm::uvec3 coords; // (ix, iy, iz) coordinates
	double absorption;
	double emittance;
	Tissue* tissue;

	// Exit direction tracking for energy conservation
	double emittance_transmitted;  // Energy exiting via transmission (no angle change)
	double emittance_reflected;    // Energy exiting via specular reflection (angle changed)
	double emittance_diffuse;      // Energy exiting via diffuse reflection (scattered)

	// Partial volume support for boundary physics
	double volume_fraction_inside;  // Fraction of voxel volume inside the geometry [0.0, 1.0]
	double volume_fraction_outside; // Fraction of voxel volume outside the geometry [0.0, 1.0]
	bool is_boundary_voxel;         // True if voxel intersects geometry boundary
	bool is_surface_voxel;          // True if voxel is at the outer surface of the voxel cluster

	Voxel() :
		coords(0), absorption(0), emittance(0), tissue(nullptr), 
		emittance_transmitted(0), emittance_reflected(0), emittance_diffuse(0),
		volume_fraction_inside(0.0), volume_fraction_outside(0.0), 
		is_boundary_voxel(false), is_surface_voxel(false) {}

	Voxel(uint32_t x, uint32_t y, uint32_t z) :
		coords(x, y, z), absorption(0), emittance(0), tissue(nullptr), 
		emittance_transmitted(0), emittance_reflected(0), emittance_diffuse(0),
		volume_fraction_inside(0.0), volume_fraction_outside(0.0), 
		is_boundary_voxel(false), is_surface_voxel(false) {}

	constexpr const uint32_t& ix() const { return coords.x; }
	constexpr const uint32_t& iy() const { return coords.y; }
	constexpr const uint32_t& iz() const { return coords.z; }

	// GLM convenience methods
	constexpr const glm::uvec3& indices() const { return coords; }
	void set_indices(const glm::uvec3& idx) { coords = idx; }

	// Volume fraction utilities
	constexpr bool is_fully_inside() const { return volume_fraction_inside >= FULLY_INSIDE_THRESHOLD; }
	constexpr bool is_fully_outside() const { return volume_fraction_outside >= FULLY_OUTSIDE_THRESHOLD; }
	bool is_partial() const {
		return is_boundary_voxel && volume_fraction_inside > PARTIAL_THRESHOLD
			   && volume_fraction_outside > PARTIAL_THRESHOLD;
	}

	// Energy tracking utilities
	double total_emittance() const { return emittance_transmitted + emittance_reflected + emittance_diffuse; }
	void reset_energy() { 
		absorption = 0; 
		emittance = 0; 
		emittance_transmitted = 0; 
		emittance_reflected = 0; 
		emittance_diffuse = 0; 
	}

		// Special comparison overloads
	constexpr bool operator==(const Voxel& other) const { return (other.coords == coords && other.tissue == tissue); }
};
