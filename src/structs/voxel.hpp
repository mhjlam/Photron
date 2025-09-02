#pragma once

#include "glm_types.hpp"
#include "tissue.hpp"

struct Voxel {
	// Use GLM internally but expose original interface
	glm::uvec3 coords;     // (ix, iy, iz) coordinates
	double absorption;
	double emittance;
	Tissue* tissue;

	Voxel() : coords(0), absorption(0), emittance(0), tissue(nullptr) {}

	Voxel(uint32_t x, uint32_t y, uint32_t z) : coords(x, y, z), absorption(0), emittance(0), tissue(nullptr) {}

	// Backwards compatibility accessors
	uint32_t& ix() { return coords.x; }
	uint32_t& iy() { return coords.y; }
	uint32_t& iz() { return coords.z; }

	const uint32_t& ix() const { return coords.x; }
	const uint32_t& iy() const { return coords.y; }
	const uint32_t& iz() const { return coords.z; }

	// GLM convenience methods
	const glm::uvec3& indices() const { return coords; }
	void set_indices(const glm::uvec3& idx) { coords = idx; }

	bool operator==(const Voxel& other) const {
		return (other.coords == coords && other.tissue == tissue);
	}
};
