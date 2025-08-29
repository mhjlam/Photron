#pragma once

#include "tissue.hpp"

#include <cstdint>

struct Voxel {
	uint32_t ix;
	uint32_t iy;
	uint32_t iz;
	double absorption;
	double emittance;
	Tissue* tissue;

	Voxel() {
		ix = 0;
		iy = 0;
		iz = 0;
		absorption = 0;
		emittance = 0;
		tissue = nullptr;
	}

	Voxel(uint32_t x, uint32_t y, uint32_t z) {
		ix = x;
		iy = y;
		iz = z;
		absorption = 0;
		emittance = 0;
		tissue = nullptr;
	}

	bool operator==(Voxel& other) const {
		return (other.ix == ix && other.iy == iy && other.iz == iz && other.tissue == tissue);
	}
};
