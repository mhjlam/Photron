#pragma once

#include "tissue.hpp"

struct Voxel {
	unsigned int ix;
	unsigned int iy;
	unsigned int iz;

	double absorption;
	double emittance;

	Tissue* tissue;

	Voxel() {
		ix = 0;
		iy = 0;
		iz = 0;

		absorption = 0;
		emittance = 0;

		tissue = NULL;
	}

	Voxel(unsigned int x, unsigned int y, unsigned int z) {
		ix = x;
		iy = y;
		iz = z;

		absorption = 0;
		emittance = 0;

		tissue = NULL;
	}

	bool operator==(Voxel& other) const {
		return (other.ix == ix && other.iy == iy && other.iz == iz && other.tissue == tissue);
	}
};
