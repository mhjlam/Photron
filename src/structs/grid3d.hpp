#pragma once

#include "voxel.hpp"

#include <cstdint>
#include <vector>

struct Grid3D {
	uint32_t nx;                // number of voxels in x direction
	uint32_t ny;                // number of voxels in y direction
	uint32_t nz;                // number of voxels in z direction

	uint64_t num_vox;           // number of voxels
	double vox_size;            // voxel size (dx = dy = dz)

	std::vector<Voxel*> tensor; // 3-dimensional array for Voxel pointers

	Grid3D() {
		nx = 0;
		ny = 0;
		nz = 0;
		num_vox = 0;
		vox_size = 0;
		tensor = std::vector<Voxel*>();
	}

	Grid3D(double vsize, uint8_t numx, uint8_t numy, uint8_t numz) {
		nx = numx;
		ny = numy;
		nz = numz;

		num_vox = nx * ny * nz;
		vox_size = vsize;

		tensor = std::vector<Voxel*>(num_vox, 0);
		for (uint8_t i = 0; i < nx; ++i) {
			for (uint8_t j = 0; j < ny; ++j) {
				for (uint8_t k = 0; k < nz; ++k) {
					tensor.at(k * nx * ny + j * nx + i) = new Voxel(i, j, k);
				}
			}
		}
	}

	~Grid3D() {
		for (uint64_t i = 0; i < tensor.size(); ++i) {
			delete tensor[i];
		}
	}

	uint32_t index(uint8_t x, uint8_t y, uint8_t z) {
		// [x][y][z] = [x + y*nx + z*nx*ny]
		return (x + y * nx + z * nx * ny);
	}

	Voxel* operator()(uint8_t x, uint8_t y, uint8_t z) { return tensor.at(x + y * nx + z * nx * ny); }
};
