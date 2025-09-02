#pragma once

#include "glm_types.hpp"
#include "voxel.hpp"

#include <vector>

struct Grid3D {
	// Use GLM internally but expose original interface
	glm::uvec3 dimensions;      // (nx, ny, nz) - number of voxels in each direction
	uint64_t num_vox;           // number of voxels
	double vox_size;            // voxel size (dx = dy = dz)

	std::vector<Voxel*> tensor; // 3-dimensional array for Voxel pointers

	Grid3D() : dimensions(0), num_vox(0), vox_size(0), tensor() {}

	Grid3D(double vsize, uint8_t numx, uint8_t numy, uint8_t numz) 
		: dimensions(numx, numy, numz), vox_size(vsize) {
		
		num_vox = dimensions.x * dimensions.y * dimensions.z;
		tensor = std::vector<Voxel*>(num_vox, 0);
		
		for (uint8_t i = 0; i < dimensions.x; ++i) {
			for (uint8_t j = 0; j < dimensions.y; ++j) {
				for (uint8_t k = 0; k < dimensions.z; ++k) {
					tensor.at(k * dimensions.x * dimensions.y + j * dimensions.x + i) = new Voxel(i, j, k);
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
		return (x + y * dimensions.x + z * dimensions.x * dimensions.y);
	}

	Voxel* operator()(uint8_t x, uint8_t y, uint8_t z) { 
		return tensor.at(x + y * dimensions.x + z * dimensions.x * dimensions.y); 
	}

	// Backwards compatibility accessors
	uint32_t& nx() { return dimensions.x; }
	uint32_t& ny() { return dimensions.y; }
	uint32_t& nz() { return dimensions.z; }

	const uint32_t& nx() const { return dimensions.x; }
	const uint32_t& ny() const { return dimensions.y; }
	const uint32_t& nz() const { return dimensions.z; }
};
