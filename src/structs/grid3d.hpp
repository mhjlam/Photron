#pragma once

#include "voxel.hpp"

struct Grid3D {
	unsigned int nx;            // number of voxels in x direction
	unsigned int ny;            // number of voxels in y direction
	unsigned int nz;            // number of voxels in z direction

	unsigned long numvox;       // number of voxels
	double voxsize;             // voxel size (dx = dy = dz)

	std::vector<Voxel*> tensor; // 3-dimensional array for Voxel pointers

	Grid3D() {
		nx = 0;
		ny = 0;
		nz = 0;

		numvox = 0;
		voxsize = 0;

		tensor = std::vector<Voxel*>();
	}

	Grid3D(double vsize, unsigned short numx, unsigned short numy, unsigned short numz) {
		nx = numx;
		ny = numy;
		nz = numz;

		numvox = nx * ny * nz;
		voxsize = vsize;

		tensor = std::vector<Voxel*>(numvox, 0);
		for (unsigned short i = 0; i < nx; ++i)
			for (unsigned short j = 0; j < ny; ++j)
				for (unsigned short k = 0; k < nz; ++k)
					tensor.at(k * nx * ny + j * nx + i) = new Voxel(i, j, k);
	}

	~Grid3D() {
		for (unsigned long i = 0; i < tensor.size(); ++i)
			delete tensor[i];
	}

	unsigned int index(unsigned short x, unsigned short y, unsigned short z) {
		// [x][y][z] = [x + y*nx + z*nx*ny]
		return (x + y * nx + z * nx * ny);
	}

	Voxel* operator()(unsigned short x, unsigned short y, unsigned short z) {
		return tensor.at(x + y * nx + z * nx * ny);
	}
};
