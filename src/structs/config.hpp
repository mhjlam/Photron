#pragma once

#include <cstdint>

struct Config {
	uint32_t nx;           // number of voxels in the x direction
	uint32_t ny;           // number of voxels in the y direction
	uint32_t nz;           // number of voxels in the z direction

	uint64_t num_layers;  // number of layers
	uint64_t num_voxels;  // number of voxels
	uint64_t num_photons; // number of photons
	uint64_t num_sources; // number of light sources

	double vox_size;      // uniform size of each voxel (dx=dy=dz)
	double ambient_eta;   // refractive index of ambient medium

	bool partial;         // partial reflection
	bool progress;        // display progress messages

	Config() {
		nx = ny = nz = 0;

		num_layers = 0;
		num_voxels = 0;
		num_photons = 0;
		num_sources = 0;

		vox_size = 0.0;
		ambient_eta = 0.0;

		partial = true;
		progress = false;
	}
};
