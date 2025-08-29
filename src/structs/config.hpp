#pragma once

struct Config {
	unsigned short nx;        // number of voxels in the x direction
	unsigned short ny;        // number of voxels in the y direction
	unsigned short nz;        // number of voxels in the z direction

	unsigned long numlayers;  // number of layers
	unsigned long numvoxels;  // number of voxels
	unsigned long numphotons; // number of photons
	unsigned long numsources; // number of light sources

	double voxsize;           // uniform size of each voxel (dx=dy=dz)
	double ambienteta;        // refractive index of ambient medium

	bool partial;             // partial reflection
	bool progress;            // display progress messages

	Config() {
		nx = ny = nz = 0;

		numlayers = 0;
		numvoxels = 0;
		numphotons = 0;
		numsources = 0;

		voxsize = 0.0;
		ambienteta = 0.0;

		partial = true;
		progress = false;
	}
};
