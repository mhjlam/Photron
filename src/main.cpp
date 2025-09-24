/*******************************************************************************
 * PHOTRON: MONTE CARLO PHOTON TRANSPORT RENDERER
 *
 * DATE: 	September 2025
 * AUTHOR:	Maurits H.J. Lam, BSc.
 *
 * DESCRIPTION:
 * 	Monte Carlo photon transport simulation of subsurface light transport
 * 	in voxelized multi-layered translucent materials.
 *
 * 	Includes a real-time renderer that shows the results of the simulation.
 *
 * NOTES:
 * 	Some of this program is based on MCVM code, published by Ting Li in 2009,
 * 	available at http://code.google.com/p/mcvm.
 *
 * 	MCVM is described in the following publication:
 * 	T. Li, H. Gong, and Q. Luo,
 * 	"MCVM: Monte Carlo Modeling of photon migration in voxelized media",
 * 	Journal of Innovative Optical Health Sciences 3(2), 91-102 (2010).
 ******************************************************************************/

#include <iostream>

#include "app.hpp"

int main(int argc, char* argv[]) {
	App app;
	if (!app.initialize(argc, argv)) {
		return 1;
	}

	app.run();
	app.shutdown();
	return 0;
}
