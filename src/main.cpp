/*******************************************************************************
 * PHOTRON: MONTE CARLO PHOTON TRANSPORT RENDERER
 *
 * DATE: 	September 2025
 * AUTHOR:	Maurits H.J. Lam, BSc.
 *
 * DESCRIPTION:
 * 	Monte Carlo photon transport simulation and renderer of
 *  subsurface light transport in voxelized multi-layered mediums.
 *
 * NOTES:
 * 	Photon transport functionality is based on MCML 3.0.0 algorithms,
 * 	which is available at: https://github.com/mhjlam/MCML
 ******************************************************************************/

#include <iostream>

#include "app.hpp"

/**
 * @brief Application entry point for Photron Monte Carlo photon transport simulation
 *
 * Initializes the application, runs the main loop (GUI or headless mode),
 * and ensures proper cleanup before termination.
 *
 * @param argc Number of command line arguments
 * @param argv Array of command line argument strings
 * @return int Exit code (0 for success, 1 for initialization failure)
 */
int main(int argc, char* argv[]) {
	// Initialize application with command line arguments
	App app;
	if (!app.initialize(argc, argv)) {
		return 1;
	}

	// Run main application loop (GUI or headless simulation)
	app.run();

	// Clean shutdown and resource deallocation
	app.shutdown();
	return 0;
}
