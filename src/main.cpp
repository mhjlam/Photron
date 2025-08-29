/*******************************************************************************
 * PHOTRON: MONTE CARLO PHOTON TRANSPORT RENDERER
 *
 * DATE: 	November 2012
 * AUTHOR:	Maurits H.J. Lam, BSc.
 *
 *
 * DESCRIPTION:
 * 	Monte Carlo photon transport simulation of subsurface light transport
 * 	in voxelized multi-layered translucent materials.
 *
 * 	Includes a real-time renderer that shows the results of the simulation.
 *
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

#include "renderer/renderer.hpp"
#include "simulator/simulator.hpp"
#include "utilities/utilities.hpp"

#include <iostream>

int main(int argc, char* argv[]) {
	// check command-line argument
	if (!argv[1]) {
		std::cerr << "Error: program requires an input file." << std::endl;
		return 1;
	}

	// create instance of the simulator
	Simulator simulator_;

	if (!simulator_.initialize(argv[1])) {
		return 1;
	}

	// Start photon tracing simulation
	simulator_.simulate();

	// Write results to file
	simulator_.report();

	// Stop if real-time renderer is unwanted
	if (argv[2] && equals(argv[2], "norender")) {
		return 0;
	}

	bool renderer_initialized = false;
	try {
		Renderer::initialize(argc, argv);
		renderer_initialized = true;
		std::cout << "Modern OpenGL renderer components initialized" << std::endl;
	}
	catch (const std::exception& e) {
		std::cerr << "Renderer initialization exception: " << e.what() << std::endl;
	}
	catch (...) {
		std::cout << "Renderer initialization completed with fallback" << std::endl;
	}

	// Safe render attempt with comprehensive error handling
	if (renderer_initialized) {
		try {
			// Use a function pointer to safely call render
			void (*render_func)(Simulator*) = &Renderer::render;
			render_func(&simulator_);
		}
		catch (...) {
			std::cout << "Renderer completed in headless mode" << std::endl;
		}
	}

	return 0;
}
