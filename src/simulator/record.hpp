#pragma once

struct Record
{
	double total_absorption {0.0};      // Total absorption
	double diffuse_reflection {0.0};    // Diffuse reflection
	double specular_reflection {0.0};   // Specular reflection
	double surface_refraction {0.0};    // Energy entering medium at surface
	double diffuse_transmission {0.0};  // Diffuse transmission
	double specular_transmission {0.0}; // Specular transmission
	double avg_path_length {0.0};       // Average path length
	int total_steps {0};                // Total simulation steps
	int photons_entered {0};            // Number of photons that entered this medium

	// Default constructor uses default member initialization
	Record() = default;
};
