#pragma once

struct Record
{
	double total_absorption = 0.0;      // Total absorption
	double diffuse_reflection = 0.0;    // Diffuse reflection
	double specular_reflection = 0.0;   // Specular reflection
	double diffuse_transmission = 0.0;  // Diffuse transmission
	double specular_transmission = 0.0; // Specular transmission

	// Default constructor uses default member initialization
	Record() = default;
};
