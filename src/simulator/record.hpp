#pragma once

struct Record
{
	double total_absorption;      // Total absorption
	double diffuse_reflection;    // Diffuse reflection
	double specular_reflection;   // Specular reflection
	double diffuse_transmission;  // Diffuse transmission
	double specular_transmission; // Specular transmission

	Record() :
		total_absorption(0), diffuse_reflection(0), specular_reflection(0), diffuse_transmission(0),
		specular_transmission(0) {}
};
