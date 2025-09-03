#pragma once

#include <cstdint>

struct Tissue
{
	uint8_t id;

	double g;    // anisotropy coefficient
	double eta;  // index of refraction
	double mu_a; // absorption coefficient
	double mu_s; // scattering coefficient

	Tissue() {
		id = 0;
		g = 0.00;
		eta = 1.37;
		mu_a = 1.00;
		mu_s = 10.00;
	}

	Tissue(uint8_t i, double g, double n, double a, double s) {
		id = i;
		g = g;
		eta = n;
		mu_a = a;
		mu_s = s;
	}

	bool operator==(const Tissue& other) const { return other.id == id; }
};
