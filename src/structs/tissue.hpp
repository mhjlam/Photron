#pragma once

#include <cstdint>

struct Tissue {
	uint8_t id;

	double ani; // anisotropy coefficient
	double eta; // index of refraction
	double mua; // absorption coefficient
	double mus; // scattering coefficient

	Tissue() {
		id = 0;
		ani = 0.00;
		eta = 1.37;
		mua = 1.00;
		mus = 10.00;
	}

	Tissue(uint8_t i, double g, double n, double a, double s) {
		id = i;
		ani = g;
		eta = n;
		mua = a;
		mus = s;
	}

	bool operator==(const Tissue& other) const { return other.id == id; }
};
