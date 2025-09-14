#pragma once

#include <cstdint>

struct Tissue
{
	uint8_t id {0};
	
	double g {0.00};       // anisotropy coefficient
	double eta {1.37};     // index of refraction
	double mu_a {1.00};    // absorption coefficient
	double mu_s {10.00};   // scattering coefficient

	// Default constructor uses default member initialization
	Tissue() = default;

	// Explicit constructor with member initializer list
	explicit Tissue(uint8_t i, double g_val, double eta_val, double mu_a_val, double mu_s_val) noexcept
		: id(i), g(g_val), eta(eta_val), mu_a(mu_a_val), mu_s(mu_s_val) {
	}

	bool operator==(const Tissue& other) const noexcept { return other.id == id; }
};
