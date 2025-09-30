/**
 * @file material.cpp
 * @brief Implementation of optical material properties for photon transport
 *
 * Implements the Material class which encapsulates optical properties used in
 * Monte Carlo photon transport simulations. Features efficient property comparison
 * through hash-based caching and provides setter methods that maintain consistency.
 */

#include "material.hpp"

Material::Material() {
	// Initialize with default optical properties and compute hash
	compute_optical_hash();
}

Material::Material(double g_val, double eta_val, double mu_a_val, double mu_s_val) noexcept :
	g_(g_val), eta_(eta_val), mu_a_(mu_a_val), mu_s_(mu_s_val) {
	// Initialize with specified properties and compute identification hash
	compute_optical_hash();
}

bool Material::has_same_optical_properties(const Material& other) const noexcept {
	// Fast optical property comparison using precomputed hashes
	return optical_hash_ == other.optical_hash_;
}

std::size_t Material::get_optical_properties_hash() const noexcept {
	// Return cached hash for external comparison and grouping
	return optical_hash_;
}

void Material::set_optical_properties(double eta_val, double mu_a_val, double mu_s_val, double g_val) noexcept {
	// Batch update all optical properties with single hash recomputation
	eta_ = eta_val;
	mu_a_ = mu_a_val;
	mu_s_ = mu_s_val;
	g_ = g_val;

	compute_optical_hash();
}

void Material::set_eta(double value) noexcept {
	// Update refractive index and refresh identification hash
	eta_ = value;
	compute_optical_hash();
}

void Material::set_mu_a(double value) noexcept {
	// Update absorption coefficient and refresh identification hash
	mu_a_ = value;
	compute_optical_hash();
}

void Material::set_mu_s(double value) noexcept {
	// Update scattering coefficient and refresh identification hash
	mu_s_ = value;
	compute_optical_hash();
}

void Material::set_g(double value) noexcept {
	// Update anisotropy factor and refresh identification hash
	g_ = value;
	compute_optical_hash();
}

void Material::compute_optical_hash() noexcept {
	// Compute FNV-1a hash for fast optical property comparison
	// Using FNV-1a hash algorithm for consistent cross-platform results
	constexpr const std::size_t FNV_PRIME = 1099511628211ULL;
	constexpr const std::size_t FNV_OFFSET_BASIS = 14695981039346656037ULL;

	optical_hash_ = FNV_OFFSET_BASIS;

	// Hash each double property as raw bytes for precision
	auto hash_double = [&](double value) {
		const auto* bytes = reinterpret_cast<const uint8_t*>(&value);
		for (size_t i = 0; i < sizeof(double); ++i) {
			optical_hash_ ^= bytes[i];
			optical_hash_ *= FNV_PRIME;
		}
	};

	// Hash properties in consistent order for deterministic results
	hash_double(eta_);  // Refractive index
	hash_double(mu_a_); // Absorption coefficient
	hash_double(mu_s_); // Scattering coefficient
	hash_double(g_);    // Anisotropy factor
}
