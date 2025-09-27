#include "material.hpp"

Material::Material() {
	compute_optical_hash();
}

Material::Material(double g_val, double eta_val, double mu_a_val, double mu_s_val) noexcept
	: g_(g_val), eta_(eta_val), mu_a_(mu_a_val), mu_s_(mu_s_val) {
	compute_optical_hash();
}

bool Material::has_same_optical_properties(const Material& other) const noexcept {
	return optical_hash_ == other.optical_hash_;
}

std::size_t Material::get_optical_properties_hash() const noexcept {
	return optical_hash_;
}

void Material::set_optical_properties(double eta_val, double mu_a_val, double mu_s_val, double g_val) noexcept {
	eta_ = eta_val;
	mu_a_ = mu_a_val;
	mu_s_ = mu_s_val;
	g_ = g_val;
	compute_optical_hash();
}

void Material::set_eta(double value) noexcept {
	eta_ = value;
	compute_optical_hash();
}

void Material::set_mu_a(double value) noexcept {
	mu_a_ = value;
	compute_optical_hash();
}

void Material::set_mu_s(double value) noexcept {
	mu_s_ = value;
	compute_optical_hash();
}

void Material::set_g(double value) noexcept {
	g_ = value;
	compute_optical_hash();
}

void Material::compute_optical_hash() noexcept {
	const std::size_t FNV_PRIME = 1099511628211ULL;
	const std::size_t FNV_OFFSET_BASIS = 14695981039346656037ULL;
	
	optical_hash_ = FNV_OFFSET_BASIS;
	
	// Hash each property as bytes for consistent results across platforms
	auto hash_double = [&](double value) {
		const auto* bytes = reinterpret_cast<const uint8_t*>(&value);
		for (size_t i = 0; i < sizeof(double); ++i) {
			optical_hash_ ^= bytes[i];
			optical_hash_ *= FNV_PRIME;
		}
	};
	
	hash_double(eta_);
	hash_double(mu_a_);  
	hash_double(mu_s_);
	hash_double(g_);
}
