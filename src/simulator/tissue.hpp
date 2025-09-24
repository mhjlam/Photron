#pragma once

#include <cstdint>
#include <cmath>
#include <functional>

struct Tissue
{
private:
	uint8_t id_ {0};
	double g_ {0.00};       // anisotropy coefficient
	double eta_ {1.37};     // index of refraction
	double mu_a_ {1.00};    // absorption coefficient
	double mu_s_ {10.00};   // scattering coefficient

	// Pre-computed hash of optical properties for fast comparison
	mutable std::size_t optical_hash_ {0};

public:
	// Getters for read access
	uint8_t id() const noexcept { return id_; }
	double g() const noexcept { return g_; }
	double eta() const noexcept { return eta_; }
	double mu_a() const noexcept { return mu_a_; }
	double mu_s() const noexcept { return mu_s_; }

	// Default constructor uses default member initialization
	Tissue() { compute_optical_hash(); }

	// Explicit constructor with member initializer list
	explicit Tissue(uint8_t i, double g_val, double eta_val, double mu_a_val, double mu_s_val) noexcept
		: id_(i), g_(g_val), eta_(eta_val), mu_a_(mu_a_val), mu_s_(mu_s_val) {
		compute_optical_hash();
	}

	bool operator==(const Tissue& other) const noexcept { return other.id_ == id_; }
	
	// Fast comparison using pre-computed optical properties hash
	bool has_same_optical_properties(const Tissue& other) const noexcept {
		return optical_hash_ == other.optical_hash_;
	}
	
	// Get the pre-computed hash of optical properties
	std::size_t get_optical_properties_hash() const noexcept {
		return optical_hash_;
	}
	
	// Update optical properties and recompute hash (controlled modification)
	void set_optical_properties(double eta_val, double mu_a_val, double mu_s_val, double g_val) noexcept {
		eta_ = eta_val;
		mu_a_ = mu_a_val;
		mu_s_ = mu_s_val;
		g_ = g_val;
		compute_optical_hash();
	}
	
	// Individual setters that maintain hash consistency
	void set_eta(double value) noexcept { eta_ = value; compute_optical_hash(); }
	void set_mu_a(double value) noexcept { mu_a_ = value; compute_optical_hash(); }
	void set_mu_s(double value) noexcept { mu_s_ = value; compute_optical_hash(); }
	void set_g(double value) noexcept { g_ = value; compute_optical_hash(); }
	void set_id(uint8_t value) noexcept { id_ = value; } // ID doesn't affect optical hash

private:
	// Compute a robust hash of all optical properties using FNV-1a algorithm
	void compute_optical_hash() noexcept {
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
};
