#pragma once

#include <cstdint>
#include <cmath>

/**
 * @class Material
 * @brief Represents optical material properties for Monte Carlo photon transport simulation
 * 
 * This class encapsulates all optical properties needed for photon interaction calculations,
 * including absorption, scattering, anisotropy, and refractive index. It provides controlled
 * access to these properties and maintains a pre-computed hash for fast boundary detection.
 */
class Material
{
private:
	double g_ {0.00};       ///< Anisotropy coefficient (Henyey-Greenberg parameter)
	double eta_ {1.37};     ///< Refractive index
	double mu_a_ {1.00};    ///< Absorption coefficient (mm^-1)
	double mu_s_ {10.00};   ///< Scattering coefficient (mm^-1)

	/// Pre-computed hash of optical properties for fast comparison
	mutable std::size_t optical_hash_ {0};

	/// Compute a robust hash of all optical properties using FNV-1a algorithm
	void compute_optical_hash() noexcept;

public:
	/// Default constructor uses default member initialization
	Material();

	/// Explicit constructor with all optical properties
	explicit Material(double g_val, double eta_val, double mu_a_val, double mu_s_val) noexcept;

	/// Copy constructor
	Material(const Material& other) = default;

	/// Copy assignment operator
	Material& operator=(const Material& other) = default;

	/// Move constructor
	Material(Material&& other) noexcept = default;

	/// Move assignment operator
	Material& operator=(Material&& other) noexcept = default;

	/// Destructor
	~Material() = default;

	// Getters for read access
	double g() const noexcept { return g_; }
	double eta() const noexcept { return eta_; }
	double mu_a() const noexcept { return mu_a_; }
	double mu_s() const noexcept { return mu_s_; }

	/// Equality operator based on optical properties hash
	bool operator==(const Material& other) const noexcept { return has_same_optical_properties(other); }

	/// Inequality operator
	bool operator!=(const Material& other) const noexcept { return !(*this == other); }

	/// Fast comparison using pre-computed optical properties hash
	bool has_same_optical_properties(const Material& other) const noexcept;

	/// Get the pre-computed hash of optical properties
	std::size_t get_optical_properties_hash() const noexcept;

	/// Update all optical properties and recompute hash (controlled modification)
	void set_optical_properties(double eta_val, double mu_a_val, double mu_s_val, double g_val) noexcept;

	/// Individual setters that maintain hash consistency
	void set_eta(double value) noexcept;
	void set_mu_a(double value) noexcept;
	void set_mu_s(double value) noexcept;
	void set_g(double value) noexcept;
};