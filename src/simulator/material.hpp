/**
 * @file material.hpp
 * @brief Optical material properties for Monte Carlo photon transport simulation
 * 
 * Defines the Material class which encapsulates all optical properties needed
 * for accurate photon-material interaction calculations in biological tissues
 * and synthetic materials.
 */

#pragma once

#include <cstdint>
#include <cmath>

/**
 * @class Material
 * @brief Optical material properties for photon transport calculations
 * 
 * The Material class represents the complete set of optical properties required
 * for Monte Carlo photon transport simulation in biological tissues or synthetic
 * materials. It provides:
 * 
 * **Optical Properties:**
 * - **Absorption Coefficient (μₐ)**: Rate of photon absorption per unit length
 * - **Scattering Coefficient (μₛ)**: Rate of photon scattering per unit length  
 * - **Anisotropy Factor (g)**: Henyey-Greenberg phase function parameter
 * - **Refractive Index (n)**: Optical density for Fresnel calculations
 * 
 * **Key Features:**
 * - **Hash-based Comparison**: Fast material boundary detection using pre-computed hashes
 * - **Physics Validation**: Ensures optical properties remain within physical bounds
 * - **Efficient Storage**: Compact representation with computed derived properties
 * - **Type Safety**: Controlled access to prevent inconsistent state
 * 
 * **Physical Interpretation:**
 * - g = 0: Isotropic scattering (equal probability in all directions)
 * - g > 0: Forward scattering (biological tissues typically 0.8-0.95)
 * - g < 0: Backward scattering (rare in biological applications)
 * - μₐ + μₛ = μₜ (total attenuation coefficient)
 * 
 * **Typical Values for Biological Tissues:**
 * - Skin: μₐ=0.1-1.0 mm⁻¹, μₛ=10-40 mm⁻¹, g=0.8-0.9, n=1.4
 * - Blood: μₐ=0.2-10 mm⁻¹, μₛ=20-100 mm⁻¹, g=0.95-0.99, n=1.35
 */
class Material
{
private:
	double g_ {0.00};       ///< Anisotropy coefficient [-1,1] (Henyey-Greenberg g parameter)
	double eta_ {1.37};     ///< Refractive index [1,inf) (typically 1.33-1.45 for biological tissues)
	double mu_a_ {1.00};    ///< Absorption coefficient [0,inf) mm⁻¹
	double mu_s_ {10.00};   ///< Scattering coefficient [0,inf) mm⁻¹

	/// Pre-computed hash for O(1) material comparison and boundary detection
	mutable std::size_t optical_hash_ {0};

	/**
	 * @brief Compute robust hash of optical properties using FNV-1a algorithm
	 * 
	 * Creates a hash that uniquely identifies material properties for fast
	 * comparison operations. Uses bit-level representation to handle floating
	 * point precision consistently.
	 */
	void compute_optical_hash() noexcept;

public:
	/**
	 * @brief Default constructor with typical biological tissue properties
	 * 
	 * Initializes material with moderate absorption and scattering properties
	 * representative of human skin tissue.
	 */
	Material();

	/**
	 * @brief Construct material with specified optical properties
	 * 
	 * @param g_val Anisotropy factor [-1,1] (0=isotropic, >0=forward scattering)
	 * @param eta_val Refractive index [1,inf) (1.33-1.45 typical for tissues)
	 * @param mu_a_val Absorption coefficient [0,inf) mm⁻¹
	 * @param mu_s_val Scattering coefficient [0,inf) mm⁻¹
	 */
	explicit Material(double g_val, double eta_val, double mu_a_val, double mu_s_val) noexcept;

	/// Standard copy/move semantics with automatic hash recomputation
	Material(const Material& other) = default;
	Material& operator=(const Material& other) = default;
	Material(Material&& other) noexcept = default;
	Material& operator=(Material&& other) noexcept = default;
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

	/// Set all optical properties and recompute hash (controlled modification)
	void set_optical_properties(double eta_val, double mu_a_val, double mu_s_val, double g_val) noexcept;

	/// Individual setters that maintain hash consistency
	void set_eta(double value) noexcept;
	void set_mu_a(double value) noexcept;
	void set_mu_s(double value) noexcept;
	void set_g(double value) noexcept;
};
