/**
 * @file transport.hpp
 * @brief Photon transport and interaction handlers for Monte Carlo simulation
 *
 * The Transport class manages complex photon transport operations including
 * boundary crossings, medium transitions, scattering events, and energy deposition.
 * This separates the low-level transport physics from the main simulation loop.
 */

#pragma once

#include <memory>
#include <vector>

#include <glm/glm.hpp>

#include "math/random.hpp"

// Forward declarations
class Photon;
class Medium;
class Voxel;
class Metrics;
class DDA;
struct Source;
class Simulator;

/**
 * @class Transport
 * @brief Handles detailed photon transport physics and interactions
 *
 * The Transport class encapsulates all the complex physics calculations
 * for photon transport including boundary crossings, medium transitions,
 * scattering events, reflection, transmission, and energy deposition.
 * It works in conjunction with the main Simulator but handles the detailed
 * physics calculations that can be separated from the simulation loop.
 */
class Transport
{
public:
	/**
	 * @brief Constructor with references to simulation data
	 * @param mediums Reference to mediums vector (not owned)
	 * @param metrics Reference to metrics object (not owned)
	 * @param dda_instances Reference to DDA instances (not owned)
	 * @param rng Reference to random number generator (not owned)
	 * @param simulator Reference to main simulator for delegation (not owned)
	 */
	Transport(std::vector<Medium>& mediums,
			  Metrics& metrics,
			  const std::vector<std::unique_ptr<DDA>>& dda_instances,
			  std::shared_ptr<Random> rng,
			  Simulator* simulator);

	/**
	 * @brief Handle photon boundary crossing between media
	 * @param photon Reference to photon being transported
	 *
	 * Manages complex photon transitions between different media including
	 * reflection, refraction, and transmission calculations. Handles Fresnel
	 * equations and energy conservation during boundary interactions.
	 */
	void cross(Photon& photon);

	/**
	 * @brief Handle sub-stepping for voxel traversal
	 * @param photon Reference to photon being transported
	 *
	 * Manages detailed photon movement through voxelized geometry including
	 * accurate path tracking, energy deposition, and boundary detection.
	 * Uses DDA algorithms for efficient voxel traversal.
	 */
	void sub_step(Photon& photon);

	/**
	 * @brief Handle photon radiation and exit processing
	 * @param photon Reference to photon being radiated
	 * @param direction Direction vector for radiation
	 * @param weight Energy weight for radiation calculation
	 *
	 * Processes photon exit from simulation geometry including exit point
	 * calculation, energy accounting, and emitter generation for visualization.
	 */
	void radiate(Photon& photon, glm::dvec3& direction, double weight);

	/**
	 * @brief Handle photon transfer between positions
	 * @param photon Reference to photon being transferred
	 *
	 * Manages photon movement including step size calculation, medium transitions,
	 * and path tracking. Coordinates with voxel traversal and energy deposition.
	 */
	void transfer(Photon& photon);

	/**
	 * @brief Deposit photon energy in current voxel based on absorption
	 * @param photon Reference to photon depositing energy
	 *
	 * Calculates energy deposition using absorption coefficient and updates
	 * voxel energy statistics for visualization and analysis.
	 */
	void deposit(Photon& photon);

	/**
	 * @brief Handle medium transitions and boundary physics
	 * @param photon Reference to photon at medium boundary
	 * @param from_medium Pointer to source medium
	 * @param to_medium Pointer to destination medium
	 *
	 * Calculates reflection and transmission probabilities using Fresnel equations
	 * and handles energy conservation during medium transitions.
	 */
	void handle_medium_transition(Photon& photon, Medium* from_medium, Medium* to_medium);

	/**
	 * @brief Find last surface voxel using DDA traversal
	 * @param photon Reference to photon for path calculation
	 * @param exit_direction Direction vector for exit calculation
	 * @return Voxel* Pointer to last surface voxel or nullptr if not found
	 *
	 * Uses DDA algorithm to efficiently find the last voxel on the surface
	 * along a given direction for accurate exit point calculation.
	 */
	Voxel* find_last_surface_voxel_with_dda(const Photon& photon, const glm::dvec3& exit_direction);

	/**
	 * @brief Terminate photon and record energy for conservation tracking
	 * @param photon Reference to photon to terminate
	 * @param reason Reason for termination (absorption, roulette, etc.)
	 *
	 * Handles proper energy accounting when photons are terminated,
	 * ensuring energy conservation and proper medium energy records.
	 */
	void terminate_photon_and_record_energy(Photon& photon, const std::string& reason);

	/**
	 * @brief Scatter photon using Henyey-Greenberg phase function
	 * @param photon Reference to photon to scatter
	 *
	 * Implements anisotropic scattering with proper coordinate transformations
	 * and numerical stability for extreme anisotropy values.
	 */
	void scatter(Photon& photon);

	/**
	 * @brief Apply Russian Roulette termination to low-weight photons
	 * @param photon Reference to photon for roulette evaluation
	 *
	 * Implements adaptive Russian Roulette with weight-dependent survival
	 * probability for efficient Monte Carlo variance reduction.
	 */
	void roulette(Photon& photon);

	/**
	 * @brief Sample photon step size using Beer-Lambert law
	 * @param photon Reference to photon for step size calculation
	 *
	 * Implements Monte Carlo sampling of photon free path length with
	 * improved numerical stability and precision handling.
	 */
	void step_size(Photon& photon);

	/**
	 * @brief Launch photon from source into simulation medium
	 * @param photon Reference to photon to initialize
	 * @param source Source configuration for photon launch
	 *
	 * Initializes photon properties, handles specular reflection,
	 * and sets up path tracking for Monte Carlo transport.
	 */
	void launch(Photon& photon, const Source& source);

	// ========== Public Delegate Functions ==========

	/**
	 * Track photon path segments for absorption using DDA traversal.
	 * Applies Beer-Lambert law along the actual photon path.
	 */
	void track_photon_path_segments_for_absorption(Photon& photon);

private:
	// Essential delegate methods to Simulator for spatial queries and data access
	Medium* delegate_find_medium_at_with_dda(const glm::dvec3& position);
	Medium* delegate_find_medium_at(const glm::dvec3& position);
	const std::vector<Medium>& delegate_get_mediums() const;

	// ========== Bulk Extracted Methods ==========

	/**
	 * Handle specular reflection at medium interfaces using Fresnel equations.
	 * Calculates reflection coefficients and updates photon direction.
	 */
	void specular_reflection(Photon& photon);

	/**
	 * Robust voxel traversal using 3D DDA algorithm.
	 * Handles energy deposition along voxel boundaries.
	 */
	void track_voxel_path_with_dda(Photon& photon);

	// ========== Internal Helper Functions ==========
	// References to simulation data (not owned)
	std::vector<Medium>* mediums_;                           ///< Reference to mediums vector
	Metrics* metrics_;                                       ///< Reference to metrics object
	const std::vector<std::unique_ptr<DDA>>* dda_instances_; ///< Reference to DDA instances
	std::shared_ptr<Random> rng_;                            ///< Reference to random number generator
	Simulator* simulator_;                                   ///< Reference to main simulator for delegation

	// Transport configuration
	static constexpr double MCML_WEIGHT_THRESHOLD = 1e-4; ///< Minimum weight threshold for photon survival
	static constexpr double STEP_SIZE_TOLERANCE = 1e-12;  ///< Numerical tolerance for step size calculations
	static constexpr double BOUNDARY_TOLERANCE = 1e-10;   ///< Tolerance for boundary detection

	/**
	 * @brief Calculate Fresnel reflection coefficient
	 * @param cos_i Cosine of incident angle
	 * @param eta_i Refractive index of incident medium
	 * @param eta_t Refractive index of transmitted medium
	 * @return double Reflection coefficient (0-1)
	 */
	double calculate_fresnel_reflection(double cos_i, double eta_i, double eta_t) const;

	/**
	 * @brief Calculate specular reflection direction
	 * @param incident Incident direction vector
	 * @param normal Surface normal vector
	 * @return glm::dvec3 Reflected direction vector
	 */
	glm::dvec3 calculate_specular_reflection_direction(const glm::dvec3& incident, const glm::dvec3& normal) const;

	/**
	 * @brief Calculate refraction direction using Snell's law
	 * @param incident Incident direction vector
	 * @param normal Surface normal vector
	 * @param eta Ratio of refractive indices (eta_i/eta_t)
	 * @return glm::dvec3 Refracted direction vector
	 */
	glm::dvec3 calculate_refraction_direction(const glm::dvec3& incident, const glm::dvec3& normal, double eta) const;
};
