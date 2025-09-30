/**
 * @file simulator.hpp
 * @brief Monte Carlo photon transport simulation engine
 * 
 * The Simulator class implements the core Monte Carlo photon transport algorithm
 * based on MCML (Monte Carlo Multi-Layered) methodology for simulating photon
 * migration in voxelized multi-layered biological tissues and materials.
 */

#pragma once

// Standard library includes
#include <functional>
#include <list>
#include <map>
#include <memory>
#include <string>
#include <vector>

// Forward declarations
class Cuboid;
class Triangle;
class DDA;
class Config;
class Logger;
class Layer;
class Photon;
class Material;
class Voxel;  // Make consistent with voxel.hpp
class Volume;
class Medium;
class Random;
struct Source;
struct Emitter;

// Include essential headers for interface
#include "math/range.hpp"  // Needed for Range3 in interface
#include "simulator/metrics.hpp"  // Needed for member variable
#include "common/result.hpp"
#include "common/error_types.hpp"

/**
 * @class Simulator
 * @brief Monte Carlo photon transport simulation engine for multi-layered materials
 * 
 * The Simulator class implements the MCML (Monte Carlo Multi-Layered) algorithm
 * for simulating photon migration through voxelized biological tissues and materials.
 * It supports:
 * - Multi-layered geometry with varying optical properties
 * - Voxel-based spatial discretization for accurate energy deposition tracking
 * - 3D DDA (Digital Differential Analyzer) for robust voxel traversal
 * - Real-time progress tracking and energy conservation monitoring
 * - Multiple light source configurations (point, beam, etc.)
 * - Comprehensive energy statistics (absorption, reflection, transmission)
 * 
 * The simulation follows these key steps for each photon:
 * 1. Launch photon from configured light source
 * 2. Calculate step size based on optical properties
 * 3. Move photon and handle medium transitions
 * 4. Deposit energy in traversed voxels
 * 5. Handle scattering events using Monte Carlo sampling
 * 6. Process interface crossings with Fresnel reflection/transmission
 * 7. Apply Russian roulette for low-weight photons
 * 8. Terminate when photon exits geometry or weight becomes negligible
 */
class Simulator
{
public:
	/**
	 * @brief Construct a new Simulator object with default parameters
	 * 
	 * Initializes empty data structures and sets default MCML parameters.
	 * Call initialize() to load configuration and prepare for simulation.
	 */
	Simulator();
	
	/**
	 * @brief Destroy the Simulator object and clean up resources
	 */
	~Simulator();

	/**
	 * @brief Initialize simulator with configuration file
	 * 
	 * Loads configuration from TOML file, initializes geometry, materials,
	 * light sources, and prepares data structures for simulation.
	 * 
	 * @param file Path to configuration file (TOML format)
	 * @return Result<void, SimulationError> Success or structured error information
	 */
	Result<void, SimulationError> initialize(std::string file);
	
	/**
	 * @brief Run complete Monte Carlo photon transport simulation
	 * 
	 * Executes the main simulation loop, launching the configured number
	 * of photons and tracking their transport through the medium. Updates
	 * progress via callback if configured and aggregates final results.
	 * 
	 * @return Result<void, SimulationError> Success or structured error information
	 */
	Result<void, SimulationError> simulate();
	
	/**
	 * @brief Simulate transport of a single photon for interactive use
	 * 
	 * Useful for step-by-step debugging or real-time visualization.
	 * Adds one photon to the simulation results.
	 */
	void simulate_single_photon();
	
	/**
	 * @brief Aggregate voxel-level energy data into medium-level records
	 * 
	 * Consolidates detailed voxel energy deposition data into summary
	 * statistics for each medium layer for analysis and visualization.
	 */
	void aggregate_voxel_data();
	
	/**
	 * @brief Generate simulation reports and export results
	 * 
	 * Creates comprehensive output including energy conservation analysis,
	 * statistical summaries, and optionally CSV data files.
	 * 
	 * @param generate_csv If true, exports detailed data as CSV files
	 */
	void report(bool generate_csv = true);

public:
	Metrics metrics;                                      ///< Simulation statistics and performance metrics

	std::vector<Photon> photons;                         ///< Complete photon transport paths and interactions
	std::vector<Source> sources;                          ///< Light source definitions and configurations  
	std::vector<std::shared_ptr<Emitter>> emitters;      ///< Energy emission points from photon exits
	std::vector<Medium> mediums;                          ///< Multi-layered material geometry definitions

	std::shared_ptr<Random> rng;                          ///< Random number generator for Monte Carlo sampling
	double mcml_weight_threshold;                         ///< MCML photon weight termination threshold

	std::vector<std::unique_ptr<DDA>> medium_ddas_;       ///< 3D DDA traversal instances (one per medium)

	mutable uint64_t simulation_version_{0};              ///< Version counter for renderer cache invalidation
	
	/**
	 * @brief Get current simulation version for cache invalidation
	 * @return uint64_t Current version counter, increments when data changes
	 */
	uint64_t get_simulation_version() const { return simulation_version_; }
	
	/**
	 * @brief Increment simulation version when data changes (internal use)
	 */
	void increment_simulation_version() const { ++simulation_version_; }

	/**
	 * @brief Get aggregated material definitions from all mediums
	 * @return std::vector<Material> Combined tissue/material properties
	 */
	std::vector<Material> get_all_tissues() const;
	
	/**
	 * @brief Get aggregated layer definitions from all mediums
	 * @return const std::vector<Layer>& Combined layer structure (read-only)
	 */
	const std::vector<Layer>& get_all_layers() const;
	
	/**
	 * @brief Get aggregated voxel data from all mediums (mutable)
	 * @return std::vector<std::shared_ptr<Voxel>>& Combined voxel collection
	 */
	std::vector<std::shared_ptr<Voxel>>& get_all_voxels();
	
	/**
	 * @brief Get aggregated voxel data from all mediums (read-only)
	 * @return const std::vector<std::shared_ptr<Voxel>>& Combined voxel collection
	 */
	const std::vector<std::shared_ptr<Voxel>>& get_all_voxels() const;
	
	/**
	 * @brief Get total number of voxels across all mediums
	 * @return size_t Combined voxel count
	 */
	size_t get_total_voxel_count() const;
	
	/**
	 * @brief Get combined spatial bounds of all mediums
	 * @return Range3 Bounding box encompassing entire simulation geometry
	 */
	Range3 get_combined_bounds() const;
	
	/**
	 * @brief Get combined specular reflection coefficient across all mediums
	 * @return double Aggregated specular reflection value
	 */
	double get_combined_specular_reflection() const; 
	
	/**
	 * @brief Get combined surface refraction coefficient across all mediums
	 * @return double Aggregated surface refraction value
	 */
	double get_combined_surface_refraction() const;
	
	/**
	 * @brief Aggregate energy deposition data from all mediums
	 * @return Metrics::MediumEnergyData Combined energy statistics
	 */
	Metrics::MediumEnergyData aggregate_medium_energy_data() const;
	
	/**
	 * @brief Get unified energy data for display/visualization
	 * @return Metrics::EnergyDisplayData Formatted energy data for UI display
	 */
	Metrics::EnergyDisplayData get_energy_display_data() const;
	
	/**
	 * @brief Get shared metrics instance for external components
	 * @return std::shared_ptr<Metrics> Shared metrics object
	 */
	std::shared_ptr<Metrics> get_shared_metrics() const { return shared_metrics_; }
	
	// Shared metrics instance for energy statistics
	std::shared_ptr<Metrics> shared_metrics_;
	
	/**
	 * @brief Set shared metrics instance for cross-component data sharing
	 * @param shared_metrics Shared metrics object to use
	 */
	void set_shared_metrics(std::shared_ptr<Metrics> shared_metrics) { shared_metrics_ = shared_metrics; }
	
	/**
	 * @brief Check if 3D position is inside simulation geometry
	 * @param position 3D world position to test
	 * @return bool True if position is within any medium boundary
	 */
	bool is_point_inside_geometry(const glm::dvec3& position) const;
	
	/**
	 * @brief Access voxel by discrete grid coordinates
	 * @param x Grid X coordinate
	 * @param y Grid Y coordinate  
	 * @param z Grid Z coordinate
	 * @return Voxel* Pointer to voxel at coordinates, nullptr if invalid
	 */
	Voxel* voxel_grid(uint32_t x, uint32_t y, uint32_t z) const;

	// Path access for backward compatibility (creates PhotonPaths from photon data)
	// Removed get_paths() - use photons member directly
	
	/**
	 * @brief Set progress callback function for UI updates during simulation
	 * @param callback Function receiving (current_photon, total_photons) progress
	 */
	void set_progress_callback(std::function<void(uint64_t, uint64_t)> callback) {
		progress_callback_ = callback;
	}

private:
	/**
	 * @brief Initialize photon launch from source with proper direction sampling
	 * 
	 * Sets initial position, direction, and weight based on source configuration.
	 * Handles different source types (pencil beam, Gaussian, etc.) with appropriate
	 * spatial and angular distributions.
	 * 
	 * @param photon Reference to photon being launched
	 * @param source Source configuration defining launch parameters
	 */
	void launch(Photon& photon, const Source& source);

	/**
	 * @brief Calculate photon step size based on optical properties and random sampling
	 * 
	 * Uses Beer-Lambert law with exponential sampling to determine distance
	 * to next interaction event. Accounts for varying optical properties
	 * along the photon path through different materials.
	 * 
	 * @param photon Reference to photon for step size calculation
	 */
	void step_size(Photon& photon);

	/**
	 * @brief Move photon along its trajectory by the calculated step size
	 * 
	 * Updates photon position and handles boundary crossings between
	 * different media. Ensures energy conservation during transport.
	 * 
	 * @param photon Reference to photon being transported
	 */
	void transfer(Photon& photon);

	/**
	 * @brief Handle partial step when photon crosses medium boundaries
	 * 
	 * Adjusts step size when photon path intersects material interfaces,
	 * ensuring accurate position tracking and energy conservation.
	 * 
	 * @param photon Reference to photon crossing boundaries
	 */
	void sub_step(Photon& photon);

	/**
	 * @brief Deposit photon energy in current voxel based on absorption
	 * 
	 * Calculates energy deposition using absorption coefficient and
	 * updates voxel energy statistics for visualization and analysis.
	 * 
	 * @param photon Reference to photon depositing energy
	 */
	void deposit(Photon& photon);

	/**
	 * @brief Track photon path through voxels and deposit energy using 3D DDA
	 * 
	 * Uses 3D Digital Differential Analyzer for efficient voxel traversal
	 * and accurate energy deposition along the photon path.
	 * 
	 * @param photon Reference to photon being tracked
	 */
	void track_voxel_path_and_deposit(Photon& photon);

	/**
	 * @brief Perform Henyey-Greenberg scattering to change photon direction
	 * 
	 * Implements anisotropic scattering using HG phase function with
	 * proper cosine and azimuthal angle sampling for realistic light transport.
	 * 
	 * @param photon Reference to photon being scattered
	 */
	void scatter(Photon& photon);

	/**
	 * @brief Handle photon crossing between different optical media
	 * 
	 * Processes Fresnel reflection/transmission at interfaces using
	 * Snell's law and proper refractive index handling.
	 * 
	 * @param photon Reference to photon at medium interface
	 */
	void cross(Photon& photon);

	/**
	 * @brief Record photon energy exit from simulation volume
	 * 
	 * Classifies exit type (diffuse/specular reflection/transmission) and
	 * records energy for boundary condition analysis.
	 * 
	 * @param photon Reference to exiting photon
	 * @param direction Exit direction vector
	 * @param weight Energy weight being recorded
	 */
	void radiate(Photon& photon, glm::dvec3& direction, double weight);

	/**
	 * @brief Find last surface voxel using DDA traversal for accurate exit recording
	 * 
	 * Traces backward from exit point to find the last voxel containing
	 * geometry for proper energy classification and recording.
	 * 
	 * @param photon Photon that has exited the geometry
	 * @param exit_direction Direction of photon exit
	 * @return Pointer to last surface voxel, or nullptr if not found
	 */
	Voxel* find_last_surface_voxel_with_dda(const Photon& photon, const glm::dvec3& exit_direction);

	/**
	 * @brief Apply Russian roulette termination for low-weight photons
	 * 
	 * Stochastically terminates photons with low energy to improve
	 * computational efficiency while maintaining unbiased results.
	 * 
	 * @param photon Reference to photon for roulette evaluation
	 */
	void roulette(Photon& photon);

	/**
	 * @brief Normalize all accumulated energy data by total photon count
	 * 
	 * Converts raw energy accumulation to proper Monte Carlo estimates
	 * by normalizing with the number of simulated photons.
	 */
	void normalize();
	
	/**
	 * @brief Terminate photon simulation and record final energy state
	 * 
	 * Handles photon termination due to various conditions (absorption,
	 * boundary exit, roulette) and ensures proper energy accounting.
	 * 
	 * @param photon Reference to photon being terminated
	 * @param reason String description of termination cause for debugging
	 */
	void terminate_photon_and_record_energy(Photon& photon, const std::string& reason);

	/**
	 * @brief Handle specular reflection at smooth optical interfaces
	 * 
	 * Computes reflection direction using geometric optics and updates
	 * photon trajectory for mirror-like surface interactions.
	 * 
	 * @param photon Reference to photon undergoing specular reflection
	 */
	void specular_reflection(Photon& photon);

	/**
	 * @brief Calculate Fresnel reflection/transmission coefficients at interfaces
	 * 
	 * Uses Fresnel equations to compute reflection probability based on
	 * incident angle and refractive indices. Handles total internal reflection.
	 * 
	 * @param photon Reference to photon at interface
	 * @param eta_t Transmission medium refractive index ratio
	 * @param tran Computed transmission direction (output)
	 * @param refl Computed reflection direction (output)
	 * @return Fresnel reflection coefficient [0,1]
	 */
	double internal_reflection(Photon& photon, double& eta_t, glm::dvec3& tran, glm::dvec3& refl);

	/**
	 * @brief Determine if photon undergoes reflection at current interface
	 * 
	 * Uses Fresnel coefficients and random sampling to decide between
	 * reflection and transmission at optical interfaces.
	 * 
	 * @param photon Photon at interface for reflection test
	 * @return true if photon reflects, false if it transmits
	 */
	bool is_photon_reflecting(const Photon& photon) const;

	/**
	 * @brief Move point along direction vector by specified distance
	 * 
	 * Performs 3D vector translation for photon position updates.
	 * Updates position in-place for efficiency.
	 * 
	 * @param position Current position (modified in-place)
	 * @param direction Normalized direction vector
	 * @param d Distance to move along direction
	 * @return New position after movement
	 */
	glm::dvec3 move(glm::dvec3& position, glm::dvec3& direction, double d);

	/**
	 * @brief Move point by remaining step size in current direction
	 * 
	 * Convenience method that moves by the photon's current step size.
	 * Used for boundary crossing and final position updates.
	 * 
	 * @param position Current position (modified in-place)
	 * @param direction Normalized direction vector
	 * @return New position after delta movement
	 */
	glm::dvec3 move_delta(glm::dvec3& position, glm::dvec3& direction);

	/**
	 * @brief Parse TOML configuration file with structured error handling
	 * 
	 * Loads simulation parameters from TOML file including materials,
	 * geometry, sources, and simulation settings with comprehensive
	 * error reporting for invalid configurations.
	 * 
	 * @param fconfig Path to TOML configuration file
	 * @return Result containing success or detailed parsing error
	 */
	Result<void, SimulationError> parse(const std::string& fconfig);

	/**
	 * @brief Initialize photon sources from parsed configuration
	 * 
	 * Validates and sets up photon source parameters including position,
	 * direction, and emission characteristics with error checking.
	 * 
	 * @return Result containing success or initialization error details
	 */
	Result<void, SimulationError> initialize_sources();
	
	/**
	 * @brief Reset all simulation data for new simulation run
	 * 
	 * Clears accumulated energy data, resets voxel states, and initializes
	 * counters for a fresh simulation while preserving configuration.
	 */
	void reset_simulation_data();

	/**
	 * @brief Find the medium containing a given 3D position
	 * 
	 * Performs spatial query to determine which medium (if any) contains
	 * the specified point using geometry intersection tests.
	 * 
	 * @param position 3D world coordinate to test
	 * @return Pointer to containing medium, or nullptr if outside all media
	 */
	Medium* find_medium_at(const glm::dvec3& position) const;

	/**
	 * @brief Check if position is inside any simulation medium
	 * 
	 * Fast test to determine if a point is within the simulation domain
	 * without identifying the specific medium.
	 * 
	 * @param position 3D world coordinate to test
	 * @return true if inside any medium, false if in ambient space
	 */
	bool is_inside_any_medium(const glm::dvec3& position) const;

	/**
	 * @brief Handle photon transition between different optical media
	 * 
	 * Manages photon state updates when crossing medium boundaries,
	 * including optical property changes and interface physics.
	 * 
	 * @param photon Reference to photon crossing boundary
	 * @param from Source medium being exited (may be nullptr)
	 * @param to Destination medium being entered (may be nullptr)
	 */
	void handle_medium_transition(Photon& photon, Medium* from, Medium* to);

	/**
	 * @brief Validate photon state consistency after interface crossing
	 * 
	 * Debug validation to ensure proper photon state after medium
	 * transitions, checking for energy conservation and valid parameters.
	 * 
	 * @param photon Reference to photon after transition
	 * @param from_medium Previous medium (may be nullptr)
	 * @param to_medium Current medium (may be nullptr)
	 */
	void validate_photon_state_after_interface_transition(Photon& photon, Medium* from_medium, Medium* to_medium);
	
	/**
	 * @brief Initialize 3D DDA instances for efficient voxel traversal
	 * 
	 * Sets up Digital Differential Analyzer objects for fast voxel grid
	 * traversal during photon path tracking and energy deposition.
	 */
	void initialize_dda_instances();

	/**
	 * @brief Track photon path through voxel grid using 3D DDA algorithm
	 * 
	 * Efficiently traverses voxel grid along photon trajectory using
	 * Digital Differential Analyzer for accurate path tracking.
	 * 
	 * @param photon Reference to photon being tracked through voxels
	 */
	void track_voxel_path_with_dda(Photon& photon);

	/**
	 * @brief Track photon path segments for accurate absorption calculation
	 * 
	 * Computes absorption along photon path by tracking segments through
	 * different materials with varying optical properties.
	 * 
	 * @param photon Reference to photon for path segment analysis
	 */
	void track_photon_path_segments_for_absorption(Photon& photon);

	/**
	 * @brief Find medium at position using DDA-optimized spatial queries
	 * 
	 * Uses DDA acceleration structures for fast medium lookup at
	 * specified 3D coordinates with improved performance.
	 * 
	 * @param position 3D world coordinate for medium query
	 * @return Pointer to containing medium, or nullptr if outside geometry
	 */
	Medium* find_medium_at_with_dda(const glm::dvec3& position) const;
	
	/**
	 * @brief Find voxel at specified 3D position across all media
	 * 
	 * Delegates to appropriate medium to find voxel containing the
	 * given position, handling multi-medium scenarios efficiently.
	 * 
	 * @param position 3D world coordinate for voxel lookup
	 * @return Pointer to containing voxel, or nullptr if outside voxelized regions
	 */
	Voxel* voxel_at(const glm::dvec3& position) const;

	/**
	 * @brief Compute 3D bounding box corners for given voxel
	 * 
	 * Calculates world-space corner coordinates of voxel for geometric
	 * operations and visualization. Delegates to owning medium.
	 * 
	 * @param voxel Pointer to voxel for corner calculation
	 * @return Cuboid containing min and max corner coordinates
	 */
	Cuboid voxel_corners(Voxel* voxel) const;
	
	std::function<void(uint64_t, uint64_t)> progress_callback_;  ///< Callback for simulation progress reporting
};
