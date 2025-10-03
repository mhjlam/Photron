/**
 * @file simulator.hpp
 * @brief Monte Carlo photon transport simulation engine
 *
 * The Simulator class impleme core Monte Carlo photon transport algorithm
 * based on MCML (Monte Carlo Multi-Layered) methodology for simulating photon
 * migration in voxelized multi-layered biological tissues and materials.
 */

#pragma once

#include <functional>
#include <list>
#include <map>
#include <memory>
#include <string>
#include <vector>

#include "common/error_types.hpp"
#include "common/result.hpp"
#include "math/dda.hpp"
#include "math/range.hpp"
#include "simulator/aggregator.hpp"
#include "simulator/geom_lookup.hpp"
#include "simulator/medium.hpp"
#include "simulator/metrics.hpp"
#include "simulator/photon.hpp"
#include "simulator/transport.hpp"

// Forward declarations
class Cuboid;
class Layer;
class Material;
class Voxel;
class Random;

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
	friend class Transport; // Allow Transport to access private methods for delegation

public:
	/**
	 * @brief Construct a new Simulator object with default parameters
	 *
	 * Initializes empty data structures and sets default MCML parameters.
	 * Call initialize() to load configuration and prepare for simulation.
	 */
	Simulator();

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
	Metrics metrics;                                ///< Simulation statistics and performance metrics

	std::vector<Photon> photons;                    ///< Complete photon transport paths and interactions
	std::vector<Source> sources;                    ///< Light source definitions and configurations
	std::vector<std::shared_ptr<Emitter>> emitters; ///< Energy emission points from photon exits
	std::vector<Medium> mediums;                    ///< Multi-layered material geometry definitions

	std::shared_ptr<Random> rng;                    ///< Random number generator for Monte Carlo sampling
	double mcml_weight_threshold;                   ///< MCML photon weight termination threshold

	std::vector<std::unique_ptr<DDA>> medium_ddas_; ///< 3D DDA traversal instances (one per medium)

	mutable uint64_t simulation_version_ {0};       ///< Version counter for renderer cache invalidation

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
	 * @brief Access to mediums for Transport delegation
	 * @return const std::vector<Medium>& Reference to mediums vector
	 */
	const std::vector<Medium>& get_mediums() const { return mediums; }

	/**
	 * @brief Access to emitters for Transport delegation
	 * @return std::vector<std::shared_ptr<Emitter>>& Reference to emitters vector
	 */
	std::vector<std::shared_ptr<Emitter>>& get_emitters() { return emitters; }

	/**
	 * @brief Access to DDA instances for Transport delegation
	 * @return const std::vector<std::unique_ptr<DDA>>& Reference to DDA instances
	 */
	const std::vector<std::unique_ptr<DDA>>& get_medium_ddas() const { return medium_ddas_; }

	/**
	 * @brief Get shared metrics instance for external components
	 * @return std::shared_ptr<Metrics> Shared metrics object
	 */
	std::shared_ptr<Metrics> get_shared_metrics() const { return shared_metrics_; }

	/**
	 * @brief Get direct access to metrics instance
	 * @return const Metrics& Reference to metrics for direct access
	 */
	const Metrics& get_metrics() const { return metrics; }

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
	void set_progress_callback(std::function<void(uint64_t, uint64_t)> callback) { progress_callback_ = callback; }

private:
	/**
	 * @brief Initialize 3D DDA instances for voxel traversal
	 *
	 * Creates DDA (Digital Differential Analyzer) objects for each medium
	 * to enable efficient ray-voxel intersection calculations.
	 */
	void initialize_dda_instances();

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
	 * @brief Initialize component instances for modular functionality
	 *
	 * Creates instances of GeomLookup and Aggregator with proper
	 * references to simulation data. Called during initialization after
	 * mediums and DDA instances are fully configured.
	 */
	void initialize_components();

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

	std::function<void(uint64_t, uint64_t)> progress_callback_; ///< Callback for simulation progress reporting

	// Component instances for modular functionality
	std::unique_ptr<GeomLookup> geom_lookup_;     ///< Handles spatial queries and geometry operations
	std::unique_ptr<Aggregator> aggregator_;      ///< Handles results aggregation and export
	std::unique_ptr<Transport> photon_transport_; ///< Handles detailed photon transport physics
};
