/**
 * @file metrics.hpp
 * @brief Comprehensive simulation metrics and energy conservation tracking
 *
 * The Metrics class provides detailed tracking and analysis of Monte Carlo
 * photon transport simulation results, including energy conservation validation,
 * statistical analysis, and performance monitoring.
 */

#pragma once

#include <chrono>
#include <fstream>
#include <iostream>
#include <optional>
#include <string>
#include <vector>

#include <glm/glm.hpp>

#include "common/config.hpp"

// Forward declarations
class Simulator;
class Medium;

/**
 * @class Metrics
 * @brief Comprehensive simulation metrics and energy conservation analysis
 *
 * The Metrics class serves as the central repository for all simulation
 * statistics and provides energy conservation validation. It tracks:
 *
 * - **Energy Statistics**: Absorption, reflection, and transmission measurements
 * - **Performance Metrics**: Execution time, photon processing rates
 * - **Conservation Analysis**: Validates energy conservation laws
 * - **Statistical Summaries**: Mean path lengths, scattering events, etc.
 * - **Progress Tracking**: Real-time simulation progress monitoring
 *
 * The class ensures that the fundamental principle of energy conservation
 * is maintained throughout the simulation by tracking all energy pathways
 * and providing detailed analysis of any discrepancies.
 */
class Metrics
{
public:
	/**
	 * @struct MediumEnergyData
	 * @brief Energy statistics for a single medium layer
	 *
	 * Tracks all forms of energy interaction within a specific medium
	 * including absorption, scattering, and boundary interactions.
	 */
	struct MediumEnergyData
	{
		double total_absorption = 0.0;      ///< Total energy absorbed within the medium
		double specular_reflection = 0.0;   ///< Energy lost to specular reflection at boundaries
		double diffuse_reflection = 0.0;    ///< Energy lost to diffuse reflection at boundaries
		double surface_refraction = 0.0;    ///< Energy transmitted through boundaries via refraction
		double specular_transmission = 0.0; ///< Energy transmitted via specular transmission
		double diffuse_transmission = 0.0;  ///< Energy transmitted via diffuse transmission
		double scatter_events = 0.0;        ///< Number of scattering events within medium
	};

	/**
	 * @struct EnergyConservation
	 * @brief Global energy conservation tracking across all media
	 *
	 * Provides comprehensive energy accounting to validate that the
	 * fundamental principle of energy conservation is maintained.
	 */
	struct EnergyConservation
	{
		double total_absorption = 0.0;   ///< Total energy absorbed across all media
		double total_reflection = 0.0;   ///< Total energy reflected from all surfaces
		double total_transmission = 0.0; ///< Total energy transmitted through all boundaries
		double surface_reflection = 0.0; ///< Energy reflected at medium interfaces
		double surface_refraction = 0.0; ///< Energy refracted at medium interfaces
		double total_energy = 0.0;       ///< Total input energy (should equal sum of outputs)
		double total_diffusion = 0.0;    ///< Total diffusive energy transport
	};

	/**
	 * @struct EnergyConservationPercentages
	 * @brief Percentage-based energy conservation analysis
	 *
	 * Converts absolute energy values to percentages for easy validation
	 * of energy conservation laws and identification of energy leaks.
	 */
	struct EnergyConservationPercentages
	{
		double baseline_energy = 0.0;            ///< Reference energy level (100%)
		double surface_reflection_percent = 0.0; ///< Surface reflection as percentage of input
		double absorption_percent = 0.0;         ///< Absorption as percentage of input
		double reflection_percent = 0.0;         ///< Total reflection as percentage of input
		double transmission_percent = 0.0;       ///< Total transmission as percentage of input
		double total_percent = 0.0;              ///< Sum of all energy pathways as percentage
		double total_accounted = 0.0;            ///< Total accounted energy
		bool is_conserved {true};                ///< True if total_percent is approximately 100%
	};

	/**
	 * @struct EnergyDisplayData
	 * @brief Unified energy data structure for GUI display and export
	 *
	 * Combines conservation data and percentage analysis with caching
	 * support for efficient real-time display updates.
	 */
	struct EnergyDisplayData
	{
		EnergyConservation conservation;           ///< Absolute energy conservation data
		EnergyConservationPercentages percentages; ///< Percentage-based analysis
		bool is_valid {false};                     ///< True if data is current and valid
		uint64_t cached_version {0};               ///< Simulation version when data was cached
	};

public:
	/**
	 * @brief Construct a new Metrics object with default values
	 */
	Metrics() = default;

	/**
	 * @brief Destroy the Metrics object
	 */
	~Metrics() = default;

	/**
	 * @brief Add a vertex position for path length calculation
	 * @param x X coordinate of photon position
	 * @param y Y coordinate of photon position
	 * @param z Z coordinate of photon position
	 */
	void add_vertex(double x, double y, double z);

	/**
	 * @brief Add a step size measurement for average calculation
	 * @param s Step size in millimeters
	 */
	void add_step_size(double s);

	/**
	 * @brief Increment the scatter event counter
	 */
	void increment_scatters();

	/**
	 * @brief Collect energy distribution data from simulation
	 * @param at Absorption coefficient
	 * @param rs Specular reflection coefficient
	 * @param rd Diffuse reflection coefficient
	 * @param sr Surface refraction coefficient
	 * @param ts Specular transmission coefficient
	 * @param td Diffuse transmission coefficient
	 */
	void collect_data(double at, double rs, double rd, double sr, double ts, double td);

	/**
	 * @brief Calculate total diffusion distance from collected vertices
	 * @return double Total path length in millimeters
	 */
	double compute_diffusion_distance() const;

	/**
	 * @brief Calculate average step size from collected measurements
	 * @return double Average step size in millimeters
	 */
	double compute_average_step_size() const;

	/**
	 * @brief Calculate total path length from collected data
	 * @return double Total path length in millimeters
	 */
	double compute_path_length() const;

	/**
	 * @brief Start timing measurement for performance analysis
	 */
	void start_clock();

	/**
	 * @brief Stop timing measurement and record elapsed time
	 */
	void stop_clock();

	/**
	 * @brief Write metrics data to output file
	 */
	void write_to_file();

	/**
	 * @brief Print comprehensive simulation report with energy conservation analysis
	 * @param simulator Reference to simulator for energy data access
	 */
	void print_report(const class Simulator& simulator);

	/**
	 * @brief Reset all metrics to initial state
	 */
	void reset();

	/**
	 * @brief Normalize raw energy values by photon count
	 * @param divisor Number of photons for normalization
	 */
	void normalize_raw_values(double divisor);

	/**
	 * @brief Reset absorption and diffuse reflection accumulators
	 */
	void reset_raw_absorption_and_diffuse();

	/**
	 * @brief Aggregate energy data across all mediums in simulation
	 * @param simulator Reference to simulator for medium access
	 * @return MediumEnergyData Combined energy statistics
	 */
	MediumEnergyData aggregate_medium_energy_data(const Simulator& simulator) const;

	/**
	 * @brief Calculate energy conservation for validation
	 * @param simulator Reference to simulator for energy data access
	 * @return EnergyConservation Energy conservation analysis
	 */
	EnergyConservation calculate_energy_conservation(const Simulator& simulator) const;

	/**
	 * @brief Calculate energy distribution percentages
	 * @param simulator Reference to simulator for energy data access
	 * @return EnergyConservationPercentages Percentage breakdown of energy pathways
	 */
	EnergyConservationPercentages calculate_energy_percentages(const Simulator& simulator) const;

	/**
	 * @brief Get unified energy data for GUI display and export
	 * @param simulator Reference to simulator for energy data access
	 * @return EnergyDisplayData Formatted energy data for visualization
	 */
	EnergyDisplayData get_energy_display_data(const Simulator& simulator) const;

	/**
	 * @brief Validate energy conservation within tolerance
	 * @param simulator Reference to simulator for validation
	 * @param tolerance_percent Acceptable error percentage (default 2.0%)
	 * @return bool True if energy is conserved within tolerance
	 */
	bool is_energy_conserved(const Simulator& simulator, double tolerance_percent = 2.0) const;

	/**
	 * @brief Calculate energy conservation error percentage
	 * @param simulator Reference to simulator for error analysis
	 * @return double Conservation error as percentage
	 */
	double get_conservation_error_percent(const Simulator& simulator) const;

	// Energy reporting methods
	void export_energy_conservation_log(const Simulator& simulator, std::ofstream& ofs) const;
	void print_energy_conservation_console(const Simulator& simulator) const;
	std::string get_energy_summary_text(const Simulator& simulator) const;

	// Results export methods (moved from ResultsExporter)
	void export_results(const Simulator& simulator, bool generate_csv = true) const;
	void export_simulation_log(const Simulator& simulator, std::ofstream& ofs) const;
	void export_medium_statistics(const Simulator& simulator, std::ofstream& ofs) const;
	void export_csv_files(const Simulator& simulator) const;

	// CSV export methods (moved from ResultsExporter)
	void export_absorption_csv(const Simulator& simulator, std::ofstream& ofs) const;
	void export_photon_paths_csv(const Simulator& simulator, std::ofstream& ofs) const;
	void export_emittance_data_csv(const Simulator& simulator, std::ofstream& ofs) const;

	// Direct energy data access - use aggregate_medium_energy_data() for all energy calculations

	// Raw data getters (for internal calculations)
	double get_path_length() const { return path_length_; }
	double get_scatter_events() const { return scatter_events_; }
	double get_average_step_size() const { return average_step_size_; }
	double get_diffusion_distance() const { return diffusion_distance_; }
	double get_total_absorption_raw() const { return total_absorption_; }
	double get_diffuse_reflection_raw() const { return diffuse_reflection_; }
	double get_total_reflection_raw() const { return total_reflection_; }
	double get_diffuse_transmission_raw() const { return diffuse_transmission_; }
	double get_specular_transmission_raw() const { return specular_transmission_; }
	double get_total_transmission_raw() const { return total_transmission_; }
	double get_total_diffusion_raw() const { return total_diffusion_; }
	double get_surface_reflection_raw() const { return surface_reflection_; }
	double get_surface_refraction_raw() const { return surface_refraction_; }
	int get_total_steps() const { return total_steps_; }
	int get_photons_entered() const { return photons_entered_; }

	// Display-normalized getters (for UI) - return per-photon normalized values
	double get_total_absorption() const { return normalized_absorption_; }
	double get_diffuse_reflection() const { return normalized_diffuse_reflection_; }
	double get_total_reflection() const { return normalized_total_reflection_; }
	double get_diffuse_transmission() const { return normalized_diffuse_transmission_; }
	double get_specular_transmission() const { return normalized_specular_transmission_; }
	double get_total_transmission() const { return normalized_total_transmission_; }
	double get_total_diffusion() const { return normalized_total_diffusion_; }
	double get_surface_reflection() const { return normalized_surface_reflection_; }
	double get_surface_refraction() const { return normalized_surface_refraction_; }

	// Accumulator methods for energy tracking
	void add_total_absorption(double weight) { total_absorption_ += weight; }
	void add_diffuse_reflection(double weight) {
		diffuse_reflection_ += weight;
		total_reflection_ = diffuse_reflection_ + surface_reflection_;        // keep absolute total in sync
		total_diffusion_ = diffuse_reflection_ + total_transmission_;         // keep absolute total in sync
	}
	void add_surface_reflection(double weight) {
		surface_reflection_ += weight;
		total_reflection_ = diffuse_reflection_ + surface_reflection_;        // keep absolute total in sync
	}
	void add_surface_refraction(double weight) { surface_refraction_ += weight; }
	void add_diffuse_transmission(double weight) {
		diffuse_transmission_ += weight;
		total_transmission_ = diffuse_transmission_ + specular_transmission_; // keep absolute total in sync
		total_diffusion_ = diffuse_reflection_ + total_transmission_;         // keep absolute total in sync
	}
	void add_specular_transmission(double weight) {
		specular_transmission_ += weight;
		total_transmission_ = diffuse_transmission_ + specular_transmission_; // keep absolute total in sync
		total_diffusion_ = diffuse_reflection_ + total_transmission_;         // keep absolute total in sync
	}
	void set_path_length(double length) { path_length_ = length; }
	void increment_total_steps() { total_steps_++; }
	void increment_photons_entered() { photons_entered_++; }

	// Getters for transport data (for aggregation)
	const std::vector<double>& get_step_sizes() const { return step_sizes_; }
	const std::vector<glm::dvec3>& get_path_vertices() const { return path_vertices_; }

	// Setter for aggregated scatter events
	void set_scatter_events(double scatter_events) { scatter_events_ = scatter_events; }

	// Get elapsed time in milliseconds
	double get_elapsed_time_ms() const { return elapsed_time_ms_; }

	// Set simulation completion data for accurate reporting
	void set_simulation_completion(size_t photons_simulated) {
		simulation_timestamp_ = std::chrono::system_clock::now();
		actual_photons_simulated_ = photons_simulated;
	}

private:
	// Modern C++ chrono timing
	std::chrono::high_resolution_clock::time_point start_time_;
	double elapsed_time_ms_ {0.0};                               // elapsed time in milliseconds
	std::chrono::system_clock::time_point simulation_timestamp_; // when simulation completed
	size_t actual_photons_simulated_ {0};                        // actual number of photons simulated

	// Absolute values (never modified - used for energy conservation)
	double total_absorption_ {0.0};                              // total absorption (absolute)
	double diffuse_reflection_ {0.0};                            // diffuse reflection (absolute)
	double total_reflection_ {0.0};                              // total reflection (absolute)
	double diffuse_transmission_ {0.0};                          // diffuse transmission (absolute)
	double specular_transmission_ {0.0};                         // specular transmission (absolute)
	double total_transmission_ {0.0};                            // total transmission (absolute)
	double surface_reflection_ {0.0};                            // surface reflection (absolute)
	double surface_refraction_ {0.0};                            // surface refraction (absolute)
	double total_diffusion_ {0.0};                               // total diffusion (absolute)
	// Normalized values (for display only)
	double normalized_absorption_ {0.0};                         // normalized absorption for display
	double normalized_diffuse_reflection_ {0.0};                 // normalized diffuse reflection for display
	double normalized_total_reflection_ {0.0};                   // normalized total reflection for display
	double normalized_diffuse_transmission_ {0.0};               // normalized diffuse transmission for display
	double normalized_specular_transmission_ {0.0};              // normalized specular transmission for display
	double normalized_total_transmission_ {0.0};                 // normalized total transmission for display
	double normalized_surface_reflection_ {0.0};                 // normalized surface reflection for display
	double normalized_surface_refraction_ {0.0};                 // normalized surface refraction for display
	double normalized_total_diffusion_ {0.0};                    // normalized total diffusion for display

	double path_length_ {0.0};              // path length of photon migration
	double scatter_events_ {0.0};           // number of scattering events
	int total_steps_ {0};                   // total simulation steps
	int photons_entered_ {0};               // number of photons that entered this medium
	double average_step_size_ {0.0};        // average distance between interaction sites
	double diffusion_distance_ {0.0};       // maximum distance between two points on the path

	std::vector<double> step_sizes_;        // step sizes
	std::vector<glm::dvec3> path_vertices_; // path vertices

	// Caching for expensive calculations
	mutable std::optional<MediumEnergyData> cached_medium_energy_;
	mutable std::optional<EnergyConservation> cached_conservation_;
	mutable std::optional<EnergyConservationPercentages> cached_percentages_;
	mutable uint64_t medium_energy_cache_version_ {0};
	mutable uint64_t conservation_cache_version_ {0};
	mutable uint64_t percentages_cache_version_ {0};

	// Helper methods for energy statistics (moved from EnergyStatisticsManager)
	void write_percentage_line(std::ofstream& ofs, const std::string& label, double percent) const;
	void print_percentage_line(const std::string& label, double percent) const;

	// Helper methods for results export (moved from ResultsExporter)
	std::string get_output_filename(const std::string& base_name, const std::string& extension) const;
	void write_header(std::ofstream& ofs, const std::string& title) const;
};
