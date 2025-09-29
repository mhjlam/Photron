#pragma once

#include <vector>
#include <chrono>
#include <iostream>
#include <string>
#include <fstream>
#include <optional>

#include <glm/glm.hpp>
#include "config.hpp"

// Forward declarations
class Simulator;
class Medium;

class Metrics
{
public:
	// Energy data structures
	struct MediumEnergyData {
		double total_absorption = 0.0;
		double specular_reflection = 0.0;
		double diffuse_reflection = 0.0;
		double surface_refraction = 0.0;
		double specular_transmission = 0.0;
		double diffuse_transmission = 0.0;
		double scatter_events = 0.0;
	};

	struct EnergyConservation {
		double total_absorption = 0.0;
		double total_reflection = 0.0; 
		double total_transmission = 0.0;
		double surface_reflection = 0.0;
		double surface_refraction = 0.0;
		double total_energy = 0.0;
		double total_diffusion = 0.0;
	};

	struct EnergyConservationPercentages {
		double baseline_energy = 0.0;
		double surface_reflection_percent = 0.0;
		double absorption_percent = 0.0;
		double reflection_percent = 0.0;
		double transmission_percent = 0.0;
		double total_percent = 0.0;
		double total_accounted = 0.0;
		bool is_conserved = true; // True if total_percent is close to 100%
	};

	// Unified energy display data structure for GUI and export
	struct EnergyDisplayData {
		EnergyConservation conservation;
		EnergyConservationPercentages percentages;
		bool is_valid = false;
		uint64_t cached_version = 0;  // Track when data was cached
	};

public:
	Metrics() = default;
	~Metrics() = default;

	// Data collection methods
	void add_vertex(double x, double y, double z);
	void add_step_size(double s);
	void increment_scatters();

	void collect_data(double at, double rs, double rd, double sr, double ts, double td);
	double compute_diffusion_distance() const;
	double compute_average_step_size() const;
	double compute_path_length() const;

	void start_clock();
	void stop_clock();

	void write_to_file();
	void print_report(const class Simulator& simulator); // Unified energy conservation reporting
	void reset();
	void normalize_raw_values(double divisor);
	void reset_raw_absorption_and_diffuse();

	// Energy statistics methods (moved from EnergyStatisticsManager)
	MediumEnergyData aggregate_medium_energy_data(const Simulator& simulator) const;
	EnergyConservation calculate_energy_conservation(const Simulator& simulator) const;
	EnergyConservationPercentages calculate_energy_percentages(const Simulator& simulator) const;
	
	// Unified energy display data (single call for GUI and export)
	EnergyDisplayData get_energy_display_data(const Simulator& simulator) const;
	
	// Energy validation methods
	bool is_energy_conserved(const Simulator& simulator, double tolerance_percent = 2.0) const;
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

	// Getters for UI display
	double get_path_length() const { return path_length_; }
	double get_scatter_events() const { return scatter_events_; }
	double get_average_step_size() const { return average_step_size_; }
	double get_diffusion_distance() const { return diffusion_distance_; }
	double get_total_absorption() const { return total_absorption_; }
	double get_diffuse_reflection() const { return diffuse_reflection_; }
	double get_total_reflection() const { return total_reflection_; }
	double get_diffuse_transmission() const { return diffuse_transmission_; }
	double get_specular_transmission() const { return specular_transmission_; }
	double get_total_transmission() const { return total_transmission_; }
	double get_total_diffusion() const { return total_diffusion_; }
	double get_surface_reflection() const { return surface_reflection_; }
	double get_surface_refraction() const { return surface_refraction_; }
	int get_total_steps() const { return total_steps_; }
	int get_photons_entered() const { return photons_entered_; }

	// Accumulator methods for energy tracking
	void add_total_absorption(double weight) { 
		total_absorption_ += weight; 
	}
	void add_diffuse_reflection(double weight) { 
		diffuse_reflection_ += weight; 
		total_reflection_ = diffuse_reflection_ + surface_reflection_; // keep total in sync
		total_diffusion_ = diffuse_reflection_ + total_transmission_; // keep total in sync
	}
	void add_surface_reflection(double weight) { 
		surface_reflection_ += weight;
		total_reflection_ = diffuse_reflection_ + surface_reflection_; // keep total in sync
	}
	void add_surface_refraction(double weight) { surface_refraction_ += weight; }
	void add_diffuse_transmission(double weight) { 
		diffuse_transmission_ += weight; 
		total_transmission_ = diffuse_transmission_ + specular_transmission_; // keep total in sync
		total_diffusion_ = diffuse_reflection_ + total_transmission_; // keep total in sync
	}
	void add_specular_transmission(double weight) { 
		specular_transmission_ += weight;
		total_transmission_ = diffuse_transmission_ + specular_transmission_; // keep total in sync
		total_diffusion_ = diffuse_reflection_ + total_transmission_; // keep total in sync
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
	double elapsed_time_ms_ {0.0};          // elapsed time in milliseconds
	std::chrono::system_clock::time_point simulation_timestamp_; // when simulation completed
	size_t actual_photons_simulated_ {0};   // actual number of photons simulated

	double total_absorption_ {0.0};         // total absorption
	double diffuse_reflection_ {0.0};       // diffuse reflection (separate from specular)
	double total_reflection_ {0.0};         // total reflection (diffuse + specular) 
	double diffuse_transmission_ {0.0};     // diffuse transmission (separate from specular)
	double specular_transmission_ {0.0};    // specular transmission (separate from diffuse)
	double total_transmission_ {0.0};       // total transmission (diffuse + specular)
	double total_diffusion_ {0.0};          // total diffusion (diffuse reflection + total transmission)

	double surface_reflection_ {0.0};       // amount of weight that reflects at the surface (specular reflection)
	double surface_refraction_ {0.0};       // amount of weight that refracts at the surface

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
	mutable uint64_t medium_energy_cache_version_ = 0;
	mutable uint64_t conservation_cache_version_ = 0;
	mutable uint64_t percentages_cache_version_ = 0;
	
	// Helper methods for energy statistics (moved from EnergyStatisticsManager)
	void write_percentage_line(std::ofstream& ofs, const std::string& label, double percent) const;
	void print_percentage_line(const std::string& label, double percent) const;
	
	// Helper methods for results export (moved from ResultsExporter)
	std::string get_output_filename(const std::string& base_name, const std::string& extension) const;
	void write_header(std::ofstream& ofs, const std::string& title) const;
};
