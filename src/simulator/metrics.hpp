#pragma once

#include <vector>
#include <chrono>
#include <iostream>

#include <glm/glm.hpp>
#include "config.hpp"

class Metrics
{
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
	void reset(); // Reset all accumulated data
	
	// Normalize raw accumulator values by given factor
	void normalize_raw_values(double factor);
	
	// Reset specific raw accumulator values (used during aggregation)
	void reset_raw_absorption_and_diffuse();

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
		if (Config::get().verbose()) {
			std::cout << "ENERGY DEBUG: add_total_absorption(" << weight << ") to " << total_absorption_ << " = " << (total_absorption_ + weight) << std::endl;
		}
		total_absorption_ += weight; 
	}
	void add_diffuse_reflection(double weight) { 
		if (Config::get().verbose()) {
			std::cout << "ENERGY DEBUG: add_diffuse_reflection(" << weight << ") to " << diffuse_reflection_ << " = " << (diffuse_reflection_ + weight) << std::endl;
		}
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
		if (Config::get().verbose()) {
			std::cout << "ENERGY DEBUG: add_diffuse_transmission(" << weight << ") to " << diffuse_transmission_ << " = " << (diffuse_transmission_ + weight) << std::endl;
		}
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

private:
	// Modern C++ chrono timing
	std::chrono::high_resolution_clock::time_point start_time_;
	double elapsed_time_ms_ {0.0};          // elapsed time in milliseconds

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
};
