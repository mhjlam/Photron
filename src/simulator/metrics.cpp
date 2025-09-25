#include "metrics.hpp"

#include <algorithm>
#include <ranges>
#include <cmath>
#include <cstdint>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <string>

// Include for the unified energy conservation method
#include "simulator.hpp"
#include "config.hpp"
#include <chrono>
#include "../app.hpp" // For output path utilities

void Metrics::increment_scatters() {
	scatter_events_++;
}

void Metrics::add_vertex(double x, double y, double z) {
	path_vertices_.push_back(glm::dvec3(x, y, z));
}

void Metrics::add_step_size(double s) {
	step_sizes_.push_back(s);
}

void Metrics::collect_data(double at, double rs, double rd, double sr, double ts, double td) {
	(void)sr; // Suppress unused parameter warning
	(void)ts; // Suppress unused parameter warning
	
	total_absorption_ = at;
	total_reflection_ = rd;      // Diffuse reflection (energy that entered medium and exited back)
	total_transmission_ = td;    // Diffuse transmission (energy that entered medium and exited forward)
	total_diffusion_ = rd + td;  // Total diffuse emittance (reflection + transmission)

	surface_reflection_ = rs;    // Specular reflection (energy reflected at surface, never entered)
	(void)sr; // Unused parameter - suppress warning
	(void)ts; // Unused parameter - suppress warning
	surface_refraction_ = sr;    // Energy entering medium at surface

	path_length_ = compute_path_length();
	average_step_size_ = compute_average_step_size();
	diffusion_distance_ = compute_diffusion_distance();
}

double Metrics::compute_diffusion_distance() const {
	if (path_vertices_.empty()) {
		return 0.0;
	}

	// Modern C++20: Use ranges::minmax_element for better performance
	const auto [min_x, max_x] = std::ranges::minmax_element(path_vertices_, 
		[](const glm::dvec3& a, const glm::dvec3& b) noexcept { return a.x < b.x; });
	const auto [min_y, max_y] = std::ranges::minmax_element(path_vertices_, 
		[](const glm::dvec3& a, const glm::dvec3& b) noexcept { return a.y < b.y; });
	const auto [min_z, max_z] = std::ranges::minmax_element(path_vertices_, 
		[](const glm::dvec3& a, const glm::dvec3& b) noexcept { return a.z < b.z; });

	const double dx = max_x->x - min_x->x;
	const double dy = max_y->y - min_y->y;
	const double dz = max_z->z - min_z->z;

	return std::sqrt(dx * dx + dy * dy + dz * dz);
}

double Metrics::compute_average_step_size() const {
	if (step_sizes_.empty()) {
		return 0;
	}

	double total_step_size = 0;
	size_t num_steps = step_sizes_.size();

	for (size_t i = 0; i < num_steps; ++i) {
		total_step_size += step_sizes_[i];
	}
	return total_step_size / static_cast<double>(num_steps);
}

double Metrics::compute_path_length() const {
	if (path_vertices_.empty()) {
		return 0;
	}

	double path_len = 0;
	for (std::size_t i = 1; i < path_vertices_.size(); ++i) {
		glm::dvec3 diff = path_vertices_[i] - path_vertices_[i - 1];
		path_len += glm::length(diff);
	}

	return path_len;
}

void Metrics::start_clock() {
	// Modern C++: Cross-platform high-resolution timing
	start_time_ = std::chrono::high_resolution_clock::now();
}

void Metrics::stop_clock() {
	// Calculate elapsed time using modern C++ chrono
	auto end_time = std::chrono::high_resolution_clock::now();
	auto duration = std::chrono::duration_cast<std::chrono::microseconds>(end_time - start_time_);
	elapsed_time_ms_ = duration.count() / 1000.0; // Convert microseconds to milliseconds
}

void Metrics::write_to_file() {
	std::ofstream out(App::get_output_path("experiment.out").c_str(), std::ios_base::app);
	if (!out.good()) {
		return;
	}

	out << "Path length        " << path_length_ << std::endl;
	out << "Scatter events     " << scatter_events_ << std::endl;
	out << "Average step size  " << average_step_size_ << std::endl;
	out << "Diffusion distance " << diffusion_distance_ << std::endl;

	out << "Total absorption   " << total_absorption_ << std::endl;
	out << "Total reflection   " << total_reflection_ << std::endl;
	out << "Total transmission " << total_transmission_ << std::endl;
	out << "Total diffusion    " << total_diffusion_ << std::endl;

	out << "Surface reflection " << surface_reflection_ << std::endl;
	out << "Surface refraction " << surface_refraction_ << std::endl;

	out << "Total time taken   " << elapsed_time_ms_ << " ms" << std::endl;

	out << std::endl;
}

void Metrics::print_report(const class Simulator& simulator) {
	std::cout << std::endl << "=== Monte Carlo Simulation Results ===" << std::endl << std::endl;

	// Per-medium statistics (matching overlay format)
	auto& mediums = simulator.mediums;
	for (size_t i = 0; i < mediums.size(); ++i) {
		auto& medium = mediums[i];
		const auto& metrics = medium.get_metrics();  // const reference
		
		// Medium header with ID
		std::cout << "Medium " << (i + 1) << std::endl;
		
		// Volume Statistics
		std::cout << "Volume Statistics" << std::endl;
		const auto& config = Config::get();
		std::cout << "  Volume Grid:         " << config.nx() << "x" 
		          << config.ny() << "x" << config.nz() << std::endl;
		std::cout << "  Total Voxels:        " << medium.get_volume().size() << std::endl;
		
		// Count surface voxels manually
		size_t surface_count = 0;
		const auto& volume = medium.get_volume();
		for (uint64_t idx = 0; idx < volume.size(); ++idx) {
			auto voxel = volume.at(static_cast<uint32_t>(idx));
			if (voxel && voxel->is_surface_voxel) {
				surface_count++;
			}
		}
		std::cout << "  Surface Voxels:      " << surface_count << std::endl;
		
		std::cout << std::endl;
		
		// Transport Statistics
		std::cout << "Transport Statistics" << std::endl;
		std::cout << "  Total photons:       " << simulator.get_paths().size() << std::endl;
		std::cout << "  Photons entered:     " << metrics.get_photons_entered() << std::endl;
		std::cout << "  Scatter events:      " << std::fixed << std::setprecision(0) 
		          << metrics.get_scatter_events() << std::endl;
		std::cout << "  Total path length:   " << std::fixed << std::setprecision(6) 
		          << metrics.compute_path_length() << std::endl;
		std::cout << "  Average step size:   " << std::fixed << std::setprecision(6) 
		          << metrics.compute_average_step_size() << std::endl;
		std::cout << "  Diffusion distance:  " << std::fixed << std::setprecision(6) 
		          << metrics.compute_diffusion_distance() << std::endl;
		
		std::cout << std::endl;
		
		// Use unified energy conservation calculation (same as overlay)
		auto energy = simulator.calculate_energy_conservation();

		std::cout << "Radiance Properties" << std::endl;
		std::cout << "  Total absorption:    " << std::fixed << std::setprecision(6) 
		          << energy.total_absorption << std::endl;
		std::cout << "  Total diffusion:     " << std::fixed << std::setprecision(6) 
		          << energy.total_diffusion << std::endl;
		std::cout << "    Reflection:        " << std::fixed << std::setprecision(6) 
		          << energy.total_reflection << std::endl;
		std::cout << "    Transmission:      " << std::fixed << std::setprecision(6) 
		          << energy.total_transmission << std::endl;
		
		std::cout << std::endl;
		
		// Energy Conservation
		std::cout << "Energy Conservation" << std::endl;
		if (energy.surface_refraction > 0) {
			// Use TOTAL initial energy as baseline (specular reflection + refracted energy)
			double baseline_energy = energy.surface_reflection + energy.surface_refraction;
			
			// Calculate percentages relative to total initial energy
			double surface_reflection_percent = (energy.surface_reflection / baseline_energy) * 100.0;
			double absorption_percent = (energy.total_absorption / baseline_energy) * 100.0;
			double reflection_percent = (energy.total_reflection / baseline_energy) * 100.0;
			double transmission_percent = (energy.total_transmission / baseline_energy) * 100.0;
			
			// Total should equal baseline_energy for perfect conservation
			// Include surface_reflection in total for accurate energy conservation
			double total_accounted = energy.surface_reflection + energy.total_absorption + 
			                         energy.total_reflection + energy.total_transmission;
			double total_percent = (total_accounted / baseline_energy) * 100.0;

			std::cout << "  Specular reflection: " << std::fixed << std::setprecision(1) 
			          << surface_reflection_percent << "%" << std::endl;
			std::cout << "  Absorption:          " << std::fixed << std::setprecision(1) 
			          << absorption_percent << "%" << std::endl;
			std::cout << "  Reflection:          " << std::fixed << std::setprecision(1) 
			          << reflection_percent << "%" << std::endl;
			std::cout << "  Transmission:        " << std::fixed << std::setprecision(1) 
			          << transmission_percent << "%" << std::endl;
			std::cout << "  Total:               " << std::fixed << std::setprecision(1) 
			          << total_percent << "%" << std::endl;

			// Check if energy is conserved
			if (std::abs(total_percent - 100.0) > 2.0) {
				std::cout << "  WARNING: Energy not conserved!" << std::endl;
			}
		} else {
			std::cout << "  Specular reflection: 0.0%" << std::endl;
			std::cout << "  Absorption:          0.0%" << std::endl;
			std::cout << "  Reflection:          0.0%" << std::endl;
			std::cout << "  Transmission:        0.0%" << std::endl;
			std::cout << "  Total:               0.0%" << std::endl;
		}
		
		std::cout << std::endl;
	}

	std::cout << "Execution time:        " << elapsed_time_ms_ << " ms" << std::endl;
	std::cout << "======================================" << std::endl << std::endl;
}

void Metrics::reset() {
	// Reset timing
	elapsed_time_ms_ = 0.0;

	total_absorption_ = 0.0;
	total_reflection_ = 0.0;
	total_transmission_ = 0.0;
	total_diffusion_ = 0.0;

	surface_reflection_ = 0.0;
	surface_refraction_ = 0.0;

	path_length_ = 0.0;
	scatter_events_ = 0.0;
	diffusion_distance_ = 0.0;
	average_step_size_ = 0.0;

	step_sizes_.clear();
	path_vertices_.clear();

	// Reset all accumulated data
	total_absorption_ = 0.0;
	diffuse_reflection_ = 0.0;
	total_reflection_ = 0.0;
	diffuse_transmission_ = 0.0;
	specular_transmission_ = 0.0;
	total_transmission_ = 0.0;
	total_diffusion_ = 0.0;
	surface_reflection_ = 0.0;
	surface_refraction_ = 0.0;
	path_length_ = 0.0;
	scatter_events_ = 0.0;
	total_steps_ = 0;
	photons_entered_ = 0;
}

void Metrics::normalize_raw_values(double divisor) {
	total_absorption_ /= divisor;
	diffuse_reflection_ /= divisor;
	surface_reflection_ /= divisor;
	surface_refraction_ /= divisor;
	diffuse_transmission_ /= divisor;
	specular_transmission_ /= divisor;
	// Update totals after normalization
	total_reflection_ = diffuse_reflection_ + surface_reflection_;
	total_transmission_ = diffuse_transmission_ + specular_transmission_;
	total_diffusion_ = diffuse_reflection_ + total_transmission_;
	// Note: path_length_, total_steps_, and photons_entered_ are not normalized
}

void Metrics::reset_raw_absorption_and_diffuse() {
	total_absorption_ = 0.0;
	diffuse_reflection_ = 0.0;
	diffuse_transmission_ = 0.0;
	specular_transmission_ = 0.0;
	// Update totals after reset
	total_reflection_ = diffuse_reflection_ + surface_reflection_;
	total_transmission_ = diffuse_transmission_ + specular_transmission_;
	total_diffusion_ = diffuse_reflection_ + total_transmission_;
	// Preserve surface_reflection_ and surface_refraction_ as per comment
}
