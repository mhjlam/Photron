/**
 * @file metrics.cpp
 * @brief Implementation of simulation metrics collection and analysis
 * 
 * Implements the Metrics class which tracks simulation performance, energy
 * conservation, path statistics, and provides export functionality for
 * results analysis. Features modern C++20 algorithms and efficient
 * statistical calculations.
 */

#include "metrics.hpp"

// Add includes for complete type definitions
#include "simulator/photon.hpp"
#include "simulator/voxel.hpp"

#include "../common/file_utils.hpp"
#include "../common/result.hpp"
#include "../common/error_types.hpp"

#include <algorithm>
#include <ranges>
#include <cmath>
#include <cstdint>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <string>
#include <sstream>

// Include for the unified energy conservation method
#include "simulator.hpp"
#include "config.hpp"
#include "logger.hpp"
#include "medium.hpp"
#include <chrono>
#include "../app.hpp" // For output path utilities

void Metrics::increment_scatters() {
	// Track scattering event statistics
	scatter_events_++;
}

void Metrics::add_vertex(double x, double y, double z) {
	// Collect photon path vertices for analysis
	path_vertices_.push_back(glm::dvec3(x, y, z));
}

void Metrics::add_step_size(double s) {
	// Collect step size data for statistical analysis
	step_sizes_.push_back(s);
}

void Metrics::collect_data(double at, double rs, double rd, double sr, double ts, double td) {
	// Collect energy conservation data from simulation
	(void)sr; // Suppress unused parameter warning
	(void)ts; // Suppress unused parameter warning
	
	// Store primary energy components
	total_absorption_ = at;
	total_reflection_ = rd;      // Diffuse reflection (energy that entered medium and exited back)
	total_transmission_ = td;    // Diffuse transmission (energy that entered medium and exited forward)
	total_diffusion_ = rd + td;  // Total diffuse emittance (reflection + transmission)

	// Store surface interaction components
	surface_reflection_ = rs;    // Specular reflection (energy reflected at surface, never entered)
	(void)sr; // Unused parameter - suppress warning
	(void)ts; // Unused parameter - suppress warning
	surface_refraction_ = sr;    // Energy entering medium at surface

	// Calculate derived metrics from collected data
	path_length_ = compute_path_length();
	average_step_size_ = compute_average_step_size();
	diffusion_distance_ = compute_diffusion_distance();
}

double Metrics::compute_diffusion_distance() const {
	// Calculate spatial extent of photon diffusion
	if (path_vertices_.empty()) {
		return 0.0;
	}

	// Use modern C++20 ranges for efficient min/max calculation
	const auto [min_x, max_x] = std::ranges::minmax_element(path_vertices_, 
		[](const glm::dvec3& a, const glm::dvec3& b) noexcept { return a.x < b.x; });
	const auto [min_y, max_y] = std::ranges::minmax_element(path_vertices_, 
		[](const glm::dvec3& a, const glm::dvec3& b) noexcept { return a.y < b.y; });
	const auto [min_z, max_z] = std::ranges::minmax_element(path_vertices_, 
		[](const glm::dvec3& a, const glm::dvec3& b) noexcept { return a.z < b.z; });

	// Calculate 3D bounding box diagonal as diffusion distance
	const double dx = max_x->x - min_x->x;
	const double dy = max_y->y - min_y->y;
	const double dz = max_z->z - min_z->z;

	return std::sqrt(dx * dx + dy * dy + dz * dz);
}

double Metrics::compute_average_step_size() const {
	// Calculate mean step size from collected samples
	if (step_sizes_.empty()) {
		return 0;
	}

	// Use accumulation for numerical stability
	double total_step_size = 0;
	size_t num_steps = step_sizes_.size();

	for (size_t i = 0; i < num_steps; ++i) {
		total_step_size += step_sizes_[i];
	}
	return total_step_size / static_cast<double>(num_steps);
}

double Metrics::compute_path_length() const {
	// Calculate total photon path length from vertex sequence
	if (path_vertices_.empty()) {
		return 0;
	}

	// Sum distances between consecutive vertices
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

void Metrics::print_report(const class Simulator& simulator) {
	// Always show simulation results header - this is core output users expect
	std::cout << std::endl << "=== Monte Carlo Simulation Results ===" << std::endl << std::endl;
	
	if (Config::get().log()) {
		Logger::instance().log_info("=== Monte Carlo Simulation Results ===");
	}

	// General Simulation Statistics (not per-medium)
	const auto& config = Config::get();
    std::cout << "Total photons:         " << simulator.photons.size() << std::endl;

	std::cout << "Volume Grid:           " << config.nx() << "x" << config.ny() << "x" << config.nz() << std::endl;
	std::cout << "Voxel Size:            " << config.vox_size() << std::endl;
	
	// Calculate overall voxel statistics (using first medium)
	auto& mediums = simulator.mediums;
	if (!mediums.empty()) {
		const auto& volume = mediums[0].get_volume();
		size_t total_voxels = volume.size();
		size_t material_voxels = 0;
		
		for (uint64_t idx = 0; idx < volume.size(); ++idx) {
			auto voxel = volume.at(static_cast<uint32_t>(idx));
			if (voxel && voxel->material) {
				material_voxels++;
			}
		}
		
		std::cout << "Total Voxels:          " << total_voxels << std::endl;
		std::cout << "Empty Voxels:          " << (total_voxels - material_voxels) << std::endl;
	}
	
	std::cout << std::endl;

	// Per-medium statistics (transport only)
	for (size_t i = 0; i < mediums.size(); ++i) {
		auto& medium = mediums[i];
		const auto& metrics = medium.get_metrics();  // const reference
		
		// Medium header with ID
		std::cout << "Medium " << (i + 1) << std::endl;
		
		// Medium Voxel Statistics
		const auto& volume = medium.get_volume();
		size_t material_voxels = 0;
		size_t surface_voxels = 0;
		
		for (uint64_t idx = 0; idx < volume.size(); ++idx) {
			auto voxel = volume.at(static_cast<uint32_t>(idx));
			if (voxel && voxel->material) {
				material_voxels++;
				if (voxel->is_surface_voxel) {
					surface_voxels++;
				}
			}
		}
		
        // Volume Statistics
        std::cout << "Volume Statistics" << std::endl;
		std::cout << "  Medium Voxels:       " << material_voxels << std::endl;
		std::cout << "  Surface Voxels:      " << surface_voxels << std::endl;
		std::cout << std::endl;
		
		// Transport Statistics
		std::cout << "Transport Statistics" << std::endl;
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
		
		// Use unified energy display data (single call for both conservation and percentages)
		auto energy_data = get_energy_display_data(simulator);
		auto& energy = energy_data.conservation;
		auto& percentages = energy_data.percentages;

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
		
		if (energy_data.is_valid) {
			std::cout << "  Specular reflection: " << std::fixed << std::setprecision(1) 
			          << percentages.surface_reflection_percent << "%" << std::endl;
			std::cout << "  Absorption:          " << std::fixed << std::setprecision(1) 
			          << percentages.absorption_percent << "%" << std::endl;
			std::cout << "  Reflection:          " << std::fixed << std::setprecision(1) 
			          << percentages.reflection_percent << "%" << std::endl;
			std::cout << "  Transmission:        " << std::fixed << std::setprecision(1) 
			          << percentages.transmission_percent << "%" << std::endl;
			std::cout << "  Total:               " << std::fixed << std::setprecision(1) 
			          << percentages.total_percent << "%" << std::endl;

			// Check if energy is conserved
			if (!percentages.is_conserved) {
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
	
	// Always log summary to file (regardless of console output setting)
	if (Config::get().log()) {
		Logger::instance().log_info("Simulation completed successfully");
		Logger::instance().log_info("Execution time: " + std::to_string(elapsed_time_ms_) + " ms");
		Logger::instance().log_info("Total photons processed: " + std::to_string(simulator.photons.size()));
		Logger::instance().log_info("======================================");
	}
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
	
	// Invalidate cache
	cached_medium_energy_.reset();
	cached_conservation_.reset();
	cached_percentages_.reset();
	medium_energy_cache_version_ = 0;
	conservation_cache_version_ = 0;
	percentages_cache_version_ = 0;
}

void Metrics::normalize_raw_values(double divisor) {
	total_absorption_ /= divisor;
	diffuse_reflection_ /= divisor;
	surface_reflection_ /= divisor;
	surface_refraction_ /= divisor;
	diffuse_transmission_ /= divisor;
	specular_transmission_ /= divisor;
	// Recalculate totals after normalization
	total_reflection_ = diffuse_reflection_ + surface_reflection_;
	total_transmission_ = diffuse_transmission_ + specular_transmission_;
	total_diffusion_ = diffuse_reflection_ + total_transmission_;
	// path_length_, total_steps_, and photons_entered_ are not normalized
}

void Metrics::reset_raw_absorption_and_diffuse() {
	total_absorption_ = 0.0;
	diffuse_reflection_ = 0.0;
	diffuse_transmission_ = 0.0;
	specular_transmission_ = 0.0;
	// Recalculate totals after reset
	total_reflection_ = diffuse_reflection_ + surface_reflection_;
	total_transmission_ = diffuse_transmission_ + specular_transmission_;
	total_diffusion_ = diffuse_reflection_ + total_transmission_;
	// Preserve surface_reflection_ and surface_refraction_ as per comment
}

/***********************************************************
 * ENERGY STATISTICS METHODS (moved from EnergyStatisticsManager)
 * Consolidates all energy conservation functionality in Metrics
 ***********************************************************/

/***********************************************************
 * AGGREGATE MEDIUM ENERGY DATA
 * Consolidates energy data from all media
 ***********************************************************/
Metrics::MediumEnergyData Metrics::aggregate_medium_energy_data(const Simulator& simulator) const {
    // Check cache first - use simulation version as the version counter
    const uint64_t current_version = simulator.get_simulation_version();
    if (cached_medium_energy_ && medium_energy_cache_version_ == current_version) {
        return *cached_medium_energy_;
    }
    
    MediumEnergyData data;
    
    for (const auto& medium : simulator.mediums) {
        const auto& medium_metrics = medium.get_metrics();
        data.total_absorption += medium_metrics.get_total_absorption();
        data.specular_reflection += medium_metrics.get_surface_reflection();
        data.diffuse_reflection += medium_metrics.get_diffuse_reflection();
        data.surface_refraction += medium_metrics.get_surface_refraction();
        data.specular_transmission += medium_metrics.get_specular_transmission();
        data.diffuse_transmission += medium_metrics.get_diffuse_transmission();
        data.scatter_events += medium_metrics.get_scatter_events();
    }
    
    // Cache the result
    cached_medium_energy_ = data;
    medium_energy_cache_version_ = current_version;
    
    return data;
}

/***********************************************************
 * CALCULATE ENERGY CONSERVATION
 * Computes comprehensive energy conservation data
 ***********************************************************/
Metrics::EnergyConservation Metrics::calculate_energy_conservation(const Simulator& simulator) const {
    // Check cache first - use simulation version as the version counter
    const uint64_t current_version = simulator.get_simulation_version();
    if (cached_conservation_ && conservation_cache_version_ == current_version) {
        return *cached_conservation_;
    }
    
    EnergyConservation result;
    
    // Calculate total voxel-based energy for accurate accounting
    double total_voxel_absorption = 0.0;
    double total_voxel_specular_reflection = 0.0;
    double total_voxel_diffuse_reflection = 0.0;
    double total_voxel_transmission = 0.0;
    
    for (const auto& medium : simulator.mediums) {
        const auto& volume = medium.get_volume();
        for (const auto& voxel : volume) {
            if (voxel && voxel->material != nullptr) {
                total_voxel_absorption += voxel->absorption;
                total_voxel_specular_reflection += voxel->specular_reflection;
                total_voxel_diffuse_reflection += voxel->diffuse_reflection;
                total_voxel_transmission += voxel->diffuse_transmission;
            }
        }
    }
    
    result.total_absorption = total_voxel_absorption;
    result.total_reflection = total_voxel_diffuse_reflection;
    result.total_transmission = total_voxel_transmission;
    result.total_diffusion = total_voxel_diffuse_reflection + total_voxel_transmission;
    result.surface_reflection = total_voxel_specular_reflection;
    result.surface_refraction = aggregate_medium_energy_data(simulator).surface_refraction;
    result.total_energy = result.surface_reflection + result.total_absorption + 
                         result.total_reflection + result.total_transmission;
    
    // Cache the result
    cached_conservation_ = result;
    conservation_cache_version_ = current_version;
    
    return result;
}

/***********************************************************
 * CALCULATE ENERGY CONSERVATION PERCENTAGES
 * Eliminates code duplication across all output methods
 ***********************************************************/
Metrics::EnergyConservationPercentages Metrics::calculate_energy_percentages(const Simulator& simulator) const {
    // Check cache first - use simulation version as the version counter
    const uint64_t current_version = simulator.get_simulation_version();
    if (cached_percentages_ && percentages_cache_version_ == current_version) {
        return *cached_percentages_;
    }
    
    EnergyConservationPercentages percentages;
    
    auto energy = calculate_energy_conservation(simulator);
    
    if (energy.surface_refraction > 0) {
        percentages.baseline_energy = energy.surface_reflection + energy.surface_refraction;
        
        // Calculate percentages relative to total initial energy
        percentages.surface_reflection_percent = (energy.surface_reflection / percentages.baseline_energy) * 100.0;
        percentages.absorption_percent = (energy.total_absorption / percentages.baseline_energy) * 100.0;
        percentages.reflection_percent = (energy.total_reflection / percentages.baseline_energy) * 100.0;
        percentages.transmission_percent = (energy.total_transmission / percentages.baseline_energy) * 100.0;
        
        // Calculate total accounted energy and percentage
        percentages.total_accounted = energy.surface_reflection + energy.total_absorption + 
                                    energy.total_reflection + energy.total_transmission;
        percentages.total_percent = (percentages.total_accounted / percentages.baseline_energy) * 100.0;
        
        // Check energy conservation (within tolerance)
        percentages.is_conserved = std::abs(percentages.total_percent - 100.0) <= 2.0;
    } else {
        // No energy in simulation - all percentages remain 0.0
        percentages.is_conserved = true;
    }
    
    // Cache the result
    cached_percentages_ = percentages;
    percentages_cache_version_ = current_version;
    
    return percentages;
}

/***********************************************************
 * UNIFIED ENERGY DISPLAY DATA
 * Single method for both GUI and export to eliminate redundant calls
 ***********************************************************/
Metrics::EnergyDisplayData Metrics::get_energy_display_data(const Simulator& simulator) const {
    EnergyDisplayData display_data;
    
    const uint64_t current_version = simulator.get_simulation_version();
    display_data.cached_version = current_version;
    
    // Get both conservation and percentage data (uses existing caching)
    display_data.conservation = calculate_energy_conservation(simulator);
    display_data.percentages = calculate_energy_percentages(simulator);
    display_data.is_valid = (display_data.percentages.baseline_energy > 0.0);
    
    return display_data;
}

// =============================================================================
// ENERGY CONSERVATION VALIDATION
// =============================================================================
bool Metrics::is_energy_conserved(const Simulator& simulator, double tolerance_percent) const {
    auto percentages = calculate_energy_percentages(simulator);
    return std::abs(percentages.total_percent - 100.0) <= tolerance_percent;
}

double Metrics::get_conservation_error_percent(const Simulator& simulator) const {
    auto percentages = calculate_energy_percentages(simulator);
    return std::abs(percentages.total_percent - 100.0);
}

// =============================================================================
// REPORTING METHODS
// =============================================================================
void Metrics::export_energy_conservation_log(const Simulator& simulator, std::ofstream& ofs) const {
    ofs << "Energy Conservation Percentages" << std::endl;
    ofs << "################################################################" << std::endl;
    
    auto percentages = calculate_energy_percentages(simulator);
    
    if (percentages.baseline_energy > 0) {
        write_percentage_line(ofs, "Specular reflection", percentages.surface_reflection_percent);
        write_percentage_line(ofs, "Absorption", percentages.absorption_percent);
        write_percentage_line(ofs, "Reflection", percentages.reflection_percent);
        write_percentage_line(ofs, "Transmission", percentages.transmission_percent);
        write_percentage_line(ofs, "Total", percentages.total_percent);
        
        if (!percentages.is_conserved) {
            ofs << "WARNING: Energy not conserved!" << std::endl;
        }
    } else {
        write_percentage_line(ofs, "Specular reflection", 0.0);
        write_percentage_line(ofs, "Absorption", 0.0);
        write_percentage_line(ofs, "Reflection", 0.0);
        write_percentage_line(ofs, "Transmission", 0.0);
        write_percentage_line(ofs, "Total", 0.0);
    }
    
    ofs << std::endl;
}

void Metrics::print_energy_conservation_console(const Simulator& simulator) const {
    std::cout << "Energy Conservation" << std::endl;
    
    auto percentages = calculate_energy_percentages(simulator);
    
    if (percentages.baseline_energy > 0) {
        print_percentage_line("Specular reflection", percentages.surface_reflection_percent);
        print_percentage_line("Absorption", percentages.absorption_percent);
        print_percentage_line("Reflection", percentages.reflection_percent);
        print_percentage_line("Transmission", percentages.transmission_percent);
        print_percentage_line("Total", percentages.total_percent);
        
        if (!percentages.is_conserved) {
            std::cout << "  WARNING: Energy not conserved!" << std::endl;
        }
    } else {
        print_percentage_line("Specular reflection", 0.0);
        print_percentage_line("Absorption", 0.0);
        print_percentage_line("Reflection", 0.0);
        print_percentage_line("Transmission", 0.0);
        print_percentage_line("Total", 0.0);
    }
    
    std::cout << std::endl;
}

std::string Metrics::get_energy_summary_text(const Simulator& simulator) const {
    auto percentages = calculate_energy_percentages(simulator);
    std::ostringstream oss;
    
    oss << "Energy Conservation: ";
    if (percentages.baseline_energy > 0) {
        oss << std::fixed << std::setprecision(1);
        oss << "Spec=" << percentages.surface_reflection_percent << "%, ";
        oss << "Abs=" << percentages.absorption_percent << "%, ";
        oss << "Refl=" << percentages.reflection_percent << "%, ";
        oss << "Trans=" << percentages.transmission_percent << "%, ";
        oss << "Total=" << percentages.total_percent << "%";
        if (!percentages.is_conserved) {
            oss << " (NOT CONSERVED!)";
        }
    } else {
        oss << "No simulation data";
    }
    
    return oss.str();
}

/***********************************************************
 * DIRECT ENERGY ACCESS METHODS
 * Note: get_combined_* methods removed - use aggregate_medium_energy_data() directly
 ***********************************************************/

// =============================================================================
// HELPER METHODS
// =============================================================================
void Metrics::write_percentage_line(std::ofstream& ofs, const std::string& label, double percent) const {
    ofs << std::left << std::setw(25) << (label + ":") 
        << std::fixed << std::setprecision(1) << percent << "%" << std::endl;
}

void Metrics::print_percentage_line(const std::string& label, double percent) const {
    std::cout << "  " << std::left << std::setw(18) << (label + ":") 
              << std::fixed << std::setprecision(1) << percent << "%" << std::endl;
}

/***********************************************************
 * RESULTS EXPORT METHODS (moved from ResultsExporter)
 * Consolidates all simulation result output functionality
 ***********************************************************/

/***********************************************************
 * MAIN EXPORT METHOD
 * Replaces the monolithic report() function
 ***********************************************************/
void Metrics::export_results(const Simulator& simulator, bool generate_csv) const {
    // Always generate simulation.log
    std::string log_path = get_output_filename("simulation", "log");
    std::ofstream log_file = FileUtils::create_output_file(log_path);
    
    if (!log_file.is_open()) {
        return; // Error already logged by FileUtils
    }
    
    // Export main log sections (now consolidated into single function)
    export_simulation_log(simulator, log_file);
    log_file.close();
    
    // Export CSV files if requested
    if (generate_csv) {
        export_csv_files(simulator);
    }
}

/***********************************************************
 * SIMULATION LOG EXPORT
 * Configuration and recorded parameters
 ***********************************************************/
void Metrics::export_simulation_log(const Simulator& simulator, std::ofstream& ofs) const {
    write_header(ofs, "SIMULATION REPORT");
    
    auto now = std::chrono::system_clock::now();
    auto time_t = std::chrono::system_clock::to_time_t(now);
    auto tm = *std::localtime(&time_t);
    ofs << std::left << std::setw(28) << "Generated:" << std::put_time(&tm, "%Y-%m-%d %H:%M:%S") << std::endl;
    ofs << std::left << std::setw(28) << "Configuration file:" << Config::get().config_filename() << std::endl;
    ofs << std::left << std::setw(28) << "Execution time:" << std::fixed << std::setprecision(1) << elapsed_time_ms_ << " ms" << std::endl;
    ofs << std::endl << "################################################################" << std::endl << std::endl;
    
    // General Simulation Statistics
    const auto& config = Config::get();
    ofs << std::left << std::setw(28) << "Volume Grid:" << config.nx() << "x" << config.ny() << "x" << config.nz() << std::endl;
    ofs << std::left << std::setw(28) << "Voxel Size:" << config.vox_size() << std::endl;
    
    // Use actual photons simulated if available, otherwise config value
    size_t photons_to_report = (simulator.photons.size() > 0) ? simulator.photons.size() : config.num_photons();
    ofs << std::left << std::setw(28) << "Total Photons:" << photons_to_report << std::endl;
    
    // Calculate overall voxel statistics
    auto& mediums = simulator.mediums;
    if (!mediums.empty()) {
        const auto& volume = mediums[0].get_volume();
        size_t total_voxels = volume.size();
        size_t material_voxels = 0;
        
        for (uint64_t idx = 0; idx < volume.size(); ++idx) {
            auto voxel = volume.at(static_cast<uint32_t>(idx));
            if (voxel && voxel->material) {
                material_voxels++;
            }
        }
        
        ofs << std::left << std::setw(28) << "Total Voxels:" << total_voxels << std::endl;
        ofs << std::left << std::setw(28) << "Empty Voxels:" << (total_voxels - material_voxels) << std::endl;
    }
    
    ofs << std::endl << "################################################################" << std::endl << std::endl;
    
    // Layer properties section
    for (const auto& medium : simulator.mediums) {
        const auto& layers = medium.get_layers();
        const auto& materials = medium.get_tissues();
        
        for (size_t i = 0; i < layers.size(); ++i) {
            const auto& layer = layers[i];
            ofs << "Layer " << (i + 1) << ":" << std::endl;
            ofs << std::left << std::setw(28) << "  Triangle count:" << layer.mesh.size() << std::endl;
            
            if (layer.tissue_id < materials.size()) {
                const auto& material = materials[layer.tissue_id];
                ofs << std::left << std::setw(28) << "  Refractive index (eta):" << material.eta() << std::endl;
                ofs << std::left << std::setw(28) << "  Absorption coeff (mua):" << material.mu_a() << std::endl;
                ofs << std::left << std::setw(28) << "  Scattering coeff (mus):" << material.mu_s() << std::endl;
                ofs << std::left << std::setw(28) << "  Anisotropy factor (g):" << material.g() << std::endl;
            }
            ofs << std::endl;
        }
    }

    ofs << "################################################################" << std::endl << std::endl;
    
    // Medium statistics section
    for (size_t i = 0; i < mediums.size(); ++i) {
        auto& medium = mediums[i];
        const auto& metrics = medium.get_metrics();
        
        // Medium header with ID (same as console)
        ofs << "Medium " << (i + 1) << ":" << std::endl;
        
        // Medium Voxel Statistics
        const auto& volume = medium.get_volume();
        size_t material_voxels = 0;
        size_t surface_voxels = 0;
        
        for (uint64_t idx = 0; idx < volume.size(); ++idx) {
            auto voxel = volume.at(static_cast<uint32_t>(idx));
            if (voxel && voxel->material) {
                material_voxels++;
                if (voxel->is_surface_voxel) {
                    surface_voxels++;
                }
            }
        }
        
        // Volume Statistics (aligned to 28 characters)
        ofs << "Volume Statistics" << std::endl;
        ofs << std::left << std::setw(28) << "  Medium Voxels:" << material_voxels << std::endl;
        ofs << std::left << std::setw(28) << "  Surface Voxels:" << surface_voxels << std::endl;
        ofs << std::endl;
        
        // Transport Statistics (aligned to 28 characters)
        ofs << "Transport Statistics" << std::endl;
        ofs << std::left << std::setw(28) << "  Photons entered:" << metrics.get_photons_entered() << std::endl;
        ofs << std::left << std::setw(28) << "  Scatter events:" << std::fixed << std::setprecision(0) 
            << metrics.get_scatter_events() << std::endl;
        ofs << std::left << std::setw(28) << "  Total path length:" << std::fixed << std::setprecision(6) 
            << metrics.compute_path_length() << std::endl;
        ofs << std::left << std::setw(28) << "  Average step size:" << std::fixed << std::setprecision(6) 
            << metrics.compute_average_step_size() << std::endl;
        ofs << std::left << std::setw(28) << "  Diffusion distance:" << std::fixed << std::setprecision(6) 
            << metrics.compute_diffusion_distance() << std::endl;
        
        ofs << std::endl;
        
        // Use unified energy display data (single call for both conservation and percentages)
        auto energy_data = get_energy_display_data(simulator);
        auto& energy = energy_data.conservation;
        auto& percentages = energy_data.percentages;

        ofs << "Radiance Properties" << std::endl;
        ofs << std::left << std::setw(28) << "  Total absorption:" << std::fixed << std::setprecision(6) 
            << energy.total_absorption << std::endl;
        ofs << std::left << std::setw(28) << "  Total diffusion:" << std::fixed << std::setprecision(6) 
            << energy.total_diffusion << std::endl;
        ofs << std::left << std::setw(28) << "    Reflection:" << std::fixed << std::setprecision(6) 
            << energy.total_reflection << std::endl;
        ofs << std::left << std::setw(28) << "    Transmission:" << std::fixed << std::setprecision(6) 
            << energy.total_transmission << std::endl;
        
        ofs << std::endl;
        
        // Energy Conservation (exact console format)
        ofs << "Energy Conservation" << std::endl;
        
        if (energy_data.is_valid) {
            ofs << std::left << std::setw(28) << "  Specular reflection:" << std::fixed << std::setprecision(1) 
                << percentages.surface_reflection_percent << "%" << std::endl;
            ofs << std::left << std::setw(28) << "  Absorption:" << std::fixed << std::setprecision(1) 
                << percentages.absorption_percent << "%" << std::endl;
            ofs << std::left << std::setw(28) << "  Reflection:" << std::fixed << std::setprecision(1) 
                << percentages.reflection_percent << "%" << std::endl;
            ofs << std::left << std::setw(28) << "  Transmission:" << std::fixed << std::setprecision(1) 
                << percentages.transmission_percent << "%" << std::endl;
            ofs << std::left << std::setw(28) << "  Total:" << std::fixed << std::setprecision(1) 
                << percentages.total_percent << "%" << std::endl;

            // Check if energy is conserved
            if (!percentages.is_conserved) {
                ofs << "  WARNING: Energy not conserved!" << std::endl;
            }
        } else {
            ofs << std::left << std::setw(28) << "  Specular reflection:" << "0.0%" << std::endl;
            ofs << std::left << std::setw(28) << "  Absorption:" << "0.0%" << std::endl;
            ofs << std::left << std::setw(28) << "  Reflection:" << "0.0%" << std::endl;
            ofs << std::left << std::setw(28) << "  Transmission:" << "0.0%" << std::endl;
            ofs << std::left << std::setw(28) << "  Total:" << "0.0%" << std::endl;
        }
        
        ofs << std::endl;
    }
}

/***********************************************************
 * MEDIUM STATISTICS EXPORT
 * Per-medium detailed information
 ***********************************************************/
void Metrics::export_medium_statistics(const Simulator& simulator, std::ofstream& ofs) const {
    ofs << "Per-Medium Details" << std::endl;
    ofs << "################################################################" << std::endl;
    
    // Get aggregated energy data for accurate statistics
    auto energy_data = aggregate_medium_energy_data(simulator);
    
    // Access mediums through the public member variable
    const auto& mediums = simulator.mediums;
    for (size_t i = 0; i < mediums.size(); ++i) {
        const auto& medium = mediums[i];
        const auto& medium_metrics = medium.get_metrics();
        
        ofs << "Medium " << (i + 1) << ":" << std::endl;
        
        // Medium Voxel Statistics
        const auto& volume = medium.get_volume();
        size_t material_voxels = 0;
        size_t surface_voxels = 0;
        
        for (uint64_t idx = 0; idx < volume.size(); ++idx) {
            auto voxel = volume.at(static_cast<uint32_t>(idx));
            if (voxel && voxel->material) {
                material_voxels++;
                if (voxel->is_surface_voxel) {
                    surface_voxels++;
                }
            }
        }
        
        // Use correct aggregated energy data with percentage formatting
        // Calculate percentages for single medium case
        auto percentages = calculate_energy_percentages(simulator);
        
        ofs << "  Medium voxels:          " << material_voxels << std::endl;
        ofs << "  Surface voxels:         " << surface_voxels << std::endl;
        ofs << "  Photons entered:        " << medium_metrics.get_photons_entered() << std::endl;
        ofs << "  Scatter events:         " << static_cast<int>(medium_metrics.get_scatter_events()) << std::endl;
        
        ofs << "  Total absorption:       " << std::fixed << std::setprecision(1) << percentages.absorption_percent << "%" << std::endl;
        ofs << "  Surface reflection:     " << std::fixed << std::setprecision(1) << percentages.surface_reflection_percent << "%" << std::endl;
        ofs << "  Surface refraction:     " << std::fixed << std::setprecision(1) << (100.0 - percentages.surface_reflection_percent) << "%" << std::endl;
        ofs << "  Diffuse reflection:     " << std::fixed << std::setprecision(1) << percentages.reflection_percent << "%" << std::endl;
        ofs << "  Diffuse transmission:   " << std::fixed << std::setprecision(1) << percentages.transmission_percent << "%" << std::endl;
        ofs << "  Specular transmission:  " << std::fixed << std::setprecision(1) << "0.0%" << std::endl;
        ofs << std::endl;
    }
}

/***********************************************************
 * CSV FILES EXPORT
 * All CSV output functionality
 ***********************************************************/
void Metrics::export_csv_files(const Simulator& simulator) const {
    // Export absorption CSV
    std::string absorption_path = get_output_filename("absorption", "csv");
    std::ofstream absorption_file = FileUtils::create_output_file(absorption_path);
    if (absorption_file.is_open()) {
        export_absorption_csv(simulator, absorption_file);
        absorption_file.close();
    }
    
    // Export photon paths CSV
    std::string photon_path = get_output_filename("photons", "csv");
    std::ofstream photon_file = FileUtils::create_output_file(photon_path);
    if (photon_file.is_open()) {
        export_photon_paths_csv(simulator, photon_file);
        photon_file.close();
    }
    
    // Export emittance data CSV (only voxels with emittance > 0)
    std::string emittance_path = get_output_filename("emittance", "csv");
    std::ofstream emittance_file = FileUtils::create_output_file(emittance_path);
    if (emittance_file.is_open()) {
        export_emittance_data_csv(simulator, emittance_file);
        emittance_file.close();
    }
}

// =============================================================================
// HELPER METHODS FOR RESULTS EXPORT
// =============================================================================
std::string Metrics::get_output_filename(const std::string& base_name, const std::string& extension) const {
    return App::get_output_path(base_name + "." + extension);
}

void Metrics::write_header(std::ofstream& ofs, const std::string& title) const {
    FileUtils::write_file_header(ofs, title);
}

// Placeholder implementations for CSV methods (to be filled based on existing report() logic)
void Metrics::export_absorption_csv(const Simulator& simulator, std::ofstream& ofs) const {
    ofs << "MediumID,VoxelX,VoxelY,VoxelZ,WorldX,WorldY,WorldZ,Absorption,Weight\n";
    
    // Export absorption data for each voxel that has absorbed energy
    for (const auto& medium : simulator.mediums) {
        const auto& volume = medium.get_volume();
        const auto& bounds = medium.get_bounds();
        double voxel_size = volume.voxel_size();
        
        for (const auto& voxel : volume) {
            if (voxel && voxel->material != nullptr && voxel->absorption > 0.0) {
                // Calculate world coordinates of voxel center using bounds
                glm::dvec3 world_pos = bounds.min_bounds + glm::dvec3(
                    (voxel->ix() + 0.5) * voxel_size,
                    (voxel->iy() + 0.5) * voxel_size,
                    (voxel->iz() + 0.5) * voxel_size
                );
                
                ofs << 0 << ","  // MediumID
                    << voxel->ix() << "," << voxel->iy() << "," << voxel->iz() << ","
                    << world_pos.x << "," << world_pos.y << "," << world_pos.z << ","
                    << voxel->absorption << ","
                    << 1.0 << "\n";  // Weight (normalized)
            }
        }
    }
}

void Metrics::export_photon_paths_csv(const Simulator& simulator, std::ofstream& ofs) const {
    ofs << "PhotonID,EntranceX,EntranceY,EntranceZ,ExitX,ExitY,ExitZ,TerminationX,TerminationY,TerminationZ,"
        << "PathLength,ScatterCount,AbsorptionDeposited,ExitedMedium,TerminationReason\n";
    
    const auto& photons = simulator.photons;
    for (const auto& photon : photons) {
        // Use new Photon class methods for detailed path analysis
        auto entrance_pos = photon.get_entrance_position();
        auto exit_pos = photon.get_exit_position();
        auto termination_pos = photon.get_termination_position();
        
        ofs << photon.id << ","
            << entrance_pos.x << "," << entrance_pos.y << "," << entrance_pos.z << ","
            << exit_pos.x << "," << exit_pos.y << "," << exit_pos.z << ","
            << termination_pos.x << "," << termination_pos.y << "," << termination_pos.z << ","
            << photon.get_total_path_length() << ","
            << photon.get_path_scatter_count() << ","
            << photon.get_total_absorption_deposited() << ","
            << (photon.exited_medium() ? "TRUE" : "FALSE") << ","
            << photon.get_termination_reason() << "\n";
    }
}


void Metrics::export_emittance_data_csv(const Simulator& simulator, std::ofstream& ofs) const {
    ofs << "MediumID,VoxelX,VoxelY,VoxelZ,WorldX,WorldY,WorldZ,IsSurface,HasTissue,Absorption,Emittance,EmittanceReflected,EmittanceTransmitted,TotalEmittance,EnergyNormalized\n";
    
    // Calculate total photon energy launched (each photon starts with weight 1.0)
    double total_photon_energy = static_cast<double>(Config::get().num_photons());
    int total_voxels = 0;
    int surface_voxels = 0;
    int voxels_with_emittance = 0;
    
    // Export detailed voxel information
    for (const auto& medium : simulator.mediums) {
        const auto& volume = medium.get_volume();
        const auto& bounds = medium.get_bounds();
        double voxel_size = volume.voxel_size();
        
        for (const auto& voxel : volume) {
            if (voxel && voxel->material != nullptr) {
                total_voxels++;
                
                // Proper surface detection: check if voxel has any neighbor without material
                bool is_surface = false;
                uint32_t vx = voxel->ix();
                uint32_t vy = voxel->iy(); 
                uint32_t vz = voxel->iz();
                
                // Check all 6 neighbors (up, down, left, right, front, back)
                std::vector<std::tuple<int, int, int>> neighbors = {
                    {-1, 0, 0}, {1, 0, 0},   // left, right
                    {0, -1, 0}, {0, 1, 0},   // down, up
                    {0, 0, -1}, {0, 0, 1}    // back, front
                };
                
                for (const auto& [dx, dy, dz] : neighbors) {
                    int nx = static_cast<int>(vx) + dx;
                    int ny = static_cast<int>(vy) + dy;
                    int nz = static_cast<int>(vz) + dz;
                    
                    // Check bounds
                    if (nx < 0 || ny < 0 || nz < 0 || 
                        nx >= static_cast<int>(volume.width()) ||
                        ny >= static_cast<int>(volume.height()) ||
                        nz >= static_cast<int>(volume.depth())) {
                        is_surface = true;
                        break;
                    }
                    
                    // Check if neighbor has no material
                    const Voxel* neighbor = volume(static_cast<uint32_t>(nx), 
                                                  static_cast<uint32_t>(ny), 
                                                  static_cast<uint32_t>(nz));
                    if (!neighbor || neighbor->material == nullptr) {
                        is_surface = true;
                        break;
                    }
                }
                
                if (is_surface) surface_voxels++;
                
                double total_voxel_emittance = voxel->specular_reflection + voxel->diffuse_transmission;
                if (total_voxel_emittance > 0.0) voxels_with_emittance++;
                
                // Only export voxels with emittance > 0
                if (total_voxel_emittance > 0.0) {
                    // Calculate world coordinates of voxel center using bounds
                    glm::dvec3 world_pos = bounds.min_bounds + glm::dvec3(
                        (voxel->ix() + 0.5) * voxel_size,
                        (voxel->iy() + 0.5) * voxel_size,
                        (voxel->iz() + 0.5) * voxel_size
                    );
                    
                    // Energy normalized by total photon energy launched
                    double energy_normalized = total_photon_energy > 0.0 ? 
                        (voxel->absorption + voxel->total_emittance()) / total_photon_energy : 0.0;
                    
                    ofs << 0 << ","  // MediumID
                        << voxel->ix() << "," << voxel->iy() << "," << voxel->iz() << ","
                        << world_pos.x << "," << world_pos.y << "," << world_pos.z << ","
                        << (is_surface ? "TRUE" : "FALSE") << ","
                        << "TRUE,"  // HasTissue
                        << voxel->absorption << ","
                        << voxel->emittance << ","
                        << voxel->specular_reflection << ","
                        << voxel->diffuse_transmission << ","
                        << total_voxel_emittance << ","
                        << energy_normalized << "\n";
                }
            }
        }
    }
    
    // Log summary statistics (matching original implementation)
    Logger::instance().log_photon_event(
        -1, "SUMMARY", glm::dvec3(0), glm::dvec3(0), 0.0, glm::ivec3(0), -1, total_photon_energy,
        "Total voxels: " + std::to_string(total_voxels) + 
        ", Surface: " + std::to_string(surface_voxels) + 
        ", With emittance: " + std::to_string(voxels_with_emittance)
    );
}
