#include "metrics.hpp"

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

void Metrics::print_report(const class Simulator& simulator) {
	std::cout << std::endl << "=== Monte Carlo Simulation Results ===" << std::endl << std::endl;
	
	if (Config::get().log()) {
		Logger::instance().log_info("=== Monte Carlo Simulation Results ===");
	}

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
		auto energy = calculate_energy_conservation(simulator);

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
		
		// Use consolidated percentage calculation from Metrics
		auto percentages = calculate_energy_percentages(simulator);
		
		if (percentages.baseline_energy > 0) {
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
	
	if (Config::get().log()) {
		Logger::instance().log_info("Simulation completed successfully");
		Logger::instance().log_info("Execution time: " + std::to_string(elapsed_time_ms_) + " ms");
		Logger::instance().log_info("Total photons processed: " + std::to_string(simulator.get_paths().size()));
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

/***********************************************************
 * ENERGY STATISTICS METHODS (moved from EnergyStatisticsManager)
 * Consolidates all energy conservation functionality in Metrics
 ***********************************************************/

/***********************************************************
 * AGGREGATE MEDIUM ENERGY DATA
 * Consolidates energy data from all media
 ***********************************************************/
Metrics::MediumEnergyData Metrics::aggregate_medium_energy_data(const Simulator& simulator) const {
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
    
    return data;
}

/***********************************************************
 * CALCULATE ENERGY CONSERVATION
 * Computes comprehensive energy conservation data
 ***********************************************************/
Metrics::EnergyConservation Metrics::calculate_energy_conservation(const Simulator& simulator) const {
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
    result.surface_refraction = get_combined_surface_refraction(simulator);
    result.total_energy = result.surface_reflection + result.total_absorption + 
                         result.total_reflection + result.total_transmission;
    
    return result;
}

/***********************************************************
 * CALCULATE ENERGY CONSERVATION PERCENTAGES
 * Eliminates code duplication across all output methods
 ***********************************************************/
Metrics::EnergyConservationPercentages Metrics::calculate_energy_percentages(const Simulator& simulator) const {
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
    
    return percentages;
}

/***********************************************************
 * ENERGY CONSERVATION VALIDATION
 ***********************************************************/
bool Metrics::is_energy_conserved(const Simulator& simulator, double tolerance_percent) const {
    auto percentages = calculate_energy_percentages(simulator);
    return std::abs(percentages.total_percent - 100.0) <= tolerance_percent;
}

double Metrics::get_conservation_error_percent(const Simulator& simulator) const {
    auto percentages = calculate_energy_percentages(simulator);
    return std::abs(percentages.total_percent - 100.0);
}

/***********************************************************
 * REPORTING METHODS
 ***********************************************************/
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
 * BACKWARD COMPATIBILITY METHODS
 ***********************************************************/
double Metrics::get_combined_total_absorption(const Simulator& simulator) const {
    return aggregate_medium_energy_data(simulator).total_absorption;
}

double Metrics::get_combined_diffuse_reflection(const Simulator& simulator) const {
    return aggregate_medium_energy_data(simulator).diffuse_reflection;
}

double Metrics::get_combined_specular_reflection(const Simulator& simulator) const {
    return aggregate_medium_energy_data(simulator).specular_reflection;
}

double Metrics::get_combined_surface_refraction(const Simulator& simulator) const {
    return aggregate_medium_energy_data(simulator).surface_refraction;
}

double Metrics::get_combined_diffuse_transmission(const Simulator& simulator) const {
    return aggregate_medium_energy_data(simulator).diffuse_transmission;
}

double Metrics::get_combined_specular_transmission(const Simulator& simulator) const {
    return aggregate_medium_energy_data(simulator).specular_transmission;
}

/***********************************************************
 * HELPER METHODS
 ***********************************************************/
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
    std::ofstream log_file(log_path);
    
    if (!log_file.good()) {
        std::cerr << "Error: simulation.log file could not be opened." << std::endl;
        return;
    }
    
    // Export main log sections
    export_simulation_log(simulator, log_file);
    export_energy_conservation_log(simulator, log_file);
    export_medium_statistics(simulator, log_file);
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
    
    // Add timestamp
    auto now = std::chrono::system_clock::now();
    auto time_t = std::chrono::system_clock::to_time_t(now);
    auto tm = *std::localtime(&time_t);
    ofs << "Generated: " << std::put_time(&tm, "%Y-%m-%d %H:%M:%S") << std::endl;
    ofs << std::endl;
    
    // Configuration section
    ofs << "Configuration" << std::endl;
    ofs << "################################################################" << std::endl;
    ofs << std::endl;
    ofs << "Number of photons:        " << Config::get().num_photons() << std::endl;
    ofs << "Number of layers:         " << Config::get().num_layers() << std::endl;
    ofs << "Number of voxels:         " << Config::get().num_voxels() << std::endl;
    ofs << "Grid dimensions:          " << Config::get().nx() 
        << " x " << Config::get().ny() 
        << " x " << Config::get().nz() << std::endl;
    ofs << "Voxel dimensions:         " << Config::get().vox_size() 
        << " x " << Config::get().vox_size() 
        << " x " << Config::get().vox_size() << std::endl;
    ofs << "Deterministic mode:       " << (Config::get().deterministic() ? "Yes" : "No") << std::endl;
    ofs << std::endl << std::endl;
    
    // Recorded parameters section using consolidated method
    ofs << "Recorded parameters" << std::endl;
    ofs << "################################################################" << std::endl;
    
    auto energy_data = aggregate_medium_energy_data(simulator);
    
    ofs << "Total absorption:        " << std::fixed << energy_data.total_absorption << std::endl;
    ofs << "Surface reflection:      " << std::fixed << energy_data.specular_reflection << std::endl;
    ofs << "Surface refraction:      " << std::fixed << energy_data.surface_refraction << std::endl;
    ofs << "Diffuse reflection:      " << std::fixed << energy_data.diffuse_reflection << std::endl;
    ofs << "Diffuse transmission:    " << std::fixed << energy_data.diffuse_transmission << std::endl;
    ofs << "Specular transmission:   " << std::fixed << energy_data.specular_transmission << std::endl;
    ofs << std::endl;
}

/***********************************************************
 * MEDIUM STATISTICS EXPORT
 * Per-medium detailed information
 ***********************************************************/
void Metrics::export_medium_statistics(const Simulator& simulator, std::ofstream& ofs) const {
    ofs << "Per-Medium Details" << std::endl;
    ofs << "################################################################" << std::endl;
    
    // Access mediums through the public member variable
    const auto& mediums = simulator.mediums;
    for (size_t i = 0; i < mediums.size(); ++i) {
        const auto& medium = mediums[i];
        const auto& volume = medium.get_volume();
        const auto& medium_metrics = medium.get_metrics();
        
        ofs << "Medium " << (i + 1) << ":" << std::endl;
        ofs << "  Volume dimensions:      " 
            << volume.width() << " x " 
            << volume.height() << " x " 
            << volume.depth() << std::endl;
        ofs << "  Total voxels:           " << volume.size() << std::endl;
        
        // Count different types of voxels
        size_t material_count = 0, surface_count = 0, boundary_count = 0;
        for (const auto& voxel_ptr : volume) {
            if (voxel_ptr->material) material_count++;
            if (voxel_ptr->is_surface_voxel) surface_count++;
            if (voxel_ptr->is_boundary_voxel) boundary_count++;
        }
        
        ofs << "  Material voxels:        " << material_count << std::endl;
        ofs << "  Surface voxels:         " << surface_count << std::endl;
        ofs << "  Boundary voxels:        " << boundary_count << std::endl;
        ofs << "  Photons entered:        " << medium_metrics.get_photons_entered() << std::endl;
        ofs << "  Total absorption:       " << std::fixed << medium_metrics.get_total_absorption() << std::endl;
        ofs << "  Surface reflection:     " << std::fixed << medium_metrics.get_surface_reflection() << std::endl;
        ofs << "  Surface refraction:     " << std::fixed << medium_metrics.get_surface_refraction() << std::endl;
        ofs << "  Diffuse reflection:     " << std::fixed << medium_metrics.get_diffuse_reflection() << std::endl;
        ofs << "  Diffuse transmission:   " << std::fixed << medium_metrics.get_diffuse_transmission() << std::endl;
        ofs << "  Specular transmission:  " << std::fixed << medium_metrics.get_specular_transmission() << std::endl;
        ofs << "  Scatter events:         " << static_cast<int>(medium_metrics.get_scatter_events()) << std::endl;
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
    std::ofstream absorption_file(absorption_path);
    if (absorption_file.good()) {
        export_absorption_csv(simulator, absorption_file);
        absorption_file.close();
    }
    
    // Export optical properties CSV
    std::string optical_path = get_output_filename("optics", "csv");
    std::ofstream optical_file(optical_path);
    if (optical_file.good()) {
        export_optical_properties_csv(simulator, optical_file);
        optical_file.close();
    }
    
    // Export photon paths CSV
    std::string photon_path = get_output_filename("photons", "csv");
    std::ofstream photon_file(photon_path);
    if (photon_file.good()) {
        export_photon_paths_csv(simulator, photon_file);
        photon_file.close();
    }
    
    // Export voxel statistics CSV
    std::string voxel_path = get_output_filename("voxels", "csv");
    std::ofstream voxel_file(voxel_path);
    if (voxel_file.good()) {
        export_voxel_data_csv(simulator, voxel_file);
        voxel_file.close();
    }
}

/***********************************************************
 * HELPER METHODS FOR RESULTS EXPORT
 ***********************************************************/
std::string Metrics::get_output_filename(const std::string& base_name, const std::string& extension) const {
    return App::get_output_path(base_name + "." + extension);
}

void Metrics::write_header(std::ofstream& ofs, const std::string& title) const {
    ofs.precision(8);
    ofs << "################################################################" << std::endl;
    ofs << "# " << title << std::endl;
    ofs << "################################################################" << std::endl;
    ofs << std::endl;
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
    ofs << "PhotonID,X,Y,Z,DirX,DirY,DirZ,Weight,ScatterCount,VoxelX,VoxelY,VoxelZ\n";
    
    const auto& photons = simulator.get_paths();
    for (const auto& photon : photons) {
        // Use current position and direction
        ofs << photon.id << ","
            << photon.position.x << "," << photon.position.y << "," << photon.position.z << ","
            << photon.direction.x << "," << photon.direction.y << "," << photon.direction.z << ","
            << photon.weight << ","
            << photon.scatter_count << ",";
        
        // Add voxel coordinates if voxel exists
        if (photon.voxel) {
            ofs << photon.voxel->ix() << "," << photon.voxel->iy() << "," << photon.voxel->iz() << "\n";
        } else {
            ofs << "-1,-1,-1\n";
        }
    }
}

void Metrics::export_optical_properties_csv(const Simulator& simulator, std::ofstream& ofs) const {
    ofs << "MediumID,VoxelX,VoxelY,VoxelZ,WorldX,WorldY,WorldZ,Eta,MuA,MuS,G\n";
    
    // Export optical properties for each voxel with material
    for (const auto& medium : simulator.mediums) {
        const auto& volume = medium.get_volume();
        const auto& bounds = medium.get_bounds();
        double voxel_size = volume.voxel_size();
        
        for (const auto& voxel : volume) {
            if (voxel && voxel->material != nullptr) {
                // Calculate world coordinates of voxel center using bounds
                glm::dvec3 world_pos = bounds.min_bounds + glm::dvec3(
                    (voxel->ix() + 0.5) * voxel_size,
                    (voxel->iy() + 0.5) * voxel_size,
                    (voxel->iz() + 0.5) * voxel_size
                );
                
                const auto* material = voxel->material;
                ofs << 0 << ","  // MediumID
                    << voxel->ix() << "," << voxel->iy() << "," << voxel->iz() << ","
                    << world_pos.x << "," << world_pos.y << "," << world_pos.z << ","
                    << material->eta() << ","
                    << material->mu_a() << ","
                    << material->mu_s() << ","
                    << material->g() << "\n";
            }
        }
    }
}

void Metrics::export_voxel_data_csv(const Simulator& simulator, std::ofstream& ofs) const {
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
    
    // Log summary statistics (matching original implementation)
    Logger::instance().log_photon_event(
        -1, "SUMMARY", glm::dvec3(0), glm::dvec3(0), 0.0, glm::ivec3(0), -1, total_photon_energy,
        "Total voxels: " + std::to_string(total_voxels) + 
        ", Surface: " + std::to_string(surface_voxels) + 
        ", With emittance: " + std::to_string(voxels_with_emittance)
    );
}
