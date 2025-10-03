/**
 * @file results_exporter.cpp
 * @brief Monte Carlo simulation results export in multiple formats
 *
 * Implements comprehensive results export functionality for Photron simulations:
 * - JSON format export for machine-readable analysis
 * - Human-readable text format for reports
 * - Per-medium energy conservation statistics
 * - Photon path data and scattering metrics
 * - Timestamped output with configuration metadata
 *
 * Essential for post-simulation analysis, validation, and scientific
 * publication of Monte Carlo photon transport results.
 */

#include "results_exporter.hpp"

// Include complete type definitions for export operations
#include <algorithm>
#include <chrono>
#include <filesystem>
#include <iomanip>
#include <sstream>

#include "app.hpp"
#include "common/config.hpp"
#include "common/error_handler.hpp"
#include "common/file_utils.hpp"
#include "simulator/medium.hpp"
#include "simulator/metrics.hpp"
#include "simulator/photon.hpp"
#include "simulator/simulator.hpp"

/**
 * @brief Get singleton ResultsExporter instance
 * @return Reference to global ResultsExporter instance
 *
 * Thread-safe singleton providing centralized results export functionality
 * across the application with consistent formatting and error handling.
 */
ResultsExporter& ResultsExporter::instance() {
	static ResultsExporter exporter;
	return exporter;
}

/**
 * @brief Export simulation results in JSON format
 * @param simulator Simulator instance containing results data
 * @param filepath Target file path for JSON output
 *
 * Creates machine-readable JSON export containing complete simulation results:
 * - Configuration metadata and timestamps
 * - Per-medium energy conservation statistics
 * - Photon path data and scattering metrics
 * - Structured data for automated analysis tools
 */
void ResultsExporter::export_json(const Simulator& simulator, const std::string& filepath) {
	// Construct standardized output path within results directory
	std::filesystem::path results_dir = std::filesystem::path(App::get_executable_directory()) / "out" / "results";
	std::filesystem::path full_path = results_dir / filepath;

	// Create output file with automatic directory structure creation
	std::ofstream file = FileUtils::create_output_file(full_path.string());
	if (!file.is_open()) {
		return; // Error already logged by FileUtils
	}

	// Begin JSON structure with metadata header
	file << "{\n";

	// Include configuration and execution metadata
	std::string sanitized_config = "unknown";
	if (Config::is_initialized()) {
		// Extract config file path if available for traceability
		sanitized_config = sanitize_path_for_json("config_file_path");
	}

	file << "  \"config_file\": \"" << sanitized_config << "\",\n";
	file << "  \"timestamp\": \"" << get_timestamp() << "\",\n";
	file << "  \"medium_statistics\": [\n";

	// Export detailed results for each material medium
	auto& mediums = simulator.mediums;
	for (size_t i = 0; i < mediums.size(); ++i) {
		write_medium_json(file, simulator, i, i == mediums.size() - 1);
	}

	// Complete JSON structure and finalize file
	file << "  ]\n";
	file << "}\n";
	file.close();

	REPORT_INFO("Results exported to " + filepath);
}

/**
 * @brief Export simulation results in human-readable text format
 * @param simulator Simulator instance containing results data
 * @param filepath Target file path for text output
 *
 * Creates human-readable text report containing simulation results:
 * - Formatted tables with energy conservation data
 * - Statistical summaries and analysis
 * - Readable format for reports and documentation
 * - Consistent formatting with headers and sections
 */
void ResultsExporter::export_text(const Simulator& simulator, const std::string& filepath) {
	// Use consistent output directory structure
	std::filesystem::path results_dir = std::filesystem::path(App::get_executable_directory()) / "out" / "results";
	std::filesystem::path full_path = results_dir / filepath;

	// Create formatted text output file
	std::ofstream file = FileUtils::create_output_file(full_path.string());
	if (!file.is_open()) {
		return; // Error already logged by FileUtils
	}

	// Write report header with generation metadata
	write_header(file, "PHOTRON SIMULATION RESULTS");
	file << "Generated: " << get_timestamp() << "\n\n";

	// Export formatted statistics for all simulation mediums
	auto& mediums = simulator.mediums;
	for (size_t i = 0; i < mediums.size(); ++i) {
		file << "Medium " << (i + 1) << "\n";
		file << "================================================================\n";
		write_medium_statistics_text(file, simulator, i);
		write_energy_conservation_text(file, simulator, i);
		if (i < mediums.size() - 1) {
			file << "\n";
		}
	}

	file.close();
	REPORT_INFO("Text results exported to " + filepath);
}

void ResultsExporter::export_csv(const Simulator& simulator, const std::string& base_path) {
	// Export all CSV files with unified naming (matching existing Metrics pattern)
	std::filesystem::path base_dir = std::filesystem::path(App::get_executable_directory()) / "out";

	export_absorption_csv(simulator, (base_dir / (base_path + "_absorption.csv")).string());
	export_emittance_csv(simulator, (base_dir / (base_path + "_emittance.csv")).string());
	export_paths_csv(simulator, (base_dir / (base_path + "_paths.csv")).string());

	REPORT_INFO("CSV files exported with base path: " + base_path);
}

void ResultsExporter::export_simulation_log(const Simulator& simulator, const std::string& filepath) {
	std::filesystem::path full_path = std::filesystem::path(App::get_executable_directory()) / filepath;

	std::ofstream file = FileUtils::create_output_file(full_path.string());
	if (!file.is_open()) {
		return; // Error already logged by FileUtils
	}

	write_header(file, "SIMULATION REPORT");
	file << "Generated: " << get_timestamp() << "\n\n";

	// Configuration section
	if (Config::is_initialized()) {
		file << "Configuration\n";
		file << "################################################################\n\n";
		file << "Number of photons:        " << Config::get().num_photons() << "\n";
		file << "Number of layers:         " << Config::get().num_layers() << "\n";
		file << "Voxel size:               " << Config::get().vox_size() << "\n";
		file << "Grid dimensions:          " << Config::get().nx() << "x" << Config::get().ny() << "x"
			 << Config::get().nz() << "\n\n";
	}

	// Medium statistics and energy conservation
	auto& mediums = simulator.mediums;
	for (size_t i = 0; i < mediums.size(); ++i) {
		file << "Medium " << (i + 1) << " Statistics\n";
		file << "################################################################\n";
		write_medium_statistics_text(file, simulator, i);
		write_energy_conservation_text(file, simulator, i);
		if (i < mediums.size() - 1) {
			file << "\n";
		}
	}

	file.close();
	REPORT_INFO("Simulation log exported to " + filepath);
}

// Absorption CSV export
void ResultsExporter::export_absorption_csv(const Simulator& simulator, const std::string& filepath) {
	std::ofstream file = FileUtils::create_output_file(filepath);
	if (!file.is_open()) {
		return;
	}

	file << "MediumID,VoxelX,VoxelY,VoxelZ,Absorption\n";

	auto& mediums = simulator.mediums;
	for (size_t med_idx = 0; med_idx < mediums.size(); ++med_idx) {
		const auto& medium = mediums[med_idx];
		const auto& volume = medium.get_volume();

		for (const auto& voxel : volume) {
			if (voxel && voxel->material && voxel->absorption > 0.0) {
				file << med_idx << "," << voxel->ix() << "," << voxel->iy() << "," << voxel->iz() << ","
					 << std::scientific << std::setprecision(6) << voxel->absorption << "\n";
			}
		}
	}

	file.close();
}

// Emittance CSV export
void ResultsExporter::export_emittance_csv(const Simulator& simulator, const std::string& filepath) {
	std::ofstream file = FileUtils::create_output_file(filepath);
	if (!file.is_open()) {
		return;
	}

	file << "MediumID,VoxelX,VoxelY,VoxelZ,Emittance\n";

	auto& mediums = simulator.mediums;
	for (size_t med_idx = 0; med_idx < mediums.size(); ++med_idx) {
		const auto& medium = mediums[med_idx];
		const auto& volume = medium.get_volume();

		for (const auto& voxel : volume) {
			if (voxel && voxel->material) {
				double emittance = voxel->total_emittance();
				if (emittance > 0.0) {
					file << med_idx << "," << voxel->ix() << "," << voxel->iy() << "," << voxel->iz() << ","
						 << std::scientific << std::setprecision(6) << emittance << "\n";
				}
			}
		}
	}

	file.close();
}

// Paths CSV export
void ResultsExporter::export_paths_csv(const Simulator& simulator, const std::string& filepath) {
	std::ofstream file = FileUtils::create_output_file(filepath);
	if (!file.is_open()) {
		return;
	}

	file << "PhotonID,EventType,PosX,PosY,PosZ,DirX,DirY,DirZ,Weight\n";

	// Export photon paths using new Photon class methods
	const auto& photons = simulator.photons;
	for (const auto& photon : photons) {
		auto entrance_pos = photon.get_entrance_position();
		auto entrance_dir = photon.get_entrance_direction();

		file << photon.id << ",entrance," << std::scientific << std::setprecision(6) << entrance_pos.x << ","
			 << entrance_pos.y << "," << entrance_pos.z << "," << entrance_dir.x << "," << entrance_dir.y << ","
			 << entrance_dir.z << "," << photon.weight << "\n";

		// Add exit position if photon has exited
		if (photon.has_exit()) {
			auto exit_pos = photon.get_exit_position();
			auto exit_dir = photon.get_exit_direction();

			file << photon.id << ",exit," << std::scientific << std::setprecision(6) << exit_pos.x << "," << exit_pos.y
				 << "," << exit_pos.z << "," << exit_dir.x << "," << exit_dir.y << "," << exit_dir.z << ","
				 << photon.weight << "\n";
		}
	}

	file.close();
}

// Utility methods
std::string ResultsExporter::get_timestamp() {
	auto now = std::chrono::system_clock::now();
	auto time_t = std::chrono::system_clock::to_time_t(now);
	auto tm = *std::localtime(&time_t);

	std::ostringstream oss;
	oss << std::put_time(&tm, "%Y-%m-%d %H:%M:%S");
	return oss.str();
}

void ResultsExporter::write_header(std::ostream& out, const std::string& title) {
	out << "################################################################\n";
	out << "# " << title << "\n";
	out << "################################################################\n";
}

std::string ResultsExporter::sanitize_path_for_json(const std::string& path) {
	std::string sanitized = path;
	std::replace(sanitized.begin(), sanitized.end(), '\\', '/');
	return sanitized;
}

// JSON helpers
void ResultsExporter::write_medium_json(std::ostream& out, const Simulator& sim, size_t medium_index, bool last) {
	if (medium_index >= sim.mediums.size()) {
		return;
	}

	const auto& medium = sim.mediums[medium_index];
	const auto& metrics = medium.get_metrics();
	const auto& volume = medium.get_volume();
	const auto& dimensions = volume.dimensions();
	auto& materials = medium.get_tissues();

	// Count surface voxels (matching App::save_results_as_json)
	size_t surface_count = 0;
	for (uint64_t idx = 0; idx < volume.size(); ++idx) {
		auto voxel = volume.at(static_cast<uint32_t>(idx));
		if (voxel && voxel->is_surface_voxel) {
			surface_count++;
		}
	}

	out << "    {\n";
	out << "      \"medium_id\": " << (medium_index + 1) << ",\n";
	out << "      \"volume_statistics\": {\n";
	out << "        \"grid_size\": [" << dimensions.x << ", " << dimensions.y << ", " << dimensions.z << "],\n";
	out << "        \"voxel_size\": " << volume.voxel_size() << ",\n";
	out << "        \"total_voxels\": " << volume.size() << ",\n";
	out << "        \"surface_voxels\": " << surface_count << "\n";
	out << "      },\n";
	out << "      \"transport_statistics\": {\n";
	out << "        \"total_photons\": " << sim.photons.size() << ",\n";
	out << "        \"photons_entered\": " << metrics.get_photons_entered() << ",\n";
	out << "        \"scatter_events\": " << metrics.get_scatter_events() << ",\n";
	out << "        \"path_length\": " << metrics.compute_path_length() << ",\n";
	out << "        \"average_step_size\": " << metrics.compute_average_step_size() << ",\n";
	out << "        \"diffusion_distance\": " << metrics.compute_diffusion_distance() << "\n";
	out << "      },\n";
	out << "      \"energy_conservation\": {\n";
	out << "        \"total_absorption\": " << metrics.get_total_absorption() << ",\n";
	out << "        \"diffuse_reflection\": " << metrics.get_diffuse_reflection() << ",\n";
	out << "        \"specular_reflection\": " << metrics.get_surface_reflection() << ",\n";
	out << "        \"diffuse_transmission\": " << metrics.get_diffuse_transmission() << ",\n";
	out << "        \"specular_transmission\": " << metrics.get_specular_transmission() << ",\n";
	out << "        \"surface_refraction\": " << metrics.get_surface_refraction() << "\n";
	out << "      },\n";
	out << "      \"tissue_properties\": [\n";

	// Export materials for this medium (matching App::save_results_as_json)
	for (size_t j = 0; j < materials.size(); ++j) {
		const auto& material = materials[j];
		out << "        {\n";
		out << "          \"index\": " << j << ",\n";
		out << "          \"hash\": " << material.get_optical_properties_hash() << ",\n";
		out << "          \"eta\": " << material.eta() << ",\n";
		out << "          \"mua\": " << material.mu_a() << ",\n";
		out << "          \"mus\": " << material.mu_s() << ",\n";
		out << "          \"ani\": " << material.g() << "\n";
		out << "        }" << (j < materials.size() - 1 ? "," : "") << "\n";
	}
	out << "      ]\n";

	out << "    }" << (last ? "" : ",") << "\n";
}

// Text helpers
void ResultsExporter::write_medium_statistics_text(std::ostream& out, const Simulator& sim, size_t medium_index) {
	if (medium_index >= sim.mediums.size())
		return;

	const auto& medium = sim.mediums[medium_index];
	const auto& metrics = medium.get_metrics();
	const auto& volume = medium.get_volume();
	const auto& dimensions = volume.dimensions();

	out << "Volume Statistics\n";
	out << "  Volume Grid:         " << dimensions.x << "x" << dimensions.y << "x" << dimensions.z << "\n";

	// Count different types of voxels
	size_t total_voxels = volume.size();
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

	out << "  Total Volume Voxels: " << total_voxels << "\n";
	out << "  Material Voxels:     " << material_voxels << "\n";
	out << "  Empty Voxels:        " << (total_voxels - material_voxels) << "\n";
	out << "  Surface Voxels:      " << surface_voxels << "\n";

	out << "\nTransport Statistics\n";
	out << "  Total photons:       " << sim.photons.size() << "\n";
	out << "  Photons entered:     " << metrics.get_photons_entered() << "\n";
	out << "  Scatter events:      " << static_cast<int>(metrics.get_scatter_events()) << "\n";
	out << "  Total path length:   " << std::fixed << std::setprecision(6) << metrics.compute_path_length() << "\n";
	out << "  Average step size:   " << std::fixed << std::setprecision(6) << metrics.compute_average_step_size()
		<< "\n";
	out << "  Diffusion distance:  " << std::fixed << std::setprecision(6) << metrics.compute_diffusion_distance()
		<< "\n";

	// Use unified energy display data (single call for both conservation and percentages)
	auto energy_data = sim.get_metrics().get_energy_display_data(sim);
	const auto& energy = energy_data.conservation;

	out << "\nRadiance Properties\n";
	out << "  Total absorption:    " << std::fixed << std::setprecision(6) << energy.total_absorption << "\n";
	out << "  Total diffusion:     " << std::fixed << std::setprecision(6) << energy.total_diffusion << "\n";
	out << "    Reflection:        " << std::fixed << std::setprecision(6) << energy.total_reflection << "\n";
	out << "    Transmission:      " << std::fixed << std::setprecision(6) << energy.total_transmission << "\n";
}

void ResultsExporter::write_energy_conservation_text(std::ostream& out, const Simulator& sim, size_t medium_index) {
	if (medium_index >= sim.mediums.size()) {
		return;
	}

	const auto& medium = sim.mediums[medium_index];
	auto& materials = medium.get_tissues();

	// Use unified energy display data (same data as used for radiance properties)
	auto energy_data = sim.get_metrics().get_energy_display_data(sim);
	const auto& percentages = energy_data.percentages;

	out << "\nEnergy Conservation\n";
	if (energy_data.is_valid) {
		out << "  Surface reflection:  " << std::fixed << std::setprecision(1) << percentages.surface_reflection_percent
			<< "%\n";
		out << "  Absorption:          " << std::fixed << std::setprecision(1) << percentages.absorption_percent
			<< "%\n";
		out << "  Reflection:          " << std::fixed << std::setprecision(1) << percentages.reflection_percent
			<< "%\n";
		out << "  Transmission:        " << std::fixed << std::setprecision(1) << percentages.transmission_percent
			<< "%\n";
		out << "  Total:               " << std::fixed << std::setprecision(1) << percentages.total_percent << "%\n";

		if (!percentages.is_conserved) {
			out << "  WARNING: Energy not conserved!\n";
		}
	}
	else {
		out << "  Surface reflection:  0.0%\n";
		out << "  Absorption:          0.0%\n";
		out << "  Reflection:          0.0%\n";
		out << "  Transmission:        0.0%\n";
		out << "  Total:               0.0%\n";
	}

	// Add material properties for this medium (matching App::save_results_as_text)
	out << "\nTissue Properties\n";
	for (size_t idx = 0; idx < materials.size(); ++idx) {
		const auto& material = materials[idx];
		out << "  material " << idx << " (hash: " << material.get_optical_properties_hash() << "):\n";
		out << "    Refractive Index (eta): " << material.eta() << "\n";
		out << "    Absorption Coefficient (mua): " << material.mu_a() << " cm^-1\n";
		out << "    Scattering Coefficient (mus): " << material.mu_s() << " cm^-1\n";
		out << "    Anisotropy Factor (g): " << material.g() << "\n";
	}
}
