/**
 * @file results_exporter.hpp
 * @brief Centralized results export system for simulation data
 *
 * Consolidates all file output operations from various system components
 * into a unified interface. Supports multiple export formats (JSON, text, CSV)
 * with consistent formatting and error handling.
 */

#pragma once

#include <fstream>
#include <memory>
#include <string>

class Simulator;

/**
 * @class ResultsExporter
 * @brief Unified results export manager with multiple format support
 *
 * The ResultsExporter class consolidates all file output operations from App,
 * Metrics, and Medium classes into a single, consistent interface. Key features:
 *
 * **Export Formats:**
 * - JSON: Structured data for programmatic analysis
 * - Text: Human-readable summary reports
 * - CSV: Spreadsheet-compatible data tables
 *
 * **Data Categories:**
 * - Complete simulation results with all metrics
 * - Individual absorption/emittance data
 * - Photon path traces for visualization
 * - Energy conservation summaries
 *
 * **Design Pattern:**
 * - Singleton pattern for global access
 * - Consistent file naming and directory structure
 * - Automatic timestamp generation
 * - Standardized headers and formatting
 *
 * **Usage Examples:**
 * ```cpp
 * auto& exporter = ResultsExporter::instance();
 * exporter.export_json(simulator, "results/simulation_001.json");
 * exporter.export_csv(simulator, "results/simulation_001");  // Creates multiple CSV files
 * ```
 */
class ResultsExporter
{
public:
	/**
	 * @brief Get singleton instance of ResultsExporter
	 * @return ResultsExporter& Global results exporter instance
	 */
	static ResultsExporter& instance();

	// Main export methods

	/**
	 * @brief Export complete simulation results as structured JSON
	 *
	 * Creates comprehensive JSON export containing all simulation metrics,
	 * energy conservation data, medium properties, and configuration parameters.
	 *
	 * @param simulator Simulator instance containing results to export
	 * @param filepath Output JSON file path
	 */
	void export_json(const Simulator& simulator, const std::string& filepath);

	/**
	 * @brief Export human-readable text summary report
	 *
	 * Creates formatted text report with simulation overview, energy conservation
	 * analysis, and key metrics summary. Ideal for quick review and documentation.
	 *
	 * @param simulator Simulator instance containing results to export
	 * @param filepath Output text file path
	 */
	void export_text(const Simulator& simulator, const std::string& filepath);

	/**
	 * @brief Export simulation data as multiple CSV files
	 *
	 * Creates separate CSV files for different data categories:
	 * - Absorption data per voxel
	 * - Emittance data per surface
	 * - Photon path traces
	 *
	 * @param simulator Simulator instance containing results to export
	 * @param base_path Base path for CSV files (extensions added automatically)
	 */
	void export_csv(const Simulator& simulator, const std::string& base_path);

	/**
	 * @brief Export detailed simulation execution log
	 *
	 * Replaces Metrics::export_simulation_log with unified interface.
	 * Contains timing information, iteration counts, and performance metrics.
	 *
	 * @param simulator Simulator instance containing execution metrics
	 * @param filepath Output log file path
	 */
	void export_simulation_log(const Simulator& simulator, const std::string& filepath);

	// Individual CSV export methods

	/**
	 * @brief Export voxel absorption data as CSV
	 * @param simulator Simulator instance with absorption data
	 * @param filepath Output CSV file path
	 */
	void export_absorption_csv(const Simulator& simulator, const std::string& filepath);

	/**
	 * @brief Export surface emittance data as CSV
	 * @param simulator Simulator instance with emittance data
	 * @param filepath Output CSV file path
	 */
	void export_emittance_csv(const Simulator& simulator, const std::string& filepath);

	/**
	 * @brief Export photon path traces as CSV
	 * @param simulator Simulator instance with path data
	 * @param filepath Output CSV file path
	 */
	void export_paths_csv(const Simulator& simulator, const std::string& filepath);

private:
	ResultsExporter() = default;

	// Utility methods
	static std::string get_timestamp();
	static void write_header(std::ostream& out, const std::string& title);
	static std::string sanitize_path_for_json(const std::string& path);

	// JSON export helpers
	void write_medium_json(std::ostream& out, const Simulator& simulator, size_t medium_index, bool is_last);

	// Text export helpers
	void write_medium_statistics_text(std::ostream& out, const Simulator& simulator, size_t medium_index);
	void write_energy_conservation_text(std::ostream& out, const Simulator& simulator, size_t medium_index);

	// CSV export helpers
	void write_csv_header_absorption(std::ostream& out);
	void write_csv_header_emittance(std::ostream& out);
	void write_csv_header_paths(std::ostream& out);
};