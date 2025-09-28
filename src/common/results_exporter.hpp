#pragma once

#include <string>
#include <fstream>
#include <memory>

class Simulator;

/**
 * @brief Unified results export manager
 * 
 * Consolidates all file output operations from App, Metrics, and Medium classes
 * into a single, consistent interface. Handles JSON, text, and CSV export formats.
 */
class ResultsExporter {
public:
    // Singleton access
    static ResultsExporter& instance();
    
    // Main export methods
    void export_json(const Simulator& simulator, const std::string& filepath);
    void export_text(const Simulator& simulator, const std::string& filepath);
    void export_csv(const Simulator& simulator, const std::string& base_path);
    
    // Unified simulation log export (replaces Metrics::export_simulation_log)
    void export_simulation_log(const Simulator& simulator, const std::string& filepath);
    
    // Individual CSV exports
    void export_absorption_csv(const Simulator& simulator, const std::string& filepath);
    void export_emittance_csv(const Simulator& simulator, const std::string& filepath);
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