#pragma once

#include <fstream>
#include <string>
#include <sstream>
#include <chrono>
#include <iomanip>
#include <glm/glm.hpp>

class Logger {
public:
    static Logger& instance();

    void initialize(const std::string& csv_filepath, const std::string& log_filepath, bool enable_logging = true);

    void log_photon_event(int photon_id, const std::string& event, 
                         const glm::dvec3& position, const glm::dvec3& direction, 
                         double weight, const glm::ivec3& voxel_coords = glm::ivec3(-1),
                         int medium_id = -1, double energy = 0.0,
                         const std::string& description = "");

    void log_voxel_emittance(int photon_id, const glm::dvec3& position, const glm::dvec3& direction,
                           double weight, const glm::ivec3& voxel_coords, double emittance, 
                           const std::string& surface_type = "");

    // General logging methods for replacing log console output
    void log_info(const std::string& message);
    void log_debug(const std::string& message);
    void log_warning(const std::string& message);
    void log_error(const std::string& message);

    // Simple formatted logging using stringstream
    void log_info_fmt(const std::string& message);
    void log_debug_fmt(const std::string& message);
    void log_warning_fmt(const std::string& message);
    void log_error_fmt(const std::string& message);

    ~Logger();

private:
    std::ofstream csv_file_;
    std::ofstream log_file_;
    bool logging_enabled_ = false;
    
    Logger() = default;

    std::string get_timestamp() const;
};
