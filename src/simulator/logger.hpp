#pragma once

#include <fstream>
#include <string>
#include <sstream>
#include <chrono>
#include <iomanip>
#include <glm/glm.hpp>

// Forward declaration for Config (needed by FAST_LOG macros)
class Config;

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

    ~Logger();

    // Convenience macros that eliminate Config::get().log() checks
#ifdef _DEBUG
    #define FAST_LOG_INFO(msg) do { if (Config::get().log()) { Logger::instance().log_info(msg); } } while(0)
    #define FAST_LOG_DEBUG(msg) do { if (Config::get().log()) { Logger::instance().log_debug(msg); } } while(0)  
    #define FAST_LOG_WARNING(msg) do { if (Config::get().log()) { Logger::instance().log_warning(msg); } } while(0)
    #define FAST_LOG_ERROR(msg) do { if (Config::get().log()) { Logger::instance().log_error(msg); } } while(0)
    #define FAST_LOG_PHOTON_EVENT(id, event, pos, dir, weight, voxel, med, energy, desc) \
        do { if (Config::get().log()) { Logger::instance().log_photon_event(id, event, pos, dir, weight, voxel, med, energy, desc); } } while(0)
#else
    // Complete no-ops in release build - compiler eliminates entirely
    #define FAST_LOG_INFO(msg) ((void)0)
    #define FAST_LOG_DEBUG(msg) ((void)0)
    #define FAST_LOG_WARNING(msg) ((void)0)
    #define FAST_LOG_ERROR(msg) ((void)0)
    #define FAST_LOG_PHOTON_EVENT(id, event, pos, dir, weight, voxel, med, energy, desc) ((void)0)
#endif

private:
    std::ofstream csv_file_;
    std::ofstream log_file_;
    bool logging_enabled_ = false;
    
    Logger() = default;

    std::string get_timestamp() const;
};
