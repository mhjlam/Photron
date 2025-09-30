/**
 * @file logger.hpp
 * @brief Debug logging system for Monte Carlo photon transport simulation
 * 
 * Provides comprehensive logging capabilities for debugging photon paths,
 * energy conservation validation, and performance analysis. The logging system
 * is designed to be lightweight with compile-time optimization for release builds.
 */

#pragma once

// Standard library includes
#include <fstream>
#include <string>
#include <sstream>
#include <chrono>
#include <iomanip>

// Third-party includes for interface types
#define GLM_ENABLE_EXPERIMENTAL
#include <glm/glm.hpp>

// Forward declarations
class Config;

/**
 * @class Logger
 * @brief Singleton logging system for simulation debugging and analysis
 * 
 * The Logger class provides comprehensive logging capabilities specifically
 * designed for Monte Carlo photon transport simulation debugging:
 * 
 * **Key Features:**
 * - **Photon Event Tracking**: Detailed logging of photon interactions and state changes
 * - **Energy Conservation**: Validation logging for energy conservation analysis
 * - **Performance Monitoring**: Execution time and resource usage tracking
 * - **Conditional Compilation**: Zero-overhead logging in release builds
 * 
 * **Output Formats:**
 * - **CSV Files**: Structured data for analysis in spreadsheet applications
 * - **Debug Logs**: Human-readable text output for development debugging
 * - **Console Output**: Real-time feedback during simulation execution
 * 
 * **Optimization:**
 * - Compile-time macros eliminate logging overhead in release builds
 * - Singleton pattern ensures minimal memory footprint
 * - Buffered file I/O for optimal performance
 * 
 * **Usage Pattern:**
 * ```cpp
 * FAST_LOG_PHOTON_EVENT(photon_id, "scatter", position, direction, weight, voxel_coords, medium_id, energy, "Rayleigh scattering");
 * FAST_LOG_INFO("Simulation completed with " + std::to_string(photon_count) + " photons");
 * ```
 */
class Logger {
public:
    /**
     * @brief Get singleton Logger instance
     * 
     * Thread-safe singleton access using Meyer's singleton pattern.
     * 
     * @return Logger& Reference to global Logger instance
     */
    static Logger& instance();

    /**
     * @brief Initialize logging system with output file paths
     * 
     * Sets up file streams for CSV data export and debug text logging.
     * Must be called before any logging operations.
     * 
     * @param csv_filepath Path for CSV format photon event data
     * @param log_filepath Path for human-readable debug log file
     * @param enable_logging True to enable file output, false for console only
     */
    void initialize(const std::string& csv_filepath, const std::string& log_filepath, bool enable_logging = true);

    /**
     * @brief Log detailed photon interaction event
     * 
     * Records comprehensive information about photon state changes during
     * Monte Carlo simulation for detailed analysis and debugging.
     * 
     * @param photon_id Unique identifier for the photon being tracked
     * @param event Event type ("launch", "scatter", "absorb", "transmit", etc.)
     * @param position Current photon position in world coordinates
     * @param direction Current photon direction (normalized vector)
     * @param weight Current photon statistical weight [0,1]
     * @param voxel_coords Grid coordinates of current voxel (optional)
     * @param medium_id ID of current medium/material (optional)
     * @param energy Energy value associated with event (optional)
     * @param description Human-readable event description (optional)
     */
    void log_photon_event(int photon_id, const std::string& event, 
                          const glm::dvec3& position, const glm::dvec3& direction, 
                          double weight, const glm::ivec3& voxel_coords = glm::ivec3(-1),
                          int medium_id = -1, double energy = 0.0,
                          const std::string& description = "");

    /**
     * @brief Log voxel energy emittance for conservation validation
     * 
     * Records energy emission events from voxels for energy conservation
     * analysis and validation of simulation accuracy.
     * 
     * @param photon_id ID of photon contributing to emittance
     * @param position World position of emittance event
     * @param direction Direction of energy emission
     * @param weight Photon weight contributing to emittance
     * @param voxel_coords Grid coordinates of emitting voxel
     * @param emittance Amount of energy being emitted
     * @param surface_type Type of surface interaction ("boundary", "interface", etc.)
     */
    void log_voxel_emittance(int photon_id, const glm::dvec3& position, const glm::dvec3& direction,
                             double weight, const glm::ivec3& voxel_coords, double emittance, 
                             const std::string& surface_type = "");

    // General purpose logging methods
    
    /**
     * @brief Log informational message
     * @param message Information message text
     */
    void log_info(const std::string& message);
    
    /**
     * @brief Log debug message (development builds only)
     * @param message Debug message text
     */
    void log_debug(const std::string& message);
    
    /**
     * @brief Log warning message
     * @param message Warning message text
     */
    void log_warning(const std::string& message);
    
    /**
     * @brief Log error message
     * @param message Error message text
     */
    void log_error(const std::string& message);

    /**
     * @brief Destructor ensures proper file stream closure
     */
    ~Logger();

    // Compile-time optimization macros for zero-overhead logging
    
    /**
     * @def FAST_LOG_INFO
     * @brief High-performance info logging with compile-time optimization
     * 
     * In debug builds, checks logging state and calls Logger::log_info().
     * In release builds, expands to complete no-op for zero overhead.
     */
    
    /**
     * @def FAST_LOG_PHOTON_EVENT
     * @brief High-performance photon event logging
     * 
     * Optimized macro for frequent photon event logging during simulation.
     * Completely eliminated in release builds for maximum performance.
     */
     
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
    std::ofstream csv_file_;        ///< CSV format output stream for structured data
    std::ofstream log_file_;        ///< Text format output stream for log messages
    bool logging_enabled_ {false};  ///< Global logging enable/disable flag
    
    /**
     * @brief Private constructor for singleton pattern
     */
    Logger() = default;

    /**
     * @brief Generate ISO timestamp string for log entries
     * @return std::string Current timestamp in ISO format
     */
    std::string get_timestamp() const;
};
