#pragma once

#include <string>
#include <fstream>
#include <sstream>
#include <chrono>
#include <iomanip>
#include <glm/glm.hpp>

/**
 * @brief Compile-time optimized logging system
 * 
 * Eliminates runtime overhead in release builds by using constexpr templates.
 * Logging calls are completely eliminated by the compiler when logging is disabled.
 */
template<bool EnableLogging = true>
class OptimizedLogger {
public:
    static OptimizedLogger& instance() {
        static OptimizedLogger logger;
        return logger;
    }
    
    void initialize(const std::string& csv_filepath, const std::string& log_filepath, bool enable_runtime_logging = true) {
        if constexpr (EnableLogging) {
            runtime_logging_enabled_ = enable_runtime_logging;
            
            if (!enable_runtime_logging) {
                return; // Don't initialize if runtime logging disabled
            }
            
            // Initialize CSV file for photon traces
            if (csv_file_.is_open()) {
                csv_file_.close();
            }
            csv_file_.open(csv_filepath);
            if (csv_file_.is_open()) {
                csv_file_ << "PhotonID,Event,PosX,PosY,PosZ,DirX,DirY,DirZ,Weight,VoxelX,VoxelY,VoxelZ,MediumID,Energy,Description\n";
            }
            
            // Initialize log file for general messages
            if (log_file_.is_open()) {
                log_file_.close();
            }
            log_file_.open(log_filepath);
            if (log_file_.is_open()) {
                log_file_ << "=== Photron Debug Log ===\n";
                log_file_ << "Started at: " << get_timestamp() << "\n\n";
            }
        }
    }
    
    void log_info(const std::string& message) {
        if constexpr (EnableLogging) {
            if (runtime_logging_enabled_ && log_file_.is_open()) {
                log_file_ << "[" << get_timestamp() << "] INFO: " << message << "\n";
                log_file_.flush();
            }
        }
        // Compiler completely eliminates this method call when EnableLogging=false
    }
    
    void log_debug(const std::string& message) {
        if constexpr (EnableLogging) {
            if (runtime_logging_enabled_ && log_file_.is_open()) {
                log_file_ << "[" << get_timestamp() << "] DEBUG: " << message << "\n";
                log_file_.flush();
            }
        }
    }
    
    void log_warning(const std::string& message) {
        if constexpr (EnableLogging) {
            if (runtime_logging_enabled_ && log_file_.is_open()) {
                log_file_ << "[" << get_timestamp() << "] WARN: " << message << "\n";
                log_file_.flush();
            }
        }
    }
    
    void log_error(const std::string& message) {
        if constexpr (EnableLogging) {
            if (runtime_logging_enabled_ && log_file_.is_open()) {
                log_file_ << "[" << get_timestamp() << "] ERROR: " << message << "\n";
                log_file_.flush();
            }
        }
    }
    
    void log_photon_event(int photon_id, const std::string& event, 
                         const glm::dvec3& position, const glm::dvec3& direction, 
                         double weight, const glm::ivec3& voxel_coords = glm::ivec3(-1),
                         int medium_id = -1, double energy = 0.0,
                         const std::string& description = "") {
        if constexpr (EnableLogging) {
            if (runtime_logging_enabled_ && csv_file_.is_open()) {
                csv_file_ << photon_id << "," << event << ","
                          << position.x << "," << position.y << "," << position.z << ","
                          << direction.x << "," << direction.y << "," << direction.z << ","
                          << weight << ","
                          << voxel_coords.x << "," << voxel_coords.y << "," << voxel_coords.z << ","
                          << medium_id << "," << energy << "," << description << "\n";
                csv_file_.flush();
            }
        }
    }
    
    ~OptimizedLogger() {
        if constexpr (EnableLogging) {
            if (csv_file_.is_open()) {
                csv_file_.close();
            }
            if (log_file_.is_open()) {
                log_file_ << "\n=== Log ended at: " << get_timestamp() << " ===\n";
                log_file_.close();
            }
        }
    }
    
private:
    OptimizedLogger() = default;
    
    std::string get_timestamp() const {
        if constexpr (EnableLogging) {
            auto now = std::chrono::system_clock::now();
            auto time_t = std::chrono::system_clock::to_time_t(now);
            auto ms = std::chrono::duration_cast<std::chrono::milliseconds>(now.time_since_epoch()) % 1000;
            
            std::ostringstream oss;
            oss << std::put_time(std::localtime(&time_t), "%H:%M:%S");
            oss << '.' << std::setfill('0') << std::setw(3) << ms.count();
            return oss.str();
        }
        else {
            return "";
        }
    }
    
    std::ofstream csv_file_;
    std::ofstream log_file_;
    bool runtime_logging_enabled_ = false;
};

// Type aliases for different build configurations
#ifdef _DEBUG
    using FastLogger = OptimizedLogger<true>;
#else
    using FastLogger = OptimizedLogger<false>; // Zero runtime overhead in release
#endif

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
