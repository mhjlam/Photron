#include "logger.hpp"

Logger& Logger::instance() {
    static Logger logger;
    return logger;
}

void Logger::initialize(const std::string& csv_filepath, const std::string& log_filepath, bool enable_logging) {
#ifdef _DEBUG
    logging_enabled_ = enable_logging;
    
    if (!enable_logging) {
        return;  // Don't initialize if logging is disabled
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
#endif
}

void Logger::log_photon_event(int photon_id, const std::string& event, 
                             const glm::dvec3& position, const glm::dvec3& direction, 
                             double weight, const glm::ivec3& voxel_coords,
                             int medium_id, double energy,
                             const std::string& description) {
#ifdef _DEBUG
    if (logging_enabled_ && csv_file_.is_open()) {
        csv_file_ << photon_id << "," << event << ","
                  << position.x << "," << position.y << "," << position.z << ","
                  << direction.x << "," << direction.y << "," << direction.z << ","
                  << weight << ","
                  << voxel_coords.x << "," << voxel_coords.y << "," << voxel_coords.z << ","
                  << medium_id << "," << energy << "," << description << "\n";
        csv_file_.flush();
    }
#endif
}

void Logger::log_voxel_emittance(int photon_id, const glm::dvec3& position, const glm::dvec3& direction,
                               double weight, const glm::ivec3& voxel_coords, double emittance, 
                               const std::string& surface_type) {
#ifdef _DEBUG
    if (logging_enabled_ && csv_file_.is_open()) {
        csv_file_ << photon_id << ",VOXEL_EMITTANCE,"
                  << position.x << "," << position.y << "," << position.z << ","
                  << direction.x << "," << direction.y << "," << direction.z << ","
                  << weight << ","
                  << voxel_coords.x << "," << voxel_coords.y << "," << voxel_coords.z << ","
                  << "-1," << emittance << "," << surface_type << "\n";
        csv_file_.flush();
    }
#endif
}

void Logger::log_info(const std::string& message) {
#ifdef _DEBUG
    if (logging_enabled_ && log_file_.is_open()) {
        log_file_ << "[" << get_timestamp() << "] INFO: " << message << "\n";
        log_file_.flush();
    }
#endif
}

void Logger::log_debug(const std::string& message) {
#ifdef _DEBUG
    if (logging_enabled_ && log_file_.is_open()) {
        log_file_ << "[" << get_timestamp() << "] DEBUG: " << message << "\n";
        log_file_.flush();
    }
#endif
}

void Logger::log_warning(const std::string& message) {
#ifdef _DEBUG
    if (logging_enabled_ && log_file_.is_open()) {
        log_file_ << "[" << get_timestamp() << "] WARN: " << message << "\n";
        log_file_.flush();
    }
#endif
}

void Logger::log_error(const std::string& message) {
#ifdef _DEBUG
    if (logging_enabled_ && log_file_.is_open()) {
        log_file_ << "[" << get_timestamp() << "] ERROR: " << message << "\n";
        log_file_.flush();
    }
#endif
}

void Logger::log_info_fmt(const std::string& message) {
    log_info(message);
}

void Logger::log_debug_fmt(const std::string& message) {
    log_debug(message);
}

void Logger::log_warning_fmt(const std::string& message) {
    log_warning(message);
}

void Logger::log_error_fmt(const std::string& message) {
    log_error(message);
}

Logger::~Logger() {
#ifdef _DEBUG
    if (csv_file_.is_open()) {
        csv_file_.close();
    }
    if (log_file_.is_open()) {
        log_file_ << "\n=== Log ended at: " << get_timestamp() << " ===\n";
        log_file_.close();
    }
#endif
}

std::string Logger::get_timestamp() const {
    auto now = std::chrono::system_clock::now();
    auto time_t = std::chrono::system_clock::to_time_t(now);
    auto ms = std::chrono::duration_cast<std::chrono::milliseconds>(
        now.time_since_epoch()) % 1000;
    
    std::ostringstream oss;
    oss << std::put_time(std::localtime(&time_t), "%H:%M:%S");
    oss << '.' << std::setfill('0') << std::setw(3) << ms.count();
    return oss.str();
}