/**
 * @file logger.cpp
 * @brief Implementation of debug logging system with compile-time optimization
 *
 * Implements the Logger singleton providing structured debug output for
 * Monte Carlo photon transport analysis. Features conditional compilation
 * for zero-overhead release builds and dual-format output (CSV for data,
 * text for messages).
 */

#include "logger.hpp"

#include <glm/glm.hpp>

#include "common/config.hpp"

Logger& Logger::instance() {
	// Singleton pattern with thread-safe initialization
	static Logger logger;
	return logger;
}

void Logger::initialize(const std::string& csv_filepath, const std::string& log_filepath, bool enable_logging) {
#ifdef _DEBUG
	// Set up debug logging (compile-time optimized for release builds)
	logging_enabled_ = enable_logging;

	if (!enable_logging) {
		return; // Early exit if logging disabled
	}

	// Set up CSV file for structured photon event data
	if (csv_file_.is_open()) {
		csv_file_.close();
	}

	csv_file_.open(csv_filepath);
	if (csv_file_.is_open()) {
		csv_file_
			<< "PhotonID,Event,PosX,PosY,PosZ,DirX,DirY,DirZ,Weight,VoxelX,VoxelY,VoxelZ,MediumID,Energy,Description\n";
	}

	// Initialize text log file for general debug messages
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

void Logger::log_photon_event(int photon_id,
							  const std::string& event,
							  const glm::dvec3& position,
							  const glm::dvec3& direction,
							  double weight,
							  const glm::ivec3& voxel_coords,
							  int medium_id,
							  double energy,
							  const std::string& description) {
#ifdef _DEBUG
	// Write event data to CSV for post-processing analysis
	if (logging_enabled_ && csv_file_.is_open()) {
		csv_file_ << photon_id << "," << event << "," << position.x << "," << position.y << "," << position.z << ","
				  << direction.x << "," << direction.y << "," << direction.z << "," << weight << "," << voxel_coords.x
				  << "," << voxel_coords.y << "," << voxel_coords.z << "," << medium_id << "," << energy << ","
				  << description << "\n";
		csv_file_.flush();
	}
#endif
}

void Logger::log_voxel_emittance(int photon_id,
								 const glm::dvec3& position,
								 const glm::dvec3& direction,
								 double weight,
								 const glm::ivec3& voxel_coords,
								 double emittance,
								 const std::string& surface_type) {
#ifdef _DEBUG
	// Record voxel emittance for energy conservation tracking
	if (logging_enabled_ && csv_file_.is_open()) {
		csv_file_ << photon_id << ",VOXEL_EMITTANCE," << position.x << "," << position.y << "," << position.z << ","
				  << direction.x << "," << direction.y << "," << direction.z << "," << weight << "," << voxel_coords.x
				  << "," << voxel_coords.y << "," << voxel_coords.z << ","
				  << "-1," << emittance << "," << surface_type << "\n";
		csv_file_.flush();
	}
#endif
}

void Logger::log_info(const std::string& message) {
#ifdef _DEBUG
	// Log informational messages with timestamp
	if (logging_enabled_ && log_file_.is_open()) {
		log_file_ << "[" << get_timestamp() << "] INFO: " << message << "\n";
		log_file_.flush();
	}
#endif
}

void Logger::log_debug(const std::string& message) {
#ifdef _DEBUG
	// Log detailed debug information for troubleshooting
	if (logging_enabled_ && log_file_.is_open()) {
		log_file_ << "[" << get_timestamp() << "] DEBUG: " << message << "\n";
		log_file_.flush();
	}
#endif
}

void Logger::log_warning(const std::string& message) {
#ifdef _DEBUG
	// Log warning messages for potential issues
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
	auto ms = std::chrono::duration_cast<std::chrono::milliseconds>(now.time_since_epoch()) % 1000;

	std::ostringstream oss;
	oss << std::put_time(std::localtime(&time_t), "%H:%M:%S");
	oss << '.' << std::setfill('0') << std::setw(3) << ms.count();
	return oss.str();
}
