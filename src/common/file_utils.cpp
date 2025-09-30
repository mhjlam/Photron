/**
 * @file file_utils.cpp
 * @brief Implementation of standardized file operations with robust error handling
 *
 * Provides consistent file I/O operations used throughout the Photron application
 * for simulation results export, debug logging, and configuration file management.
 */

#include "file_utils.hpp"

namespace FileUtils
{

std::ofstream create_output_file(const std::string& filepath) {
	// Ensure directory structure exists before file creation
	ensure_output_directory(filepath);

	// Create file stream with error checking
	std::ofstream ofs(filepath);
	if (!ofs.is_open()) {
		std::cerr << "Error: Could not open file for writing: " << filepath << std::endl;
	}

	return ofs;
}

void ensure_output_directory(const std::string& filepath) {
	// Extract directory path from complete file path
	std::filesystem::path file_path(filepath);
	std::filesystem::path dir_path = file_path.parent_path();

	// Create directory hierarchy if it doesn't exist
	if (!dir_path.empty() && !std::filesystem::exists(dir_path)) {
		std::filesystem::create_directories(dir_path);
	}
}

void write_file_header(std::ofstream& ofs, const std::string& title) {
	// Set numerical precision for consistent output formatting
	ofs.precision(8);

	// Write standardized header with decorative border
	ofs << "################################################################" << std::endl;
	ofs << "# " << title << std::endl;
	ofs << "################################################################" << std::endl;
	ofs << std::endl;
}
} // namespace FileUtils