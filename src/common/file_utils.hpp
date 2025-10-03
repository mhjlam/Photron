/**
 * @file file_utils.hpp
 * @brief File system utilities for robust output file management
 *
 * Provides standardized file I/O operations with automatic directory creation,
 * error handling, and consistent formatting. Used throughout the application
 * for simulation results export and debug output.
 */

#pragma once

#include <filesystem>
#include <fstream>
#include <iostream>
#include <string>

/**
 * @namespace FileUtils
 * @brief Standardized file operations with robust error handling
 *
 * The FileUtils namespace provides a consistent interface for file operations
 * throughout the Photron application. All functions handle common edge cases
 * such as missing directories, permission issues, and provide standardized
 * error reporting.
 *
 * **Key Features:**
 * - Automatic directory creation for output paths
 * - Consistent file header formatting
 * - Robust error handling with informative messages
 * - Cross-platform compatibility
 *
 * **Usage Pattern:**
 * ```cpp
 * auto ofs = FileUtils::create_output_file("results/simulation_data.txt");
 * if (ofs.is_open()) {
 *     FileUtils::write_file_header(ofs, "Simulation Results");
 *     // ... write data
 * }
 * ```
 */
namespace FileUtils
{
/**
 * @brief Create output file with automatic directory creation
 *
 * Creates a standardized output file stream with automatic parent directory
 * creation. Handles path validation and provides error reporting if file
 * creation fails.
 *
 * @param filepath Complete file path for the output file
 * @return std::ofstream Opened file stream (caller should check is_open())
 *
 * @note If parent directories don't exist, they will be created automatically.
 *       If file creation fails, error message is written to stderr.
 */
std::ofstream create_output_file(const std::string& filepath);

/**
 * @brief Ensure output directory structure exists
 *
 * Creates all parent directories for the specified file path if they don't
 * already exist. Safe to call multiple times - will not fail if directories
 * already exist.
 *
 * @param filepath File path whose parent directories should be created
 *
 * @note Uses std::filesystem for cross-platform directory creation.
 */
void ensure_output_directory(const std::string& filepath);

/**
 * @brief Write standardized file header with decorative border
 *
 * Creates a consistent header format used across all output files,
 * including simulation results, debug logs, and CSV exports.
 *
 * @param ofs Output file stream to write header to
 * @param title Header title text to display
 *
 * @note Sets precision to 8 decimal places for numerical output consistency.
 */
void write_file_header(std::ofstream& ofs, const std::string& title);
} // namespace FileUtils
