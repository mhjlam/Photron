#pragma once
#include <string>
#include <fstream>
#include <filesystem>
#include <iostream>

namespace FileUtils {
    /**
     * Create a standardized output file with directory creation and error handling
     * @param filepath The full file path
     * @return An opened ofstream, caller should check is_open()
     */
    std::ofstream create_output_file(const std::string& filepath);
    
    /**
     * Ensure output directory exists for the given file path
     * @param filepath The file path to create directories for
     */
    void ensure_output_directory(const std::string& filepath);
    
    /**
     * Write a standardized file header with title and borders
     * @param ofs The output stream to write to
     * @param title The header title
     */
    void write_file_header(std::ofstream& ofs, const std::string& title);
}