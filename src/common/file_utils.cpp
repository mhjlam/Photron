#include "file_utils.hpp"

namespace FileUtils {
    
    std::ofstream create_output_file(const std::string& filepath) {
        // Ensure the output directory exists
        ensure_output_directory(filepath);
        
        // Create and open the file
        std::ofstream ofs(filepath);
        if (!ofs.is_open()) {
            std::cerr << "Error: Could not open file for writing: " << filepath << std::endl;
        }
        
        return ofs;
    }
    
    void ensure_output_directory(const std::string& filepath) {
        std::filesystem::path file_path(filepath);
        std::filesystem::path dir_path = file_path.parent_path();
        
        if (!dir_path.empty() && !std::filesystem::exists(dir_path)) {
            std::filesystem::create_directories(dir_path);
        }
    }
    
    void write_file_header(std::ofstream& ofs, const std::string& title) {
        ofs.precision(8);
        ofs << "################################################################" << std::endl;
        ofs << "# " << title << std::endl;
        ofs << "################################################################" << std::endl;
        ofs << std::endl;
    }
}