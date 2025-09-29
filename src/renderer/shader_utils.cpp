#include "shader_utils.hpp"

#include <iostream>
#include <fstream>
#include <sstream>
#include <array>

namespace ShaderUtils {
    
    GLuint compile_shader(const std::string& source, GLenum shader_type) {
        GLuint shader = glCreateShader(shader_type);
        const char* source_cstr = source.c_str();
        glShaderSource(shader, 1, &source_cstr, nullptr);
        glCompileShader(shader);

        GLint success;
        glGetShaderiv(shader, GL_COMPILE_STATUS, &success);
        if (!success) {
            std::array<char, 512> info_log{};
            glGetShaderInfoLog(shader, static_cast<GLsizei>(info_log.size()), nullptr, info_log.data());
            std::cerr << "Shader compilation failed: " << info_log.data() << std::endl;
            glDeleteShader(shader);
            return 0;
        }

        return shader;
    }

    std::string load_shader_source(const std::string& file_path) {
        std::ifstream file(file_path);
        if (!file.is_open()) {
            std::cerr << "Failed to open shader file: " << file_path << std::endl;
            return "";
        }

        std::stringstream buffer;
        buffer << file.rdbuf();
        return buffer.str();
    }
    
}