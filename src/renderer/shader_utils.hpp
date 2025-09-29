#pragma once

#include <string>
#include <GL/glew.h>

namespace ShaderUtils {
    /**
     * Compile a shader from source code
     * @param source The GLSL source code
     * @param shader_type GL_VERTEX_SHADER, GL_FRAGMENT_SHADER, etc.
     * @return OpenGL shader handle, or 0 on failure
     */
    GLuint compile_shader(const std::string& source, GLenum shader_type);

    /**
     * Load shader source from a file
     * @param file_path Path to shader file
     * @return Shader source code, or empty string on failure
     */
    std::string load_shader_source(const std::string& file_path);
}