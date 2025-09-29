#pragma once

#include <string>
#include <vector>
#include <functional>
#include <iostream>
#include <GL/glew.h>
#include "shader_utils.hpp"

namespace RenderSetup {
    /**
     * Configuration structure for basic shader setup
     */
    struct ShaderConfig {
        std::string vertex_shader_path;
        std::string fragment_shader_path;
        std::string error_name; // For error messages
    };

    /**
     * Configuration structure for uniform caching
     */
    struct UniformConfig {
        std::string name;
        GLint* location_ptr; // Pointer to store the uniform location
    };

    /**
     * Configuration structure for vertex attribute setup
     */
    struct VertexAttributeConfig {
        GLuint location;
        GLint size;
        GLenum type;
        GLboolean normalized;
        GLsizei stride;
        const void* pointer;
        GLuint divisor{0}; // 0 = per vertex, 1 = per instance
    };

    /**
     * Base geometry setup function type for instanced rendering
     */
    using GeometrySetupFunction = std::function<bool(GLuint vbo)>;

    /**
     * Template function for basic rendering setup (lines, points, triangles)
     * Handles: shader loading, program creation, uniform caching, VAO/VBO setup, vertex attributes
     */
    template<typename VertexType>
    inline bool setup_basic_rendering(
        const ShaderConfig& shader_config,
        const std::vector<UniformConfig>& uniforms,
        const std::vector<VertexAttributeConfig>& vertex_attributes,
        GLuint& shader_program,
        GLuint& vao,
        GLuint& vbo
    ) {
        // Load shaders
        std::string vertex_source = ShaderUtils::load_shader_source(shader_config.vertex_shader_path);
        std::string fragment_source = ShaderUtils::load_shader_source(shader_config.fragment_shader_path);

        if (vertex_source.empty() || fragment_source.empty()) {
            std::cerr << "Failed to load " << shader_config.error_name << " shaders" << std::endl;
            return false;
        }

        // Create shader program
        GLuint vertex_shader = ShaderUtils::compile_shader(vertex_source, GL_VERTEX_SHADER);
        GLuint fragment_shader = ShaderUtils::compile_shader(fragment_source, GL_FRAGMENT_SHADER);

        if (vertex_shader == 0 || fragment_shader == 0) {
            if (vertex_shader) glDeleteShader(vertex_shader);
            if (fragment_shader) glDeleteShader(fragment_shader);
            return false;
        }

        shader_program = glCreateProgram();
        glAttachShader(shader_program, vertex_shader);
        glAttachShader(shader_program, fragment_shader);
        glLinkProgram(shader_program);

        GLint success;
        glGetProgramiv(shader_program, GL_LINK_STATUS, &success);
        if (!success) {
            char info_log[512];
            glGetProgramInfoLog(shader_program, 512, nullptr, info_log);
            std::cerr << "Program linking failed for " << shader_config.error_name << ": " << info_log << std::endl;
            glDeleteProgram(shader_program);
            shader_program = 0;
            glDeleteShader(vertex_shader);
            glDeleteShader(fragment_shader);
            return false;
        }

        glDeleteShader(vertex_shader);
        glDeleteShader(fragment_shader);

        // Cache uniform locations
        for (const auto& uniform : uniforms) {
            if (uniform.location_ptr) {
                *uniform.location_ptr = glGetUniformLocation(shader_program, uniform.name.c_str());
            }
        }

        // Create VAO and VBO
        glGenVertexArrays(1, &vao);
        glGenBuffers(1, &vbo);

        glBindVertexArray(vao);
        glBindBuffer(GL_ARRAY_BUFFER, vbo);

        // Setup vertex attributes
        for (const auto& attr : vertex_attributes) {
            glVertexAttribPointer(attr.location, attr.size, attr.type, attr.normalized, attr.stride, attr.pointer);
            glEnableVertexAttribArray(attr.location);
            if (attr.divisor > 0) {
                glVertexAttribDivisor(attr.location, attr.divisor);
            }
        }

        glBindVertexArray(0);
        return true;
    }

    /**
     * Template function for instanced rendering setup (point instanced, line instanced)
     * Handles: shader loading, base geometry, instance attributes
     */
    template<typename InstanceType>
    inline bool setup_instanced_rendering(
        const ShaderConfig& shader_config,
        const std::vector<UniformConfig>& uniforms,
        const GeometrySetupFunction& geometry_setup,
        const std::vector<VertexAttributeConfig>& instance_attributes,
        GLuint& shader_program,
        GLuint& vao,
        GLuint& base_vbo,
        GLuint& instance_vbo
    ) {
        // Load shaders
        std::string vertex_source = ShaderUtils::load_shader_source(shader_config.vertex_shader_path);
        std::string fragment_source = ShaderUtils::load_shader_source(shader_config.fragment_shader_path);

        if (vertex_source.empty() || fragment_source.empty()) {
            std::cerr << "Failed to load " << shader_config.error_name << " shaders" << std::endl;
            return false;
        }

        // Create shader program
        GLuint vertex_shader = ShaderUtils::compile_shader(vertex_source, GL_VERTEX_SHADER);
        GLuint fragment_shader = ShaderUtils::compile_shader(fragment_source, GL_FRAGMENT_SHADER);

        if (vertex_shader == 0 || fragment_shader == 0) {
            if (vertex_shader) glDeleteShader(vertex_shader);
            if (fragment_shader) glDeleteShader(fragment_shader);
            return false;
        }

        shader_program = glCreateProgram();
        glAttachShader(shader_program, vertex_shader);
        glAttachShader(shader_program, fragment_shader);
        glLinkProgram(shader_program);

        GLint success;
        glGetProgramiv(shader_program, GL_LINK_STATUS, &success);
        if (!success) {
            char info_log[512];
            glGetProgramInfoLog(shader_program, 512, nullptr, info_log);
            std::cerr << "Program linking failed for " << shader_config.error_name << ": " << info_log << std::endl;
            glDeleteProgram(shader_program);
            shader_program = 0;
            glDeleteShader(vertex_shader);
            glDeleteShader(fragment_shader);
            return false;
        }

        glDeleteShader(vertex_shader);
        glDeleteShader(fragment_shader);

        // Cache uniform locations
        for (const auto& uniform : uniforms) {
            if (uniform.location_ptr) {
                *uniform.location_ptr = glGetUniformLocation(shader_program, uniform.name.c_str());
            }
        }

        // Create VAO and VBOs
        glGenVertexArrays(1, &vao);
        glGenBuffers(1, &base_vbo);
        glGenBuffers(1, &instance_vbo);

        glBindVertexArray(vao);

        // Setup base geometry using the provided function
        if (!geometry_setup(base_vbo)) {
            std::cerr << "Failed to setup base geometry for " << shader_config.error_name << std::endl;
            glBindVertexArray(0);
            return false;
        }

        // Setup instance data buffer
        glBindBuffer(GL_ARRAY_BUFFER, instance_vbo);

        // Setup instance attributes
        for (const auto& attr : instance_attributes) {
            glVertexAttribPointer(attr.location, attr.size, attr.type, attr.normalized, attr.stride, attr.pointer);
            glEnableVertexAttribArray(attr.location);
            if (attr.divisor > 0) {
                glVertexAttribDivisor(attr.location, attr.divisor);
            }
        }

        glBindVertexArray(0);
        return true;
    }
}