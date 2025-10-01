/**
 * @file shader.cpp
 * @brief Implementation of OpenGL shader program management
 *
 * Implements the Shader class for loading, compiling, and managing OpenGL
 * shader programs. Provides type-safe uniform variable setting and robust
 * error handling for shader compilation and linking.
 */

#include "shader.hpp"

#include <array>
#include <fstream>
#include <iostream>
#include <sstream>

#include <glm/gtc/type_ptr.hpp>

Shader::~Shader() {
	if (program_id_ != 0) {
		glDeleteProgram(program_id_);
	}
}

bool Shader::load_from_file(const std::string& vertex_path, const std::string& fragment_path) {
	std::string vertex_source = load_shader_source(vertex_path);
	std::string fragment_source = load_shader_source(fragment_path);

	if (vertex_source.empty() || fragment_source.empty()) {
		return false;
	}

	return load_from_source(vertex_source, fragment_source);
}

bool Shader::load_from_source(const std::string& vertex_source, const std::string& fragment_source) {
	GLuint vertex_shader = compile_shader(vertex_source, GL_VERTEX_SHADER);
	GLuint fragment_shader = compile_shader(fragment_source, GL_FRAGMENT_SHADER);

	if (vertex_shader == 0 || fragment_shader == 0) {
		if (vertex_shader) {
			glDeleteShader(vertex_shader);
		}
		if (fragment_shader) {
			glDeleteShader(fragment_shader);
		}
		return false;
	}

	bool success = link_program(vertex_shader, fragment_shader);

	glDeleteShader(vertex_shader);
	glDeleteShader(fragment_shader);

	return success;
}

void Shader::use() const {
	if (program_id_ != 0) {
		glUseProgram(program_id_);
	}
}

void Shader::unuse() const {
	glUseProgram(0);
}

void Shader::set_uniform(const std::string& name, bool value) {
	glUniform1i(get_uniform_location(name), static_cast<int>(value));
}

void Shader::set_uniform(const std::string& name, int value) {
	glUniform1i(get_uniform_location(name), value);
}

void Shader::set_uniform(const std::string& name, float value) {
	glUniform1f(get_uniform_location(name), value);
}

void Shader::set_uniform(const std::string& name, const glm::vec2& value) {
	glUniform2fv(get_uniform_location(name), 1, &value[0]);
}

void Shader::set_uniform(const std::string& name, const glm::vec3& value) {
	glUniform3fv(get_uniform_location(name), 1, &value[0]);
}

void Shader::set_uniform(const std::string& name, const glm::vec4& value) {
	glUniform4fv(get_uniform_location(name), 1, &value[0]);
}

void Shader::set_uniform(const std::string& name, const glm::mat3& value) {
	glUniformMatrix3fv(get_uniform_location(name), 1, GL_FALSE, glm::value_ptr(value));
}

void Shader::set_uniform(const std::string& name, const glm::mat4& value) {
	glUniformMatrix4fv(get_uniform_location(name), 1, GL_FALSE, glm::value_ptr(value));
}

std::string Shader::load_shader_source(const std::string& file_path) {
	// Open shader file for reading
	std::ifstream file(file_path);
	if (!file.is_open()) {
		std::cerr << "Failed to open shader file: " << file_path << std::endl;
		return "";
	}

	// Read entire file content into string
	std::stringstream buffer;
	buffer << file.rdbuf();
	return buffer.str();
}

std::string Shader::read_file(const std::string& path) {
	return load_shader_source(path);
}

GLuint Shader::compile_shader(const std::string& source, GLenum shader_type) {
	// Create shader object for specified type
	GLuint shader = glCreateShader(shader_type);
	const char* source_cstr = source.c_str();

	// Load source code and compile shader
	glShaderSource(shader, 1, &source_cstr, nullptr);
	glCompileShader(shader);

	// Check compilation status and handle errors
	GLint success;
	glGetShaderiv(shader, GL_COMPILE_STATUS, &success);
	if (!success) {
		// Extract and report compilation error details
		std::array<char, 512> info_log {};
		glGetShaderInfoLog(shader, static_cast<GLsizei>(info_log.size()), nullptr, info_log.data());
		std::cerr << "Shader compilation failed: " << info_log.data() << std::endl;

		// Clean up failed shader object
		glDeleteShader(shader);
		return 0;
	}

	return shader;
}

bool Shader::link_program(GLuint vertex_shader, GLuint fragment_shader) {
	program_id_ = glCreateProgram();

	glAttachShader(program_id_, vertex_shader);
	glAttachShader(program_id_, fragment_shader);
	glLinkProgram(program_id_);

	GLint success;
	glGetProgramiv(program_id_, GL_LINK_STATUS, &success);

	if (!success) {
		std::array<GLchar, 512> info_log {};
		glGetProgramInfoLog(program_id_, static_cast<GLsizei>(info_log.size()), nullptr, info_log.data());
		std::cerr << "Error: Program linking failed: " << info_log.data() << std::endl;
		glDeleteProgram(program_id_);
		program_id_ = 0;
		return false;
	}

	return true;
}

GLint Shader::get_uniform_location(const std::string& name) {
	if (program_id_ == 0) {
		return -1;
	}

	// Check cache first to avoid repeated glGetUniformLocation calls
	auto it = uniform_location_cache_.find(name);
	if (it != uniform_location_cache_.end()) {
		return it->second;
	}

	// Cache miss - query OpenGL and store result
	GLint location = glGetUniformLocation(program_id_, name.c_str());
	uniform_location_cache_[name] = location;
	return location;
}
