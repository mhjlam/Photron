#include "shader.hpp"

#include <glm/gtc/type_ptr.hpp>

#include <fstream>
#include <iostream>
#include <sstream>

Shader::~Shader() {
	if (program_id_ != 0) {
		glDeleteProgram(program_id_);
	}
}

bool Shader::load_from_file(const std::string& vertex_path, const std::string& fragment_path) {
	std::string vertex_source = read_file(vertex_path);
	std::string fragment_source = read_file(fragment_path);

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

std::string Shader::read_file(const std::string& path) {
	std::ifstream file(path);
	if (!file.is_open()) {
		std::cerr << "Error: Cannot open shader file: " << path << std::endl;
		return "";
	}

	std::stringstream stream;
	stream << file.rdbuf();
	file.close();

	return stream.str();
}

GLuint Shader::compile_shader(const std::string& source, GLenum type) {
	GLuint shader = glCreateShader(type);
	const char* source_cstr = source.c_str();

	glShaderSource(shader, 1, &source_cstr, nullptr);
	glCompileShader(shader);

	GLint success;
	glGetShaderiv(shader, GL_COMPILE_STATUS, &success);

	if (!success) {
		GLchar info_log[512];
		glGetShaderInfoLog(shader, sizeof(info_log), nullptr, info_log);
		std::cerr << "Error: Shader compilation failed: " << info_log << std::endl;
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
		GLchar info_log[512];
		glGetProgramInfoLog(program_id_, sizeof(info_log), nullptr, info_log);
		std::cerr << "Error: Program linking failed: " << info_log << std::endl;
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
	return glGetUniformLocation(program_id_, name.c_str());
}
