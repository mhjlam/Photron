#include "shader_manager.hpp"

#include <fstream>
#include <iostream>
#include <sstream>

ShaderManager::ShaderManager() = default;

ShaderManager::~ShaderManager() {
	shutdown();
}

bool ShaderManager::initialize() {
	// Load basic shaders for wireframe rendering
	if (!load_shader("lines", "shaders/lines.vert", "shaders/lines.frag")) {
		std::cerr << "Failed to load lines shader" << std::endl;
		return false;
	}

	if (!load_shader("points", "shaders/points.vert", "shaders/points.frag")) {
		std::cerr << "Failed to load points shader" << std::endl;
		return false;
	}

	return true;
}

void ShaderManager::shutdown() {
	for (auto& [name, program] : shader_programs_) {
		glDeleteProgram(program);
	}
	shader_programs_.clear();
}

bool ShaderManager::load_shader(const std::string& name, const std::string& vertex_path,
								const std::string& fragment_path) {
	std::string vertex_source = load_shader_source(vertex_path);
	std::string fragment_source = load_shader_source(fragment_path);

	if (vertex_source.empty() || fragment_source.empty()) {
		return false;
	}

	GLuint program = create_program(vertex_source, fragment_source);
	if (program == 0) {
		return false;
	}

	shader_programs_[name] = program;
	return true;
}

void ShaderManager::use_shader(const std::string& name) {
	auto it = shader_programs_.find(name);
	if (it != shader_programs_.end()) {
		glUseProgram(it->second);
	}
	else {
		std::cerr << "Shader not found: " << name << std::endl;
	}
}

GLuint ShaderManager::get_shader_program(const std::string& name) {
	auto it = shader_programs_.find(name);
	return (it != shader_programs_.end()) ? it->second : 0;
}

void ShaderManager::set_uniform_mat4(const std::string& shader_name, const std::string& uniform_name,
									 const glm::mat4& matrix) {
	GLint location = get_uniform_location(shader_name, uniform_name);
	if (location != -1) {
		use_shader(shader_name);
		glUniformMatrix4fv(location, 1, GL_FALSE, &matrix[0][0]);
	}
}

void ShaderManager::set_uniform_vec3(const std::string& shader_name, const std::string& uniform_name,
									 const glm::vec3& vector) {
	GLint location = get_uniform_location(shader_name, uniform_name);
	if (location != -1) {
		use_shader(shader_name);
		glUniform3fv(location, 1, &vector[0]);
	}
}

void ShaderManager::set_uniform_vec4(const std::string& shader_name, const std::string& uniform_name,
									 const glm::vec4& vector) {
	GLint location = get_uniform_location(shader_name, uniform_name);
	if (location != -1) {
		use_shader(shader_name);
		glUniform4fv(location, 1, &vector[0]);
	}
}

void ShaderManager::set_uniform_float(const std::string& shader_name, const std::string& uniform_name, float value) {
	GLint location = get_uniform_location(shader_name, uniform_name);
	if (location != -1) {
		use_shader(shader_name);
		glUniform1f(location, value);
	}
}

void ShaderManager::set_uniform_int(const std::string& shader_name, const std::string& uniform_name, int value) {
	GLint location = get_uniform_location(shader_name, uniform_name);
	if (location != -1) {
		use_shader(shader_name);
		glUniform1i(location, value);
	}
}

GLuint ShaderManager::compile_shader(const std::string& source, GLenum shader_type) {
	GLuint shader = glCreateShader(shader_type);
	const char* source_cstr = source.c_str();
	glShaderSource(shader, 1, &source_cstr, nullptr);
	glCompileShader(shader);

	GLint success;
	glGetShaderiv(shader, GL_COMPILE_STATUS, &success);
	if (!success) {
		char info_log[512];
		glGetShaderInfoLog(shader, 512, nullptr, info_log);
		std::cerr << "Shader compilation failed: " << info_log << std::endl;
		glDeleteShader(shader);
		return 0;
	}

	return shader;
}

GLuint ShaderManager::create_program(const std::string& vertex_source, const std::string& fragment_source) {
	GLuint vertex_shader = compile_shader(vertex_source, GL_VERTEX_SHADER);
	GLuint fragment_shader = compile_shader(fragment_source, GL_FRAGMENT_SHADER);

	if (vertex_shader == 0 || fragment_shader == 0) {
		if (vertex_shader)
			glDeleteShader(vertex_shader);
		if (fragment_shader)
			glDeleteShader(fragment_shader);
		return 0;
	}

	GLuint program = glCreateProgram();
	glAttachShader(program, vertex_shader);
	glAttachShader(program, fragment_shader);
	glLinkProgram(program);

	GLint success;
	glGetProgramiv(program, GL_LINK_STATUS, &success);
	if (!success) {
		char info_log[512];
		glGetProgramInfoLog(program, 512, nullptr, info_log);
		std::cerr << "Program linking failed: " << info_log << std::endl;
		glDeleteProgram(program);
		program = 0;
	}

	glDeleteShader(vertex_shader);
	glDeleteShader(fragment_shader);

	return program;
}

std::string ShaderManager::load_shader_source(const std::string& file_path) {
	std::ifstream file(file_path);
	if (!file.is_open()) {
		std::cerr << "Failed to open shader file: " << file_path << std::endl;
		return "";
	}

	std::stringstream buffer;
	buffer << file.rdbuf();
	return buffer.str();
}

GLint ShaderManager::get_uniform_location(const std::string& shader_name, const std::string& uniform_name) {
	auto it = shader_programs_.find(shader_name);
	if (it == shader_programs_.end()) {
		std::cerr << "Shader not found: " << shader_name << std::endl;
		return -1;
	}

	return glGetUniformLocation(it->second, uniform_name.c_str());
}
