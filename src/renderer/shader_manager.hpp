#pragma once

#include <string>
#include <unordered_map>

#include <GL/glew.h>
#include <glm/glm.hpp>

class ShaderManager
{
public:
	ShaderManager();
	~ShaderManager();

	bool initialize();
	void shutdown();

	bool load_shader(const std::string& name, const std::string& vertex_path, const std::string& fragment_path);
	void use_shader(const std::string& name);
	GLuint get_shader_program(const std::string& name);

	// Uniform setters
	void set_uniform_mat4(const std::string& shader_name, const std::string& uniform_name, const glm::mat4& matrix);
	void set_uniform_vec3(const std::string& shader_name, const std::string& uniform_name, const glm::vec3& vector);
	void set_uniform_vec4(const std::string& shader_name, const std::string& uniform_name, const glm::vec4& vector);
	void set_uniform_float(const std::string& shader_name, const std::string& uniform_name, float value);
	void set_uniform_int(const std::string& shader_name, const std::string& uniform_name, int value);

private:
	std::unordered_map<std::string, GLuint> shader_programs_;

	GLuint compile_shader(const std::string& source, GLenum shader_type);
	GLuint create_program(const std::string& vertex_source, const std::string& fragment_source);
	std::string load_shader_source(const std::string& file_path);
	GLint get_uniform_location(const std::string& shader_name, const std::string& uniform_name);
};
