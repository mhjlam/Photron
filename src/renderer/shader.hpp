#pragma once

#include <string>

#include <GL/glew.h>
#include <glm/glm.hpp>

class Shader
{
public:
	Shader() = default;
	~Shader();

	// Load shaders from files
	bool load_from_file(const std::string& vertex_path, const std::string& fragment_path);

	// Load shaders from source strings
	bool load_from_source(const std::string& vertex_source, const std::string& fragment_source);

	// Use/unuse shader
	void use() const;
	void unuse() const;

	// Uniform setters
	void set_uniform(const std::string& name, bool value);
	void set_uniform(const std::string& name, int value);
	void set_uniform(const std::string& name, float value);
	void set_uniform(const std::string& name, const glm::vec2& value);
	void set_uniform(const std::string& name, const glm::vec3& value);
	void set_uniform(const std::string& name, const glm::vec4& value);
	void set_uniform(const std::string& name, const glm::mat3& value);
	void set_uniform(const std::string& name, const glm::mat4& value);

	// Get program ID
	GLuint id() const { return program_id_; }
	bool is_valid() const { return program_id_ != 0; }

private:
	GLuint program_id_ = 0;

	// Helper functions
	std::string read_file(const std::string& path);
	GLuint compile_shader(const std::string& source, GLenum type);
	bool link_program(GLuint vertex_shader, GLuint fragment_shader);
	GLint get_uniform_location(const std::string& name);
};
