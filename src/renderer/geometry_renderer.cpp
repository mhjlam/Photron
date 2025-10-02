/**
 * @file geometry_renderer.cpp
 * @brief Implementation of specialized medium geometry visualization system
 */

#include "geometry_renderer.hpp"

#include <algorithm>
#include <format>
#include <fstream>
#include <iostream>
#include <map>

#include <glm/gtc/matrix_transform.hpp>
#include <glm/gtc/type_ptr.hpp>

#include "math/math.hpp"
#include "renderer/camera.hpp"
#include "simulator/layer.hpp"
#include "simulator/simulator.hpp"

GeometryRenderer::GeometryRenderer() = default;

GeometryRenderer::~GeometryRenderer() {
	// Clean up OpenGL resources
	if (triangles_vao_)
		glDeleteVertexArrays(1, &triangles_vao_);
	if (triangles_vbo_)
		glDeleteBuffers(1, &triangles_vbo_);
	if (triangles_instance_vbo_)
		glDeleteBuffers(1, &triangles_instance_vbo_);
	if (triangles_shader_)
		glDeleteProgram(triangles_shader_);

	if (lines_vao_)
		glDeleteVertexArrays(1, &lines_vao_);
	if (lines_vbo_)
		glDeleteBuffers(1, &lines_vbo_);
	if (lines_instance_vbo_)
		glDeleteBuffers(1, &lines_instance_vbo_);
	if (lines_shader_)
		glDeleteProgram(lines_shader_);
}

bool GeometryRenderer::initialize() {
	return setup_triangle_rendering() && setup_line_rendering();
}

void GeometryRenderer::render(const Simulator& simulator, const Settings& settings, const Camera& camera) {
	if (!settings.draw_volume) {
		return;
	}

	// Check if geometry cache needs rebuilding
	if (!geometry_cached_) {
		collect_geometry_instances(simulator);
		geometry_cached_ = true;
		triangle_buffer_uploaded_ = false;
		line_buffer_uploaded_ = false;
	}

	// Enable blending for transparent faces
	glEnable(GL_BLEND);
	glBlendFunc(GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);

	// Render triangles first, then wireframes
	draw_triangle_instances(camera);
	draw_line_instances(camera);
}

void GeometryRenderer::invalidate_cache() {
	geometry_cached_ = false;
	triangle_instances_.clear();
	line_instances_.clear();
	triangle_buffer_uploaded_ = false;
	line_buffer_uploaded_ = false;
}

bool GeometryRenderer::setup_triangle_rendering() {
	// Load and compile triangle instance shaders
	std::string vertex_source = load_shader_source("shaders/triangles.vert");
	std::string fragment_source = load_shader_source("shaders/triangles.frag");

	if (vertex_source.empty() || fragment_source.empty()) {
		std::cerr << "GeometryRenderer: Failed to load triangle shader sources" << std::endl;
		return false;
	}

	triangles_shader_ = create_shader_program(vertex_source, fragment_source);
	if (!triangles_shader_) {
		std::cerr << "GeometryRenderer: Failed to create triangle shader program" << std::endl;
		return false;
	}

	// Create VAO for triangle instances
	glGenVertexArrays(1, &triangles_vao_);
	glBindVertexArray(triangles_vao_);

	// Create base triangle geometry (single triangle)
	glGenBuffers(1, &triangles_vbo_);
	glBindBuffer(GL_ARRAY_BUFFER, triangles_vbo_);

	// Simple triangle vertices (will be transformed per instance)
	float triangle_vertices[] = {0.0f, 0.0f, 0.0f, 1.0f, 0.0f, 0.0f, 0.0f, 1.0f, 0.0f};
	glBufferData(GL_ARRAY_BUFFER, sizeof(triangle_vertices), triangle_vertices, GL_STATIC_DRAW);

	// Setup vertex positions (layout location 0)
	glVertexAttribPointer(0, 3, GL_FLOAT, GL_FALSE, 3 * sizeof(float), (void*)0);
	glEnableVertexAttribArray(0);

	// Create instance data buffer
	glGenBuffers(1, &triangles_instance_vbo_);
	glBindBuffer(GL_ARRAY_BUFFER, triangles_instance_vbo_);

	// Setup instance attributes (per-triangle data)
	// v0 (layout location 1)
	glVertexAttribPointer(1, 3, GL_FLOAT, GL_FALSE, sizeof(TriangleInstance), (void*)offsetof(TriangleInstance, v0));
	glEnableVertexAttribArray(1);
	glVertexAttribDivisor(1, 1);

	// v1 (layout location 2)
	glVertexAttribPointer(2, 3, GL_FLOAT, GL_FALSE, sizeof(TriangleInstance), (void*)offsetof(TriangleInstance, v1));
	glEnableVertexAttribArray(2);
	glVertexAttribDivisor(2, 1);

	// v2 (layout location 3)
	glVertexAttribPointer(3, 3, GL_FLOAT, GL_FALSE, sizeof(TriangleInstance), (void*)offsetof(TriangleInstance, v2));
	glEnableVertexAttribArray(3);
	glVertexAttribDivisor(3, 1);

	// Color (layout location 4)
	glVertexAttribPointer(4, 4, GL_FLOAT, GL_FALSE, sizeof(TriangleInstance), (void*)offsetof(TriangleInstance, color));
	glEnableVertexAttribArray(4);
	glVertexAttribDivisor(4, 1);

	// Normal (layout location 5)
	glVertexAttribPointer(
		5, 3, GL_FLOAT, GL_FALSE, sizeof(TriangleInstance), (void*)offsetof(TriangleInstance, normal));
	glEnableVertexAttribArray(5);
	glVertexAttribDivisor(5, 1);

	glBindVertexArray(0);
	return true;
}

bool GeometryRenderer::setup_line_rendering() {
	// Load and compile line shaders for wireframe rendering
	std::string vertex_source = load_shader_source("shaders/lines.vert");
	std::string fragment_source = load_shader_source("shaders/lines.frag");

	if (vertex_source.empty() || fragment_source.empty()) {
		std::cerr << "GeometryRenderer: Failed to load line shader sources" << std::endl;
		return false;
	}

	lines_shader_ = create_shader_program(vertex_source, fragment_source);
	if (!lines_shader_) {
		std::cerr << "GeometryRenderer: Failed to create line shader program" << std::endl;
		return false;
	}

	// Get uniform location for MVP matrix
	line_instanced_mvp_uniform_location_ = glGetUniformLocation(lines_shader_, "uMVP");
	if (line_instanced_mvp_uniform_location_ == -1) {
		std::cerr << "GeometryRenderer: Warning - Could not find uMVP uniform in line shader" << std::endl;
	}

	// Create VAO for line instances
	glGenVertexArrays(1, &lines_vao_);
	glBindVertexArray(lines_vao_);

	// Create base line geometry VBO
	glGenBuffers(1, &lines_vbo_);
	glBindBuffer(GL_ARRAY_BUFFER, lines_vbo_);

	// Setup line geometry inline (unit line from origin to (1,0,0))
	float line_vertices[] = {
		0.0f,
		0.0f,
		0.0f, // Start point
		1.0f,
		0.0f,
		0.0f  // End point
	};
	glBufferData(GL_ARRAY_BUFFER, sizeof(line_vertices), line_vertices, GL_STATIC_DRAW);

	// Setup vertex positions (layout location 0)
	glVertexAttribPointer(0, 3, GL_FLOAT, GL_FALSE, 3 * sizeof(float), (void*)0);
	glEnableVertexAttribArray(0);

	// Create separate instance data buffer for lines
	glGenBuffers(1, &lines_instance_vbo_);
	glBindBuffer(GL_ARRAY_BUFFER, lines_instance_vbo_);

	// Setup instance attributes (per-line data)
	// Start position (layout location 1)
	glVertexAttribPointer(
		1, 3, GL_FLOAT, GL_FALSE, sizeof(MediumLineInstance), (void*)offsetof(MediumLineInstance, start));
	glEnableVertexAttribArray(1);
	glVertexAttribDivisor(1, 1);

	// End position (layout location 2)
	glVertexAttribPointer(
		2, 3, GL_FLOAT, GL_FALSE, sizeof(MediumLineInstance), (void*)offsetof(MediumLineInstance, end));
	glEnableVertexAttribArray(2);
	glVertexAttribDivisor(2, 1);

	// Start color (layout location 3)
	glVertexAttribPointer(
		3, 4, GL_FLOAT, GL_FALSE, sizeof(MediumLineInstance), (void*)offsetof(MediumLineInstance, start_color));
	glEnableVertexAttribArray(3);
	glVertexAttribDivisor(3, 1);

	// End color (layout location 4)
	glVertexAttribPointer(
		4, 4, GL_FLOAT, GL_FALSE, sizeof(MediumLineInstance), (void*)offsetof(MediumLineInstance, end_color));
	glEnableVertexAttribArray(4);
	glVertexAttribDivisor(4, 1);

	glBindVertexArray(0);
	return true;
}

void GeometryRenderer::collect_geometry_instances(const Simulator& simulator) {
	triangle_instances_.clear();
	line_instances_.clear();

	glm::vec4 wireframe_color(0.7f, 0.7f, 0.7f, 0.8f); // Bright wireframe
	glm::vec4 face_color(0.2f, 0.2f, 0.2f, 0.1f);      // Subtle transparent faces

	// Process each layer's geometry
	const auto& layers = simulator.get_all_layers();
	for (const auto& layer : layers) {
		// For each triangle in the layer's mesh, we'll determine if it's part of a planar face
		// and group triangles that share the same plane into faces

		std::map<std::string, std::vector<glm::vec3>> planar_faces;

		for (const auto& triangle : layer.mesh) {
			glm::vec3 v0 = to_float(triangle.v0());
			glm::vec3 v1 = to_float(triangle.v1());
			glm::vec3 v2 = to_float(triangle.v2());

			// Calculate the triangle normal
			glm::vec3 normal = normalize(cross(v1 - v0, v2 - v0));

			// Calculate plane equation: ax + by + cz = d
			float d = dot(normal, v0);

			// Create a key for this plane (quantized to handle floating point precision)
			std::string plane_key = std::format(
				"{},{},{},{}", int(normal.x * 1000), int(normal.y * 1000), int(normal.z * 1000), int(d * 1000));

			// Add vertices to this plane's face
			auto& face_vertices = planar_faces[plane_key];
			face_vertices.push_back(v0);
			face_vertices.push_back(v1);
			face_vertices.push_back(v2);
		}

		// Now process each planar face
		for (const auto& face_pair : planar_faces) {
			const std::vector<glm::vec3>& vertices = face_pair.second;

			if (vertices.size() == 3) {
				// Single triangle - render as triangle
				glm::vec3 normal = normalize(cross(vertices[1] - vertices[0], vertices[2] - vertices[0]));
				add_triangle_instance(vertices[0], vertices[1], vertices[2], face_color, normal);

				// Add wireframe edges
				add_line_instance(vertices[0], vertices[1], wireframe_color);
				add_line_instance(vertices[1], vertices[2], wireframe_color);
				add_line_instance(vertices[2], vertices[0], wireframe_color);
			}
			else if (vertices.size() >= 6 && vertices.size() % 3 == 0) {
				// Multiple triangles forming a face - check if they form a quad

				// For now, let's find unique vertices and try to form a quadrilateral
				std::vector<glm::vec3> unique_vertices;
				const float epsilon = MathConstants::UV_MATCH_EPSILON;

				// Collect unique vertices
				for (const auto& v : vertices) {
					bool is_unique = true;
					for (const auto& uv : unique_vertices) {
						if (length(v - uv) < epsilon) {
							is_unique = false;
							break;
						}
					}
					if (is_unique) {
						unique_vertices.push_back(v);
					}
				}

				if (unique_vertices.size() == 4) {
					// We have a quadrilateral! Order the vertices properly
					// Find the centroid
					glm::vec3 center(0.0f);
					for (const auto& v : unique_vertices) {
						center += v;
					}
					center /= static_cast<float>(unique_vertices.size());

					// Calculate the face normal from first triangle
					glm::vec3 normal = normalize(cross(vertices[1] - vertices[0], vertices[2] - vertices[0]));

					// Create a reference vector perpendicular to normal
					glm::vec3 ref_vec =
						abs(normal.x) < 0.9f ? glm::vec3(1.0f, 0.0f, 0.0f) : glm::vec3(0.0f, 1.0f, 0.0f);
					glm::vec3 tangent = normalize(cross(normal, ref_vec));
					glm::vec3 bitangent = normalize(cross(normal, tangent));

					// Sort vertices by angle around the center
					std::ranges::sort(unique_vertices,
									  [&center, &tangent, &bitangent](const glm::vec3& a, const glm::vec3& b) noexcept {
										  const glm::vec3 dir_a = a - center;
										  const glm::vec3 dir_b = b - center;

										  const float angle_a = atan2f(dot(dir_a, bitangent), dot(dir_a, tangent));
										  const float angle_b = atan2f(dot(dir_b, bitangent), dot(dir_b, tangent));

										  return angle_a < angle_b;
									  });

					// Render as two triangles forming a quad
					add_triangle_instance(
						unique_vertices[0], unique_vertices[1], unique_vertices[2], face_color, normal);
					add_triangle_instance(
						unique_vertices[0], unique_vertices[2], unique_vertices[3], face_color, normal);

					// Add wireframe edges for the quad
					add_line_instance(unique_vertices[0], unique_vertices[1], wireframe_color);
					add_line_instance(unique_vertices[1], unique_vertices[2], wireframe_color);
					add_line_instance(unique_vertices[2], unique_vertices[3], wireframe_color);
					add_line_instance(unique_vertices[3], unique_vertices[0], wireframe_color);
				}
				else {
					// Fallback: render all triangles individually
					for (size_t i = 0; i < vertices.size(); i += 3) {
						glm::vec3 normal =
							normalize(cross(vertices[i + 1] - vertices[i], vertices[i + 2] - vertices[i]));
						add_triangle_instance(vertices[i], vertices[i + 1], vertices[i + 2], face_color, normal);

						// Add wireframe edges
						add_line_instance(vertices[i], vertices[i + 1], wireframe_color);
						add_line_instance(vertices[i + 1], vertices[i + 2], wireframe_color);
						add_line_instance(vertices[i + 2], vertices[i], wireframe_color);
					}
				}
			}
		}
	}
}

void GeometryRenderer::draw_triangle_instances(const Camera& camera) {
	if (triangle_instances_.empty() || !triangles_shader_) {
		return;
	}

	if (!triangle_buffer_uploaded_) {
		glBindBuffer(GL_ARRAY_BUFFER, triangles_instance_vbo_);
		glBufferData(GL_ARRAY_BUFFER,
					 triangle_instances_.size() * sizeof(TriangleInstance),
					 triangle_instances_.data(),
					 GL_DYNAMIC_DRAW);
		triangle_buffer_uploaded_ = true;
	}

	glUseProgram(triangles_shader_);

	glm::mat4 mvp = camera.get_projection_matrix() * camera.get_view_matrix();
	GLint mvp_location = glGetUniformLocation(triangles_shader_, "mvp");
	if (mvp_location >= 0) {
		glUniformMatrix4fv(mvp_location, 1, GL_FALSE, &mvp[0][0]);
	}

	glBindVertexArray(triangles_vao_);
	glDrawArraysInstanced(GL_TRIANGLES, 0, 3, static_cast<GLsizei>(triangle_instances_.size()));
	glBindVertexArray(0);
}

void GeometryRenderer::draw_line_instances(const Camera& camera) {
	if (line_instances_.empty() || !lines_shader_) {
		return;
	}

	if (!line_buffer_uploaded_) {
		glBindBuffer(GL_ARRAY_BUFFER, lines_instance_vbo_);
		glBufferData(GL_ARRAY_BUFFER,
					 line_instances_.size() * sizeof(MediumLineInstance),
					 line_instances_.data(),
					 GL_DYNAMIC_DRAW);
		line_buffer_uploaded_ = true;
	}

	glUseProgram(lines_shader_);

	glm::mat4 mvp = camera.get_projection_matrix() * camera.get_view_matrix();
	glUniformMatrix4fv(line_instanced_mvp_uniform_location_, 1, GL_FALSE, &mvp[0][0]);

	glBindVertexArray(lines_vao_);
	glDrawArraysInstanced(GL_LINES, 0, 2, static_cast<GLsizei>(line_instances_.size()));
	glBindVertexArray(0);
}

void GeometryRenderer::add_triangle_instance(
	const glm::vec3& v0, const glm::vec3& v1, const glm::vec3& v2, const glm::vec4& color, const glm::vec3& normal) {
	TriangleInstance instance;
	instance.v0 = v0;
	instance.v1 = v1;
	instance.v2 = v2;
	instance.color = color;
	instance.normal = normal;
	triangle_instances_.push_back(instance);
}

void GeometryRenderer::add_line_instance(const glm::vec3& start, const glm::vec3& end, const glm::vec4& color) {
	MediumLineInstance instance;
	instance.start = start;
	instance.end = end;
	instance.start_color = color;
	instance.end_color = color;
	line_instances_.push_back(instance);
}

bool GeometryRenderer::is_point_inside_mesh(const glm::vec3& point, const Simulator& simulator) const {
	// Use simulator's built-in geometry testing
	glm::dvec3 dpoint(point.x, point.y, point.z);
	return simulator.is_point_inside_geometry(dpoint);
}

GLuint GeometryRenderer::create_shader_program(const std::string& vertex_source, const std::string& fragment_source) {
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
		std::cerr << "GeometryRenderer: Program linking failed: " << info_log << std::endl;
		glDeleteProgram(program);
		program = 0;
	}

	glDeleteShader(vertex_shader);
	glDeleteShader(fragment_shader);
	return program;
}

GLuint GeometryRenderer::compile_shader(const std::string& source, GLenum type) {
	GLuint shader = glCreateShader(type);
	const char* source_cstr = source.c_str();
	glShaderSource(shader, 1, &source_cstr, nullptr);
	glCompileShader(shader);

	GLint success;
	glGetShaderiv(shader, GL_COMPILE_STATUS, &success);

	if (!success) {
		char info_log[512];
		glGetShaderInfoLog(shader, 512, nullptr, info_log);
		std::cerr << "GeometryRenderer: Shader compilation failed: " << info_log << std::endl;
		glDeleteShader(shader);
		return 0;
	}

	return shader;
}

std::string GeometryRenderer::load_shader_source(const std::string& file_path) {
	std::ifstream file(file_path);
	if (!file.is_open()) {
		std::cerr << "GeometryRenderer: Failed to open shader file: " << file_path << std::endl;
		return "";
	}

	std::string content;
	std::string line;
	while (std::getline(file, line)) {
		content += line + "\n";
	}

	return content;
}