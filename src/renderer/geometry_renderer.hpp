/**
 * @file geometry_renderer.hpp
 * @brief Specialized medium geometry visualization system
 *
 * Handles rendering of medium boundaries, wireframes, and 3D geometry faces.
 * Provides both solid face rendering and wireframe visualization for
 * understanding simulation geometry.
 */

#pragma once

#include <map>
#include <string>
#include <vector>

#include <GL/glew.h>
#include <glm/glm.hpp>

#include "renderer/settings.hpp"

// Forward declarations
class Simulator;
class Camera;

/**
 * @class GeometryRenderer
 * @brief Specialized renderer for medium geometry visualization
 *
 * Handles all aspects of geometry rendering including:
 * - Medium boundary triangles with transparency
 * - Wireframe edges for structure clarity
 * - Planar face detection and optimization
 * - Instanced rendering for performance
 */
class GeometryRenderer
{
public:
	/**
	 * @brief Construct a new GeometryRenderer
	 */
	GeometryRenderer();

	/**
	 * @brief Destroy the GeometryRenderer and clean up resources
	 */
	~GeometryRenderer();

	/**
	 * @brief Initialize OpenGL resources for geometry rendering
	 * @return true if initialization succeeded, false otherwise
	 */
	bool initialize();

	/**
	 * @brief Render medium geometry with current settings
	 * @param simulator Simulator containing geometry data
	 * @param settings Current rendering settings
	 * @param camera Camera for view/projection matrices
	 */
	void render(const Simulator& simulator, const Settings& settings, const Camera& camera);

	/**
	 * @brief Invalidate cached geometry data
	 */
	void invalidate_cache();

	/**
	 * @brief Check if point is inside rendered mesh geometry
	 * @param point 3D point to test
	 * @param simulator Simulator containing geometry data
	 * @return true if point is inside mesh, false otherwise
	 */
	bool is_point_inside_mesh(const glm::vec3& point, const Simulator& simulator) const;

private:
	/**
	 * @brief Instance data structure for triangle rendering
	 */
	struct TriangleInstance
	{
		glm::vec3 v0 {};
		glm::vec3 v1 {};
		glm::vec3 v2 {};
		glm::vec4 color {1.0f};
		glm::vec3 normal {}; // For lighting calculations
	};

	/**
	 * @brief Instance data structure for medium line rendering
	 */
	struct MediumLineInstance
	{
		glm::vec3 start {};
		glm::vec3 end {};
		glm::vec4 start_color {1.0f};
		glm::vec4 end_color {1.0f};
	};

	// Setup methods
	bool setup_triangle_rendering();
	bool setup_line_rendering();

	// Geometry processing
	void collect_geometry_instances(const Simulator& simulator);
	void process_layer_geometry(const auto& layer);
	void detect_planar_faces(const std::vector<glm::vec3>& vertices,
							 const glm::vec4& face_color,
							 const glm::vec4& wireframe_color);

	// Rendering
	void draw_triangle_instances(const Camera& camera);
	void draw_line_instances(const Camera& camera);

	// Utility
	void add_triangle_instance(
		const glm::vec3& v0, const glm::vec3& v1, const glm::vec3& v2, const glm::vec4& color, const glm::vec3& normal);
	void add_line_instance(const glm::vec3& start, const glm::vec3& end, const glm::vec4& color);

	// Shader utilities
	GLuint create_shader_program(const std::string& vertex_source, const std::string& fragment_source);
	GLuint compile_shader(const std::string& source, GLenum type);
	std::string load_shader_source(const std::string& file_path);

	// OpenGL resources for triangles
	GLuint triangles_vao_ {0};
	GLuint triangles_vbo_ {0};
	GLuint triangles_instance_vbo_ {0};
	GLuint triangles_shader_ {0};

	// OpenGL resources for lines
	GLuint lines_vao_ {0};
	GLuint lines_vbo_ {0};
	GLuint lines_instance_vbo_ {0};
	GLuint lines_shader_ {0};
	GLint line_instanced_mvp_uniform_location_ {-1};

	// Instance data
	std::vector<TriangleInstance> triangle_instances_;
	std::vector<MediumLineInstance> line_instances_;

	// Performance caching
	mutable bool geometry_cached_ {false};
	mutable size_t last_geometry_version_ {0};
	mutable bool triangle_buffer_uploaded_ {false};
	mutable bool line_buffer_uploaded_ {false};
};