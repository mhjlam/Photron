#pragma once

#include <glm/glm.hpp>

struct RenderVertex
{
	glm::vec3 position{0.0f};
	glm::vec3 normal{0.0f, 1.0f, 0.0f};
	glm::vec2 tex_coords{0.0f};
	glm::vec4 color{1.0f};

	// Default constructor uses default member initialization
	RenderVertex() = default;

	explicit RenderVertex(const glm::vec3& pos) noexcept
		: position(pos) {}

	RenderVertex(const glm::vec3& pos, const glm::vec4& col) noexcept
		: position(pos), color(col) {}

	RenderVertex(const glm::vec3& pos, const glm::vec3& norm, const glm::vec2& tex, const glm::vec4& col) noexcept
		: position(pos), normal(norm), tex_coords(tex), color(col) {}
};
