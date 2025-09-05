#pragma once

#include <cstdint>
#include <vector>

#include "math/mesh_geometry.hpp"
#include "math/triangle.hpp"

struct Layer
{
	uint8_t id = 0;
	uint8_t tissue_id = 0;
	std::vector<Triangle> mesh;
	MeshGeometry geometry;

	// Default constructor uses default member initialization
	Layer() = default;

	// Update the mesh geometry based on the mesh triangles
	void update_geometry() noexcept {
		geometry.set_triangles(mesh);
	}

	// Check if a point is inside this layer's mesh
	bool contains_point(const glm::dvec3& point) const noexcept {
		return geometry.contains(point);
	}

	bool operator==(const Layer& other) const noexcept { return other.id == id; }
};
