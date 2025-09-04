#pragma once

#include <cstdint>
#include <vector>

#include "math/mesh_geometry.hpp"
#include "math/triangle.hpp"

struct Layer
{
	uint8_t id;
	uint8_t tissue_id;
	std::vector<Triangle> mesh;
	MeshGeometry geometry;

	Layer() {
		id = 0;
		tissue_id = 0;
		mesh = std::vector<Triangle>();
	}

	// Update the mesh geometry based on the mesh triangles
	void update_geometry() {
		geometry.set_triangles(mesh);
	}

	// Check if a point is inside this layer's mesh
	bool contains_point(const glm::dvec3& point) const {
		return geometry.contains(point);
	}

	bool operator==(const Layer& other) const { return other.id == id; }
};
