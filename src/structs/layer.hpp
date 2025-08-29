#pragma once

#include "triangle.hpp"

#include <vector>

struct Layer {
	unsigned short id;
	unsigned short tissueid;
	std::vector<Triangle> mesh;

	Layer() {
		id = 0;
		tissueid = 0;

		mesh = std::vector<Triangle>();
	}

	bool operator==(const Layer& other) const { return other.id == id; }
};
