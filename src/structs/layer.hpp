#pragma once

#include "triangle.hpp"

#include <cstdint>
#include <vector>

struct Layer {
	uint8_t id;
	uint8_t tissue_id;
	std::vector<Triangle> mesh;

	Layer() {
		id = 0;
		tissue_id = 0;

		mesh = std::vector<Triangle>();
	}

	bool operator==(const Layer& other) const { return other.id == id; }
};
