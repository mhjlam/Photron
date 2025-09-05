#pragma once

#include <cstdint>

#include <glm/glm.hpp>

#include "math/triangle.hpp"

struct Source
{
	uint64_t id {0};                     // identifier
	glm::dvec3 origin {0.0};             // origin
	glm::dvec3 direction {0.0};          // direction of incidence
	glm::dvec3 specular_direction {0.0}; // direction of specular reflectance

	glm::dvec3 intersect {0.0};          // intersection point
	Triangle triangle;                   // triangle at intersection point

	// Default constructor uses default member initialization
	Source() = default;

	explicit Source(uint64_t i, const glm::dvec3& p, const glm::dvec3& v) noexcept : id(i), origin(p), direction(v) {}
};
