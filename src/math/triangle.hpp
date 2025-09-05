#pragma once

#include "glm_types.hpp"

class Triangle
{
private:
	glm::dvec3 v0_{0.0};
	glm::dvec3 v1_{0.0};
	glm::dvec3 v2_{0.0};
	glm::dvec3 normal_{0.0, 1.0, 0.0};

public:
	Triangle() = default;
	Triangle(const glm::dvec3& v0, const glm::dvec3& v1, const glm::dvec3& v2) noexcept;

	// Getters
	const glm::dvec3& v0() const noexcept { return v0_; }
	const glm::dvec3& v1() const noexcept { return v1_; }
	const glm::dvec3& v2() const noexcept { return v2_; }
	const glm::dvec3& normal() const noexcept { return normal_; }

	// Setters
	void set_vertices(const glm::dvec3& v0, const glm::dvec3& v1, const glm::dvec3& v2) noexcept;
	void set_normal(const glm::dvec3& normal) noexcept { normal_ = normal; }

	// Utility methods
	bool is_invalid() const noexcept;
	glm::dvec3 center() const noexcept;
	double area() const noexcept;

private:
	void compute_normal() noexcept;
};
