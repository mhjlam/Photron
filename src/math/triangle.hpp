#pragma once

#include "glm_types.hpp"

class Triangle
{
private:
	glm::dvec3 v0_;
	glm::dvec3 v1_;
	glm::dvec3 v2_;
	glm::dvec3 normal_;

public:
	Triangle();
	Triangle(const glm::dvec3& v0, const glm::dvec3& v1, const glm::dvec3& v2);

	// Getters
	const glm::dvec3& v0() const { return v0_; }
	const glm::dvec3& v1() const { return v1_; }
	const glm::dvec3& v2() const { return v2_; }
	const glm::dvec3& normal() const { return normal_; }

	// Setters
	void set_vertices(const glm::dvec3& v0, const glm::dvec3& v1, const glm::dvec3& v2);
	void set_normal(const glm::dvec3& normal) { normal_ = normal; }

	// Utility methods
	bool is_invalid() const;
	glm::dvec3 center() const;
	double area() const;

private:
	void compute_normal();
};
