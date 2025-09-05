#include "triangle.hpp"

Triangle::Triangle(const glm::dvec3& v0, const glm::dvec3& v1, const glm::dvec3& v2) noexcept 
	: v0_(v0), v1_(v1), v2_(v2) {
	compute_normal();
}

void Triangle::set_vertices(const glm::dvec3& v0, const glm::dvec3& v1, const glm::dvec3& v2) noexcept {
	v0_ = v0;
	v1_ = v1;
	v2_ = v2;
	compute_normal();
}

void Triangle::compute_normal() noexcept {
	// Compute normal using cross product
	glm::dvec3 edge1 = v1_ - v0_;
	glm::dvec3 edge2 = v2_ - v0_;
	normal_ = glm::normalize(glm::cross(edge1, edge2));
}

bool Triangle::is_invalid() const noexcept {
	return (v0_ == v1_ || v0_ == v2_ || v1_ == v2_);
}

glm::dvec3 Triangle::center() const noexcept {
	return (v0_ + v1_ + v2_) / 3.0;
}

double Triangle::area() const noexcept {
	glm::dvec3 edge1 = v1_ - v0_;
	glm::dvec3 edge2 = v2_ - v0_;
	return 0.5 * glm::length(glm::cross(edge1, edge2));
}
