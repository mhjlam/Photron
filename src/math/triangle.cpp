/**
 * @file triangle.cpp
 * @brief Triangle geometry operations and mesh processing utilities
 *
 * Implements comprehensive triangle geometry operations for 3D mesh processing:
 * - Triangle construction with automatic normal computation
 * - Geometric property calculations (area, normal, centroid)
 * - Validation and degenerate case detection
 * - Efficient vertex management with cached normals
 *
 * All algorithms use GLM vectorized operations for optimal performance
 * in Monte Carlo simulations and 3D rendering applications.
 */

#include "triangle.hpp"

/**
 * @brief Construct triangle from three vertices with automatic normal computation
 * @param v0,v1,v2 Triangle vertices in counter-clockwise order
 *
 * Creates triangle with explicit vertex specification and computes surface normal.
 * Vertex ordering determines normal direction via right-hand rule.
 */
Triangle::Triangle(const glm::dvec3& v0, const glm::dvec3& v1, const glm::dvec3& v2) noexcept :
	v0_(v0), v1_(v1), v2_(v2) {
	compute_normal();
}

/**
 * @brief Update all triangle vertices and recompute normal
 * @param v0,v1,v2 New vertex positions
 *
 * Efficiently updates triangle geometry and maintains cached normal vector.
 */
void Triangle::set_vertices(const glm::dvec3& v0, const glm::dvec3& v1, const glm::dvec3& v2) noexcept {
	v0_ = v0;
	v1_ = v1;
	v2_ = v2;
	compute_normal();
}

/**
 * @brief Compute and cache triangle surface normal
 *
 * Internal method that calculates normal via cross product of edge vectors.
 * Uses right-hand rule based on vertex ordering (v0->v1->v2).
 */
void Triangle::compute_normal() noexcept {
	// Calculate edge vectors for cross product computation
	glm::dvec3 edge1 = v1_ - v0_;
	glm::dvec3 edge2 = v2_ - v0_;

	// Compute and store normalized surface normal
	normal_ = glm::normalize(glm::cross(edge1, edge2));
}

/**
 * @brief Test if triangle is degenerate (has duplicate vertices)
 * @return True if any two vertices are identical
 *
 * Detects degenerate triangles that have zero area and can cause
 * numerical issues in geometric calculations.
 */
bool Triangle::is_invalid() const noexcept {
	return (v0_ == v1_ || v0_ == v2_ || v1_ == v2_);
}

/**
 * @brief Calculate triangle centroid
 * @return Geometric center point
 *
 * Computes centroid as arithmetic mean of three vertices.
 */
glm::dvec3 Triangle::center() const noexcept {
	return (v0_ + v1_ + v2_) / 3.0;
}

/**
 * @brief Calculate triangle area
 * @return Triangle area as positive scalar
 *
 * Uses cross product magnitude for numerically stable area calculation.
 */
double Triangle::area() const noexcept {
	glm::dvec3 edge1 = v1_ - v0_;
	glm::dvec3 edge2 = v2_ - v0_;
	return 0.5 * glm::length(glm::cross(edge1, edge2));
}
