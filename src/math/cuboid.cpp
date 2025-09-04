#include "cuboid.hpp"

#include <algorithm>
#include <iostream>
#include <limits>

#include "ray.hpp"

Cuboid::Cuboid(double x1, double y1, double z1, double x2, double y2, double z2) :
	min_point_(x1, y1, z1), max_point_(x2, y2, z2) {
}

Cuboid::Cuboid(const glm::dvec3& min_pt, const glm::dvec3& max_pt) : min_point_(min_pt), max_point_(max_pt) {
}

void Cuboid::set_bounds(const glm::dvec3& min_pt, const glm::dvec3& max_pt) {
	min_point_ = min_pt;
	max_point_ = max_pt;
}

glm::dvec3 Cuboid::center() const {
	return (min_point_ + max_point_) * 0.5;
}

glm::dvec3 Cuboid::size() const {
	return max_point_ - min_point_;
}

bool Cuboid::contains(const glm::dvec3& point) const {
	return point.x >= min_point_.x && point.x <= max_point_.x && point.y >= min_point_.y && point.y <= max_point_.y
		   && point.z >= min_point_.z && point.z <= max_point_.z;
}

double Cuboid::volume() const {
	glm::dvec3 size_vec = size();
	return size_vec.x * size_vec.y * size_vec.z;
}

double Cuboid::intersect_ray_internal(const Ray& ray, glm::dvec3& intersection, glm::dvec3& normal) const {
	// Standard ray-AABB intersection algorithm
	// Based on "Real-Time Rendering" and other graphics literature
	
	glm::dvec3 ray_origin = ray.origin();
	glm::dvec3 ray_dir = glm::normalize(ray.direction());
	
	// Calculate intersection distances with each slab
	glm::dvec3 t_min = (min_point_ - ray_origin) / ray_dir;
	glm::dvec3 t_max = (max_point_ - ray_origin) / ray_dir;
	
	// Handle potential division by zero (ray parallel to slab)
	// When ray_dir component is very small, use a large value
	const double epsilon = 1e-10;
	if (std::abs(ray_dir.x) < epsilon) {
		if (ray_origin.x < min_point_.x || ray_origin.x > max_point_.x) {
			return std::numeric_limits<double>::max(); // No intersection
		}
		t_min.x = -std::numeric_limits<double>::max();
		t_max.x = std::numeric_limits<double>::max();
	}
	if (std::abs(ray_dir.y) < epsilon) {
		if (ray_origin.y < min_point_.y || ray_origin.y > max_point_.y) {
			return std::numeric_limits<double>::max(); // No intersection
		}
		t_min.y = -std::numeric_limits<double>::max();
		t_max.y = std::numeric_limits<double>::max();
	}
	if (std::abs(ray_dir.z) < epsilon) {
		if (ray_origin.z < min_point_.z || ray_origin.z > max_point_.z) {
			return std::numeric_limits<double>::max(); // No intersection
		}
		t_min.z = -std::numeric_limits<double>::max();
		t_max.z = std::numeric_limits<double>::max();
	}
	
	// Ensure t_min <= t_max for each axis
	if (t_min.x > t_max.x) std::swap(t_min.x, t_max.x);
	if (t_min.y > t_max.y) std::swap(t_min.y, t_max.y);
	if (t_min.z > t_max.z) std::swap(t_min.z, t_max.z);
	
	// Find the largest t_min and smallest t_max
	double t_near = std::max({t_min.x, t_min.y, t_min.z});
	double t_far = std::min({t_max.x, t_max.y, t_max.z});
	
	// Use a small tolerance for floating point precision issues
	const double tolerance = 1e-6;
	
	// Check if there's an intersection
	if (t_near > t_far || t_far < -tolerance) {
		return std::numeric_limits<double>::max(); // No intersection
	}
	
	// For internal intersection, we want the exit point (t_far)
	// if we're inside the box (t_near < 0), otherwise t_near
	double t;
	if (t_near < 0) {
		// Ray starts inside the box - use exit point
		t = t_far;
	} else {
		// Ray starts outside the box - use entry point
		t = t_near;
	}
	
	// Accept small positive t values and handle floating-point precision issues
	// For rays exactly on boundaries
	if (t < -tolerance) {
		return std::numeric_limits<double>::max(); // No valid forward intersection
	}
	
	// Clamp very small negative t to 0 (floating-point precision issue)
	if (t < 0) {
		t = 0;
	}
	
	// Calculate intersection point
	intersection = ray_origin + t * ray_dir;
	
	// Determine which face was hit by checking which axis achieved t
	if (std::abs(t - t_min.x) < epsilon) {
		normal = glm::dvec3(-1, 0, 0); // Left face (inward normal)
	} else if (std::abs(t - t_max.x) < epsilon) {
		normal = glm::dvec3(1, 0, 0);  // Right face (inward normal)
	} else if (std::abs(t - t_min.y) < epsilon) {
		normal = glm::dvec3(0, -1, 0); // Bottom face (inward normal)
	} else if (std::abs(t - t_max.y) < epsilon) {
		normal = glm::dvec3(0, 1, 0);  // Top face (inward normal)
	} else if (std::abs(t - t_min.z) < epsilon) {
		normal = glm::dvec3(0, 0, -1); // Front face (inward normal)
	} else {
		normal = glm::dvec3(0, 0, 1);  // Back face (inward normal)
	}
	
	return t;
}
