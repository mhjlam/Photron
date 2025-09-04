#pragma once

#include <vector>
#include <optional>
#include "math/ray.hpp"
#include "math/triangle.hpp"
#include "math/glm_types.hpp"
#include "math/range.hpp"
#include "math/bvh.hpp"

// Forward declarations
class BoundingVolume;

/**
 * Represents a triangular mesh geometry for simulation.
 * Provides methods for ray-mesh intersection and point-in-mesh testing.
 * Now uses BVH acceleration for fast intersection testing.
 */
class MeshGeometry
{
private:
    std::vector<Triangle> triangles_;
    Range3 bounds_;
    bool has_bounds_;
    BVH bvh_; // Spatial acceleration structure
    
public:
    MeshGeometry() : has_bounds_(false) {}
    
    MeshGeometry(const std::vector<Triangle>& triangles)
        : triangles_(triangles), has_bounds_(false)
    {
        calculate_bounds();
        bvh_.build(triangles_);
    }
    
    // Delete copy constructor and copy assignment (non-copyable due to BVH)
    MeshGeometry(const MeshGeometry&) = delete;
    MeshGeometry& operator=(const MeshGeometry&) = delete;
    
    // Enable move constructor and move assignment
    MeshGeometry(MeshGeometry&&) = default;
    MeshGeometry& operator=(MeshGeometry&&) = default;
    
    // Add a triangle to the mesh
    void add_triangle(const Triangle& triangle) {
        triangles_.push_back(triangle);
        if (has_bounds_) {
            // Update bounds
            extend_bounds(triangle);
        }
        // Rebuild BVH when triangles are added
        bvh_.build(triangles_);
    }
    
    // Set all triangles at once
    void set_triangles(const std::vector<Triangle>& triangles) {
        triangles_ = triangles;
        calculate_bounds();
        bvh_.build(triangles_);
    }
    
    // Get the triangles in the mesh
    const std::vector<Triangle>& triangles() const { return triangles_; }
    
    // Get the bounding box of the mesh
    const Range3& bounds() const { return bounds_; }
    
    // Fast check if a point might be inside (bounding box test)
    bool might_contain(const glm::dvec3& point) const {
        return has_bounds_ && bounds_.includes(point.x, point.y, point.z);
    }
    
    // Check if a point is inside the mesh using BVH-accelerated ray casting
    bool contains(const glm::dvec3& point) const;
    
    // Find the first intersection of a ray with the mesh using BVH
    // Returns the distance to intersection and the triangle hit
    std::pair<double, Triangle> find_first_intersection(const Ray& ray) const;
    
    // Find all intersections of a ray with the mesh (uses BVH for initial filtering)
    std::vector<std::pair<double, Triangle>> find_all_intersections(const Ray& ray) const;
    
    // Check if the mesh has any triangles
    bool is_empty() const { return triangles_.empty(); }
    
    // Check if the mesh bounds have been calculated
    bool has_bounds() const { return has_bounds_; }
    
private:
    // Calculate the bounding box of the mesh
    void calculate_bounds();
    
    // Extend the current bounds to include a triangle
    void extend_bounds(const Triangle& triangle);
};
