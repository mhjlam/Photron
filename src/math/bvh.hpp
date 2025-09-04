#pragma once

#include "glm_types.hpp"
#include "triangle.hpp"
#include "ray.hpp"
#include <vector>
#include <memory>

/**
 * Axis-Aligned Bounding Box for BVH nodes
 */
struct AABB {
    glm::dvec3 min_point;
    glm::dvec3 max_point;
    
    AABB() : min_point(DBL_MAX), max_point(-DBL_MAX) {}
    AABB(const glm::dvec3& min_pt, const glm::dvec3& max_pt) 
        : min_point(min_pt), max_point(max_pt) {}
    
    // Expand to include another AABB
    void expand(const AABB& other) {
        min_point = glm::min(min_point, other.min_point);
        max_point = glm::max(max_point, other.max_point);
    }
    
    // Expand to include a point
    void expand(const glm::dvec3& point) {
        min_point = glm::min(min_point, point);
        max_point = glm::max(max_point, point);
    }
    
    // Check if ray intersects this AABB
    bool intersect_ray(const Ray& ray, double& t_min, double& t_max) const;
    
    // Check if point is inside AABB
    bool contains(const glm::dvec3& point) const {
        return point.x >= min_point.x && point.x <= max_point.x &&
               point.y >= min_point.y && point.y <= max_point.y &&
               point.z >= min_point.z && point.z <= max_point.z;
    }
    
    // Get the center of the AABB
    glm::dvec3 center() const {
        return (min_point + max_point) * 0.5;
    }
    
    // Get the surface area (for SAH heuristic)
    double surface_area() const {
        glm::dvec3 extent = max_point - min_point;
        return 2.0 * (extent.x * extent.y + extent.y * extent.z + extent.z * extent.x);
    }
};

/**
 * BVH Node - can be either internal or leaf
 */
struct BVHNode {
    AABB bounds;
    std::unique_ptr<BVHNode> left;
    std::unique_ptr<BVHNode> right;
    std::vector<int> triangle_indices; // Only used for leaf nodes
    
    bool is_leaf() const {
        return left == nullptr && right == nullptr;
    }
};

/**
 * Bounding Volume Hierarchy for fast ray-triangle intersection
 * and point-in-mesh testing
 */
class BVH {
private:
    std::unique_ptr<BVHNode> root_;
    std::vector<Triangle> triangles_;
    
    // Build the BVH recursively
    std::unique_ptr<BVHNode> build_recursive(std::vector<int>& triangle_indices, int depth = 0);
    
    // Calculate AABB for a set of triangles
    AABB calculate_bounds(const std::vector<int>& triangle_indices) const;
    
    // Find best split using Surface Area Heuristic (SAH)
    int find_best_split(std::vector<int>& triangle_indices, int& best_axis) const;
    
    // Ray intersection traversal
    bool intersect_recursive(const BVHNode* node, const Ray& ray, double& closest_t, 
                           glm::dvec3& intersection, Triangle& hit_triangle) const;
    
    // Point-in-mesh testing traversal
    int count_intersections_recursive(const BVHNode* node, const Ray& ray, const glm::dvec3& point) const;
    
public:
    // Maximum triangles per leaf node
    static const int MAX_TRIANGLES_PER_LEAF = 4;
    
    // Maximum tree depth
    static const int MAX_DEPTH = 20;
    
    BVH() = default;
    BVH(const std::vector<Triangle>& triangles);
    
    // Delete copy constructor and copy assignment (non-copyable due to unique_ptr)
    BVH(const BVH&) = delete;
    BVH& operator=(const BVH&) = delete;
    
    // Enable move constructor and move assignment
    BVH(BVH&&) = default;
    BVH& operator=(BVH&&) = default;
    
    // Build the BVH from triangles
    void build(const std::vector<Triangle>& triangles);
    
    // Find closest ray-triangle intersection
    bool intersect(const Ray& ray, double& distance, glm::dvec3& intersection, Triangle& hit_triangle) const;
    
    // Test if point is inside the mesh using ray casting
    bool contains_point(const glm::dvec3& point) const;
    
    // Get all triangles (for debugging)
    const std::vector<Triangle>& triangles() const { return triangles_; }
    
    // Check if BVH is empty
    bool empty() const { return triangles_.empty() || !root_; }
};
