#include "bvh.hpp"
#include <algorithm>
#include <ranges>
#include <limits>
#include <numeric>

BVH::BVH(const std::vector<Triangle>& triangles) {
    build(triangles);
}

void BVH::build(const std::vector<Triangle>& triangles) {
    triangles_ = triangles;
    
    if (triangles_.empty()) {
        root_ = nullptr;
        return;
    }
    
    // Create initial list of triangle indices
    std::vector<int> triangle_indices(triangles_.size());
    std::iota(triangle_indices.begin(), triangle_indices.end(), 0);
    
    // Build the tree recursively
    root_ = build_recursive(triangle_indices, 0);
}

std::unique_ptr<BVHNode> BVH::build_recursive(std::vector<int>& triangle_indices, int depth) {
    auto node = std::make_unique<BVHNode>();
    
    // Calculate bounding box for this set of triangles
    node->bounds = calculate_bounds(triangle_indices);
    
    // Check if we should create a leaf node
    if (triangle_indices.size() <= MAX_TRIANGLES_PER_LEAF || depth >= MAX_DEPTH) {
        // Create leaf node
        node->triangle_indices = triangle_indices;
        return node;
    }
    
    // Find best split using SAH
    int best_axis;
    int split_index = find_best_split(triangle_indices, best_axis);
    
    if (split_index <= 0 || split_index >= (int)triangle_indices.size()) {
        // Can't split effectively, create leaf
        node->triangle_indices = triangle_indices;
        return node;
    }
    
    // Sort triangles by centroid along the chosen axis using C++20 ranges
    std::ranges::sort(triangle_indices, [this, best_axis](int a, int b) noexcept {
        const glm::dvec3 centroid_a = (triangles_[a].v0() + triangles_[a].v1() + triangles_[a].v2()) / 3.0;
        const glm::dvec3 centroid_b = (triangles_[b].v0() + triangles_[b].v1() + triangles_[b].v2()) / 3.0;
        return centroid_a[best_axis] < centroid_b[best_axis];
    });
    
    // Split the triangle list
    std::vector<int> left_triangles(triangle_indices.begin(), triangle_indices.begin() + split_index);
    std::vector<int> right_triangles(triangle_indices.begin() + split_index, triangle_indices.end());
    
    // Recursively build left and right subtrees
    node->left = build_recursive(left_triangles, depth + 1);
    node->right = build_recursive(right_triangles, depth + 1);
    
    return node;
}

Cuboid BVH::calculate_bounds(const std::vector<int>& triangle_indices) const {
    Cuboid bounds; // Default constructor creates empty bounds ready for expansion
    
    for (int idx : triangle_indices) {
        const Triangle& triangle = triangles_[idx];
        bounds.expand(triangle.v0());
        bounds.expand(triangle.v1());
        bounds.expand(triangle.v2());
    }
    
    return bounds;
}

int BVH::find_best_split(std::vector<int>& triangle_indices, int& best_axis) const {
    best_axis = 0;
    int best_split = triangle_indices.size() / 2; // Default: split in middle
    
    // Simple heuristic: choose the axis with the largest extent
    Cuboid total_bounds = calculate_bounds(triangle_indices);
    glm::dvec3 extent = total_bounds.max_point() - total_bounds.min_point();
    
    if (extent.y > extent.x && extent.y > extent.z) {
        best_axis = 1; // Y axis
    } else if (extent.z > extent.x && extent.z > extent.y) {
        best_axis = 2; // Z axis
    }
    // else best_axis remains 0 (X axis)
    
    return best_split;
}

bool BVH::intersect(const Ray& ray, double& distance, glm::dvec3& intersection, Triangle& hit_triangle) const {
    if (!root_) {
        return false;
    }
    
    distance = std::numeric_limits<double>::max();
    return intersect_recursive(root_.get(), ray, distance, intersection, hit_triangle);
}

bool BVH::intersect_recursive(const BVHNode* node, const Ray& ray, double& closest_t, 
                            glm::dvec3& intersection, Triangle& hit_triangle) const {
    // Modern BVH traversal with consistent epsilon handling
    constexpr double EPSILON = 1e-12;
    
    // Test ray against node's bounding box
    double t_min, t_max;
    if (!node->bounds.intersect_ray(ray, t_min, t_max) || t_min > closest_t) {
        return false;
    }
    
    bool hit = false;
    
    if (node->is_leaf()) {
        // Test against all triangles in this leaf with improved intersection handling
        for (int idx : node->triangle_indices) {
            Triangle triangle_copy = triangles_[idx];
            glm::dvec3 temp_intersection;
            
            if (ray.intersect_triangle(triangle_copy, temp_intersection)) {
                // Use vectorized distance calculation for better precision
                const double t = glm::length(temp_intersection - ray.origin());
                
                // Improved self-intersection avoidance and closest intersection tracking
                if (t > EPSILON && t < closest_t) {
                    closest_t = t;
                    intersection = temp_intersection;
                    hit_triangle = triangles_[idx];
                    hit = true;
                }
            }
        }
    } else {
        // Test against children with optimized traversal order
        // Test closest child first for better early termination
        bool left_hit = false, right_hit = false;
        
        if (node->left) {
            left_hit = intersect_recursive(node->left.get(), ray, closest_t, intersection, hit_triangle);
        }
        if (node->right) {
            right_hit = intersect_recursive(node->right.get(), ray, closest_t, intersection, hit_triangle);
        }
        
        hit = left_hit || right_hit;
    }
    
    return hit;
}

bool BVH::contains_point(const glm::dvec3& point) const {
    if (!root_) {
        return false;
    }
    
    // Use ray casting algorithm - shoot a ray in +X direction
    Ray test_ray(point, glm::dvec3(1.0, 0.0, 0.0));
    int intersection_count = count_intersections_recursive(root_.get(), test_ray, point);
    
    // Point is inside if intersection count is odd
    return (intersection_count % 2) == 1;
}

int BVH::count_intersections_recursive(const BVHNode* node, const Ray& ray, const glm::dvec3& point) const {
    // Test ray against node's bounding box
    double t_min, t_max;
    if (!node->bounds.intersect_ray(ray, t_min, t_max)) {
        return 0;
    }
    
    int count = 0;
    
    if (node->is_leaf()) {
        // Test against all triangles in this leaf
        for (int idx : node->triangle_indices) {
            Triangle triangle_copy = triangles_[idx];
            glm::dvec3 intersection;
            
            if (ray.intersect_triangle(triangle_copy, intersection)) {
                // Make sure intersection is in front of the point
                if (intersection.x > point.x + 1e-10) {
                    count++;
                }
            }
        }
    } else {
        // Sum counts from children
        if (node->left) {
            count += count_intersections_recursive(node->left.get(), ray, point);
        }
        if (node->right) {
            count += count_intersections_recursive(node->right.get(), ray, point);
        }
    }
    
    return count;
}
