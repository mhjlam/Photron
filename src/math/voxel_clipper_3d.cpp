#include "voxel_clipper_3d.hpp"
#include "triangle.hpp"
#include <algorithm>
#include <cmath>

const float EPSILON = 1e-6f;

ClippedVoxel VoxelClipper3D::clip_voxel_against_mesh(
    const glm::vec3& voxel_center,
    float voxel_half_size,
    const std::vector<Triangle>& triangles,
    const glm::vec4& voxel_color
) {
    ClippedVoxel result;
    result.color = voxel_color;
    
    // Start with the original cube vertices
    std::vector<glm::vec3> current_vertices = create_cube_vertices(voxel_center, voxel_half_size);
    
    if (current_vertices.empty()) {
        return result; // Empty result
    }
    
    // For each triangle, clip the current polyhedron against the triangle's plane
    for (const auto& triangle : triangles) {
        if (current_vertices.empty()) {
            break; // Nothing left to clip
        }
        
        // Get triangle vertices and compute plane normal
        glm::vec3 v0 = glm::vec3(triangle.v0().x, triangle.v0().y, triangle.v0().z);
        glm::vec3 v1 = glm::vec3(triangle.v1().x, triangle.v1().y, triangle.v1().z);
        glm::vec3 v2 = glm::vec3(triangle.v2().x, triangle.v2().y, triangle.v2().z);
        
        glm::vec3 edge1 = v1 - v0;
        glm::vec3 edge2 = v2 - v0;
        glm::vec3 normal = glm::normalize(glm::cross(edge1, edge2));
        
        // Clip current vertices against this plane
        current_vertices = clip_polygon_by_plane(current_vertices, v0, normal);
        
        if (current_vertices.size() < 3) {
            // Not enough vertices for a valid polyhedron
            current_vertices.clear();
            break;
        }
    }
    
    if (current_vertices.size() >= 4) { // Minimum for a tetrahedron
        result.vertices = current_vertices;
        
        // Create triangles using convex hull approach
        // For simplicity, use fan triangulation from first vertex
        result.triangles = triangulate_convex_polygon(current_vertices.size(), 0);
    }
    
    return result;
}

std::vector<glm::vec3> VoxelClipper3D::clip_polygon_by_plane(
    const std::vector<glm::vec3>& polygon,
    const glm::vec3& plane_point,
    const glm::vec3& plane_normal
) {
    if (polygon.size() < 3) {
        return {}; // Invalid polygon
    }
    
    std::vector<glm::vec3> clipped;
    
    for (size_t i = 0; i < polygon.size(); i++) {
        const glm::vec3& current = polygon[i];
        const glm::vec3& next = polygon[(i + 1) % polygon.size()];
        
        bool current_inside = is_point_inside_plane(current, plane_point, plane_normal);
        bool next_inside = is_point_inside_plane(next, plane_point, plane_normal);
        
        if (current_inside) {
            // Current vertex is inside
            clipped.push_back(current);
            
            if (!next_inside) {
                // Leaving the inside - add intersection point
                glm::vec3 intersection = intersect_line_plane(current, next, plane_point, plane_normal);
                clipped.push_back(intersection);
            }
        } else if (next_inside) {
            // Entering the inside - add intersection point
            glm::vec3 intersection = intersect_line_plane(current, next, plane_point, plane_normal);
            clipped.push_back(intersection);
        }
        // If both outside, add nothing
    }
    
    return clipped;
}

bool VoxelClipper3D::is_point_inside_plane(
    const glm::vec3& point,
    const glm::vec3& plane_point,
    const glm::vec3& plane_normal
) {
    glm::vec3 to_point = point - plane_point;
    float dot_product = glm::dot(to_point, plane_normal);
    return dot_product <= EPSILON; // Inside if on or behind the plane
}

glm::vec3 VoxelClipper3D::intersect_line_plane(
    const glm::vec3& p1,
    const glm::vec3& p2,
    const glm::vec3& plane_point,
    const glm::vec3& plane_normal
) {
    glm::vec3 line_dir = p2 - p1;
    float denom = glm::dot(line_dir, plane_normal);
    
    if (std::abs(denom) < EPSILON) {
        // Line is parallel to plane
        return p1; // Return arbitrary point
    }
    
    glm::vec3 to_plane = plane_point - p1;
    float t = glm::dot(to_plane, plane_normal) / denom;
    
    // Clamp t to [0, 1] to stay within line segment
    t = std::clamp(t, 0.0f, 1.0f);
    
    return p1 + t * line_dir;
}

std::vector<glm::vec3> VoxelClipper3D::create_cube_vertices(
    const glm::vec3& center,
    float half_size
) {
    return {
        glm::vec3(center.x - half_size, center.y - half_size, center.z - half_size), // 0
        glm::vec3(center.x + half_size, center.y - half_size, center.z - half_size), // 1
        glm::vec3(center.x + half_size, center.y + half_size, center.z - half_size), // 2
        glm::vec3(center.x - half_size, center.y + half_size, center.z - half_size), // 3
        glm::vec3(center.x - half_size, center.y - half_size, center.z + half_size), // 4
        glm::vec3(center.x + half_size, center.y - half_size, center.z + half_size), // 5
        glm::vec3(center.x + half_size, center.y + half_size, center.z + half_size), // 6
        glm::vec3(center.x - half_size, center.y + half_size, center.z + half_size)  // 7
    };
}

std::vector<glm::ivec3> VoxelClipper3D::triangulate_convex_polygon(int vertex_count, int start_index) {
    std::vector<glm::ivec3> triangles;
    
    if (vertex_count < 3) {
        return triangles;
    }
    
    // Simple fan triangulation from the first vertex
    for (int i = 1; i < vertex_count - 1; i++) {
        triangles.emplace_back(
            start_index,
            start_index + i,
            start_index + i + 1
        );
    }
    
    return triangles;
}
