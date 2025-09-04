#pragma once

#include "glm_types.hpp"
#include "triangle.hpp"
#include <vector>

struct ClippedVoxel {
    std::vector<glm::vec3> vertices;
    std::vector<glm::ivec3> triangles; // Indices into vertices array
    glm::vec4 color;
    
    bool is_empty() const {
        return vertices.empty() || triangles.empty();
    }
};

class VoxelClipper3D {
public:
    // Clip a voxel cube against a set of triangular faces
    static ClippedVoxel clip_voxel_against_mesh(
        const glm::vec3& voxel_center,
        float voxel_half_size,
        const std::vector<Triangle>& triangles,
        const glm::vec4& voxel_color
    );

private:
    // Clip a convex polyhedron against a plane
    static std::vector<glm::vec3> clip_polygon_by_plane(
        const std::vector<glm::vec3>& polygon,
        const glm::vec3& plane_point,
        const glm::vec3& plane_normal
    );
    
    // Check if a point is inside (behind) a plane
    static bool is_point_inside_plane(
        const glm::vec3& point,
        const glm::vec3& plane_point,
        const glm::vec3& plane_normal
    );
    
    // Compute intersection of line segment with plane
    static glm::vec3 intersect_line_plane(
        const glm::vec3& p1,
        const glm::vec3& p2,
        const glm::vec3& plane_point,
        const glm::vec3& plane_normal
    );
    
    // Create initial cube vertices
    static std::vector<glm::vec3> create_cube_vertices(
        const glm::vec3& center,
        float half_size
    );
    
    // Triangulate a convex polygon (fan triangulation)
    static std::vector<glm::ivec3> triangulate_convex_polygon(int vertex_count, int start_index);
};
