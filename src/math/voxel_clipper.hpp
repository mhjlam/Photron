#pragma once

#include <vector>
#include <functional>
#include <glm/glm.hpp>
#include "math/triangle.hpp"

namespace VoxelClipping {

struct ClippedVertex {
    glm::vec3 position;
    bool is_intersection = false; // true if this vertex is from triangle intersection
};

struct ClippedFace {
    std::vector<ClippedVertex> vertices;
};

// Result of clipping a voxel against geometry
struct ClipResult {
    std::vector<ClippedFace> faces;        // Triangulated faces of clipped volume
    std::vector<ClippedFace> boundary_caps; // Surface caps at geometry boundary
    bool is_fully_inside = false;          // true if voxel is completely inside
    bool is_fully_outside = false;         // true if voxel is completely outside
    float volume_ratio = 0.0f;             // Approximate ratio of clipped volume to original
};

class VoxelClipper {
public:
    VoxelClipper() = default;
    ~VoxelClipper() = default;

    // Main clipping function - clips a voxel cube against a set of triangles
    ClipResult clip_voxel_against_triangles(
        const glm::vec3& voxel_center,
        float voxel_half_size,
        const std::vector<Triangle>& triangles
    );

    // Alternative interface using mesh geometry point-in-mesh testing
    ClipResult clip_voxel_against_mesh(
        const glm::vec3& voxel_center,
        float voxel_half_size,
        std::function<bool(const glm::dvec3&)> is_point_inside_mesh
    );

private:
    // Core geometric algorithms
    bool cube_triangle_intersect(
        const glm::vec3& cube_min,
        const glm::vec3& cube_max,
        const Triangle& triangle
    );

    std::vector<glm::vec3> clip_cube_by_plane(
        const std::vector<glm::vec3>& cube_vertices,
        const glm::vec3& plane_point,
        const glm::vec3& plane_normal
    );

    // Sutherland-Hodgman polygon clipping adapted for 3D
    std::vector<glm::vec3> sutherland_hodgman_clip_3d(
        const std::vector<glm::vec3>& input_vertices,
        const glm::vec3& plane_point,
        const glm::vec3& plane_normal
    );

    // Line-plane intersection
    bool line_plane_intersect(
        const glm::vec3& line_start,
        const glm::vec3& line_end,
        const glm::vec3& plane_point,
        const glm::vec3& plane_normal,
        glm::vec3& intersection,
        float& t
    );

    // Point-plane classification
    float point_plane_distance(
        const glm::vec3& point,
        const glm::vec3& plane_point,
        const glm::vec3& plane_normal
    );

    // Triangulation helpers
    std::vector<ClippedFace> triangulate_convex_polygon(
        const std::vector<glm::vec3>& vertices,
        const glm::vec3& normal
    );

    // Generate standard cube vertices and faces
    std::vector<glm::vec3> generate_cube_vertices(
        const glm::vec3& center,
        float half_size
    );

    std::vector<std::vector<int>> get_cube_faces();

    // Volume estimation
    float estimate_polyhedron_volume(const std::vector<ClippedFace>& faces);
};

} // namespace VoxelClipping
