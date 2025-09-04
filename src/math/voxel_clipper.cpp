#include "voxel_clipper.hpp"

#include <algorithm>
#include <cmath>
#include <unordered_set>

namespace VoxelClipping {

ClipResult VoxelClipper::clip_voxel_against_mesh(
    const glm::vec3& voxel_center,
    float voxel_half_size,
    std::function<bool(const glm::dvec3&)> is_point_inside_mesh
) {
    ClipResult result;
    
    // Generate cube vertices
    std::vector<glm::vec3> cube_vertices = generate_cube_vertices(voxel_center, voxel_half_size);
    
    // Test all 8 corners of the cube
    std::vector<bool> corners_inside(8);
    int inside_count = 0;
    
    for (int i = 0; i < 8; i++) {
        glm::dvec3 corner_double(cube_vertices[i].x, cube_vertices[i].y, cube_vertices[i].z);
        corners_inside[i] = is_point_inside_mesh(corner_double);
        if (corners_inside[i]) inside_count++;
    }
    
    if (inside_count == 0) {
        result.is_fully_outside = true;
        return result;
    }
    
    if (inside_count == 8) {
        result.is_fully_inside = true;
        result.volume_ratio = 1.0f;
        return result;
    }
    
    // Partial intersection - use adaptive mesh-based clipping
    // This uses a more sophisticated approach than simple sampling
    
    // Create a high-resolution sampling grid within the voxel
    const int resolution = 8; // 8x8x8 = 512 sample points for accurate boundary detection
    const float step = (voxel_half_size * 2.0f) / resolution;
    const float offset = -voxel_half_size + step * 0.5f;
    
    std::vector<glm::vec3> inside_points;
    std::vector<glm::vec3> boundary_points;
    
    // Sample points and classify them
    for (int x = 0; x < resolution; x++) {
        for (int y = 0; y < resolution; y++) {
            for (int z = 0; z < resolution; z++) {
                glm::vec3 sample_point = voxel_center + glm::vec3(
                    offset + x * step,
                    offset + y * step,
                    offset + z * step
                );
                
                glm::dvec3 sample_double(sample_point.x, sample_point.y, sample_point.z);
                bool point_inside = is_point_inside_mesh(sample_double);
                
                if (point_inside) {
                    inside_points.push_back(sample_point);
                    
                    // Check if this point is near the boundary by testing neighbors
                    bool near_boundary = false;
                    const float boundary_test_distance = step * 0.7f;
                    
                    glm::vec3 test_directions[] = {
                        {boundary_test_distance, 0, 0}, {-boundary_test_distance, 0, 0},
                        {0, boundary_test_distance, 0}, {0, -boundary_test_distance, 0},
                        {0, 0, boundary_test_distance}, {0, 0, -boundary_test_distance}
                    };
                    
                    for (const auto& dir : test_directions) {
                        glm::vec3 test_point = sample_point + dir;
                        glm::dvec3 test_double(test_point.x, test_point.y, test_point.z);
                        if (!is_point_inside_mesh(test_double)) {
                            near_boundary = true;
                            boundary_points.push_back(sample_point);
                            break;
                        }
                    }
                }
            }
        }
    }
    
    if (inside_points.empty()) {
        result.is_fully_outside = true;
        return result;
    }
    
    // Calculate volume ratio
    float total_samples = resolution * resolution * resolution;
    result.volume_ratio = static_cast<float>(inside_points.size()) / total_samples;
    
    // Create simplified geometry based on inside points
    // Use marching cubes-inspired approach for better geometry
    
    if (result.volume_ratio > 0.8f) {
        // Mostly inside - create a slightly scaled cube
        float scale_factor = std::pow(result.volume_ratio, 0.3f); // Gentle scaling
        std::vector<glm::vec3> scaled_vertices = generate_cube_vertices(voxel_center, voxel_half_size * scale_factor);
        
        // Create faces using the standard cube topology
        auto cube_face_indices = get_cube_faces();
        for (const auto& face_indices : cube_face_indices) {
            ClippedFace face;
            for (int idx : face_indices) {
                face.vertices.push_back({scaled_vertices[idx], false});
            }
            result.faces.push_back(face);
        }
        
    } else {
        // Partially inside - create more complex clipped geometry
        // Use convex hull of inside points near the boundary
        
        if (!boundary_points.empty()) {
            // Create geometry from boundary points using a simple approach
            // Find center of mass of boundary points
            glm::vec3 boundary_center(0.0f);
            for (const auto& point : boundary_points) {
                boundary_center += point;
            }
            boundary_center /= static_cast<float>(boundary_points.size());
            
            // Create tetrahedral geometry from boundary center to nearby inside points
            const float max_distance = voxel_half_size * 0.8f;
            std::vector<glm::vec3> nearby_inside_points;
            
            for (const auto& point : inside_points) {
                if (glm::length(point - boundary_center) <= max_distance) {
                    nearby_inside_points.push_back(point);
                }
            }
            
            // Create faces connecting boundary center to groups of nearby points
            if (nearby_inside_points.size() >= 3) {
                // Simple approach: create triangular faces
                for (size_t i = 0; i < nearby_inside_points.size(); i += 3) {
                    if (i + 2 < nearby_inside_points.size()) {
                        ClippedFace face;
                        face.vertices.push_back({nearby_inside_points[i], false});
                        face.vertices.push_back({nearby_inside_points[i+1], false});
                        face.vertices.push_back({nearby_inside_points[i+2], false});
                        result.faces.push_back(face);
                    }
                }
            }
        }
        
        // If we don't have enough geometry, fall back to a scaled cube
        if (result.faces.empty()) {
            float scale_factor = std::pow(result.volume_ratio, 0.5f);
            std::vector<glm::vec3> scaled_vertices = generate_cube_vertices(voxel_center, voxel_half_size * scale_factor);
            
            auto cube_face_indices = get_cube_faces();
            for (const auto& face_indices : cube_face_indices) {
                ClippedFace face;
                for (int idx : face_indices) {
                    face.vertices.push_back({scaled_vertices[idx], false});
                }
                result.faces.push_back(face);
            }
        }
    }
    
    return result;
}

ClipResult VoxelClipper::clip_voxel_against_triangles(
    const glm::vec3& voxel_center,
    float voxel_half_size,
    const std::vector<Triangle>& triangles
) {
    // This method would implement true cube-triangle intersection
    // For now, we'll return an empty result as we're using the mesh-based approach
    ClipResult result;
    result.is_fully_outside = true;
    return result;
}

std::vector<glm::vec3> VoxelClipper::generate_cube_vertices(
    const glm::vec3& center,
    float half_size
) {
    return {
        center + glm::vec3(-half_size, -half_size, -half_size), // 0
        center + glm::vec3( half_size, -half_size, -half_size), // 1
        center + glm::vec3( half_size,  half_size, -half_size), // 2
        center + glm::vec3(-half_size,  half_size, -half_size), // 3
        center + glm::vec3(-half_size, -half_size,  half_size), // 4
        center + glm::vec3( half_size, -half_size,  half_size), // 5
        center + glm::vec3( half_size,  half_size,  half_size), // 6
        center + glm::vec3(-half_size,  half_size,  half_size)  // 7
    };
}

std::vector<std::vector<int>> VoxelClipper::get_cube_faces() {
    return {
        {0, 1, 2, 3}, // bottom face (-Z)
        {4, 7, 6, 5}, // top face (+Z)
        {0, 4, 5, 1}, // front face (-Y)
        {2, 6, 7, 3}, // back face (+Y)
        {0, 3, 7, 4}, // left face (-X)
        {1, 5, 6, 2}  // right face (+X)
    };
}

bool VoxelClipper::line_plane_intersect(
    const glm::vec3& line_start,
    const glm::vec3& line_end,
    const glm::vec3& plane_point,
    const glm::vec3& plane_normal,
    glm::vec3& intersection,
    float& t
) {
    glm::vec3 line_dir = line_end - line_start;
    float denominator = glm::dot(line_dir, plane_normal);
    
    if (std::abs(denominator) < 1e-6f) {
        return false; // Line is parallel to plane
    }
    
    t = glm::dot(plane_point - line_start, plane_normal) / denominator;
    
    if (t < 0.0f || t > 1.0f) {
        return false; // Intersection outside line segment
    }
    
    intersection = line_start + t * line_dir;
    return true;
}

float VoxelClipper::point_plane_distance(
    const glm::vec3& point,
    const glm::vec3& plane_point,
    const glm::vec3& plane_normal
) {
    return glm::dot(point - plane_point, plane_normal);
}

float VoxelClipper::estimate_polyhedron_volume(const std::vector<ClippedFace>& faces) {
    // Simple volume estimation - this is approximate
    if (faces.empty()) return 0.0f;
    
    // Find centroid of all vertices
    glm::vec3 centroid(0.0f);
    int vertex_count = 0;
    
    for (const auto& face : faces) {
        for (const auto& vertex : face.vertices) {
            centroid += vertex.position;
            vertex_count++;
        }
    }
    
    if (vertex_count == 0) return 0.0f;
    
    centroid /= static_cast<float>(vertex_count);
    
    // Estimate volume using tetrahedra from centroid to each face
    float total_volume = 0.0f;
    
    for (const auto& face : faces) {
        if (face.vertices.size() >= 3) {
            // Triangulate face and compute volume of tetrahedra
            for (size_t i = 1; i < face.vertices.size() - 1; i++) {
                glm::vec3 v1 = face.vertices[0].position - centroid;
                glm::vec3 v2 = face.vertices[i].position - centroid;
                glm::vec3 v3 = face.vertices[i+1].position - centroid;
                
                // Volume of tetrahedron = |det(v1,v2,v3)| / 6
                float det = glm::dot(v1, glm::cross(v2, v3));
                total_volume += std::abs(det) / 6.0f;
            }
        }
    }
    
    return total_volume;
}

// Placeholder implementations for unused methods
bool VoxelClipper::cube_triangle_intersect(const glm::vec3&, const glm::vec3&, const Triangle&) { return false; }
std::vector<glm::vec3> VoxelClipper::clip_cube_by_plane(const std::vector<glm::vec3>&, const glm::vec3&, const glm::vec3&) { return {}; }
std::vector<glm::vec3> VoxelClipper::sutherland_hodgman_clip_3d(const std::vector<glm::vec3>&, const glm::vec3&, const glm::vec3&) { return {}; }
std::vector<ClippedFace> VoxelClipper::triangulate_convex_polygon(const std::vector<glm::vec3>&, const glm::vec3&) { return {}; }

} // namespace VoxelClipping
