#include "voxel_volume_calculator.hpp"
#include "random.hpp"
#include "simulator/simulator.hpp"
#include <random>

double VoxelVolumeCalculator::compute_volume_fraction_inside(
    const glm::dvec3& voxel_min,
    const glm::dvec3& voxel_max,
    const Simulator* simulator,
    int samples
) {
    if (samples <= 0 || !simulator) return 0.0;
    
    // Monte Carlo sampling
    std::random_device rd;
    std::mt19937 gen(rd());
    std::uniform_real_distribution<double> dis(0.0, 1.0);
    
    int inside_count = 0;
    
    for (int i = 0; i < samples; i++) {
        // Generate random point within voxel
        glm::dvec3 point = glm::dvec3(
            voxel_min.x + dis(gen) * (voxel_max.x - voxel_min.x),
            voxel_min.y + dis(gen) * (voxel_max.y - voxel_min.y),
            voxel_min.z + dis(gen) * (voxel_max.z - voxel_min.z)
        );
        
        if (simulator->is_point_inside_geometry(point)) {
            inside_count++;
        }
    }
    
    return static_cast<double>(inside_count) / static_cast<double>(samples);
}

double VoxelVolumeCalculator::compute_volume_fraction_inside_fast(
    const glm::dvec3& voxel_min,
    const glm::dvec3& voxel_max,
    const Simulator* simulator,
    int max_subdivisions
) {
    if (!simulator) return 0.0;
    return subdivide_and_test(voxel_min, voxel_max, simulator, 0, max_subdivisions);
}

double VoxelVolumeCalculator::subdivide_and_test(
    const glm::dvec3& min_corner,
    const glm::dvec3& max_corner,
    const Simulator* simulator,
    int depth,
    int max_depth
) {
    // Test corner points
    std::vector<glm::dvec3> corners = {
        min_corner,
        {max_corner.x, min_corner.y, min_corner.z},
        {min_corner.x, max_corner.y, min_corner.z},
        {max_corner.x, max_corner.y, min_corner.z},
        {min_corner.x, min_corner.y, max_corner.z},
        {max_corner.x, min_corner.y, max_corner.z},
        {min_corner.x, max_corner.y, max_corner.z},
        max_corner
    };
    
    int inside_corners = 0;
    for (const auto& corner : corners) {
        if (simulator->is_point_inside_geometry(corner)) {
            inside_corners++;
        }
    }
    
    // If all corners have same state and we're not at max depth, we're done
    if (inside_corners == 0) {
        return 0.0; // Fully outside
    } else if (inside_corners == 8) {
        return 1.0; // Fully inside
    } else if (depth >= max_depth) {
        // At max depth, estimate based on corner ratio
        return static_cast<double>(inside_corners) / 8.0;
    }
    
    // Mixed state - subdivide
    glm::dvec3 center = (min_corner + max_corner) * 0.5;
    
    // 8 sub-cubes
    std::vector<std::pair<glm::dvec3, glm::dvec3>> sub_cubes = {
        {min_corner, center},
        {{center.x, min_corner.y, min_corner.z}, {max_corner.x, center.y, center.z}},
        {{min_corner.x, center.y, min_corner.z}, {center.x, max_corner.y, center.z}},
        {{center.x, center.y, min_corner.z}, {max_corner.x, max_corner.y, center.z}},
        {{min_corner.x, min_corner.y, center.z}, {center.x, center.y, max_corner.z}},
        {{center.x, min_corner.y, center.z}, {max_corner.x, center.y, max_corner.z}},
        {{min_corner.x, center.y, center.z}, {center.x, max_corner.y, max_corner.z}},
        {center, max_corner}
    };
    
    double total_fraction = 0.0;
    for (const auto& sub_cube : sub_cubes) {
        total_fraction += subdivide_and_test(sub_cube.first, sub_cube.second, simulator, depth + 1, max_depth);
    }
    
    return total_fraction / 8.0; // Average of 8 sub-cubes
}
