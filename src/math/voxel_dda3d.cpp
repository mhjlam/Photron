#include "voxel_dda3d.hpp"
#include <limits>
#include <cmath>

VoxelDDA3D::VoxelDDA3D(const glm::ivec3& grid_dimensions, 
                       const glm::dvec3& grid_origin, 
                       double voxel_size)
    : grid_dimensions_(grid_dimensions)
    , grid_origin_(grid_origin)
    , voxel_size_(voxel_size)
{
    // Calculate grid bounds
    grid_bounds_max_ = grid_origin_ + glm::dvec3(grid_dimensions_) * voxel_size_;
}

void VoxelDDA3D::initialize_ray(const glm::dvec3& origin, const glm::dvec3& direction) {
    ray_origin_ = origin;
    ray_direction_ = glm::normalize(direction);
    last_distance_ = 0.0; // Initialize distance tracking
    
    // Convert ray origin to grid coordinates
    glm::dvec3 ray_pos = (ray_origin_ - grid_origin_) / voxel_size_;
    
    // Current voxel coordinates (integer part)
    current_voxel_ = glm::ivec3(floor(ray_pos.x), floor(ray_pos.y), floor(ray_pos.z));
    
    // Calculate delta distances - how far along the ray we must travel 
    // for the ray to cross one voxel boundary in each axis
    for (int i = 0; i < 3; ++i) {
        if (std::abs(ray_direction_[i]) < 1e-15) {
            delta_dist_[i] = std::numeric_limits<double>::max();
        } else {
            delta_dist_[i] = std::abs(1.0 / ray_direction_[i]) * voxel_size_;
        }
    }
    
    // Calculate step direction and initial side_dist for each axis
    for (int i = 0; i < 3; ++i) {
        if (ray_direction_[i] < 0) {
            step_[i] = -1;
            side_dist_[i] = (ray_pos[i] - current_voxel_[i]) * delta_dist_[i];
        } else {
            step_[i] = 1;
            side_dist_[i] = (current_voxel_[i] + 1.0 - ray_pos[i]) * delta_dist_[i];
        }
    }
}

VoxelDDA3D::StepResult VoxelDDA3D::step() {
    StepResult result;
    
    // Check if current voxel is within bounds
    if (!is_valid_voxel(current_voxel_)) {
        result.valid = false;
        return result;
    }
    
    // Record current voxel information
    result.voxel_coords = current_voxel_;
    result.world_position = calculate_world_position();
    result.valid = true;
    
    // Find which axis has the shortest distance to next grid line
    int min_axis = 0;
    double min_dist = side_dist_[0];
    
    for (int i = 1; i < 3; ++i) {
        if (side_dist_[i] < min_dist) {
            min_dist = side_dist_[i];
            min_axis = i;
        }
    }
    
    // Calculate entry normal (face we're coming from)
    result.entry_normal = calculate_entry_normal(min_axis, -step_[min_axis]);
    
    // Distance traveled IN THIS VOXEL (not cumulative)
    double previous_distance = last_distance_;
    result.distance_traveled = min_dist - previous_distance;
    last_distance_ = min_dist;
    
    // Perform DDA step - move to next voxel
    side_dist_[min_axis] += delta_dist_[min_axis];
    current_voxel_[min_axis] += step_[min_axis];
    side_[min_axis] = 1; // Mark which side we stepped through
    
    return result;
}

VoxelDDA3D::TraversalResult VoxelDDA3D::traverse(double max_distance) {
    TraversalResult result;
    result.total_distance = 0.0;
    result.hit_boundary = false;
    
    while (result.total_distance < max_distance) {
        StepResult step_result = step();
        
        if (!step_result.valid) {
            // Ray has exited the grid
            result.hit_boundary = true;
            if (!result.voxels.empty()) {
                const auto& last_voxel = result.voxels.back();
                result.exit_voxel = last_voxel.voxel_coords;
                result.exit_position = last_voxel.world_position;
                
                // Calculate exit normal (opposite of last entry normal)
                result.exit_normal = -last_voxel.entry_normal;
            }
            break;
        }
        
        // Check if adding this voxel would exceed max_distance
        if (result.total_distance + step_result.distance_traveled > max_distance) {
            // Trim the distance to not exceed max_distance
            step_result.distance_traveled = max_distance - result.total_distance;
            result.total_distance = max_distance;
            result.voxels.push_back(step_result);
            break;
        }
        
        result.total_distance += step_result.distance_traveled;
        result.voxels.push_back(step_result);
    }
    
    return result;
}

glm::ivec3 VoxelDDA3D::world_to_voxel(const glm::dvec3& position) const {
    glm::dvec3 grid_pos = (position - grid_origin_) / voxel_size_;
    
    glm::ivec3 voxel_coords(
        static_cast<int>(floor(grid_pos.x)),
        static_cast<int>(floor(grid_pos.y)),
        static_cast<int>(floor(grid_pos.z))
    );
    
    // Check bounds
    if (!is_valid_voxel(voxel_coords)) {
        return glm::ivec3(-1, -1, -1); // Invalid voxel
    }
    
    return voxel_coords;
}

glm::dvec3 VoxelDDA3D::voxel_to_world(const glm::ivec3& voxel_coords) const {
    // Return center of voxel
    return grid_origin_ + (glm::dvec3(voxel_coords) + 0.5) * voxel_size_;
}

bool VoxelDDA3D::is_valid_voxel(const glm::ivec3& voxel_coords) const {
    return (voxel_coords.x >= 0 && voxel_coords.x < grid_dimensions_.x &&
            voxel_coords.y >= 0 && voxel_coords.y < grid_dimensions_.y &&
            voxel_coords.z >= 0 && voxel_coords.z < grid_dimensions_.z);
}

std::pair<glm::dvec3, glm::dvec3> VoxelDDA3D::get_voxel_bounds(const glm::ivec3& voxel_coords) const {
    glm::dvec3 min_bound = grid_origin_ + glm::dvec3(voxel_coords) * voxel_size_;
    glm::dvec3 max_bound = min_bound + glm::dvec3(voxel_size_);
    
    return std::make_pair(min_bound, max_bound);
}

glm::dvec3 VoxelDDA3D::calculate_entry_normal(int face_axis, int step_direction) const {
    glm::dvec3 normal(0.0);
    normal[face_axis] = static_cast<double>(step_direction);
    return normal;
}

glm::dvec3 VoxelDDA3D::calculate_world_position() const {
    // Calculate current position along the ray
    double min_dist = std::min({side_dist_[0] - delta_dist_[0], 
                                side_dist_[1] - delta_dist_[1], 
                                side_dist_[2] - delta_dist_[2]});
    
    return ray_origin_ + ray_direction_ * min_dist;
}
