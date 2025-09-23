#pragma once

#include <vector>
#include <cmath>
#include <algorithm>
#include <glm/glm.hpp>

/**
 * @brief 3D Digital Differential Analyzer (DDA) for robust voxel traversal
 * 
 * This implementation provides accurate and efficient ray-voxel intersection
 * using the 3D DDA algorithm, which avoids floating-point precision issues
 * that plague traditional geometric intersection methods.
 * 
 * The algorithm works by stepping through voxels along a ray path using
 * integer arithmetic and precise boundary calculations.
 */
class VoxelDDA3D {
public:
    /**
     * @brief Result of voxel traversal step
     */
    struct StepResult {
        glm::ivec3 voxel_coords;    // Current voxel coordinates (ix, iy, iz)
        glm::dvec3 world_position;  // World position at voxel entry
        glm::dvec3 entry_normal;    // Normal of the face we entered through
        double distance_traveled;   // Distance from ray origin to this point
        bool valid;                 // True if this is a valid voxel step
        
        StepResult() : voxel_coords(0), world_position(0.0), entry_normal(0.0), 
                      distance_traveled(0.0), valid(false) {}
    };
    
    /**
     * @brief Complete voxel traversal result
     */
    struct TraversalResult {
        std::vector<StepResult> voxels;  // All voxels traversed
        double total_distance;           // Total distance traveled
        bool hit_boundary;               // True if ray exited the grid
        glm::ivec3 exit_voxel;          // Last valid voxel before exit
        glm::dvec3 exit_position;       // Exit position from grid
        glm::dvec3 exit_normal;         // Normal of exit face
    };

private:
    // Grid parameters
    glm::ivec3 grid_dimensions_;    // Grid size (nx, ny, nz)
    glm::dvec3 grid_origin_;        // World coordinates of grid origin
    double voxel_size_;             // Size of each voxel
    glm::dvec3 grid_bounds_max_;    // Maximum bounds of the grid
    
    // Current state
    glm::dvec3 ray_origin_;         // Ray starting point
    glm::dvec3 ray_direction_;      // Ray direction (normalized)
    
    // DDA state variables
    glm::ivec3 current_voxel_;      // Current voxel coordinates
    glm::dvec3 delta_dist_;         // Distance ray travels for each unit step in each axis
    glm::dvec3 side_dist_;          // Distance from current position to next grid line
    glm::ivec3 step_;               // Direction to step (+1 or -1 for each axis)
    glm::ivec3 side_;               // Which axis had the shortest distance to next grid line
    double last_distance_;          // Track cumulative distance for per-voxel calculations
    
public:
    /**
     * @brief Construct DDA with grid parameters
     */
    VoxelDDA3D(const glm::ivec3& grid_dimensions, 
               const glm::dvec3& grid_origin, 
               double voxel_size);
    
    /**
     * @brief Initialize DDA for a specific ray
     */
    void initialize_ray(const glm::dvec3& origin, const glm::dvec3& direction);
    
    /**
     * @brief Perform single DDA step
     * @return StepResult for the next voxel, or invalid result if ray exits grid
     */
    StepResult step();
    
    /**
     * @brief Traverse entire ray path through voxel grid
     * @param max_distance Maximum distance to traverse (optional)
     * @return Complete traversal result
     */
    TraversalResult traverse(double max_distance = std::numeric_limits<double>::max());
    
    /**
     * @brief Find voxel containing a specific world position
     * @param position World coordinates
     * @return Voxel coordinates, or (-1,-1,-1) if outside grid
     */
    glm::ivec3 world_to_voxel(const glm::dvec3& position) const;
    
    /**
     * @brief Convert voxel coordinates to world position (center of voxel)
     * @param voxel_coords Voxel coordinates
     * @return World coordinates of voxel center
     */
    glm::dvec3 voxel_to_world(const glm::ivec3& voxel_coords) const;
    
    /**
     * @brief Check if voxel coordinates are within grid bounds
     */
    bool is_valid_voxel(const glm::ivec3& voxel_coords) const;
    
    /**
     * @brief Get the world bounds of a specific voxel
     */
    std::pair<glm::dvec3, glm::dvec3> get_voxel_bounds(const glm::ivec3& voxel_coords) const;

private:
    /**
     * @brief Calculate entry normal based on which face was crossed
     */
    glm::dvec3 calculate_entry_normal(int face_axis, int step_direction) const;
    
    /**
     * @brief Calculate world position from current DDA state
     */
    glm::dvec3 calculate_world_position() const;
};