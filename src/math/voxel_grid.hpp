#pragma once

#include <vector>
#include <cstdint>

#include "glm_types.hpp"

// Forward declarations
struct Voxel;

/**
 * @brief Modern 3D voxel grid for spatial discretization in Monte Carlo simulations
 * 
 * This class provides a uniform 3D grid structure for organizing voxels in space.
 * It uses proper encapsulation with private data members and provides efficient
 * indexing and access operations for photon transport simulations.
 */
class VoxelGrid
{
public:
    // Constructors and destructor
    VoxelGrid();
    VoxelGrid(double voxel_size, uint32_t num_x, uint32_t num_y, uint32_t num_z);
    ~VoxelGrid();
    
    // Delete copy constructor and copy assignment (non-copyable due to raw pointers)
    VoxelGrid(const VoxelGrid&) = delete;
    VoxelGrid& operator=(const VoxelGrid&) = delete;
    
    // Enable move constructor and move assignment
    VoxelGrid(VoxelGrid&& other) noexcept;
    VoxelGrid& operator=(VoxelGrid&& other) noexcept;
    
    // Grid access methods
    Voxel* operator()(uint32_t x, uint32_t y, uint32_t z);
    const Voxel* operator()(uint32_t x, uint32_t y, uint32_t z) const;
    
    Voxel* at(uint32_t x, uint32_t y, uint32_t z);
    const Voxel* at(uint32_t x, uint32_t y, uint32_t z) const;
    
    // Linear indexing
    uint32_t calculate_index(uint32_t x, uint32_t y, uint32_t z) const;
    Voxel* at(uint32_t linear_index);
    const Voxel* at(uint32_t linear_index) const;
    
    // Grid properties (const accessors)
    uint32_t width() const { return dimensions_.x; }
    uint32_t height() const { return dimensions_.y; }
    uint32_t depth() const { return dimensions_.z; }
    uint64_t size() const { return total_voxels_; }
    
    uint64_t total_voxels() const { return total_voxels_; }
    double voxel_size() const { return voxel_size_; }
    
    const glm::uvec3& dimensions() const { return dimensions_; }
    
    // Grid management
    void clear();
    bool is_valid_coordinate(uint32_t x, uint32_t y, uint32_t z) const;
    bool is_empty() const;
    
    // Iterator support for range-based loops
    std::vector<Voxel*>::iterator begin() { return voxels_.begin(); }
    std::vector<Voxel*>::iterator end() { return voxels_.end(); }
    std::vector<Voxel*>::const_iterator begin() const { return voxels_.begin(); }
    std::vector<Voxel*>::const_iterator end() const { return voxels_.end(); }
    std::vector<Voxel*>::const_iterator cbegin() const { return voxels_.cbegin(); }
    std::vector<Voxel*>::const_iterator cend() const { return voxels_.cend(); }

private:
    // Private data members
    glm::uvec3 dimensions_;          ///< Grid dimensions (nx, ny, nz)
    uint64_t total_voxels_;          ///< Total number of voxels (nx * ny * nz)
    double voxel_size_;              ///< Size of each voxel (uniform cubic voxels)
    std::vector<Voxel*> voxels_;     ///< Linear storage for voxel pointers
    
    // Private helper methods
    void initialize_voxels();
    void cleanup_voxels();
    bool is_valid_index(uint32_t linear_index) const;
};
