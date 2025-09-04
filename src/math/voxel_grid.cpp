#include "voxel_grid.hpp"
#include "../simulator/voxel.hpp"
#include <stdexcept>
#include <algorithm>

// Default constructor
VoxelGrid::VoxelGrid() 
    : dimensions_(0, 0, 0)
    , total_voxels_(0)
    , voxel_size_(0.0)
    , voxels_()
{
}

// Parameterized constructor
VoxelGrid::VoxelGrid(double voxel_size, uint32_t num_x, uint32_t num_y, uint32_t num_z)
    : dimensions_(num_x, num_y, num_z)
    , total_voxels_(static_cast<uint64_t>(num_x) * num_y * num_z)
    , voxel_size_(voxel_size)
    , voxels_()
{
    if (voxel_size <= 0.0) {
        throw std::invalid_argument("Voxel size must be positive");
    }
    
    if (num_x == 0 || num_y == 0 || num_z == 0) {
        throw std::invalid_argument("Grid dimensions must be positive");
    }
    
    // Check for potential overflow
    if (total_voxels_ / num_x / num_y != num_z) {
        throw std::overflow_error("Grid dimensions would cause integer overflow");
    }
    
    initialize_voxels();
}

// Destructor
VoxelGrid::~VoxelGrid() {
    cleanup_voxels();
}

// Move constructor
VoxelGrid::VoxelGrid(VoxelGrid&& other) noexcept
    : dimensions_(other.dimensions_)
    , total_voxels_(other.total_voxels_)
    , voxel_size_(other.voxel_size_)
    , voxels_(std::move(other.voxels_))
{
    // Reset the moved-from object
    other.dimensions_ = glm::uvec3(0, 0, 0);
    other.total_voxels_ = 0;
    other.voxel_size_ = 0.0;
}

// Move assignment operator
VoxelGrid& VoxelGrid::operator=(VoxelGrid&& other) noexcept {
    if (this != &other) {
        // Clean up current resources
        cleanup_voxels();
        
        // Move data from other
        dimensions_ = other.dimensions_;
        total_voxels_ = other.total_voxels_;
        voxel_size_ = other.voxel_size_;
        voxels_ = std::move(other.voxels_);
        
        // Reset the moved-from object
        other.dimensions_ = glm::uvec3(0, 0, 0);
        other.total_voxels_ = 0;
        other.voxel_size_ = 0.0;
    }
    return *this;
}

// Grid access operators
Voxel* VoxelGrid::operator()(uint32_t x, uint32_t y, uint32_t z) {
    if (!is_valid_coordinate(x, y, z)) {
        throw std::out_of_range("Voxel coordinate out of bounds");
    }
    return voxels_[calculate_index(x, y, z)];
}

const Voxel* VoxelGrid::operator()(uint32_t x, uint32_t y, uint32_t z) const {
    if (!is_valid_coordinate(x, y, z)) {
        throw std::out_of_range("Voxel coordinate out of bounds");
    }
    return voxels_[calculate_index(x, y, z)];
}

// Safe access methods with bounds checking
Voxel* VoxelGrid::at(uint32_t x, uint32_t y, uint32_t z) {
    if (!is_valid_coordinate(x, y, z)) {
        throw std::out_of_range("Voxel coordinate out of bounds");
    }
    return voxels_.at(calculate_index(x, y, z));
}

const Voxel* VoxelGrid::at(uint32_t x, uint32_t y, uint32_t z) const {
    if (!is_valid_coordinate(x, y, z)) {
        throw std::out_of_range("Voxel coordinate out of bounds");
    }
    return voxels_.at(calculate_index(x, y, z));
}

// Linear indexing
uint32_t VoxelGrid::calculate_index(uint32_t x, uint32_t y, uint32_t z) const {
    // Standard 3D to 1D mapping: index = x + y*nx + z*nx*ny
    return x + y * dimensions_.x + z * dimensions_.x * dimensions_.y;
}

Voxel* VoxelGrid::at(uint32_t linear_index) {
    if (!is_valid_index(linear_index)) {
        throw std::out_of_range("Linear index out of bounds");
    }
    return voxels_.at(linear_index);
}

const Voxel* VoxelGrid::at(uint32_t linear_index) const {
    if (!is_valid_index(linear_index)) {
        throw std::out_of_range("Linear index out of bounds");
    }
    return voxels_.at(linear_index);
}

// Grid management
void VoxelGrid::clear() {
    cleanup_voxels();
    dimensions_ = glm::uvec3(0, 0, 0);
    total_voxels_ = 0;
    voxel_size_ = 0.0;
}

bool VoxelGrid::is_valid_coordinate(uint32_t x, uint32_t y, uint32_t z) const {
    return x < dimensions_.x && y < dimensions_.y && z < dimensions_.z;
}

bool VoxelGrid::is_empty() const {
    return total_voxels_ == 0;
}

// Private helper methods
void VoxelGrid::initialize_voxels() {
    voxels_.reserve(total_voxels_);
    voxels_.clear();
    
    // Initialize voxels in z-y-x order for better cache locality
    for (uint32_t z = 0; z < dimensions_.z; ++z) {
        for (uint32_t y = 0; y < dimensions_.y; ++y) {
            for (uint32_t x = 0; x < dimensions_.x; ++x) {
                voxels_.push_back(new Voxel(x, y, z));
            }
        }
    }
}

void VoxelGrid::cleanup_voxels() {
    for (Voxel* voxel : voxels_) {
        delete voxel;
    }
    voxels_.clear();
}

bool VoxelGrid::is_valid_index(uint32_t linear_index) const {
    return linear_index < total_voxels_;
}
