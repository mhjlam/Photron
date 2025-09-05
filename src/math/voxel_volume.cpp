#include "voxel_volume.hpp"
#include "../simulator/voxel.hpp"
#include "random.hpp"
#include "simulator/simulator.hpp"
#include <stdexcept>
#include <algorithm>
#include <memory>
#include <random>

// Default constructor
VoxelVolume::VoxelVolume() = default;

// Parameterized constructor
VoxelVolume::VoxelVolume(double voxel_size, uint32_t num_x, uint32_t num_y, uint32_t num_z)
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
VoxelVolume::~VoxelVolume() {
    cleanup_voxels();
}

// Move constructor
VoxelVolume::VoxelVolume(VoxelVolume&& other) noexcept
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
VoxelVolume& VoxelVolume::operator=(VoxelVolume&& other) noexcept {
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
Voxel* VoxelVolume::operator()(uint32_t x, uint32_t y, uint32_t z) {
    if (!is_valid_coordinate(x, y, z)) {
        throw std::out_of_range("Voxel coordinate out of bounds");
    }
    return voxels_[calculate_index(x, y, z)].get();
}

const Voxel* VoxelVolume::operator()(uint32_t x, uint32_t y, uint32_t z) const {
    if (!is_valid_coordinate(x, y, z)) {
        throw std::out_of_range("Voxel coordinate out of bounds");
    }
    return voxels_[calculate_index(x, y, z)].get();
}

// Safe access methods with bounds checking
Voxel* VoxelVolume::at(uint32_t x, uint32_t y, uint32_t z) {
    if (!is_valid_coordinate(x, y, z)) {
        throw std::out_of_range("Voxel coordinate out of bounds");
    }
    return voxels_.at(calculate_index(x, y, z)).get();
}

const Voxel* VoxelVolume::at(uint32_t x, uint32_t y, uint32_t z) const {
    if (!is_valid_coordinate(x, y, z)) {
        throw std::out_of_range("Voxel coordinate out of bounds");
    }
    return voxels_.at(calculate_index(x, y, z)).get();
}

// Linear indexing
uint32_t VoxelVolume::calculate_index(uint32_t x, uint32_t y, uint32_t z) const {
    // Standard 3D to 1D mapping: index = x + y*nx + z*nx*ny
    return x + y * dimensions_.x + z * dimensions_.x * dimensions_.y;
}

Voxel* VoxelVolume::at(uint32_t linear_index) {
    if (!is_valid_index(linear_index)) {
        throw std::out_of_range("Linear index out of bounds");
    }
    return voxels_.at(linear_index).get();
}

const Voxel* VoxelVolume::at(uint32_t linear_index) const {
    if (!is_valid_index(linear_index)) {
        throw std::out_of_range("Linear index out of bounds");
    }
    return voxels_.at(linear_index).get();
}

// Grid management
void VoxelVolume::clear() {
    cleanup_voxels();
    dimensions_ = glm::uvec3(0, 0, 0);
    total_voxels_ = 0;
    voxel_size_ = 0.0;
}

bool VoxelVolume::is_valid_coordinate(uint32_t x, uint32_t y, uint32_t z) const {
    return x < dimensions_.x && y < dimensions_.y && z < dimensions_.z;
}

bool VoxelVolume::is_empty() const {
    return total_voxels_ == 0;
}

// Private helper methods for grid management
void VoxelVolume::initialize_voxels() {
    // Modern C++20: Pre-allocate and use emplace_back for better performance
    voxels_.clear();
    voxels_.reserve(total_voxels_);
    
    // Initialize voxels in z-y-x order for better cache locality
    // Use emplace_back to construct in-place and avoid unnecessary copies
    for (uint32_t z = 0; z < dimensions_.z; ++z) {
        for (uint32_t y = 0; y < dimensions_.y; ++y) {
            for (uint32_t x = 0; x < dimensions_.x; ++x) {
                voxels_.emplace_back(std::make_unique<Voxel>(x, y, z));
            }
        }
    }
    
    // Shrink to fit for optimal memory usage
    voxels_.shrink_to_fit();
}

void VoxelVolume::cleanup_voxels() {
    // Smart pointers automatically clean up when the vector is cleared
    voxels_.clear();
}

bool VoxelVolume::is_valid_index(uint32_t linear_index) const {
    return linear_index < total_voxels_;
}

// Volume calculation methods (integrated from VoxelVolumeCalculator)
double VoxelVolume::compute_volume_fraction_inside(
    const glm::dvec3& voxel_min,
    const glm::dvec3& voxel_max,
    const Simulator* simulator,
    int samples
) const {
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

double VoxelVolume::compute_volume_fraction_inside_fast(
    const glm::dvec3& voxel_min,
    const glm::dvec3& voxel_max,
    const Simulator* simulator,
    int max_subdivisions
) const {
    if (!simulator) return 0.0;
    return subdivide_and_test(voxel_min, voxel_max, simulator, 0, max_subdivisions);
}

double VoxelVolume::compute_voxel_volume_fraction(
    uint32_t x, uint32_t y, uint32_t z,
    const glm::dvec3& grid_min,
    const Simulator* simulator,
    int samples
) const {
    if (!is_valid_coordinate(x, y, z) || !simulator) return 0.0;
    
    // Calculate voxel bounds
    glm::dvec3 voxel_min = grid_min + glm::dvec3(
        x * voxel_size_,
        y * voxel_size_,
        z * voxel_size_
    );
    
    glm::dvec3 voxel_max = voxel_min + glm::dvec3(voxel_size_);
    
    return compute_volume_fraction_inside(voxel_min, voxel_max, simulator, samples);
}

double VoxelVolume::subdivide_and_test(
    const glm::dvec3& min_corner,
    const glm::dvec3& max_corner,
    const Simulator* simulator,
    int depth,
    int max_depth
) const {
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
