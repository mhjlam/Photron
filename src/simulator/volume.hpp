#pragma once

#include <cstdint>
#include <memory>
#include <vector>

#include <glm/glm.hpp>

#include "math/concepts.hpp"

// Forward declarations
struct Voxel;
class Simulator;

/**
 * @brief Modern C++20 unified voxel grid with volume calculation capabilities
 *
 * This class combines spatial discretization and volume intersection computation
 * for Monte Carlo photon transport simulations. It provides:
 * - Efficient 3D voxel grid storage and access
 * - Volume fraction calculations for partial voxel physics
 * - Monte Carlo and subdivision-based volume sampling
 * - Modern C++20 features with RAII and smart pointer management
 */
class Volume
{
public:
	// Constructors and destructor
	Volume();
	Volume(double voxel_size, uint32_t num_x, uint32_t num_y, uint32_t num_z);
	~Volume();

	// Delete copy constructor and copy assignment (non-copyable due to unique ownership)
	Volume(const Volume&) = delete;
	Volume& operator=(const Volume&) = delete;

	// Enable move constructor and move assignment
	Volume(Volume&& other) noexcept;
	Volume& operator=(Volume&& other) noexcept;

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
	auto begin() { return voxels_.begin(); }
	auto end() { return voxels_.end(); }
	auto begin() const { return voxels_.begin(); }
	auto end() const { return voxels_.end(); }
	auto cbegin() const { return voxels_.cbegin(); }
	auto cend() const { return voxels_.cend(); }

	// Volume calculation methods (integrated from VoxelVolumeCalculator)
	/**
	 * @brief Compute the volume fraction of a voxel that lies inside the mesh geometry.
	 * @note Uses Monte Carlo sampling for complex geometries with C++20 concepts
	 *
	 * @param voxel_min Minimum corner of the voxel
	 * @param voxel_max Maximum corner of the voxel
	 * @param simulator Reference to simulator for point-in-mesh testing
	 * @param samples Number of Monte Carlo samples to use (higher = more accurate but slower)
	 * @return Fraction of voxel volume inside geometry [0.0, 1.0]
	 */
	double compute_volume_fraction_inside(const glm::dvec3& voxel_min, const glm::dvec3& voxel_max,
										  const Simulator& simulator, int samples = 1000) const;

	/**
	 * @brief Fast approximation using corner testing and adaptive subdivision.
	 * @note Less accurate but much faster than Monte Carlo for most cases
	 */
	double compute_volume_fraction_inside_fast(const glm::dvec3& voxel_min, const glm::dvec3& voxel_max,
											   const Simulator& simulator, int max_subdivisions = 3) const;

	/**
	 * @brief Convenience method to compute volume fraction for a specific voxel by coordinates
	 */
	double compute_voxel_volume_fraction(uint32_t x, uint32_t y, uint32_t z, const glm::dvec3& grid_min,
										 const Simulator& simulator, int samples = 1000) const;

private:
	// Private data members
	glm::uvec3 dimensions_ {0, 0, 0};            ///< Grid dimensions (nx, ny, nz)
	uint64_t total_voxels_ {0};                  ///< Total number of voxels (nx * ny * nz)
	double voxel_size_ {0.0};                    ///< Size of each voxel (uniform cubic voxels)
	std::vector<std::unique_ptr<Voxel>> voxels_; ///< Linear storage for voxel smart pointers

	// Private helper methods for grid management
	void initialize_voxels();
	void cleanup_voxels();
	bool is_valid_index(uint32_t linear_index) const;

	// Private helper methods for volume calculation
	/**
	 * @brief Recursive subdivision method for volume calculation.
	 */
	double subdivide_and_test(const glm::dvec3& min_corner, const glm::dvec3& max_corner, const Simulator& simulator,
							  int depth, int max_depth) const;
};
