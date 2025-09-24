#pragma once

#include <cstdint>
#include <memory>
#include <vector>

#include <glm/glm.hpp>

#include "math/concepts.hpp"
#include "simulator/layer.hpp"  // Full include instead of forward declaration

// Forward declarations
struct Voxel;
class Simulator;

/**
 * @brief Classification result from Distance Field Based voxelization
 */
struct VoxelClassification {
	bool is_inside_geometry = false;    // Whether voxel intersects any geometry
	bool is_boundary_voxel = false;     // Whether voxel intersects geometry boundary (partial volume)
	bool is_surface_voxel = false;      // Whether voxel is at external surface (to ambient)
	double volume_fraction = 0.0;       // Fraction of voxel inside geometry [0,1]
	uint32_t dominant_tissue_id = 0;    // Tissue ID of dominant layer
};

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
	constexpr uint32_t width() const { return dimensions_.x; }
	constexpr uint32_t height() const { return dimensions_.y; }
	constexpr uint32_t depth() const { return dimensions_.z; }
	constexpr uint64_t size() const { return total_voxels_; }

	constexpr uint64_t total_voxels() const { return total_voxels_; }
	constexpr double voxel_size() const { return voxel_size_; }

	constexpr const glm::uvec3& dimensions() const { return dimensions_; }

	// Grid management
	void clear();
	inline bool is_valid_coordinate(uint32_t x, uint32_t y, uint32_t z) const {
		return x < dimensions_.x && y < dimensions_.y && z < dimensions_.z;
	}
	inline bool is_empty() const { return total_voxels_ == 0; }

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
	 * @details Uses Monte Carlo sampling to estimate the fraction of the voxel volume
	 *          that intersects with the provided mesh layers. Higher sample counts
	 *          provide more accurate results at the cost of computation time.
	 * @param voxel_min Lower corner of the voxel in world coordinates
	 * @param voxel_max Upper corner of the voxel in world coordinates
	 * @param layers Vector of mesh layers to test intersection against
	 * @param samples Number of random samples to use (default: 1000)
	 * @return Volume fraction [0.0, 1.0] where 0 = no intersection, 1 = fully inside
	 */
	double fraction_inside(const glm::dvec3& voxel_min, const glm::dvec3& voxel_max,
		const std::vector<Layer>& layers, int samples = 1000) const;

	/**
	 * @brief Fast approximation using corner testing and adaptive subdivision.
	 * @note Less accurate but much faster than Monte Carlo for most cases
	 */
	double fraction_inside_fast(const glm::dvec3& voxel_min, const glm::dvec3& voxel_max,
		const std::vector<Layer>& layers, int max_subdivisions = 3) const;

	/**
	 * @brief Convenience method to compute volume fraction for a specific voxel by coordinates
	 */
	double fraction(uint32_t x, uint32_t y, uint32_t z, const glm::dvec3& grid_min,
		const std::vector<Layer>& layers, int samples = 1000) const;

	/**
	 * @brief Distance field based voxelization algorithm.
	 * Most robust approach for complex mesh voxelization using signed distance fields.
	 */
	VoxelClassification distance_field_voxelization(const glm::dvec3& voxel_min, const glm::dvec3& voxel_max,
												   const std::vector<Layer>& layers) const;
	
	/**
	 * @brief Compute accurate volume fraction for surface voxels using SDF sampling.
	 */
	double compute_volume_fraction_sdf(const glm::dvec3& voxel_min, const glm::dvec3& voxel_max,
									   const std::vector<Layer>& layers, int subdivisions) const;
	
	/**
	 * @brief Compute signed distance from a point to the closest point on a triangle mesh.
	 */
	double compute_signed_distance_to_mesh(const glm::dvec3& point, const std::vector<Triangle>& mesh) const;
	
	/**
	 * @brief Find the closest point on a triangle to a given point.
	 */
	glm::dvec3 closest_point_on_triangle(const glm::dvec3& point, const Triangle& triangle) const;

private:
	// Private data members
	glm::uvec3 dimensions_ {0, 0, 0};            ///< Grid dimensions (nx, ny, nz)
	uint64_t total_voxels_ {0};                  ///< Total number of voxels (nx * ny * nz)
	double voxel_size_ {0.0};                    ///< Size of each voxel (uniform cubic voxels)
	std::vector<std::unique_ptr<Voxel>> voxels_; ///< Linear storage for voxel smart pointers

	// Private helper methods for grid management
	void initialize_voxels();
	void cleanup_voxels();
	inline bool is_valid_index(uint32_t linear_index) const {
		return linear_index < total_voxels_;
	}

	// Private helper methods for volume calculation
	/**
	 * @brief Recursive subdivision method for volume calculation.
	 */
	double subdivide_test(const glm::dvec3& min_corner, const glm::dvec3& max_corner, 
		const std::vector<Layer>& layers, int depth, int max_depth) const;
};