#pragma once

#include "math/glm_types.hpp"
#include "math/concepts.hpp"
#include <vector>

// Forward declarations
class Simulator;
struct Layer;

/**
 * @brief Modern C++20 utility class for computing volume intersections
 * 
 * Utility class for computing the volume intersection between a voxel cuboid and mesh geometry.
 * This is needed for partial voxel physics where absorption should only occur in the volume
 * inside the geometry, and emission should only occur in the volume outside the geometry.
 * 
 * Now uses C++20 concepts for compile-time type safety.
 */
class VoxelVolumeCalculator {
public:
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
    static double compute_volume_fraction_inside(
        const glm::dvec3& voxel_min,
        const glm::dvec3& voxel_max,
        const Simulator* simulator,
        int samples = 1000
    );
    
    /**
     * @brief Fast approximation using corner testing and adaptive subdivision.
     * @note Less accurate but much faster than Monte Carlo for most cases
     */
    static double compute_volume_fraction_inside_fast(
        const glm::dvec3& voxel_min,
        const glm::dvec3& voxel_max,
        const Simulator* simulator,
        int max_subdivisions = 3
    );

private:
    /**
     * @brief Recursive subdivision method for volume calculation.
     */
    static double subdivide_and_test(
        const glm::dvec3& min_corner,
        const glm::dvec3& max_corner,
        const Simulator* simulator,
        int depth,
        int max_depth
    );
};
