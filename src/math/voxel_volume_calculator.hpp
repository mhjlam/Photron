#pragma once

#include "math/glm_types.hpp"
#include <vector>

// Forward declarations
class Simulator;
struct Layer;

/**
 * Utility class for computing the volume intersection between a voxel cuboid and mesh geometry.
 * This is needed for partial voxel physics where absorption should only occur in the volume
 * inside the geometry, and emission should only occur in the volume outside the geometry.
 */
class VoxelVolumeCalculator {
public:
    /**
     * Compute the volume fraction of a voxel that lies inside the mesh geometry.
     * Uses Monte Carlo sampling for complex geometries.
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
     * Fast approximation using corner testing and adaptive subdivision.
     * Less accurate but much faster than Monte Carlo for most cases.
     */
    static double compute_volume_fraction_inside_fast(
        const glm::dvec3& voxel_min,
        const glm::dvec3& voxel_max,
        const Simulator* simulator,
        int max_subdivisions = 3
    );

private:
    /**
     * Recursive subdivision method for volume calculation.
     */
    static double subdivide_and_test(
        const glm::dvec3& min_corner,
        const glm::dvec3& max_corner,
        const Simulator* simulator,
        int depth,
        int max_depth
    );
};
