/**
 * @file voxel.hpp
 * @brief Discrete voxel representation for Monte Carlo photon transport
 *
 * Defines the Voxel class which represents individual volume elements in the
 * discretized simulation space. Each voxel stores material properties, energy
 * deposition data, and boundary information essential for accurate photon transport.
 */

#pragma once

#include <glm/glm.hpp>

#include "simulator/material.hpp"

/**
 * @class Voxel
 * @brief Discrete volume element with energy tracking and material properties
 *
 * The Voxel class represents an individual cubic volume element within the
 * discretized simulation space. Each voxel maintains:
 *
 * **Spatial Information:**
 * - **Grid Coordinates**: Discrete (i,j,k) position in voxel grid
 * - **Boundary Classification**: Interior, boundary, or surface voxel designation
 * - **Volume Fractions**: Partial volume support for accurate boundary physics
 *
 * **Material Properties:**
 * - **Material Assignment**: Pointer to associated material properties
 * - **Layer Identification**: ID of the geometric layer containing this voxel
 *
 * **Energy Tracking:**
 * - **Absorption**: Cumulative photon energy absorbed within voxel
 * - **Emittance**: Energy exiting voxel (for validation and visualization)
 * - **Reflection/Transmission**: Directional energy flow for conservation analysis
 *
 * **Boundary Physics:**
 * - **Volume Fractions**: Support for voxels partially inside/outside geometry
 * - **Surface Detection**: Identification of external surface voxels
 * - **Interface Handling**: Special treatment for material boundary voxels
 *
 * The voxel system enables accurate energy deposition tracking while maintaining
 * computational efficiency through discrete spatial representation.
 */
class Voxel
{
private:
	// Volume fraction thresholds for voxel classification
	static constexpr double FULLY_INSIDE_THRESHOLD = 0.999;  ///< Threshold for fully interior voxels
	static constexpr double FULLY_OUTSIDE_THRESHOLD = 0.999; ///< Threshold for fully exterior voxels
	static constexpr double PARTIAL_THRESHOLD = 0.001;       ///< Minimum fraction for partial occupancy
	static constexpr double EMITTANCE_THRESHOLD = 0.01;      ///< Minimum 1% volume overlap for emittance recording

public:
	// Spatial coordinates and identification
	glm::uvec3 coords;    ///< Voxel grid coordinates (ix, iy, iz)
	uint8_t layer_id {0}; ///< ID of geometric layer containing this voxel

	// Material properties
	Material* material; ///< Pointer to associated material properties (null for empty voxels)

	// Energy deposition and transport
	double absorption;           ///< Cumulative photon energy absorbed in this voxel
	double emittance;            ///< Energy exiting this voxel (conservation validation)
	double specular_reflection;  ///< Energy exiting via specular reflection
	double diffuse_reflection;   ///< Energy exiting via diffuse reflection
	double diffuse_transmission; ///< Energy exiting via transmission

	// Voxel classification flags
	bool is_boundary_voxel; ///< True if voxel intersects layer boundaries
	bool is_surface_voxel;  ///< True if voxel is at external surface of geometry

	// Partial volume support for accurate boundary physics
	double volume_fraction_inside;  ///< Fraction [0,1] of voxel volume inside geometry
	double volume_fraction_outside; ///< Fraction [0,1] of voxel volume outside geometry

	/**
	 * @brief Default constructor creates empty voxel at origin
	 *
	 * Initializes voxel with zero energy, null material, and default classification.
	 */
	Voxel() :
		coords(0), absorption(0), emittance(0), material(nullptr), diffuse_transmission(0), specular_reflection(0),
		diffuse_reflection(0), volume_fraction_inside(0.0), volume_fraction_outside(0.0), is_boundary_voxel(false),
		is_surface_voxel(false) {}

	/**
	 * @brief Construct voxel at specified grid coordinates
	 *
	 * @param x Grid X coordinate (ix)
	 * @param y Grid Y coordinate (iy)
	 * @param z Grid Z coordinate (iz)
	 */
	Voxel(uint32_t x, uint32_t y, uint32_t z) :
		coords(x, y, z), absorption(0), emittance(0), material(nullptr), diffuse_transmission(0),
		specular_reflection(0), diffuse_reflection(0), volume_fraction_inside(0.0), volume_fraction_outside(0.0),
		is_boundary_voxel(false), is_surface_voxel(false) {}

	// Coordinate access methods for backward compatibility

	/**
	 * @brief Get X grid coordinate
	 * @return const uint32_t& X coordinate (ix)
	 */
	constexpr const uint32_t& ix() const { return coords.x; }

	/**
	 * @brief Get Y grid coordinate
	 * @return const uint32_t& Y coordinate (iy)
	 */
	constexpr const uint32_t& iy() const { return coords.y; }

	/**
	 * @brief Get Z grid coordinate
	 * @return const uint32_t& Z coordinate (iz)
	 */
	constexpr const uint32_t& iz() const { return coords.z; }

	// GLM convenience methods
	constexpr const glm::uvec3& indices() const { return coords; }
	void set_indices(const glm::uvec3& idx) { coords = idx; }

	// Volume fraction utilities for geometric classification
	/**
	 * @brief Test if voxel is completely inside mesh geometry
	 * @return True if volume fraction inside exceeds threshold
	 */
	constexpr bool is_fully_inside() const { return volume_fraction_inside >= FULLY_INSIDE_THRESHOLD; }

	/**
	 * @brief Test if voxel is completely outside mesh geometry
	 * @return True if volume fraction outside exceeds threshold
	 */
	constexpr bool is_fully_outside() const { return volume_fraction_outside >= FULLY_OUTSIDE_THRESHOLD; }

	/**
	 * @brief Test if voxel is partially intersected by mesh boundaries
	 * @return True if voxel has significant volume both inside and outside mesh
	 */
	bool is_partial() const {
		return is_boundary_voxel && volume_fraction_inside > PARTIAL_THRESHOLD
			   && volume_fraction_outside > PARTIAL_THRESHOLD;
	}

	// Energy tracking utilities for Monte Carlo analysis
	/**
	 * @brief Calculate total emitted energy from all sources
	 * @return Sum of all emittance components for this voxel
	 */
	double total_emittance() const { return diffuse_transmission + specular_reflection + diffuse_reflection; }

	/**
	 * @brief Reset all energy accumulation counters to zero
	 *
	 * Used between simulation runs to clear previous energy data
	 * while preserving geometric properties and material assignments.
	 */
	void reset_energy() {
		absorption = 0;
		emittance = 0;
		diffuse_transmission = 0;
		specular_reflection = 0;
		diffuse_reflection = 0;
	}

	// Special comparison overloads
	constexpr bool operator==(const Voxel& other) const {
		return (other.coords == coords && other.material == material);
	}
};
