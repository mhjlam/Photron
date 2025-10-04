/**
 * @file medium.hpp
 * @brief Multi-layered medium representation for Monte Carlo photon transport
 *
 * The Medium class represents a complete multi-layered material system with
 * heterogeneous optical properties. It manages voxelized geometry, layer
 * definitions, material properties, and energy statistics for photon transport
 * simulation in complex biological or synthetic materials.
 */

#pragma once

#include <memory>
#include <vector>

#include "common/config.hpp"
#include "math/range.hpp"
#include "simulator/layer.hpp"
#include "simulator/material.hpp"
#include "simulator/metrics.hpp"
#include "simulator/volume.hpp"

/**
 * @class Medium
 * @brief Multi-layered medium with voxelized geometry for photon transport simulation
 *
 * The Medium class represents a complete heterogeneous material system consisting
 * of multiple layers with varying optical properties. Key responsibilities include:
 *
 * - **Layer Management**: Defines and manages multiple material layers with different properties
 * - **Voxelization**: Converts geometric layer definitions into discrete voxel grids
 * - **Material Properties**: Associates optical properties (absorption, scattering) with spatial regions
 * - **Energy Tracking**: Monitors energy deposition and transport statistics
 * - **Boundary Handling**: Manages interfaces between different material layers
 *
 * The class supports complex geometries including:
 * - Planar layers (typical biological tissues)
 * - Curved interfaces
 * - Embedded objects and inclusions
 * - Varying voxel resolutions
 *
 * Each medium maintains its own coordinate system and can be combined with
 * other media for more complex simulation scenarios.
 */
class Medium
{
public:
	/**
	 * @brief Construct a new Medium with configuration reference
	 *
	 * Initializes the medium structure based on the provided configuration.
	 * The configuration defines layer geometry, material properties, and
	 * voxelization parameters.
	 *
	 * @param config Reference to global configuration object
	 */
	Medium(Config& config);

	/**
	 * @brief Default destructor
	 */
	~Medium() = default;

	// Move semantics (copy disabled due to Volume complexity)
	Medium(Medium&&) = default;                ///< Move constructor
	Medium& operator=(Medium&&) = default;     ///< Move assignment

	Medium(const Medium&) = delete;            ///< Copy constructor disabled
	Medium& operator=(const Medium&) = delete; ///< Copy assignment disabled

	/**
	 * @brief Initialize medium geometry and data structures
	 *
	 * Performs complete initialization of the medium including layer setup,
	 * volume allocation, and material property assignment. Must be called
	 * before using the medium for simulation.
	 *
	 * @return true if initialization succeeded, false otherwise
	 */
	bool initialize();

	/**
	 * @brief Convert layer definitions into discrete voxel representation
	 *
	 * Discretizes the continuous layer geometry into a regular voxel grid
	 * with appropriate material property assignment. This is a critical step
	 * that affects both simulation accuracy and performance.
	 *
	 * @return true if voxelization succeeded, false otherwise
	 */
	bool voxelize_layers();

	/**
	 * @brief Reset all simulation data for new run
	 *
	 * Clears energy deposition data, photon counts, and other simulation
	 * state while preserving geometry and material properties.
	 */
	void reset_simulation_data();

	/**
	 * @brief Normalize medium-level metrics for display purposes
	 *
	 * Normalizes only the medium-level energy metrics by photon count.
	 * Voxel-level data remains as absolute values for correct runtime
	 * accumulation when individual photons are added.
	 */
	void normalize();

	/**
	 * @brief Export simulation results to output files
	 *
	 * Writes energy deposition maps, statistical summaries, and
	 * visualization data to configured output files.
	 */
	void write_results() const;

	// Getters
	const std::vector<Layer>& get_layers() const { return layers_; }
	std::vector<Layer>& get_layers() { return layers_; }
	const std::vector<Material>& get_tissues() const { return materials_; }
	std::vector<Material>& get_tissues() { return materials_; }
	const Volume& get_volume() const { return volume_; }
	Volume& get_volume() { return volume_; }
	const Metrics& get_metrics() const { return metrics_; }
	Metrics& get_metrics() { return metrics_; }

	// Reset methods for aggregation
	void reset_record_absorption_and_diffuse() { metrics_.reset_raw_absorption_and_diffuse(); }
	const Range3& get_bounds() const { return bounds_; }
	Range3& get_bounds() { return bounds_; }

	// Setters
	void set_layers(std::vector<Layer>&& layers) { layers_ = std::move(layers); }

	// Photon tracking
	void increment_photons_entered() { metrics_.increment_photons_entered(); }

	// voxel computation
	Voxel* voxel_at(glm::dvec3& position);
	Cuboid voxel_corners(Voxel* voxel) const;

	// Geometry queries
	double intersection(Source& source) const;
	bool contains_point(const glm::dvec3& point) const;

private:
	std::vector<Layer> layers_;
	std::vector<Material> materials_;

	Config& config_;
	Range3 bounds_;
	Metrics metrics_;
	Volume volume_;

	// Initialize the voxel grid based on layer bounds
	bool initialize_volume();
	bool initialize_layers();

	// Surface detection methods
	void detect_external_surface_voxels();
};
