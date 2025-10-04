/**
 * @file aggregator.hpp
 * @brief Results aggregation and export for Monte Carlo simulation
 *
 * The Aggregator class handles all post-simulation processing including
 * data aggregation across mediums, result normalization, reporting, and
 * providing accessor methods for external visualization components.
 */

#pragma once

#include <memory>
#include <string>
#include <vector>

#include "math/range.hpp"
#include "simulator/metrics.hpp"

// Forward declarations
class Medium;
class Photon;
class Material;
class Layer;
class Voxel;
class Simulator;

/**
 * @class Aggregator
 * @brief Handles simulation results processing and data aggregation
 *
 * The Aggregator manages all operations related to processing simulation
 * results after photon transport is complete. It aggregates data from voxels
 * across all mediums, normalizes results, generates reports, and provides
 * unified data access methods for external visualization components.
 *
 * Key responsibilities:
 * - Voxel data aggregation across all mediums
 * - Result normalization and energy conservation validation
 * - Report generation and export to various formats
 * - Unified data accessor methods for renderers and UI components
 * - Metrics processing and statistical analysis
 *
 * The class operates on references to simulation data without owning it,
 * ensuring proper lifetime management and efficient data access.
 */
class Aggregator
{
public:
	/**
	 * @brief Construct Aggregator with references to simulation data
	 *
	 * @param metrics Reference to metrics object (not owned)
	 * @param mediums Reference to vector of mediums (not owned)
	 * @param photons Reference to vector of photons (not owned)
	 * @param simulator Reference to simulator instance for metrics delegation (not owned)
	 */
	explicit Aggregator(Metrics& metrics,
						std::vector<Medium>& mediums,
						std::vector<Photon>& photons,
						Simulator& simulator);

	/**
	 * @brief Aggregate voxel-level energy data into medium-level records
	 *
	 * Consolidates detailed voxel energy deposition data into summary
	 * statistics for each medium layer for analysis and visualization.
	 * Prevents double-counting by using only voxel-level data.
	 */
	void aggregate_voxel_data();

	/**
	 * @brief Normalize medium-level metrics for display purposes
	 *
	 * Normalizes only the medium-level energy metrics by total photon count
	 * for display in UI. Voxel-level data remains as absolute values to allow
	 * proper accumulation when individual photons are added during runtime.
	 * The voxel renderer handles display normalization separately.
	 */
	void normalize();

	/**
	 * @brief Get aggregated material definitions from all mediums
	 *
	 * Combines material properties from all mediums into a single
	 * vector for unified access by visualization components.
	 *
	 * @return std::vector<Material> Combined tissue/material properties
	 */
	std::vector<Material> get_all_tissues() const;

	/**
	 * @brief Get aggregated layer definitions from all mediums
	 *
	 * Provides unified access to layer structure across all mediums
	 * for visualization and analysis purposes.
	 *
	 * @return const std::vector<Layer>& Combined layer structure (read-only)
	 */
	const std::vector<Layer>& get_all_layers() const;

	/**
	 * @brief Get aggregated voxel data from all mediums (mutable)
	 *
	 * Combines voxel collections from all mediums for unified access.
	 * Uses shared_ptr with custom deleter to avoid ownership issues.
	 *
	 * @return std::vector<std::shared_ptr<Voxel>>& Combined voxel collection
	 */
	std::vector<std::shared_ptr<Voxel>>& get_all_voxels();

	/**
	 * @brief Get aggregated voxel data from all mediums (read-only)
	 *
	 * Combines voxel collections from all mediums for unified access.
	 * Uses shared_ptr with custom deleter to avoid ownership issues.
	 *
	 * @return const std::vector<std::shared_ptr<Voxel>>& Combined voxel collection
	 */
	const std::vector<std::shared_ptr<Voxel>>& get_all_voxels() const;

	/**
	 * @brief Get total number of voxels across all mediums
	 *
	 * Calculates the combined voxel count from all mediums
	 * for memory usage estimation and progress tracking.
	 *
	 * @return size_t Combined voxel count
	 */
	size_t get_total_voxel_count() const;

	/**
	 * @brief Get combined spatial bounds of all mediums
	 *
	 * Calculates the bounding box that encompasses the entire
	 * simulation geometry across all mediums.
	 *
	 * @return Range3 Bounding box encompassing entire simulation geometry
	 */
	Range3 get_combined_bounds() const;

private:
	Metrics* metrics_;             ///< Reference to metrics object (not owned)
	std::vector<Medium>* mediums_; ///< Reference to mediums vector (not owned)
	std::vector<Photon>* photons_; ///< Reference to photons vector (not owned)
	Simulator* simulator_;         ///< Reference to simulator for metrics delegation (not owned)

	// Cache for combined voxel collections to avoid repeated allocations
	mutable std::vector<std::shared_ptr<Voxel>> combined_voxels_cache_;
	mutable bool voxels_cache_valid_ = false;

	/**
	 * @brief Update the combined voxels cache if invalid
	 *
	 * Rebuilds the unified voxel collection cache when simulation
	 * data changes, ensuring efficient access for repeated queries.
	 */
	void update_voxels_cache() const;
};
