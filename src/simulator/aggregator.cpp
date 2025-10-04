/**
 * @file aggregator.cpp
 * @brief Implementation of results aggregation and export functionality
 *
 * Contains complete implementations of all results processing operations
 * including data aggregation, normalization, reporting, and data accessor
 * methods for Monte Carlo photon transport simulation results.
 */

#include "aggregator.hpp"

#include <algorithm>
#include <iostream>

#include "common/config.hpp"
#include "math/range.hpp"
#include "simulator/simulator.hpp"
#include "simulator/layer.hpp"
#include "simulator/material.hpp"
#include "simulator/medium.hpp"
#include "simulator/metrics.hpp"
#include "simulator/photon.hpp"
#include "simulator/voxel.hpp"

// Initialize aggregator with simulation components
Aggregator::Aggregator(Metrics& metrics, std::vector<Medium>& mediums, std::vector<Photon>& photons, Simulator& sim) :
	metrics_(&metrics), mediums_(&mediums), photons_(&photons), simulator_(&sim) {
}

// Aggregate simulation results from all voxels across all mediums
void Aggregator::aggregate_voxel_data() {
	// Reset medium records (preserve surface interaction data)
	for (auto& medium : *mediums_) {
		medium.reset_record_absorption_and_diffuse();
	}

	// Process voxel data for each medium
	for (auto& medium : *mediums_) {
		auto& volume = medium.get_volume();

		// Traverse all voxels in volume
		for (uint32_t z = 0; z < volume.depth(); ++z) {
			for (uint32_t y = 0; y < volume.height(); ++y) {
				for (uint32_t x = 0; x < volume.width(); ++x) {
					Voxel* voxel = volume.at(x, y, z);
					if (voxel && voxel->material) {
						// Accumulate absorption energy
						medium.get_metrics().add_total_absorption(voxel->absorption);

						// Classify and aggregate emittance by direction
						medium.get_metrics().add_diffuse_reflection(voxel->specular_reflection);
						medium.get_metrics().add_diffuse_transmission(voxel->diffuse_transmission);

						// Note: surface_refraction and specular_reflection handled separately
						// Note: specular_transmission currently unused
					}
				}
			}
		}
	}

	// Invalidate voxel cache since data has changed
	voxels_cache_valid_ = false;
}

// Keep energy values as absolutes - no normalization of raw data
void Aggregator::normalize() {
	// Calculate normalization factor for metrics only
	double normalization_factor =
		static_cast<double>(Config::get().num_photons()) * static_cast<double>(Config::get().num_sources());
	
	// Only normalize medium-level metrics for display purposes
	for (auto& medium : *mediums_) {
		medium.get_metrics().normalize_raw_values(normalization_factor);
	}

	// DO NOT normalize voxel data - keep as absolute values
	// Voxel data normalization happens only for display in voxel renderer

	// Invalidate voxel cache since metrics have changed
	voxels_cache_valid_ = false;
}

// Collect all unique materials across mediums
std::vector<Material> Aggregator::get_all_tissues() const {
	std::vector<Material> all_tissues;
	for (const auto& medium : *mediums_) {
		const auto& medium_tissues = medium.get_tissues();
		for (const auto& material : medium_tissues) {
			// Add unique materials from each medium based on optical properties
			bool found = false;
			for (const auto& existing : all_tissues) {
				if (existing.has_same_optical_properties(material)) {
					found = true;
					break;
				}
			}
			if (!found) {
				all_tissues.push_back(material);
			}
		}
	}
	return all_tissues;
}

// Get layers from first medium (simple multi-medium implementation)
const std::vector<Layer>& Aggregator::get_all_layers() const {
	// Return empty vector if no mediums available
	static const std::vector<Layer> empty_layers;
	if (mediums_->empty()) {
		return empty_layers;
	}
	return (*mediums_)[0].get_layers();
}

// Get cached voxel collection (mutable)
std::vector<std::shared_ptr<Voxel>>& Aggregator::get_all_voxels() {
	update_voxels_cache();
	return combined_voxels_cache_;
}

// Get cached voxel collection (const)
const std::vector<std::shared_ptr<Voxel>>& Aggregator::get_all_voxels() const {
	update_voxels_cache();
	return combined_voxels_cache_;
}

// Rebuild voxel cache if invalidated
void Aggregator::update_voxels_cache() const {
	if (voxels_cache_valid_) {
		return; // Current cache is valid
	}

	combined_voxels_cache_.clear();

	// Combine voxels from all mediums
	for (auto& medium : *mediums_) {
		auto& volume = medium.get_volume();
		for (auto& voxel_ptr : volume) {
			// Convert unique_ptr to shared_ptr with non-owning deleter
			combined_voxels_cache_.push_back(std::shared_ptr<Voxel>(voxel_ptr.get(), [](Voxel*) {}));
		}
	}

	voxels_cache_valid_ = true;
}

// Calculate total voxel count across all mediums
size_t Aggregator::get_total_voxel_count() const {
	size_t total = 0;
	for (const auto& medium : *mediums_) {
		total += medium.get_volume().size();
	}
	return total;
}

// Calculate combined bounding box of all mediums
Range3 Aggregator::get_combined_bounds() const {
	if (mediums_->empty()) {
		return Range3(); // Default bounds for empty collection
	}

	Range3 combined = (*mediums_)[0].get_bounds();
	for (size_t i = 1; i < mediums_->size(); ++i) {
		const Range3& medium_bounds = (*mediums_)[i].get_bounds();
		// Expand bounds to encompass this medium
		combined.min_bounds.x = std::min(combined.min_bounds.x, medium_bounds.min_bounds.x);
		combined.min_bounds.y = std::min(combined.min_bounds.y, medium_bounds.min_bounds.y);
		combined.min_bounds.z = std::min(combined.min_bounds.z, medium_bounds.min_bounds.z);
		combined.max_bounds.x = std::max(combined.max_bounds.x, medium_bounds.max_bounds.x);
		combined.max_bounds.y = std::max(combined.max_bounds.y, medium_bounds.max_bounds.y);
		combined.max_bounds.z = std::max(combined.max_bounds.z, medium_bounds.max_bounds.z);
	}
	return combined;
}
