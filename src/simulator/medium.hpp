#pragma once

#include <vector>
#include <memory>

#include "simulator/layer.hpp"
#include "simulator/volume.hpp"
#include "simulator/metrics.hpp"
#include "simulator/config.hpp"
#include "simulator/material.hpp"
#include "math/range.hpp"


class Medium
{
public:
    // Constructor
    Medium(Config& config);
    
    // Destructor
    ~Medium() = default;

    // Move semantics
    Medium(Medium&&) = default;
    Medium& operator=(Medium&&) = default;
    
    // Delete copy semantics due to Volume
    Medium(const Medium&) = delete;
    Medium& operator=(const Medium&) = delete;

    // Initialize the medium with given configuration
    bool initialize();

    // Voxelize the layers into the volume
    bool voxelize_layers();
    
    // Phase 2: Disambiguate surface voxels (remove false positives)
    void disambiguate_surface_voxels();
    
    // Identify surface voxels
    void identify_surface_voxels();
    
    // Reset simulation data
    void reset_simulation_data();
    
    // Normalize recorded values
    void normalize();
    
    // Write results to files
    void write_results() const;

    // Getters
    const std::vector<Layer>& get_layers() const { return layers_; }
    std::vector<Layer>& get_layers() { return layers_; }
    const std::vector<Material>& get_tissues() const { return tissues_; }
    std::vector<Material>& get_tissues() { return tissues_; }
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
	std::vector<Material> tissues_;

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
