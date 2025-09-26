#pragma once

#include <functional>
#include <list>
#include <map>
#include <memory>
#include <string>
#include <vector>

#include "math/cuboid.hpp"
#include "math/range.hpp"
#include "math/triangle.hpp"
#include "math/voxel_dda3d.hpp"
#include "simulator/config.hpp"
#include "simulator/debug_logger.hpp"
#include "simulator/layer.hpp"
#include "simulator/metrics.hpp"
#include "simulator/photon.hpp"
#include "simulator/material.hpp"
#include "simulator/voxel.hpp"
#include "simulator/volume.hpp"
#include "simulator/medium.hpp"

// Forward declaration for Random class
class Random;

// Detailed photon tracking structure for comprehensive reporting
struct DetailedPhotonData {
	uint64_t id {0};
	
	// Entrance information
	glm::dvec3 entrance_position {0.0};
	glm::dvec3 entrance_direction {0.0};
	
	// Exit information (if photon exits)
	bool has_exit {false};
	glm::dvec3 exit_position {0.0};
	glm::dvec3 exit_direction {0.0};
	
	// Termination information (final state when photon stops)
	glm::dvec3 termination_position {0.0};
	glm::dvec3 termination_direction {0.0};
	
	// Simulation statistics
	uint32_t scatter_count {0};
	double total_absorption_deposited {0.0};
	double remaining_weight {0.0};
	
	// Additional tracking
	double initial_weight {0.0};
	bool exited_medium {false};
	std::string termination_reason {"unknown"};
};

class Simulator
{
public:
	// constructor & destructor
	Simulator();
	~Simulator();

	// main routines
	bool initialize(std::string file);
	void simulate();
	void simulate_single_photon(); 		// Add single photon for interactive use
	void aggregate_voxel_data();		// Aggregate voxel energy data into medium records
	void report();

public:
	Metrics metrics;

	std::vector<Photon> photons;  // Unified photon storage (was PhotonPath)
	std::vector<Source> sources;
	std::vector<std::shared_ptr<Emitter>> emitters;
	std::vector<Medium> mediums;
	
	// Detailed photon tracking data for comprehensive reporting
	std::vector<DetailedPhotonData> detailed_photon_data;

	std::shared_ptr<Random> rng;
	double mcml_weight_threshold;

	// 3D DDA instances for robust voxel traversal (one per medium)
	std::vector<std::unique_ptr<VoxelDDA3D>> medium_ddas_;

	// Version tracking for renderer cache optimization
	mutable uint64_t simulation_version_ = 0;
	
	// Get current simulation version (increments when data changes)
	uint64_t get_simulation_version() const { return simulation_version_; }
	
	// Increment version when simulation data changes (called internally)
	void increment_simulation_version() const { ++simulation_version_; }

	// Accessor methods for aggregating data across all mediums
	std::vector<Material> get_all_tissues() const;
	const std::vector<Layer>& get_all_layers() const;
	std::vector<std::shared_ptr<Voxel>>& get_all_voxels();
	const std::vector<std::shared_ptr<Voxel>>& get_all_voxels() const;
	size_t get_total_voxel_count() const;
	Range3 get_combined_bounds() const;
	
	// Combined metrics access methods (replacement for get_combined_record)
	double get_combined_total_absorption() const;
	double get_combined_diffuse_reflection() const;
	double get_combined_specular_reflection() const; 
	double get_combined_surface_refraction() const;
	double get_combined_diffuse_transmission() const;
	double get_combined_specular_transmission() const;
	double get_combined_avg_path_length() const;
	int get_combined_total_steps() const;
	int get_combined_photons_entered() const;
	
	// Energy conservation calculation (shared between console and overlay)
	struct EnergyConservation {
		double total_absorption = 0.0;
		double total_reflection = 0.0; 
		double total_transmission = 0.0;
		double surface_reflection = 0.0;
		double surface_refraction = 0.0;
		double total_energy = 0.0;
		double total_diffusion = 0.0;
	};
	EnergyConservation calculate_energy_conservation() const;
	
	// Compatibility properties for existing code
	[[deprecated("Use get_all_tissues() instead")]] 
	std::vector<Material> tissues() const { return get_all_tissues(); }
	
	[[deprecated("Use get_all_layers() instead")]]
	const std::vector<Layer>& layers() const { return get_all_layers(); }
	
	[[deprecated("Use get_combined_bounds() instead")]]
	Range3 bounds() const { return get_combined_bounds(); }
	
	// Geometry queries (needed by renderer)
	bool is_point_inside_geometry(const glm::dvec3& position) const;
	
	// Voxel access by grid coordinates (for renderer compatibility)
	Voxel* voxel_grid(uint32_t x, uint32_t y, uint32_t z) const;

	// Path access for backward compatibility (creates PhotonPaths from photon data)
	std::vector<Photon> get_paths() const;
	
	// Progress callback for UI updates
	void set_progress_callback(std::function<void(uint64_t, uint64_t)> callback) {
		progress_callback_ = callback;
	}

private:
	// simulation subroutines
	void launch(Photon& photon, const Source& source);
	void step_size(Photon& photon);
	void transfer(Photon& photon);
	void sub_step(Photon& photon);
	void deposit(Photon& photon);
	void track_voxel_path_and_deposit(Photon& photon);
	void scatter(Photon& photon);
	void cross(Photon& photon);
	void radiate(Photon& photon, glm::dvec3& direction, double weight);
	Voxel* find_last_surface_voxel_with_dda(const Photon& photon, const glm::dvec3& exit_direction);
	void roulette(Photon& photon);
	void normalize();
	
	// simplified energy recording system
	void terminate_photon_and_record_energy(Photon& photon, const std::string& reason);

	// physical computation
	void specular_reflection(Photon& photon);
	double internal_reflection(Photon& photon, double& eta_t, glm::dvec3& tran, glm::dvec3& refl);
	bool is_photon_reflecting(const Photon& photon) const;

	// point translation
	glm::dvec3 move(glm::dvec3& position, glm::dvec3& direction, double d);
	glm::dvec3 move_delta(glm::dvec3& position, glm::dvec3& direction);

	// MCML algorithms
	void generate_step_size(Photon& photon);
	void scatter_photon(Photon& photon, const Material& tissue);
	void roulette_photon(Photon& photon);

	// file parsing
	bool parse(const std::string& fconfig);

	// initialization
	bool initialize_sources();

	// data management
	void reset_simulation_data();
	void output_voxel_emittance_summary();

	// medium context management
	Medium* find_medium_at(const glm::dvec3& position) const;
	bool is_inside_any_medium(const glm::dvec3& position) const;
	void handle_medium_transition(Photon& photon, Medium* from, Medium* to);
	void validate_photon_state_after_interface_transition(Photon& photon, Medium* from_medium, Medium* to_medium);
	
	// 3D DDA voxel traversal methods
	void initialize_dda_instances();
	void track_voxel_path_with_dda(Photon& photon);
	void track_photon_path_segments_for_absorption(Photon& photon);
	Medium* find_medium_at_with_dda(const glm::dvec3& position) const;
	
	// Helper methods that delegate to appropriate medium
	Voxel* voxel_at(const glm::dvec3& position) const;
	Cuboid voxel_corners(Voxel* voxel) const;
	
	// Progress callback member
	std::function<void(uint64_t, uint64_t)> progress_callback_;
};
