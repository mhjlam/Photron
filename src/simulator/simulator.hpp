#pragma once

#include <list>
#include <map>
#include <memory>
#include <string>
#include <vector>

#include "math/cuboid.hpp"
#include "math/range.hpp"
#include "math/triangle.hpp"
#include "simulator/config.hpp"
#include "simulator/layer.hpp"
#include "simulator/metrics.hpp"
#include "simulator/photon.hpp"
#include "simulator/photon_path.hpp"
#include "simulator/record.hpp"
#include "simulator/source.hpp"
#include "simulator/tissue.hpp"
#include "simulator/voxel.hpp"
#include "simulator/volume.hpp"
#include "simulator/medium.hpp"

// Forward declaration for Random class
class Random;

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
	void report();

public:
	Config config;
	Metrics metrics;

	std::vector<PhotonPath> paths;
	std::vector<Photon> photons;
	std::vector<Source> sources;
	std::vector<Emitter> emitters;
	std::vector<Medium> mediums;

	std::shared_ptr<Random> rng;
	double mcml_weight_threshold;

	// Accessor methods for aggregating data across all mediums
	std::vector<Tissue> get_all_tissues() const;
	const std::vector<Layer>& get_all_layers() const;
	std::vector<std::shared_ptr<Voxel>>& get_all_voxels();
	const std::vector<std::shared_ptr<Voxel>>& get_all_voxels() const;
	size_t get_total_voxel_count() const;
	Range3 get_combined_bounds() const;
	Record get_combined_record() const;
	
	// Compatibility properties for existing code
	[[deprecated("Use get_all_tissues() instead")]] 
	std::vector<Tissue> tissues() const { return get_all_tissues(); }
	
	[[deprecated("Use get_all_layers() instead")]]
	const std::vector<Layer>& layers() const { return get_all_layers(); }
	
	[[deprecated("Use get_combined_bounds() instead")]]
	Range3 bounds() const { return get_combined_bounds(); }
	
	[[deprecated("Use get_combined_record() instead")]]
	Record record() const { return get_combined_record(); }
	
	// Geometry queries (needed by renderer)
	bool is_point_inside_geometry(const glm::dvec3& position) const;
	
	// Voxel access by grid coordinates (for renderer compatibility)
	Voxel* voxel_grid(uint32_t x, uint32_t y, uint32_t z) const;

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

	// MCML 3.0.0 algorithms
	void generate_step_size(Photon& photon);
	void scatter_photon(Photon& photon, const Tissue& tissue);
	void roulette_photon(Photon& photon);
	void set_rng_seed(int seed);

	// file parsing
	bool parse(const std::string& fconfig);

	// initialization
	bool initialize_sources();

	// data management
	void reset_simulation_data();

	// medium context management
	Medium* find_medium_at(const glm::dvec3& position) const;
	bool is_inside_any_medium(const glm::dvec3& position) const;
	void handle_medium_transition(Photon& photon, Medium* from, Medium* to);
	
	// Helper methods that delegate to appropriate medium
	Voxel* voxel_at(const glm::dvec3& position) const;
	Cuboid voxel_corners(Voxel* voxel) const;
};
