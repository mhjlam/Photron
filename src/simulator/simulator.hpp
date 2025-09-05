#pragma once

#include <list>
#include <map>
#include <memory>
#include <string>
#include <vector>

#include "math/cuboid.hpp"
#include "math/range.hpp"
#include "math/triangle.hpp"
#include "math/mesh_geometry.hpp"
#include "math/voxel_volume.hpp"
#include "simulator/config.hpp"
#include "simulator/graph.hpp"
#include "simulator/layer.hpp"
#include "simulator/metrics.hpp"
#include "simulator/photon.hpp"
#include "simulator/record.hpp"
#include "simulator/source.hpp"
#include "simulator/tissue.hpp"
#include "simulator/voxel.hpp"

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
	void simulate_single_photon(); // Add single photon for interactive use
	void report();

	// utility functions
	bool is_point_inside_geometry(const glm::dvec3& point) const;
	bool is_point_inside_layer_mesh(const glm::dvec3& point, const Layer& layer) const;

public:
	Config config;
	Record record;
	Range3 bounds;
	Metrics metrics;

	std::vector<Graph> paths;
	std::vector<Layer> layers;
	VoxelVolume voxel_grid;
	std::vector<Photon> photons;
	std::vector<Tissue> tissues;
	std::vector<Source> sources;
	std::vector<Emitter> emitters;

	std::shared_ptr<Random> mcml_random;
	double mcml_weight_threshold;

private:
	// simulation subroutines
	void launch(Photon& photon, Source& source);
	void step_size(Photon& photon);
	void transfer(Photon& photon);
	void sub_step(Photon& photon);
	void deposit(Photon& photon);
	void scatter(Photon& photon);
	void cross(Photon& photon);
	void radiate(Photon& photon, glm::dvec3& direction, double weight);
	void roulette(Photon& photon);
	void normalize();

	// physical computation
	void specular_reflection(Source& source);
	double internal_reflection(Photon& photon, double& eta_t, glm::dvec3& tran, glm::dvec3& refl);

	// voxel computation
	Voxel* voxel_at(glm::dvec3& position);
	glm::dvec3 voxel_center(Voxel* voxel);
	Cuboid voxel_corners(Voxel* voxel);

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
	bool initialize_grid();
	bool initialize_data();
	bool voxelize_layers();

	// data management
	void reset_simulation_data();
};
