#pragma once

#include "structs/config.hpp"
#include "structs/cuboid.hpp"
#include "structs/graph.hpp"
#include "structs/layer.hpp"
#include "structs/photon.hpp"
#include "structs/range3.hpp"
#include "structs/record.hpp"
#include "structs/source.hpp"
#include "structs/tissue.hpp"
#include "structs/triangle.hpp"
#include "structs/voxel.hpp"

#include <list>
#include <map>
#include <memory>
#include <string>
#include <vector>

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
	void report();

public:
	Config config;
	Record record;
	Range3 bounds;

	std::vector<Graph> paths;
	std::vector<Layer> layers;
	std::vector<Voxel*> voxels;
	std::vector<Photon> photons;
	std::vector<Tissue> tissues;
	std::vector<Source> sources;
	std::vector<Emitter> emitters;
	std::vector<Triangle> triangles;

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
	bool parse_general(std::list<std::string>& data);
	bool parse_source(std::list<std::string>& data);
	bool parse_tissue(std::list<std::string>& data);
	bool parse_layer(std::list<std::string>& data);

	// data extraction
	bool extract(const std::string& fconfig, std::multimap<std::string, std::list<std::string> >& datalines);
	void extract_block(std::ifstream& inconfig, std::string section,
					   std::multimap<std::string, std::list<std::string> >& datamap);

	// initialization
	bool initialize_grid();
	bool initialize_data();
	bool voxelize_layers();
	
	// data management
	void reset_simulation_data();
};
