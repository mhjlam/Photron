#pragma once

#include "../structs/config.hpp"
#include "../structs/cuboid.hpp"
#include "../structs/graph.hpp"
#include "../structs/layer.hpp"
#include "../structs/photon.hpp"
#include "../structs/range3.hpp"
#include "../structs/record.hpp"
#include "../structs/source.hpp"
#include "../structs/tissue.hpp"
#include "../structs/triangle.hpp"
#include "../structs/voxel.hpp"

#include <list>
#include <map>
#include <memory>
#include <string>
#include <vector>

// Forward declaration for Random class
class Random;

#define PI     (3.1415926535897932384626433832795)
#define HALFPI (1.5707963267948966192313216916398)

class Simulator
{
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
	
	// MCML 3.0.0 Monte Carlo functionality - integrated directly
	std::shared_ptr<Random> mcml_random;
	double mcml_weight_threshold;

private:
	// simulation subroutines
	void Launch(Photon& photon, Source& source);
	void StepSize(Photon& photon);
	void Transfer(Photon& photon);
	void Substep(Photon& photon);
	void Deposit(Photon& photon);
	void Scatter(Photon& photon);
	void Cross(Photon& photon);
	void Radiate(Photon& photon, Vector3& dir, double weight);
	void Roulette(Photon& photon);
	void Normalize();

	// physical computation
	void SpecularReflection(Source& source);
	double InternalReflection(Photon& photon, double& nt, Vector3& tran, Vector3& refl);

	// voxel computation
	Voxel* VoxelAt(Point3& pos);
	Point3 VoxelCenter(Voxel* vox);
	Cuboid VoxelCorners(Voxel* vox);
	bool NearVoxelFace(Point3& pos, Voxel* vox);

	// point translation
	Point3 Move(Point3& pos, Vector3& dir, double d);
	Point3 MoveDelta(Point3& pos, Vector3& dir);

	// MCML 3.0.0 Monte Carlo algorithms - integrated directly
	void GenerateStepSize(Photon& photon);
	void ScatterPhoton(Photon& photon, const Tissue& tissue);
	void RoulettePhoton(Photon& photon);
	void SetMcmlSeed(int seed);

	// MCML integration
	bool UseMcml() const;
	void SetUseMcml(bool use_mcml);

private:
	// file parsing
	bool Parse(const std::string& fconfig);
	bool ParseGeneral(std::list<std::string>& data);
	bool ParseSource(std::list<std::string>& data);
	bool ParseTissue(std::list<std::string>& data);
	bool ParseLayer(std::list<std::string>& data);

	// data extraction
	bool Extract(const std::string& fconfig, std::multimap<std::string, std::list<std::string> >& datalines);
	void ExtractBlock(std::ifstream& inconfig, std::string section,
					  std::multimap<std::string, std::list<std::string> >& datamap);

	// initialization
	bool InitializeGrid();
	bool InitializeData();
	bool VoxelizeLayers();

public:
	// constructor & destructor
	Simulator();
	~Simulator();

	// main routines
	bool Initialize(std::string file);
	void Simulate();
	void Report();
};
