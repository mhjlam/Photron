#include "simulator.hpp"

#include "random.hpp"
#include "structs/index3.hpp"
#include "structs/range1.hpp"
#include "structs/ray.hpp"
#include "utilities/experimenter.hpp"
#include "utilities/utilities.hpp"

#include <algorithm>
#include <cfloat>
#include <cstdlib>
#include <ctime>
#include <fstream>
#include <iostream>
#include <memory>
#include <numbers>
#include <sstream>

constexpr double HALF_PI = std::numbers::pi / 2.0;

/***********************************************************
 * Simulator constructor.
 ***********************************************************/
Simulator::Simulator() : mcml_random(std::make_shared<Random>()), mcml_weight_threshold(1e-4) {
	paths = std::vector<Graph>();
	layers = std::vector<Layer>();
	voxels = std::vector<Voxel*>();
	photons = std::vector<Photon>();
	tissues = std::vector<Tissue>();
	sources = std::vector<Source>();
	triangles = std::vector<Triangle>();

	// Initialize MCML random number generator
	mcml_random->seed(static_cast<int>(std::time(nullptr)));
}

/***********************************************************
 * Simulator destructor.
 ***********************************************************/
Simulator::~Simulator() {
	// delete voxels
	for (uint64_t i = 0; i < voxels.size(); ++i) {
		delete voxels[i];
		voxels[i] = nullptr;
	}

	// delete path vertices
	for (uint32_t i = 0; i < paths.size(); ++i) {
		Graph path = paths[i];
		Vertex* item = path.head;

		if (path.head->prev) {
			delete path.head->prev;
			path.head->prev = nullptr;
		}

		while (item) {
			Vertex* old = item;
			item = item->next;

			delete old;
			old = nullptr;
		}
	}
}

/***********************************************************
 * INITIALIZATION
 ***********************************************************/

/***********************************************************
 * Parse the input file and initializes the data structures.
 ***********************************************************/
bool Simulator::initialize(std::string file) {
	std::cout << "Initializing Photron" << std::endl;

	// read and parse input configuration file
	std::cout << "Parsing configuration file: " << file << std::endl;
	if (!parse(file)) {
		std::cerr << "An error occurred while parsing the input file." << std::endl;
		return false;
	}
	std::cout << "Configuration parsed successfully." << std::endl;

	// initialize voxel grid
	if (!initialize_grid()) {
		std::cerr << "An error occurred while initializing the voxel grid." << std::endl;
		return false;
	}
	std::cout << "Voxel grid initialized successfully." << std::endl;

	// initialize other data structures
	if (!initialize_data()) {
		std::cerr << "An error occurred while initializing the data structures." << std::endl;
		return false;
	}
	std::cout << "Data structures initialized successfully." << std::endl;

	// voxelize the geometry
	if (!voxelize_layers()) {
		std::cerr << "An error occurred during geometry voxelization." << std::endl;
		return false;
	}
	std::cout << "Geometry voxelization completed successfully." << std::endl;
	return true;
}

/***********************************************************
 *	INITIALIZATION SUBROUTINES
 ***********************************************************/

/***********************************************************
 * Read and parse the given configuration file.
 ***********************************************************/
bool Simulator::parse(const std::string& fconfig) {
	std::multimap<std::string, std::list<std::string>> datamap;

	// attempt to extract all the data from the input file
	if (!Simulator::extract(fconfig, datamap)) {
		std::cerr << "Failed to extract data from config file." << std::endl;
		return false;
	}

	// parse the individual configuration blocks
	for (const auto& it : datamap) {
		std::string name = it.first;
		std::list<std::string> data = it.second;

		if (equals(name, "general") && !Simulator::parse_general(data)) {
			std::cerr << "Failed to parse general section." << std::endl;
			return false;
		}
		else if (equals(name, "source") && !Simulator::parse_source(data)) {
			std::cerr << "Failed to parse source section." << std::endl;
			return false;
		}
		else if (equals(name, "tissue") && !Simulator::parse_tissue(data)) {
			std::cerr << "Failed to parse tissue section." << std::endl;
			return false;
		}
		else if (equals(name, "layer") && !Simulator::parse_layer(data)) {
			std::cerr << "Failed to parse layer section." << std::endl;
			return false;
		}
	}

	return true;
}

/***********************************************************
 * Set the general scene configuration.
 ***********************************************************/
bool Simulator::parse_general(std::list<std::string>& data) {
	// set general settings
	for (const auto& it : parameter_values(data)) {
		std::string param = it.first;
		std::string value = it.second;

		if (equals(param, "numphotons")) {
			config.num_photons = str2num<uint64_t>(value);
		}
		else if (equals(param, "ambienteta")) {
			config.ambient_eta = str2num<double>(value);
		}
		else if (equals(param, "voxelsize")) {
			config.vox_size = str2num<double>(value);
		}
		else if (equals(param, "partial")) {
			config.partial = bool(value != "0");
		}
		else if (equals(param, "progress")) {
			config.progress = bool(value != "0");
		}
	}

	// check faulty input
	if (config.num_photons < 1) {
		return false;
	}
	if (!isbetween(config.ambient_eta, 1.0, 3.0)) {
		return false;
	}
	if (config.vox_size < 1E-5) { // too small
		return false;
	}

	return true;
}

/***********************************************************
 * Create a new light source based on the given data.
 ***********************************************************/
bool Simulator::parse_source(std::list<std::string>& data) {
	Source source;

	// set source properties
	for (const auto& it : parameter_values(data)) {
		std::string param = it.first;
		std::string value = it.second;

		if (equals(param, "id")) {
			source.id = str2num<uint64_t>(value);
		}
		else if (equals(param, "position")) {
			std::vector<double> p = split(value, ',');
			if (p.size() < 3) {
				return false;
			}
			source.origin = Point3(p[0], p[1], p[2]);
		}
		else if (equals(param, "direction")) {
			std::vector<double> v = split(value, ',');
			if (v.size() < 3) {
				return false;
			}
			source.direction = Vector3(v[0], v[1], v[2]);
		}
	}

	// check faulty input
	if (source.direction == Vector3(0)) {
		return false;
	}

	// normalize direction
	source.direction.normalize();

	// add to collection
	sources.push_back(source);

	return true;
}

/***********************************************************
 * Create a new tissue type based on the given data.
 ***********************************************************/
bool Simulator::parse_tissue(std::list<std::string>& data) {
	Tissue tissue;

	// set tissue type properties
	for (const auto& it : parameter_values(data)) {
		std::string param = it.first;
		std::string value = it.second;

		if (equals(param, "id")) {
			tissue.id = str2num<uint8_t>(value);
		}
		else if (equals(param, "ani")) {
			tissue.ani = str2num<double>(value);
		}
		else if (equals(param, "eta")) {
			tissue.eta = str2num<double>(value);
		}
		else if (equals(param, "mua")) {
			tissue.mua = str2num<double>(value);
		}
		else if (equals(param, "mus")) {
			tissue.mus = str2num<double>(value);
		}
	}

	// check faulty input
	if (!isbetween(tissue.ani, -1.0, 1.0)) {
		return false;
	}
	if (!isbetween(tissue.eta, 1.0, 3.0)) {
		return false;
	}
	if (tissue.mua < 0) {
		return false;
	}
	if (tissue.mus < 0) {
		return false;
	}

	// add to collection
	tissues.push_back(tissue);

	return true;
}

/***********************************************************
 * Create a new layer based on the given data and initialize
 * its geometry.
 ***********************************************************/
bool Simulator::parse_layer(std::list<std::string>& data) {
	Layer layer;

	std::vector<Index3> faces;           // vertex index tuple (faces)
	std::vector<Point3> vertices;        // triangle vertices
	std::vector<Vector3> normals;        // normal vectors (one per triangle)
	std::vector<Triangle> triangle_mesh; // constructed faces

	// set layer properties
	for (const auto& it : parameter_values(data)) {
		std::string param = it.first;
		std::string value = it.second;

		// identifiers
		if (equals(param, "id")) {
			layer.id = str2num<uint8_t>(value);
		}
		else if (equals(param, "tissue")) {
			layer.tissue_id = str2num<uint8_t>(value);
		}

		// geometric properties
		// vertex
		else if (equals(param, "vert3")) {
			std::vector<double> v = split(value, ',');
			if (v.size() < 3) {
				continue;
			}

			Point3 vert = Point3(v[0], v[1], v[2]);
			vertices.push_back(vert);
		}
		// triangle face
		else if (equals(param, "face3")) {
			std::vector<double> v = split(value, ',');

			// invalid face3 detected
			if (v.size() < 3) {
				return false;
			}

			faces.push_back(
				Index3(static_cast<uint32_t>(v[0]), static_cast<uint32_t>(v[1]), static_cast<uint32_t>(v[2])));
		}
		// normal
		else if (equals(param, "norm3")) {
			std::vector<double> v = split(value, ',');
			if (v.size() < 3) {
				continue;
			}

			Vector3 norm = Vector3(v[0], v[1], v[2], true);
			normals.push_back(norm);
		}
	}

	// check faulty input
	// number of vertices, faces, and edges
	std::size_t num_verts = vertices.size();
	std::size_t num_faces = faces.size();
	std::size_t num_edges = (num_faces * 3) / 2; // also counts internal edges

	// need at least 4 triangular faces to have a polyhedron
	if (num_faces < 4) {
		return false;
	}
	// so, also 4 normals
	if (normals.size() < 4) {
		return false;
	}
	// and at least 4 vertices (tetrahedron)
	if (num_verts < 4) {
		return false;
	}

	// polehedron must be convex (Euler characteristic must be 2)
	if ((num_verts - num_edges + num_faces) != 2) {
		return false;
	}

	// each triangle/face must have a normal vector
	if (faces.size() != normals.size()) {
		return false;
	}

	// construct faces (triangles) from vertex faces
	for (size_t i = 0; i < faces.size(); ++i) {
		Index3 index = faces[i];
		size_t num_verts_size = static_cast<size_t>(num_verts);
		uint32_t index_max = static_cast<uint32_t>(num_verts_size - 1);

		// check vertex index pointer correctness
		if (!isbetween(index.a, 0, index_max) || !isbetween(index.b, 0, index_max)
			|| !isbetween(index.c, 0, index_max)) {
			return false;
		}

		// create new triangle for this face
		Triangle triangle = Triangle(vertices[index.a], vertices[index.b], vertices[index.c]);

		// associate normal with triangle
		triangle.normal = normals[i];

		// add the triangle to the triangle mesh
		triangle_mesh.push_back(triangle);
		triangles.push_back(triangle);
	}

	// set the layer mesh
	layer.mesh = triangle_mesh;

	// add to collection
	layers.push_back(layer);

	return true;
}

/***********************************************************
 * Extract all blocks of data from the input file and store
 * each one as an entry in a multimap.
 ***********************************************************/
bool Simulator::extract(const std::string& fconfig, std::multimap<std::string, std::list<std::string>>& datamap) {
	std::string line;
	std::ifstream in_config(fconfig.c_str(), std::ios_base::in);

	if (!in_config.good()) {
		std::cerr << "File \"" << fconfig << "\" could not be read or opened." << std::endl;
		return false;
	}

	// pre-parse (remove comments/whitespace) and extract data
	while (std::getline(in_config, line)) {
		trim_comment(line);
		trim_spaces(line);

		if (line.empty()) {
			continue;
		}

		if (equals(line, "general")) {
			extract_block(in_config, "general", datamap);
		}

		else if (equals(line, "source")) {
			extract_block(in_config, "source", datamap);
		}

		else if (equals(line, "tissue")) {
			extract_block(in_config, "tissue", datamap);
		}

		else if (equals(line, "layer")) {
			extract_block(in_config, "layer", datamap);
		}
	}

	in_config.close();

	return true;
}

/***********************************************************
 * Extract a block of data from the input file with the
 * given section name.
 ***********************************************************/
void Simulator::extract_block(std::ifstream& in_config, std::string section,
							  std::multimap<std::string, std::list<std::string>>& datamap) {
	std::string line;
	std::list<std::string> lines;

	// store all the lines between curly brackets
	while (std::getline(in_config, line) && line[0] != '}') {
		if (line[0] == '{') {
			std::getline(in_config, line);
		}

		trim_comment(line);
		trim_spaces(line);

		if (line.empty()) {
			continue;
		}

		lines.push_back(line);
	}

	datamap.insert(std::pair<std::string, std::list<std::string>>(section, lines));
}

/***********************************************************
 * Initialize the voxel grid.
 ***********************************************************/
bool Simulator::initialize_grid() {
	// compute grid boundary extent
	for (const auto& layer : layers) {
		// see if a vertex denotes a new boundary
		for (const auto& triangle : layer.mesh) {
			Point3 v0 = triangle.v0;
			Point3 v1 = triangle.v1;
			Point3 v2 = triangle.v2;

			// get the maximum value among the previous maximum or new vertices
			bounds.x_min = min4(bounds.x_min, v0.x, v1.x, v2.x); // left (-x)
			bounds.x_max = max4(bounds.x_max, v0.x, v1.x, v2.x); // right (+x)

			bounds.y_min = min4(bounds.y_min, v0.y, v1.y, v2.y); // top (+y)
			bounds.y_max = max4(bounds.y_max, v0.y, v1.y, v2.y); // bottom (-y)

			bounds.z_min = min4(bounds.z_min, v0.z, v1.z, v2.z); // front (+z)
			bounds.z_max = max4(bounds.z_max, v0.z, v1.z, v2.z); // rear (-z)
		}
	}

	// check for inconsistency and zero width/height/depth
	if (bounds.x_min >= bounds.x_max || bounds.y_min >= bounds.y_max || bounds.z_min >= bounds.z_max) {
		return false;
	}

	// grid dimensions
	bounds.width = fabs(bounds.x_min - bounds.x_max);
	bounds.height = fabs(bounds.y_min - bounds.y_max);
	bounds.depth = fabs(bounds.z_min - bounds.z_max);

	// number of voxels in each dimension (as float)
	float nx = static_cast<float>(bounds.width / config.vox_size);
	float ny = static_cast<float>(bounds.height / config.vox_size);
	float nz = static_cast<float>(bounds.depth / config.vox_size);

	// number of voxels in each dimension
	config.nx = static_cast<uint8_t>(nx);
	config.ny = static_cast<uint8_t>(ny);
	config.nz = static_cast<uint8_t>(nz);

	// total number of voxels
	config.num_voxels =
		static_cast<uint64_t>(config.nx) * static_cast<uint64_t>(config.ny) * static_cast<uint64_t>(config.nz);

	// check for voxel sizes that are too large
	if (config.vox_size > bounds.width || config.vox_size > bounds.height || config.vox_size > bounds.depth) {
		return false;
	}

	// check grid dimensions
	if (config.num_voxels < 1 || config.nx < 1 || config.ny < 1 || config.nz < 1) {
		return false;
	}

	// initialize voxel grid with empty voxels
	voxels = std::vector<Voxel*>(config.num_voxels, nullptr);
	for (uint32_t i = 0; i < config.nx; ++i) {
		for (uint32_t j = 0; j < config.ny; ++j) {
			for (uint32_t k = 0; k < config.nz; ++k) {
				voxels.at(k * config.nx * config.ny + j * config.nx + i) = new Voxel(i, j, k);
			}
		}
	}

	return true;
}

/***********************************************************
 * Initialize configuration properties.
 ***********************************************************/
bool Simulator::initialize_data() {
	// initialize config properties
	config.num_layers = static_cast<uint64_t>(layers.size());
	config.num_sources = static_cast<uint64_t>(sources.size());

	// error checking
	if (config.num_layers < 1 || config.num_voxels < 1) {
		return false;
	}

	// check for duplicate layers
	for (uint32_t i = 1; i < layers.size(); ++i) {
		if (layers[i - 1] == layers[i]) {
			return false;
		}
	}

	// check for duplicate tissues
	for (uint32_t i = 1; i < tissues.size(); ++i) {
		if (tissues[i - 1] == tissues[i]) {
			return false;
		}
	}

	// check if layer's tissue id is out of range
	for (uint32_t i = 0; i < layers.size(); ++i) {
		if (layers[i].tissue_id >= tissues.size()) {
			return false;
		}
	}

	// associate light sources with their geometric intersections
	for (auto& source : sources) {
		// intersection point and triangle that is hit
		Point3 intersect;
		Triangle triangle;

		// find intersection of ray from this source with geometry (point, triangle, normal)
		Ray ray = Ray(source.origin, source.direction);
		if (first_ray_triangle_intersect(ray, triangles, intersect, triangle) == std::numeric_limits<double>::max()) {
			std::cerr << "Error: Source (" << source.id << ") does not intersect the geometry." << std::endl;
			return false;
		}

		source.intersect = intersect;
		source.triangle = triangle;
	}

	// initialize photons
	for (uint64_t i = 0; i < config.num_photons; ++i) {
		photons.push_back(Photon(i));
	}

	return true;
}

/***********************************************************
 * Voxelize the geometry.
 ***********************************************************/
bool Simulator::voxelize_layers() {
	// 1. for each row of each slice: shoot ray from left to right
	// 2. create a range of voxel indices for each occurring layer on that row
	// 3. set the corresponding tissue properties for each voxel based on the ranges

	float vox_pz, vox_py;                                     // voxel centers
	uint32_t vox_iz, vox_iy, vox_ix;                          // voxel indices
	float v_size = static_cast<float>(config.vox_size);       // voxel size
	float h_size = static_cast<float>(config.vox_size / 2.0); // half voxel size
	float ray_x = static_cast<float>(bounds.x_min - 0.0001);  // ray starting x-position

	// per slice (back to front); voxels are ordered from left-bottom-rear to right-top-front
	for (vox_pz = float(bounds.z_min) + h_size, vox_iz = 0; vox_pz <= float(bounds.z_max); vox_pz += v_size, ++vox_iz) {
		// rectify position due to rounding error
		float a = fabsf(fmodf(vox_pz, h_size));
		if (a < 1e-5) {
			vox_pz -= a;
		}

		// per row (bottom to top)
		for (vox_py = float(bounds.y_min) + h_size, vox_iy = 0; vox_py <= float(bounds.y_max);
			 vox_py += v_size, ++vox_iy) {
			// rectify position due to rounding error
			float b = fabsf(fmodf(vox_py, h_size));
			if (b < 1e-5) {
				vox_py -= b;
			}

			// shoot ray from left to right
			Ray ray = Ray(Point3(ray_x, vox_py, vox_pz), Vector3(1, 0, 0));

			for (auto& layer : layers) {
				// find ray intersections with the geometry
				std::vector<Point3> intersections;
				ray_triangles_intersections(ray, layer.mesh, intersections);

				if (intersections.size() != 2) {
					continue;  // next layer
				}

				Point3 il, ir; // set left and right intersection point
				if (intersections[0].x < intersections[1].x) {
					il = intersections[0];
					ir = intersections[1];
				}
				else {
					ir = intersections[0];
					il = intersections[1];
				}

				// distances from boundaries to intersection points
				float dx1 = static_cast<float>(std::fabs(bounds.x_min - il.x));
				float dx2 = static_cast<float>(std::fabs(bounds.x_min - ir.x));

				// indices of first and last voxel in the range
				uint32_t ix1 = static_cast<uint32_t>(dx1 / config.vox_size);
				uint32_t ix2 = static_cast<uint32_t>(dx2 / config.vox_size);

				// avoid index overflow
				if (ix1 >= config.nx) {
					ix1 = config.nx - 1;
				}
				if (ix2 >= config.nx) {
					ix2 = config.nx - 1;
				}

				// set tissue type for all voxels in the range
				for (vox_ix = ix1; vox_ix <= ix2; ++vox_ix) {
					size_t voxel_index =
						static_cast<size_t>(vox_iz) * static_cast<size_t>(config.nx) * static_cast<size_t>(config.ny)
						+ static_cast<size_t>(vox_iy) * static_cast<size_t>(config.nx) + static_cast<size_t>(vox_ix);

					Voxel* voxel = voxels.at(voxel_index);
					voxel->tissue = &tissues[layer.tissue_id];
				}
			}
		}
	}

	return true;
}

/***********************************************************
 * Run the Monte Carlo photon transport simulation.
 ***********************************************************/
void Simulator::simulate() {
	std::cout << "Running Monte Carlo simulation" << std::endl;

	Experimenter::start_clock();

	// For each light source
	for (auto& source : sources) {
		// Compute specular reflection
		specular_reflection(source);

		// for each photon
		uint32_t p = 0;
		for (auto& photon : photons) {
			// progress report
			if (config.progress && ((p + 1) % 1000) == 0) {
				std::cout << "Photon " << p + 1 << "/" << config.num_photons << std::endl;
			}

			// launch the photon (create a new path)
			launch(photons[p], source);

			while (photons[p].alive) {
				step_size(photon); // Set new step size
				transfer(photon);  // Propagate photon through the medium in substeps
				roulette(photon);  // Determine photon termination
				scatter(photon);   // Scatter photon into a new direction
			}
		}
	}

	// Normalize physical quantities
	normalize();

	Experimenter::stop_clock();
	Experimenter::collect_data(record.at, record.rs, record.rd, record.ts, record.td);
	Experimenter::write_to_file();
	Experimenter::print_report();
}

/***********************************************************
 * Set up photon properties for tracing.
 ***********************************************************/
void Simulator::launch(Photon& photon, Source& source) {
	photon.alive = true;
	photon.weight = 1.0 - record.rs;
	photon.source = source;
	photon.direction = source.direction;
	photon.position = source.intersect;
	photon.voxel = voxel_at(photon.position);

	// create vertices for new light path
	Vertex* light = new Vertex(source.origin, photon.weight);
	Vertex* intersection = new Vertex(source.intersect, photon.weight);
	Vertex* reflection = new Vertex(move(source.intersect, source.specular_direction, 0.1), record.rs);

	light->next = intersection;      // intersection vertex/node
	intersection->prev = light;      // light source origin
	intersection->emit = reflection; // specular reflection

	// use intersection point as head
	Graph path = Graph(static_cast<long>(photon.id), intersection);
	paths.push_back(path);

	Experimenter::add_vertex(photon.position.x, photon.position.y, photon.position.z);
}

/***********************************************************
 * Set dimensionless step size for next segment of the
 * random walk using MCML 3.0.0 algorithm.
 ***********************************************************/
void Simulator::step_size(Photon& photon) {
	generate_step_size(photon);
	Experimenter::add_step_size(photon.step / photon.mus());
}

/***********************************************************
 * Transfer a photon through a voxelized medium with
 * individual substeps.
 ***********************************************************/
void Simulator::transfer(Photon& photon) {
	/*
	 * set substep (max = distance to voxel boundary)
	 * deposit weight in current voxel
	 * if (photon crosses voxel boundary)
	 *     if (photon goes outside medium)
	 *         record partial transmission
	 *     else if (photon moves to differing refractive indexed media)
	 *         reflect from or transmit across the voxel boundary
	 *     else if (photon moves to equal refractive indexed media)
	 *         continue normal propagation
	 * decrease step size by traveled distance
	 */

	while (photon.step >= 1E-10 && photon.alive) {
		// set substep
		sub_step(photon);

		// deposit weight
		deposit(photon);

		// possibly cross voxel boundary
		if (photon.cross) {
			cross(photon);
		}
		else {
			photon.position = move(photon.position, photon.direction, photon.sub_step);
		}

		// prevent errors due to crossing to ambient medium
		if (!photon.voxel) {
			photon.alive = false;
			photon.voxel = photon.prev_voxel; // for radiate
			radiate(photon, photon.direction, photon.weight);
		}

		// update step size
		photon.step -= (photon.sub_step * photon.mus());

		Experimenter::add_vertex(photon.position.x, photon.position.y, photon.position.z);
	}
}

/***********************************************************
 * Set the photon's next substep and initialize the given
 * intersection point and voxel normal if it crosses the
 * voxel boundary.
 ***********************************************************/
void Simulator::sub_step(Photon& photon) {
	// create ray and get voxel vertices
	Ray ray = Ray(photon.position, photon.direction);
	Cuboid box = voxel_corners(photon.voxel);

	// find first intersection with voxel faces and get the distance
	double voxdist = first_ray_cuboid_intersect_internal(ray, box, photon.intersect, photon.voxel_normal);

	// check if no intersections were found
	if (voxdist == std::numeric_limits<double>::max()) {
		std::cerr << "Critical error: ray-voxel intersection test failed." << std::endl;
		exit(EXIT_FAILURE);
	}

	// compute free path for a substep (distance to next scattering event)
	double freepath = photon.step / photon.mus();

	// see if the free path crosses the voxel
	if (voxdist == 0) {
		// already on voxel face; move just beyond
		photon.sub_step = config.vox_size * 0.001;
		photon.cross = true;
	}
	// crosses the voxel boundary
	else if (freepath > voxdist) {
		photon.sub_step = voxdist;
		photon.cross = true;
	}
	// does not cross
	else {
		photon.sub_step = freepath;
		photon.cross = false;
	}
}

/***********************************************************
 * Deposit some of the photon's weight into the geometry.
 ***********************************************************/
void Simulator::deposit(Photon& photon) {
	// cancel if photon is outside of medium
	if (!photon.voxel->tissue) {
		return;
	}

	// deposited weight
	double deltaw = photon.weight * (1 - std::exp(-photon.mua() * photon.sub_step));

	// update photon weight
	photon.weight -= deltaw;

	// assign deposited weight to voxel
	photon.voxel->absorption += deltaw;

	// update total absorption
	record.at += deltaw;
}

/***********************************************************
 * Determine the action to take when a photon is about to
 * traverse a voxel face. It can either cross to the ambient
 * medium, to another medium with a differing refractive
 * index, to another medium with the same refractive index,
 * or within the same medium.
 *
 * Reflection and transmission can be handled partially at
 * external-internal medium boundaries, but is always handled
 * as an all-or-none event at internal medium boundaries.
 *
 * If appropriate, the new photon direction, position and
 * voxel are computed.
 ***********************************************************/
void Simulator::cross(Photon& photon) {
	// directions of transmission and reflection
	Vector3 transmittance, reflectance;

	// retrieve the refractive index of the medium that is struck
	Point3 newpos = move_delta(photon.intersect, photon.direction);
	Voxel* newvox = voxel_at(newpos);
	photon.prev_voxel = photon.voxel;

	// determine refractive index
	double eta = (newvox == nullptr) ? config.ambient_eta : newvox->tissue->eta;

	// 1. crossing to ambient medium
	if (newvox == nullptr) {
		// compute reflected/transmitted fraction and reflection/transmission directions
		double reflection = internal_reflection(photon, eta, transmittance, reflectance);
		double transmission = 1.0 - reflection;

		// total transmission
		if (reflection == 0.0) {
			// photon dies
			photon.direction = transmittance;
			photon.alive = false;

			radiate(photon, transmittance, photon.weight);
		}
		// total internal reflection
		else if (reflection == 1.0) {
			// photon reflects off surface
			photon.direction = reflectance;
			photon.position = move_delta(photon.intersect, photon.direction);
			photon.voxel = voxel_at(photon.position);
		}
		else {
			// partial reflection
			if (config.partial) {
				// emit partial reflectance
				radiate(photon, transmittance, photon.weight * transmission);

				// adjust photon weight
				photon.weight *= reflection;

				// update direction/position/voxel
				photon.direction = reflectance;
				photon.position = move_delta(photon.intersect, photon.direction);
				photon.voxel = voxel_at(photon.position);
			}
			else { // all-or-none transmission/reflection
				// total transmission
				if (random_num() > reflection) {
					photon.direction = transmittance;
					photon.alive = false;

					radiate(photon, transmittance, photon.weight);
				}
				// total reflection
				else {
					photon.direction = reflectance;
					photon.position = move_delta(photon.intersect, photon.direction);
					photon.voxel = voxel_at(photon.position);
				}
			}
		}

		Experimenter::increment_scatters();
	}
	// 2. crossing to another medium
	else if (newvox != nullptr && photon.voxel->tissue->eta != newvox->tissue->eta) {
		// compute reflected/transmitted fraction and reflection/transmission directions
		double reflection = internal_reflection(photon, eta, transmittance, reflectance);

		// total transmission
		if (reflection == 0.0) {
			photon.direction = transmittance;
			photon.position = move_delta(photon.intersect, photon.direction);
			photon.voxel = voxel_at(photon.position);
		}
		// total internal reflection
		else if (reflection == 1.0) {
			photon.direction = reflectance;
			photon.position = move_delta(photon.intersect, photon.direction);
			photon.voxel = voxel_at(photon.position);
		}
		else { // all-or-none transmission/reflection
			// total transmission
			if (random_num() > reflection) {
				photon.direction = transmittance;
				photon.position = move_delta(photon.intersect, photon.direction);
				photon.voxel = voxel_at(photon.position);
			}
			// total reflection
			else {
				photon.direction = reflectance;
				photon.position = move_delta(photon.intersect, photon.direction);
				photon.voxel = voxel_at(photon.position);
			}
		}

		Experimenter::increment_scatters();
	}
	// 3. crossing within the same medium (total transmission)
	else {
		// direction is unchanged
		photon.position = move_delta(photon.intersect, photon.direction);
		photon.voxel = voxel_at(photon.position);
	}
}

/***********************************************************
 * Record the emittance from a photon (partially) leaving
 * the material.
 ***********************************************************/
void Simulator::radiate(Photon& photon, Vector3& direction, double weight) {
	// set voxel emittance
	photon.voxel->emittance += weight;

	// add regular vertex that stops at boundary
	paths.back().add_internal_vertex(new Vertex(photon.intersect, weight));
	paths.back().add_external_vertex(new Vertex(move(photon.intersect, direction, 0.1), weight));
	Experimenter::add_vertex(photon.intersect.x, photon.intersect.y, photon.intersect.z);

	// add emitter at this position
	emitters.push_back(Emitter(photon.id, photon.intersect, direction, weight));

	// specular or diffuse transmission
	if (!photon.scatters) {
		record.ts += weight;
	}
	else {
		// cos(theta) = a.b / |a||b|
		double cos_theta =
			dot(photon.source.direction, direction) / (photon.source.direction.length() * direction.length());
		double theta = acos(cos_theta);

		// count as transmission or reflection
		if (theta < HALF_PI) {
			record.td += weight;
		}
		else {
			record.rd += weight;
		}
	}
}

/***********************************************************
 * Determine the survivability of a given photon using MCML 3.0.0 algorithm.
 ***********************************************************/
void Simulator::roulette(Photon& photon) {
	// Use MCML 3.0.0 Russian roulette
	roulette_photon(photon);
}

/***********************************************************
 * Scatter the photon into a new direction based on the
 * Henyey-Greenstein phase function using MCML 3.0.0 algorithm.
 ***********************************************************/
void Simulator::scatter(Photon& photon) {
	if (!photon.alive) {
		return;
	}

	// Get tissue properties for scattering
	Tissue* tissue = photon.voxel->tissue;
	if (!tissue) {
		photon.alive = false;
		return;
	}

	// Use MCML 3.0.0 scattering algorithm
	scatter_photon(photon, *tissue);

	// normalize direction vector (safety check)
	photon.direction.normalize();

	// prevent scattering into ambient medium when close to boundaries
	Point3 newpos = move_delta(photon.position, photon.direction);
	if (!voxel_at(newpos)) {
		photon.alive = false;
		radiate(photon, photon.direction, photon.weight);
	}

	Experimenter::increment_scatters();

	// add new internal position to path
	paths.back().add_internal_vertex(new Vertex(photon.position, photon.weight));
}

/***********************************************************
 * Normalize the recorded values based on the number of
 * photons traced.
 ***********************************************************/
void Simulator::normalize() {
	// normalize globally recorded parameters
	record.at /= config.num_photons * config.num_sources;
	record.rd /= config.num_photons * config.num_sources;
	record.td /= config.num_photons * config.num_sources;

	record.rs /= config.num_sources;                      // rs is only computed once per source
	record.ts /= config.num_photons * config.num_sources; // ts is computed per photon

	// normalize voxel data
	for (const auto& voxel : voxels) {
		// skip computation for voxels outside the medium
		if (!voxel->tissue) {
			continue;
		}

		voxel->absorption /= (config.num_photons * config.num_sources);
		voxel->emittance /= (config.num_photons * config.num_sources);
	}
}

/***********************************************************
 * Compute the specular reflectance from a light source at
 * the surface.
 ***********************************************************/
void Simulator::specular_reflection(Source& source) {
	Voxel* voxel = voxel_at(source.intersect);

	// voxel should never be nullptr at this point
	if (!voxel) {
		std::cerr << "Critical error: specular reflection could not be computed." << std::endl;
		exit(EXIT_FAILURE);
	}

	// refractive indices of ambient medium and medium that is hit
	double n1 = config.ambient_eta;
	double n2 = voxel->tissue->eta;

	// set the specular reflection
	record.rs = (n2 != n1) ? sq((n1 - n2) / (n1 + n2)) : 0;

	// reflection direction: R = -2(V . N)N + V
	Vector3 normal = source.triangle.normal;
	Vector3 projection = normal * dot(source.direction, normal);
	Vector3 rsdir = Vector3((projection * -2.0) + source.direction, true);

	source.specular_direction = rsdir;
}

/***********************************************************
 * Compute the fraction of incoming light that is reflected
 * back at an interface between two media. Also compute the
 * directions of transmission and reflection.
 ***********************************************************/
double Simulator::internal_reflection(Photon& photon, double& eta_t, Vector3& transmittance, Vector3& reflectance) {
	// fraction reflected
	double reflection;

	// indices of refraction
	double eta_i = photon.voxel->tissue->eta;

	// critical angle
	double cos_c = sqrt(1 - sq(eta_t) / sq(eta_i));

	// angles of reflection
	double cos_t = dot(photon.direction, photon.voxel_normal);
	double cos_p = cos_t;

	// compute the fraction of reflected light

	// matched boundary
	if (eta_i == eta_t) {
		cos_p = cos_t;
		reflection = 0.0;
	}
	// near-perpendicular incidence: cos(0)
	else if (cos_t >= 1.0 - 1E-10) {
		reflection = sq(eta_t - eta_i) / sq(eta_t + eta_i);
		cos_p = cos_t;
	}
	// near-parallel incidence: cos(90)
	else if (cos_t < 1E-10) {
		reflection = 1.0;
		cos_p = 0;
	}
	// total internal reflection
	else if ((eta_t < eta_i) && cos_t < cos_c) {
		reflection = 1.0;
		cos_p = 0;
	}
	// general case
	else {
		// angle phi of reflection
		cos_p = sqrt(1 - sq(eta_i) / sq(eta_t) * (1 - sq(cos_t)));

		double temp1 = eta_i * cos_p;
		double temp2 = eta_t * cos_t;
		double temp3 = eta_t * cos_p;
		double temp4 = eta_i * cos_t;

		// fraction reflected
		reflection = 0.5 * (sq(temp1 - temp2) / sq(temp1 + temp2) + sq(temp3 - temp4) / sq(temp3 + temp4));
	}

	// direction of transmission
	double temp1 = eta_i / eta_t;
	double temp2 = cos_p - cos_t * temp1;

	transmittance.x = photon.direction.x * temp1 + temp2 * photon.voxel_normal.x;
	transmittance.y = photon.direction.y * temp1 + temp2 * photon.voxel_normal.y;
	transmittance.z = photon.direction.z * temp1 + temp2 * photon.voxel_normal.z;

	// direction of reflection
	reflectance.x = photon.direction.x - 2 * cos_t * photon.voxel_normal.x;
	reflectance.y = photon.direction.y - 2 * cos_t * photon.voxel_normal.y;
	reflectance.z = photon.direction.z - 2 * cos_t * photon.voxel_normal.z;

	transmittance.normalize();
	reflectance.normalize();

	return reflection;
}

/***********************************************************
 * Return a pointer to a voxel that encapsulates the given
 * position, or nullptr if the position is outside the medium.
 ***********************************************************/
Voxel* Simulator::voxel_at(Point3& position) {
	if (!bounds.includes(position.x, position.y, position.z)) {
		return nullptr;
	}

	// distances from boundaries
	double dx = std::fabs(bounds.x_min - position.x);
	double dy = std::fabs(bounds.y_min - position.y);
	double dz = std::fabs(bounds.z_min - position.z);

	// indices start at minimum boundaries of voxels
	uint32_t ix = static_cast<uint32_t>(std::floor(dx / config.vox_size));
	uint32_t iy = static_cast<uint32_t>(std::floor(dy / config.vox_size));
	uint32_t iz = static_cast<uint32_t>(std::floor(dz / config.vox_size));

	// avoid index overflow
	if (ix >= config.nx) {
		ix = config.nx - 1;
	}
	if (iy >= config.ny) {
		iy = config.ny - 1;
	}
	if (iz >= config.nz) {
		iz = config.nz - 1;
	}

	// retrieve the voxel at the position
	size_t voxel_index = static_cast<size_t>(iz) * static_cast<size_t>(config.nx) * static_cast<size_t>(config.ny)
						 + static_cast<size_t>(iy) * static_cast<size_t>(config.nx) + static_cast<size_t>(ix);
	Voxel* voxel = voxels.at(voxel_index);

	// if voxel does not have a tissue, it is outside the medium
	return voxel->tissue ? voxel : nullptr;
}

/***********************************************************
 * Return the minimum and maximum positions of the given
 * voxel as a cuboid structure.
 ***********************************************************/
Cuboid Simulator::voxel_corners(Voxel* voxel) {
	uint32_t ix_min = voxel->ix;
	uint32_t iy_min = voxel->iy;
	uint32_t iz_min = voxel->iz;

	uint32_t ix_max = voxel->ix + 1;
	uint32_t iy_max = voxel->iy + 1;
	uint32_t iz_max = voxel->iz + 1;

	// minimum voxel position
	float x_min = static_cast<float>(bounds.x_min + (config.vox_size * ix_min));
	float y_min = static_cast<float>(bounds.y_min + (config.vox_size * iy_min));
	float z_min = static_cast<float>(bounds.z_min + (config.vox_size * iz_min));

	// maximum voxel position
	float x_max = static_cast<float>(bounds.x_min + (config.vox_size * ix_max));
	float y_max = static_cast<float>(bounds.y_min + (config.vox_size * iy_max));
	float z_max = static_cast<float>(bounds.z_min + (config.vox_size * iz_max));

	// round off coordinate values around the origin
	x_min = (std::fabs(x_min) < 1E-10) ? 0 : x_min;
	y_min = (std::fabs(y_min) < 1E-10) ? 0 : y_min;
	z_min = (std::fabs(z_min) < 1E-10) ? 0 : z_min;

	x_max = (std::fabs(x_max) < 1E-10) ? 0 : x_max;
	y_max = (std::fabs(y_max) < 1E-10) ? 0 : y_max;
	z_max = (std::fabs(z_max) < 1E-10) ? 0 : z_max;

	return Cuboid(x_min, y_min, z_min, x_max, y_max, z_max);
}

/***********************************************************
 * Return the destination for a given origin, direction and
 * distance.
 ***********************************************************/
Point3 Simulator::move(Point3& position, Vector3& direction, double d) {
	Point3 point = position;

	// return the end point of a sub-step
	point.x = position.x + direction.x * d;
	point.y = position.y + direction.y * d;
	point.z = position.z + direction.z * d;

	return point;
}

/***********************************************************
 * Return the destination for a given origin and direction
 * after making a small hop.
 ***********************************************************/
Point3 Simulator::move_delta(Point3& position, Vector3& direction) {
	Point3 point = position;

	// delta distance (based on voxel size)
	double d = config.vox_size * 0.00001;

	point.x = position.x + direction.x * d;
	point.y = position.y + direction.y * d;
	point.z = position.z + direction.z * d;

	return point;
}

/***********************************************************
 *	Write the resulting physical quantities to a file.
 ***********************************************************/
void Simulator::report() {
	std::string str_sim = "simulation.out";
	std::string str_abs = "absorption.out";
	std::string str_emi = "emittance.out";
	std::string str_ptd = "photons.out";

	std::ofstream ofs_rep(str_sim.c_str(), std::ios_base::out); // simulation report
	std::ofstream ofs_abs(str_abs.c_str(), std::ios_base::out); // absorption report
	std::ofstream ofs_emi(str_emi.c_str(), std::ios_base::out); // emittance report
	std::ofstream ofs_ptn(str_ptd.c_str(), std::ios_base::out); // photon exitance report

	if (!ofs_rep.good() || !ofs_abs.good() || !ofs_emi.good() || !ofs_ptn.good()) {
		std::cerr << "Error: an output file could not be opened." << std::endl;
		return;
	}

	if (ofs_rep.good()) {
		ofs_rep.precision(8);
		ofs_rep << "################################################################" << std::endl;
		ofs_rep << "# SIMULATION REPORT" << std::endl;
		ofs_rep << "################################################################" << std::endl;
		ofs_rep << std::endl;

		// write input configuration
		ofs_rep << "Configuration" << std::endl;
		ofs_rep << "################################################################" << std::endl;
		ofs_rep << std::endl;
		ofs_rep << "Number of photons: " << tab << config.num_photons << std::endl;
		ofs_rep << "Number of layers:  " << tab << config.num_layers << std::endl;
		ofs_rep << "Number of voxels:  " << tab << config.num_voxels << std::endl;
		ofs_rep << "Grid dimensions:   " << tab << config.nx << " x " << config.ny << " x " << config.nz << std::endl;
		ofs_rep << "Voxel dimensions:  " << tab << config.vox_size << " x " << config.vox_size << " x "
				<< config.vox_size << std::endl;
		ofs_rep << std::endl << std::endl;

		// write recorded parameters: a, rs, rd, (ts, td)
		ofs_rep << "Recorded parameters" << std::endl;
		ofs_rep << "################################################################" << std::endl;
		ofs_rep << "Total absorption:      " << tab << std::fixed << record.at << std::endl;
		ofs_rep << "Diffuse reflection:    " << tab << std::fixed << record.rd << std::endl;
		ofs_rep << "Specular reflection:   " << tab << std::fixed << record.rs << std::endl;
		ofs_rep << "Diffuse transmission:  " << tab << std::fixed << record.td << std::endl;
		ofs_rep << "Specular transmission: " << tab << std::fixed << record.ts << std::endl;
		ofs_rep << std::endl;
		ofs_rep.close();
	}
	else {
		std::cerr << "Error: file " << str_sim << " could not be opened." << std::endl;
		ofs_rep.close();
	}

	// write voxel absorption (each 'block' is a slice)
	if (ofs_abs.good()) {
		ofs_abs.precision(5);
		ofs_abs << "################################################################" << std::endl;
		ofs_abs << "# ABSORPTION REPORT" << std::endl;
		ofs_abs << "################################################################" << std::endl;
		ofs_abs << std::endl;

		// from rear to front
		for (uint32_t iz = 0; iz < config.nz; ++iz) {
			ofs_abs << "Slice " << iz + 1 << "/" << config.nz;
			if (iz == 0) {
				ofs_abs << " (rear)";
			}
			if (iz == config.nz - 1) {
				ofs_abs << " (front)";
			}
			ofs_abs << std::endl;

			for (uint32_t iy = config.ny - 1; iy > 0; --iy) { // top to bottom
				for (uint32_t ix = 0; ix < config.nx; ++ix) { // left to right
					size_t voxel_index =
						static_cast<size_t>(iz) * static_cast<size_t>(config.nx) * static_cast<size_t>(config.ny)
						+ static_cast<size_t>(iy) * static_cast<size_t>(config.nx) + static_cast<size_t>(ix);
					Voxel* voxel = voxels.at(voxel_index);
					ofs_abs << std::fixed << voxel->absorption << tab;
				}
				ofs_abs << std::endl;
			}
			ofs_abs << std::endl;
		}
		ofs_abs.close();
	}
	else {
		std::cerr << "Error: file " << str_abs << " could not be opened." << std::endl;
		ofs_abs.close();
	}

	// write voxel emittance
	if (ofs_emi.good()) {
		ofs_emi.precision(5);
		ofs_emi << "################################################################" << std::endl;
		ofs_emi << "# EMITTANCE REPORT" << std::endl;
		ofs_emi << "################################################################" << std::endl;
		ofs_emi << std::endl;

		// rear to front
		for (uint32_t iz = 0; iz < config.nz; ++iz) {
			ofs_emi << "Slice " << iz + 1 << "/" << config.nz;
			if (iz == 0) {
				ofs_emi << " (rear)";
			}
			if (iz == config.nz - 1) {
				ofs_emi << " (front)";
			}
			ofs_emi << std::endl;

			for (uint32_t iy = config.ny - 1; iy > 0; --iy) { // top to bottom
				for (uint32_t ix = 0; ix < config.nx; ++ix) { // left to right
					size_t voxel_index =
						static_cast<size_t>(iz) * static_cast<size_t>(config.nx) * static_cast<size_t>(config.ny)
						+ static_cast<size_t>(iy) * static_cast<size_t>(config.nx) + static_cast<size_t>(ix);
					Voxel* voxel = voxels.at(voxel_index);
					ofs_emi << std::fixed << voxel->emittance << tab;
				}
				ofs_emi << std::endl;
			}
			ofs_emi << std::endl;
		}
		ofs_emi.close();
	}
	else {
		std::cerr << "Error: file " << str_emi << " could not be opened." << std::endl;
		ofs_emi.close();
	}

	// exiting photons (position, direction, weight)
	if (ofs_ptn.good()) {
		ofs_ptn.precision(8);
		ofs_ptn << "################################################################" << std::endl;
		ofs_ptn << "# PHOTON REPORT" << std::endl;
		ofs_ptn << "################################################################" << std::endl;
		ofs_ptn << std::endl;

		if (emitters.empty()) {
			ofs_ptn << "No photons exited the medium." << std::endl;
		}

		for (const auto& emitter : emitters) {
			ofs_ptn << "Photon" << std::endl;
			ofs_ptn << "{" << std::endl;
			ofs_ptn << std::fixed << tab << "id  = " << emitter.id << std::endl;
			ofs_ptn << std::fixed << tab << "position = " << emitter.position.x << ", " << emitter.position.y << ", "
					<< emitter.position.z << std::endl;
			ofs_ptn << std::fixed << tab << "direction = " << emitter.direction.x << ", " << emitter.direction.y << ", "
					<< emitter.direction.z << std::endl;
			ofs_ptn << std::fixed << tab << "val = " << emitter.weight << std::endl;
			ofs_ptn << "}" << std::endl;
			ofs_ptn << std::endl;
		}
		ofs_ptn.close();
	}
	else {
		std::cerr << "Error: file " << str_ptd << " could not be opened." << std::endl;
		ofs_ptn.close();
	}
}

/***********************************************************
 * MCML 3.0.0 Monte Carlo Methods - Integrated Implementation
 ***********************************************************/

void Simulator::generate_step_size(Photon& photon) {
	// MCML 3.0.0 step size generation using Beer-Lambert law
	if (photon.step < 1e-10) {
		double rnd;
		// Avoid zero random number
		while ((rnd = mcml_random->next()) <= 0.0) {
			// Keep generating until we get a non-zero value
		}
		photon.step = -std::log(rnd);
	}
}

void Simulator::scatter_photon(Photon& photon, const Tissue& tissue) {
	// MCML 3.0.0 scattering implementation using Henyey-Greenstein phase function
	double cos_theta, sin_theta, cos_phi, sin_phi;
	double g = tissue.ani; // anisotropy factor
	double rnd;

	// Sample scattering angle using Henyey-Greenstein phase function
	rnd = mcml_random->next();
	if (std::abs(g) > 1e-6) {
		double temp = (1.0 - g * g) / (1.0 - g + 2.0 * g * rnd);
		cos_theta = (1.0 + g * g - temp * temp) / (2.0 * g);

		// Ensure cos_theta is within valid range
		if (cos_theta > 1.0) {
			cos_theta = 1.0;
		}
		else if (cos_theta < -1.0) {
			cos_theta = -1.0;
		}
	}
	else {
		// Isotropic scattering
		cos_theta = 2.0 * rnd - 1.0;
	}

	sin_theta = std::sqrt(1.0 - cos_theta * cos_theta);

	// Sample azimuthal angle
	rnd = mcml_random->next();
	cos_phi = std::cos(2.0 * std::numbers::pi * rnd);
	sin_phi = std::sin(2.0 * std::numbers::pi * rnd);

	// Update direction using MCML 3.0.0 spinning algorithm
	if (std::abs(photon.direction.z) > 0.99999) {
		// Special case: nearly perpendicular to z-axis
		photon.direction.x = sin_theta * cos_phi;
		photon.direction.y = sin_theta * sin_phi;
		photon.direction.z = cos_theta * (photon.direction.z > 0 ? 1.0 : -1.0);
	}
	else {
		// General case
		double temp = std::sqrt(1.0 - photon.direction.z * photon.direction.z);
		double temp_x =
			sin_theta * (photon.direction.x * photon.direction.z * cos_phi - photon.direction.y * sin_phi) / temp
			+ photon.direction.x * cos_theta;
		double temp_y =
			sin_theta * (photon.direction.y * photon.direction.z * cos_phi + photon.direction.x * sin_phi) / temp
			+ photon.direction.y * cos_theta;
		double temp_z = -sin_theta * cos_phi * temp + photon.direction.z * cos_theta;

		photon.direction.x = temp_x;
		photon.direction.y = temp_y;
		photon.direction.z = temp_z;
	}

	// Mark that photon has scattered
	photon.scatters = true;
}

void Simulator::roulette_photon(Photon& photon) {
	// MCML 3.0.0 Russian roulette implementation
	if (photon.weight < mcml_weight_threshold) {
		if (mcml_random->next() <= 0.1) {
			// Survive with weight boost
			photon.weight *= 10.0;
		}
		else {
			// Terminate photon
			photon.alive = false;
		}
	}
}

void Simulator::set_rng_seed(int seed) {
	mcml_random->seed(seed);
}
