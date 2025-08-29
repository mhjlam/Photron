#include "../structs/index3.hpp"
#include "../structs/range1.hpp"
#include "../structs/ray.hpp"
#include "../utilities/utilities.hpp"
#include "simulator.hpp"
#include "random.hpp"

#include <algorithm>
#include <cfloat>
#include <ctime>
#include <fstream>
#include <iostream>
#include <sstream>

typedef unsigned int uint;
typedef unsigned long ulong;
typedef unsigned short ushort;
typedef std::pair<std::string, std::string> stringpair;

using namespace std;

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
	for (ulong i = 0; i < voxels.size(); ++i) {
		delete voxels[i];
		voxels[i] = NULL;
	}

	// delete path vertices
	for (uint i = 0; i < paths.size(); ++i) {
		Graph path = paths[i];

		Vertex* item = path.head;

		if (path.head->prev) {
			delete path.head->prev;
			path.head->prev = NULL;
		}

		while (item) {
			Vertex* old = item;
			item = item->next;

			delete old;
			old = NULL;
		}
	}
}

/***********************************************************
 * INITIALIZATION
 ***********************************************************/

/***********************************************************
 * Parse the input file and initializes the data structures.
 ***********************************************************/
bool Simulator::Initialize(string file) {
	cout << "Initializing Photron" << endl;

	// read and parse input configuration file
	cout << "Parsing configuration file: " << file << endl;
	if (!Parse(file)) {
		cerr << "An error occurred while parsing the input file." << endl;
		return false;
	}
	cout << "Configuration parsed successfully." << endl;

	// initialize voxel grid
	if (!InitializeGrid()) {
		cerr << "An error occurred while initializing the voxel grid." << endl;
		return false;
	}
	cout << "Voxel grid initialized successfully." << endl;

	// initialize other data structures
	if (!InitializeData()) {
		cerr << "An error occurred while initializing the data structures." << endl;
		return false;
	}
	cout << "Data structures initialized successfully." << endl;

	// voxelize the geometry
	if (!VoxelizeLayers()) {
		cerr << "An error occurred during geometry voxelization." << endl;
		return false;
	}
	cout << "Geometry voxelization completed successfully." << endl;
	return true;
}

/***********************************************************
 *	INITIALIZATION SUBROUTINES
 ***********************************************************/

/***********************************************************
 * Read and parse the given configuration file.
 ***********************************************************/
bool Simulator::Parse(const std::string& fconfig) {
	std::multimap<std::string, std::list<std::string> > datamap;

	// attempt to extract all the data from the input file
	if (!Simulator::Extract(fconfig, datamap)) {
		cerr << "Failed to extract data from config file." << endl;
		return false;
	}

	// parse the individual configuration blocks
	for (std::multimap<std::string, std::list<std::string> >::iterator it = datamap.begin(); it != datamap.end();
		 ++it) {
		std::string name = it->first;
		std::list<std::string> data = it->second;

		if (equals(name, "general") && !Simulator::ParseGeneral(data)) {
			cerr << "Failed to parse general section." << endl;
			return false;
		}
		else if (equals(name, "source") && !Simulator::ParseSource(data)) {
			cerr << "Failed to parse source section." << endl;
			return false;
		}
		else if (equals(name, "tissue") && !Simulator::ParseTissue(data)) {
			cerr << "Failed to parse tissue section." << endl;
			return false;
		}
		else if (equals(name, "layer") && !Simulator::ParseLayer(data)) {
			cerr << "Failed to parse layer section." << endl;
			return false;
		}
	}

	return true;
}

/***********************************************************
 * Set the general scene configuration.
 ***********************************************************/
bool Simulator::ParseGeneral(list<string>& data) {
	vector<stringpair> generaldata = getparvals(data);

	// set general settings
	for (ushort i = 0; i < generaldata.size(); ++i) {
		string param = generaldata[i].first;
		string value = generaldata[i].second;

		if (equals(param, "numphotons")) {
			config.numphotons = str2num<ulong>(value);
		}
		else if (equals(param, "ambienteta")) {
			config.ambienteta = str2num<double>(value);
		}
		else if (equals(param, "voxelsize")) {
			config.voxsize = str2num<double>(value);
		}
		else if (equals(param, "partial")) {
			config.partial = bool(value != "0");
		}
		else if (equals(param, "progress")) {
			config.progress = bool(value != "0");
		}
	}

	// check faulty input
	if (config.numphotons < 1) {
		return false;
	}
	if (!isbetween(config.ambienteta, 1.0, 3.0)) {
		return false;
	}
	if (config.voxsize < 1E-5) { // too small
		return false;
	}

	return true;
}

/***********************************************************
 * Create a new light source based on the given data.
 ***********************************************************/
bool Simulator::ParseSource(list<string>& data) {
	Source source;
	vector<stringpair> sourcedata = getparvals(data);

	// set source properties
	for (ushort i = 0; i < sourcedata.size(); ++i) {
		string param = sourcedata[i].first;
		string value = sourcedata[i].second;

		if (equals(param, "id")) {
			source.id = str2num<ulong>(value);
		}
		else if (equals(param, "position")) {
			vector<double> p = split(value, ',');
			if (p.size() < 3) {
				return false;
			}

			source.orig = Point3(p[0], p[1], p[2]);
		}
		else if (equals(param, "direction")) {
			vector<double> v = split(value, ',');
			if (v.size() < 3) {
				return false;
			}

			source.dir = Vector3(v[0], v[1], v[2]);
		}
	}

	// check faulty input  
	if (source.dir == Vector3(0)) {
		return false;
	}

	// normalize direction
	source.dir.normalize();

	// add to collection
	sources.push_back(source);

	return true;
}

/***********************************************************
 * Create a new tissue type based on the given data.
 ***********************************************************/
bool Simulator::ParseTissue(list<string>& data) {
	Tissue tissue;
	vector<stringpair> tissuedata = getparvals(data);

	// set tissue type properties
	for (ushort i = 0; i < tissuedata.size(); ++i) {
		string param = tissuedata[i].first;
		string value = tissuedata[i].second;

		if (equals(param, "id")) {
			tissue.id = str2num<ushort>(value);
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
bool Simulator::ParseLayer(list<string>& data) {
	Layer layer;

	vector<Index3> faces;     // vertex index tuple (faces)
	vector<Point3> vertices;  // triangle vertices
	vector<Vector3> normals;  // normal vectors (one per triangle)
	vector<Triangle> trimesh; // constructed faces

	vector<stringpair> layerdata = getparvals(data);

	// set layer properties
	for (ushort i = 0; i < layerdata.size(); ++i) {
		string param = layerdata[i].first;
		string value = layerdata[i].second;

		// identifiers
		if (equals(param, "id")) {
			layer.id = str2num<ushort>(value);
		}
		else if (equals(param, "tissue")) {
			layer.tissueid = str2num<ushort>(value);
		}

		// geometric properties
		else if (equals(param, "vert3")) // vertex
		{
			vector<double> v = split(value, ',');
			if (v.size() < 3) {
				continue;
			}

			Point3 vert = Point3(v[0], v[1], v[2]);
			vertices.push_back(vert);
		}
		else if (equals(param, "face3")) // triangle face
		{
			vector<double> v = split(value, ',');
			if (v.size() < 3) {
				return false;
			} // invalid face3 detected

			faces.push_back(Index3((uint)v[0], (uint)v[1], (uint)v[2]));
		}
		else if (equals(param, "norm3")) // normal
		{
			vector<double> v = split(value, ',');
			if (v.size() < 3) {
				continue;
			}

			Vector3 norm = Vector3(v[0], v[1], v[2], true);
			normals.push_back(norm);
		}
	}

	// check faulty input
	// number of vertices, faces, and edges
	std::size_t numverts = vertices.size();
	std::size_t numfaces = faces.size();
	std::size_t numedges = (numfaces * 3) / 2; // also counts internal edges

	// need at least 4 triangular faces to have a polyhedron
	if (numfaces < 4) {
		return false;
	}
	// so, also 4 normals
	if (normals.size() < 4) {
		return false;
	}
	// and at least 4 vertices (tetrahedron)
	if (numverts < 4) {
		return false;
	}

	// polehedron must be convex (Euler characteristic must be 2)
	if ((numverts - numedges + numfaces) != 2) {
		return false;
	}

	// each triangle/face must have a normal vector
	if (faces.size() != normals.size()) {
		return false;
	}

	// construct faces (triangles) from vertex faces
	for (size_t i = 0; i < faces.size(); ++i) {
		Index3 index = faces[i];
		size_t numverts_size = static_cast<size_t>(numverts);
		unsigned int indexmax = static_cast<unsigned int>(numverts_size - 1);

		// check vertex index pointer correctness
		if (!isbetween(index.a, 0, indexmax) || !isbetween(index.b, 0, indexmax) || !isbetween(index.c, 0, indexmax)) {
			return false;
		}

		// create new triangle for this face
		Triangle triangle = Triangle(vertices[index.a], vertices[index.b], vertices[index.c]);

		// associate normal with triangle
		triangle.normal = normals[i];

		// add the triangle to the triangle mesh
		trimesh.push_back(triangle);
		triangles.push_back(triangle);
	}

	// set the layer mesh
	layer.mesh = trimesh;

	// add to collection
	layers.push_back(layer);

	return true;
}

/***********************************************************
 * Extract all blocks of data from the input file and store
 * each one as an entry in a multimap.
 ***********************************************************/
bool Simulator::Extract(const string& fconfig, multimap<string, list<string> >& datamap) {
	string line;
	ifstream inconfig(fconfig.c_str(), ios_base::in);

	if (!inconfig.good()) {
		cerr << "File \"" << fconfig << "\" could not be read or opened." << endl;
		return false;
	}

	// pre-parse (remove comments/whitespace) and extract data
	while (std::getline(inconfig, line)) {
		trimcomment(line);
		trimspaces(line);

		if (line.empty()) {
			continue;
		}

		if (equals(line, "general")) {
			ExtractBlock(inconfig, "general", datamap);
		}

		else if (equals(line, "source")) {
			ExtractBlock(inconfig, "source", datamap);
		}

		else if (equals(line, "tissue")) {
			ExtractBlock(inconfig, "tissue", datamap);
		}

		else if (equals(line, "layer")) {
			ExtractBlock(inconfig, "layer", datamap);
		}
	}

	inconfig.close();

	return true;
}

/***********************************************************
 * Extract a block of data from the input file with the
 * given section name.
 ***********************************************************/
void Simulator::ExtractBlock(ifstream& inconfig, string section, multimap<string, list<string> >& datamap) {
	string line;
	list<string> lines;

	// store all the lines between curly brackets
	while (std::getline(inconfig, line) && line[0] != '}') {
		if (line[0] == '{') {
			std::getline(inconfig, line);
		}

		trimcomment(line);
		trimspaces(line);

		if (line.empty()) {
			continue;
		}

		lines.push_back(line);
	}

	datamap.insert(pair<string, list<string> >(section, lines));
}

/***********************************************************
 * Initialize the voxel grid.
 ***********************************************************/
bool Simulator::InitializeGrid() {
	// compute grid boundary extent
	for (vector<Layer>::iterator layer = layers.begin(); layer != layers.end(); ++layer) {
		// see if a vertex denotes a new boundary
		for (vector<Triangle>::iterator triangle = layer->mesh.begin(); triangle != layer->mesh.end(); ++triangle) {
			Point3 v0 = triangle->vertex0;
			Point3 v1 = triangle->vertex1;
			Point3 v2 = triangle->vertex2;

			// get the maximum value among the previous maximum or new vertices
			bounds.xmin = min4(bounds.xmin, v0.x, v1.x, v2.x); // left (-x)
			bounds.xmax = max4(bounds.xmax, v0.x, v1.x, v2.x); // right (+x)

			bounds.ymin = min4(bounds.ymin, v0.y, v1.y, v2.y); // top (+y)
			bounds.ymax = max4(bounds.ymax, v0.y, v1.y, v2.y); // bottom (-y)

			bounds.zmin = min4(bounds.zmin, v0.z, v1.z, v2.z); // front (+z)
			bounds.zmax = max4(bounds.zmax, v0.z, v1.z, v2.z); // rear (-z)
		}
	}

	// check for inconsistency and zero width/height/depth
	if (bounds.xmin >= bounds.xmax || bounds.ymin >= bounds.ymax || bounds.zmin >= bounds.zmax) {
		return false;
	}

	// grid dimensions
	bounds.width = fabs(bounds.xmin - bounds.xmax);
	bounds.height = fabs(bounds.ymin - bounds.ymax);
	bounds.depth = fabs(bounds.zmin - bounds.zmax);

	// number of voxels in each dimension (as float)
	float nx = static_cast<float>(bounds.width / config.voxsize);
	float ny = static_cast<float>(bounds.height / config.voxsize);
	float nz = static_cast<float>(bounds.depth / config.voxsize);

	// number of voxels in each dimension
	config.nx = static_cast<ushort>(nx);
	config.ny = static_cast<ushort>(ny);
	config.nz = static_cast<ushort>(nz);

	// total number of voxels
	config.numvoxels = static_cast<unsigned long>(config.nx) * static_cast<unsigned long>(config.ny) * static_cast<unsigned long>(config.nz);

	// check for voxel sizes that are too large
	if (config.voxsize > bounds.width || config.voxsize > bounds.height || config.voxsize > bounds.depth) {
		return false;
	}

	// check grid dimensions
	if (config.numvoxels < 1 || config.nx < 1 || config.ny < 1 || config.nz < 1) {
		return false;
	}

	// initialize voxel grid with empty voxels
	voxels = std::vector<Voxel*>(config.numvoxels, NULL);
	for (ulong i = 0; i < config.nx; ++i) {
		for (ulong j = 0; j < config.ny; ++j) {
			for (ulong k = 0; k < config.nz; ++k) {
				voxels.at(k * config.nx * config.ny + j * config.nx + i) = new Voxel(i, j, k);
			}
		}
	}

	return true;
}

/***********************************************************
 * Initialize configuration properties.
 ***********************************************************/
bool Simulator::InitializeData() {
	// initialize config properties
	config.numlayers = static_cast<ulong>(layers.size());
	config.numsources = static_cast<ulong>(sources.size());

	// error checking
	if (config.numlayers < 1 || config.numvoxels < 1) {
		return false;
	}

	// check for duplicate layers
	for (uint i = 1; i < layers.size(); ++i) {
		if (layers[i - 1] == layers[i]) {
			return false;
		}
	}

	// check for duplicate tissues
	for (uint i = 1; i < tissues.size(); ++i) {
		if (tissues[i - 1] == tissues[i]) {
			return false;
		}
	}

	// check if layer's tissue id is out of range
	for (uint i = 0; i < layers.size(); ++i) {
		if (layers[i].tissueid >= tissues.size()) {
			return false;
		}
	}

	// associate light sources with their geometric intersections
	for (vector<Source>::iterator source = sources.begin(); source != sources.end(); ++source) {
		// intersection point and triangle that is hit
		Point3 inter;
		Triangle triangle;

		// find intersection of ray from this source with geometry (point, triangle, normal)
		Ray ray = Ray(source->orig, source->dir);
		if (first_ray_triangle_intersect(ray, triangles, inter, triangle) == DBL_MAX) {
			cerr << "Error: Source (" << source->id << ") does not intersect the geometry." << endl;
			return false;
		}

		source->inter = inter;
		source->intertri = triangle;
	}

	// initialize photons
	for (ulong i = 0; i < config.numphotons; ++i) {
		photons.push_back(Photon(i));
	}

	return true;
}

/***********************************************************
 * Voxelize the geometry.
 ***********************************************************/
bool Simulator::VoxelizeLayers() {
	// 1. for each row of each slice: shoot ray from left to right
	// 2. create a range of voxel indices for each occurring layer on that row
	// 3. set the corresponding tissue properties for each voxel based on the ranges

	float voxpz, voxpy;                 // voxel centers
	ushort voxiz, voxiy, voxix;         // voxel indices

	float vsize = static_cast<float>(config.voxsize);       // voxel size
	float hsize = static_cast<float>(config.voxsize / 2.0); // half voxel size

	float rayx = static_cast<float>(bounds.xmin - 0.0001);  // ray starting x-position

	// per slice (back to front); voxels are ordered from left-bottom-rear to right-top-front
	for (voxpz = float(bounds.zmin) + hsize, voxiz = 0; voxpz <= float(bounds.zmax); voxpz += vsize, ++voxiz) {
		// rectify position due to rounding error
		float a = fabsf(fmodf(voxpz, hsize));
		if (a < 1e-5) {
			voxpz -= a;
		}

		// per row (bottom to top)
		for (voxpy = float(bounds.ymin) + hsize, voxiy = 0; voxpy <= float(bounds.ymax); voxpy += vsize, ++voxiy) {
			// rectify position due to rounding error
			float b = fabsf(fmodf(voxpy, hsize));
			if (b < 1e-5) {
				voxpy -= b;
			}

			// shoot ray from left to right
			Ray ray = Ray(Point3(rayx, voxpy, voxpz), Vector3(1, 0, 0));

			for (vector<Layer>::iterator layer = layers.begin(); layer != layers.end(); ++layer) {
				// find ray intersections with the geometry
				vector<Point3> intersections;
				ray_triangles_intersections(ray, layer->mesh, intersections);

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
				float dx1 = static_cast<float>(fabs(bounds.xmin - il.x));
				float dx2 = static_cast<float>(fabs(bounds.xmin - ir.x));

				// indices of first and last voxel in the range
				float ix1f = static_cast<float>(dx1 / config.voxsize);
				float ix2f = static_cast<float>(dx2 / config.voxsize);

				// indices start at minimum boundaries of voxels
				ushort ix1 = static_cast<ushort>(ix1f);
				ushort ix2 = static_cast<ushort>(ix2f);

				// avoid index overflow
				if (ix1 >= config.nx) {
					ix1 = config.nx - 1;
				}
				if (ix2 >= config.nx) {
					ix2 = config.nx - 1;
				}

				// set tissue type for all voxels in the range
				for (voxix = ix1; voxix <= ix2; ++voxix) {
					size_t voxel_index = static_cast<size_t>(voxiz) * static_cast<size_t>(config.nx) * static_cast<size_t>(config.ny) + 
					                     static_cast<size_t>(voxiy) * static_cast<size_t>(config.nx) + static_cast<size_t>(voxix);
					Voxel* voxel = voxels.at(voxel_index);
					voxel->tissue = &tissues[layer->tissueid];
				}
			}
		}
	}

	return true;
}