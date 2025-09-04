#include "config.hpp"

#include <fstream>
#include <iostream>
#include <utility>

#include "math/triangle.hpp"

Config::Config() :
	nx_(0), ny_(0), nz_(0), num_layers_(0), num_voxels_(0), num_photons_(0), num_sources_(0), vox_size_(0.0),
	ambient_eta_(0.0), partial_(true), progress_(false) {
}

bool Config::parse_config_file(const std::string& filename) {
	// Clear existing data
	sources_.clear();
	tissues_.clear();
	layers_.clear();

	std::multimap<std::string, std::list<std::string>> datamap;

	// Extract all configuration data from file
	if (!extract_config_data(filename, datamap)) {
		std::cerr << "Failed to extract data from config file: " << filename << std::endl;
		return false;
	}

	// Parse each configuration section
	for (const auto& [section_name, section_data] : datamap) {
		// Make a copy of the data since some parsing methods may modify it
		std::list<std::string> data = section_data;

		if (equals(section_name, "general")) {
			if (!parse_general_config(data)) {
				std::cerr << "Failed to parse general section." << std::endl;
				return false;
			}
		}
		else if (equals(section_name, "source")) {
			if (!parse_source_config(data)) {
				std::cerr << "Failed to parse source section." << std::endl;
				return false;
			}
		}
		else if (equals(section_name, "tissue")) {
			if (!parse_tissue_config(data)) {
				std::cerr << "Failed to parse tissue section." << std::endl;
				return false;
			}
		}
		else if (equals(section_name, "layer")) {
			if (!parse_layer_config(data)) {
				std::cerr << "Failed to parse layer section." << std::endl;
				return false;
			}
		}
	}

	// Update calculated values
	set_num_layers(static_cast<uint64_t>(layers_.size()));
	set_num_sources(static_cast<uint64_t>(sources_.size()));

	return true;
}

bool Config::parse_general_config(std::list<std::string>& data) {
	// Extract param/value pairs from data
	std::vector<std::pair<std::string, std::string>> param_values = parameter_values(data);

	// Parse each parameter
	for (const auto& [param, value] : param_values) {
		if (equals(param, "numphotons")) {
			num_photons_ = str2num<uint64_t>(value);
		}
		else if (equals(param, "ambienteta")) {
			ambient_eta_ = str2num<double>(value);
		}
		else if (equals(param, "voxelsize")) {
			vox_size_ = str2num<double>(value);
		}
		else if (equals(param, "partial")) {
			partial_ = str2num<bool>(value);
		}
		else if (equals(param, "progress")) {
			progress_ = str2num<bool>(value);
		}
	}

	// Validation
	if (ambient_eta_ < 1.0 || ambient_eta_ > 3.0) {
		std::cerr << "Error: Ambient eta must be [1.0, 3.0], not " << ambient_eta_ << std::endl;
		return false;
	}

	return true;
}

bool Config::parse_source_config(std::list<std::string>& data) {
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
			source.origin = glm::dvec3(p[0], p[1], p[2]);
		}
		else if (equals(param, "direction")) {
			std::vector<double> v = split(value, ',');
			if (v.size() < 3) {
				return false;
			}
			source.direction = glm::dvec3(v[0], v[1], v[2]);
		}
	}

	// check faulty input
	if (source.direction == glm::dvec3(0)) {
		return false;
	}

	// normalize direction
	source.direction = glm::normalize(source.direction);

	// add to collection
	sources_.push_back(source);

	return true;
}

bool Config::parse_tissue_config(std::list<std::string>& data) {
	Tissue tissue;

	// set tissue type properties
	for (const auto& it : parameter_values(data)) {
		std::string param = it.first;
		std::string value = it.second;

		if (equals(param, "id")) {
			tissue.id = str2num<uint8_t>(value);
		}
		else if (equals(param, "ani")) {
			tissue.g = str2num<double>(value);
		}
		else if (equals(param, "eta")) {
			tissue.eta = str2num<double>(value);
		}
		else if (equals(param, "mua")) {
			tissue.mu_a = str2num<double>(value);
		}
		else if (equals(param, "mus")) {
			tissue.mu_s = str2num<double>(value);
		}
	}

	// check faulty input
	if (tissue.g < -1.0 || tissue.g > 1.0) {
		return false;
	}
	if (tissue.eta < 1.0 || tissue.eta > 3.0) {
		return false;
	}
	if (tissue.mu_a < 0) {
		return false;
	}
	if (tissue.mu_s < 0) {
		return false;
	}

	// add to collection
	tissues_.push_back(tissue);

	return true;
}

bool Config::parse_layer_config(std::list<std::string>& data) {
	Layer layer;

	std::vector<glm::uvec3> faces;       // vertex index tuple (faces)
	std::vector<glm::dvec3> vertices;    // triangle vertices
	std::vector<glm::dvec3> normals;     // normal vectors (one per triangle)
	std::vector<Triangle> triangle_mesh; // constructed faces

	// set layer properties
	for (const auto& it : parameter_values(data)) {
		std::string param = it.first;
		std::string value = it.second;

		// identifiers
		if (equals(param, "id")) {
			layer.id = static_cast<uint8_t>(str2num<int>(value));
		}
		else if (equals(param, "tissue")) {
			layer.tissue_id = static_cast<uint8_t>(str2num<int>(value));
		}

		// geometric properties
		// vertex
		else if (equals(param, "vert3")) {
			std::vector<double> v = split(value, ',');
			if (v.size() < 3) {
				continue;
			}

			glm::dvec3 vert = glm::dvec3(v[0], v[1], v[2]);
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
				glm::uvec3(static_cast<uint32_t>(v[0]), static_cast<uint32_t>(v[1]), static_cast<uint32_t>(v[2])));
		}
		// normal
		else if (equals(param, "norm3")) {
			std::vector<double> v = split(value, ',');
			if (v.size() < 3) {
				continue;
			}

			glm::dvec3 norm = glm::normalize(glm::dvec3(v[0], v[1], v[2]));
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
		glm::uvec3 index = faces[i];
		size_t num_verts_size = static_cast<size_t>(num_verts);
		uint32_t index_max = static_cast<uint32_t>(num_verts_size - 1);

		// check vertex index pointer correctness
		if (index.x > index_max || index.y > index_max || index.z > index_max) {
			return false;
		}

		// create new triangle for this face
		Triangle triangle = Triangle(vertices[index.x], vertices[index.y], vertices[index.z]);

		// associate normal with triangle
		triangle.set_normal(normals[i]);

		// add the triangle to the triangle mesh
		triangle_mesh.push_back(triangle);
	}

	// set the layer mesh
	layer.mesh = triangle_mesh;
	
	// Initialize the layer's mesh geometry with triangles
	layer.geometry.set_triangles(triangle_mesh);
	
	// add to collection using move semantics
	layers_.push_back(std::move(layer));

	return true;
}

bool Config::extract_config_data(const std::string& filename,
								 std::multimap<std::string, std::list<std::string>>& datamap) {
	std::string line;
	std::ifstream in_config(filename.c_str(), std::ios_base::in);

	if (!in_config.good()) {
		std::cerr << "File \"" << filename << "\" could not be read or opened." << std::endl;
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
			extract_config_block(in_config, "general", datamap);
		}
		else if (equals(line, "source")) {
			extract_config_block(in_config, "source", datamap);
		}
		else if (equals(line, "tissue")) {
			extract_config_block(in_config, "tissue", datamap);
		}
		else if (equals(line, "layer")) {
			extract_config_block(in_config, "layer", datamap);
		}
	}

	in_config.close();
	return true;
}

void Config::extract_config_block(std::ifstream& in_config, const std::string& section,
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
