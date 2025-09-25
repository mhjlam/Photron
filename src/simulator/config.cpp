#include "config.hpp"

#include <fstream>
#include <iostream>
#include <utility>

#include "math/triangle.hpp"

// Static member definitions
std::unique_ptr<Config> Config::instance_ = nullptr;
bool Config::initialized_ = false;

void Config::initialize() {
	if (!initialized_) {
		instance_ = std::unique_ptr<Config>(new Config());
		initialized_ = true;
	}
}

void Config::initialize(const std::string& config_file) {
	if (!initialized_) {
		instance_ = std::unique_ptr<Config>(new Config());
		if (!instance_->parse_config_file(config_file)) {
			// Keep the instance but log the error
			std::cerr << "Warning: Failed to parse config file, using defaults" << std::endl;
		}
		initialized_ = true;
	}
}

bool Config::is_initialized() {
	return initialized_;
}

Config& Config::get() {
	if (!initialized_) {
		throw std::runtime_error("Config not initialized. Call Config::initialize() first.");
	}
	return *instance_;
}

void Config::shutdown() {
	instance_.reset();
	initialized_ = false;
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
		else if (equals(section_name, "material") || equals(section_name, "tissue")) {
			std::cerr << "Error: Standalone material/tissue blocks are no longer supported. "
			          << "Define material properties directly within layer blocks." << std::endl;
			return false;
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
		if (equals(param, "photons")) {
			num_photons_ = str2num<uint64_t>(value);
		}
		else if (equals(param, "voxel_size")) {
			vox_size_ = str2num<double>(value);
		}
		else if (equals(param, "verbose")) {
			verbose_ = str2num<bool>(value);
		}
		else if (equals(param, "deterministic")) {
			deterministic_ = str2num<bool>(value);
		}
	}

	return true;
}

bool Config::parse_source_config(std::list<std::string>& data) {
	Source source;
	
	// Automatically assign ID based on order of definition
	source.id = static_cast<uint64_t>(sources_.size());

	// set source properties
	for (const auto& it : parameter_values(data)) {
		std::string param = it.first;
		std::string value = it.second;

		if (equals(param, "id")) {
			std::cerr << "Error: Explicit source IDs are no longer supported. IDs are assigned automatically based on order of definition." << std::endl;
			return false;
		}
		else if (equals(param, "position")) {
			std::vector<double> p = split(value, ',');
			if (p.size() < 3) {
				std::cerr << "Error: Source position '" << value 
				          << "' must have 3 components (x,y,z)" << std::endl;
				return false;
			}
			source.origin = glm::dvec3(p[0], p[1], p[2]);
		}
		else if (equals(param, "direction")) {
			std::vector<double> v = split(value, ',');
			if (v.size() < 3) {
				std::cerr << "Error: Source direction '" << value 
				          << "' must have 3 components (x,y,z)" << std::endl;
				return false;
			}
			source.direction = glm::dvec3(v[0], v[1], v[2]);
		}
	}

	// check faulty input
	if (source.direction == glm::dvec3(0)) {
		std::cerr << "Error: Source " << source.id 
		          << " has zero direction vector. Direction cannot be (0,0,0)" << std::endl;
		return false;
	}

	// normalize direction
	source.direction = glm::normalize(source.direction);

	// add to collection
	sources_.push_back(source);

	return true;
}

bool Config::parse_layer_config(std::list<std::string>& data) {
	Layer layer;
	
	// Automatically assign ID based on order of definition
	layer.id = static_cast<uint8_t>(tissues_.size());

	std::vector<glm::uvec3> faces;       // vertex index tuple (faces)
	std::vector<glm::dvec3> vertices;    // triangle vertices
	std::vector<glm::dvec3> normals;     // normal vectors (one per triangle)
	std::vector<Triangle> triangle_mesh; // constructed faces

	// Material properties (required inline definition)
	bool has_eta = false, has_mua = false, has_mus = false, has_ani = false;
	double eta = 1.37;
	double mua = 1.0;
	double mus = 10.0;
	double ani = 0.0;

	// set layer properties
	for (const auto& it : parameter_values(data)) {
		std::string param = it.first;
		std::string value = it.second;

		// identifiers
		if (equals(param, "material") || equals(param, "tissue")) {
			std::cerr << "Error: Material references are no longer supported in layer " << (int)layer.id
			          << ". Define material properties (eta, mua, mus, ani) directly within the layer block." << std::endl;
			return false;
		}

		// Material properties (required inline definition)
		else if (equals(param, "eta")) {
			eta = str2num<double>(value);
			has_eta = true;
		}
		else if (equals(param, "mua")) {
			mua = str2num<double>(value);
			has_mua = true;
		}
		else if (equals(param, "mus")) {
			mus = str2num<double>(value);
			has_mus = true;
		}
		else if (equals(param, "ani") || equals(param, "g")) {
			ani = str2num<double>(value);
			has_ani = true;
		}

		// geometric properties
		// vertex
		else if (equals(param, "vertex")) {
			std::vector<double> v = split(value, ',');
			if (v.size() < 3) {
				continue;
			}

			glm::dvec3 vert = glm::dvec3(v[0], v[1], v[2]);
			vertices.push_back(vert);
		}
		// triangle face
		else if (equals(param, "face")) {
			std::vector<double> v = split(value, ',');

			// invalid face detected
			if (v.size() < 3) {
				std::cerr << "Error: Layer " << (int)layer.id 
				          << " face definition '" << value 
				          << "' must have 3 vertex indices" << std::endl;
				return false;
			}

			faces.push_back(
				glm::uvec3(static_cast<uint32_t>(v[0]), static_cast<uint32_t>(v[1]), static_cast<uint32_t>(v[2])));
		}
		// normal
		else if (equals(param, "normal")) {
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
		std::cerr << "Error: Layer " << (int)layer.id << " has only " << num_faces 
		          << " faces. Need at least 4 faces to form a polyhedron." << std::endl;
		return false;
	}
	// so, also 4 normals
	if (normals.size() < 4) {
		std::cerr << "Error: Layer " << (int)layer.id << " has only " << normals.size() 
		          << " normals. Need at least 4 normals (one per face)." << std::endl;
		return false;
	}
	// and at least 4 vertices (tetrahedron)
	if (num_verts < 4) {
		std::cerr << "Error: Layer " << (int)layer.id << " has only " << num_verts 
		          << " vertices. Need at least 4 vertices to form a polyhedron." << std::endl;
		return false;
	}

	// polehedron must be convex (Euler characteristic must be 2)
	if ((num_verts - num_edges + num_faces) != 2) {
		std::cerr << "Error: Layer " << (int)layer.id << " geometry fails Euler characteristic test. "
		          << "V=" << num_verts << ", E=" << num_edges << ", F=" << num_faces 
		          << " (V-E+F=" << (num_verts - num_edges + num_faces) << ", should be 2). "
		          << "Polyhedron must be convex." << std::endl;
		return false;
	}

	// each triangle/face must have a normal vector
	if (faces.size() != normals.size()) {
		std::cerr << "Error: Layer " << (int)layer.id << " has " << faces.size() 
		          << " faces but " << normals.size() << " normals. "
		          << "Each face must have exactly one normal vector." << std::endl;
		return false;
	}

	// construct faces (triangles) from vertex faces
	size_t num_verts_size = static_cast<size_t>(num_verts);
	uint32_t index_max = static_cast<uint32_t>(num_verts_size - 1);

	for (size_t i = 0; i < faces.size(); ++i) {
		const auto& face_indices = faces[i];
		const auto& face_normal = normals[i];
		
		// check vertex index pointer correctness
		if (face_indices.x > index_max || face_indices.y > index_max || face_indices.z > index_max) {
			std::cerr << "Error: Layer " << (int)layer.id << " face " << i 
			          << " has invalid vertex indices (" << face_indices.x << ", " 
			          << face_indices.y << ", " << face_indices.z << "). "
			          << "Valid range is 0-" << index_max << " (total vertices: " 
			          << num_verts << ")" << std::endl;
			return false;
		}

		// create new triangle for this face
		Triangle triangle = Triangle(vertices[face_indices.x], vertices[face_indices.y], vertices[face_indices.z]);

		// Validate config normal against computed normal
		glm::dvec3 computed_normal = triangle.normal();
		double dot_product = glm::dot(glm::normalize(face_normal), glm::normalize(computed_normal));
		double angle_diff = acos(std::clamp(std::abs(dot_product), 0.0, 1.0)) * 180.0 / 3.14159265359;
		
		if (angle_diff > 5.0) { // Warn if normals differ by more than 5 degrees
			std::cerr << "Warning: Config normal differs from computed normal by " << angle_diff 
					  << " degrees for face " << i << std::endl;
			std::cerr << "  Config normal: (" << face_normal.x << ", " << face_normal.y << ", " << face_normal.z << ")" << std::endl;
			std::cerr << "  Computed normal: (" << computed_normal.x << ", " << computed_normal.y << ", " << computed_normal.z << ")" << std::endl;
		}

		// add the triangle to the triangle mesh
		triangle_mesh.push_back(triangle);
	}

	// set the layer mesh
	layer.mesh = triangle_mesh;
	
	// Initialize the layer's mesh geometry with triangles
	layer.set_triangles(triangle_mesh);
	
	// Validate and fix normal orientations to ensure consistency
	layer.validate_and_fix_normals();
	
	// Validate that all required material properties are provided
	if (!has_eta || !has_mua || !has_mus || !has_ani) {
		std::cerr << "Error: Layer " << (int)layer.id << " is missing required material properties. ";
		std::cerr << "All layers must define: eta, mua, mus, ani (or g)" << std::endl;
		std::cerr << "Missing: ";
		if (!has_eta) std::cerr << "eta ";
		if (!has_mua) std::cerr << "mua ";
		if (!has_mus) std::cerr << "mus ";
		if (!has_ani) std::cerr << "ani ";
		std::cerr << std::endl;
		return false;
	}

	// Validate material property ranges
	if (ani < -1.0 || ani > 1.0) {
		std::cerr << "Error: Invalid anisotropy factor ani=" << ani 
		          << " for layer " << (int)layer.id
		          << ". Must be between -1.0 and 1.0" << std::endl;
		return false;
	}
	if (eta < 1.0 || eta > 3.0) {
		std::cerr << "Error: Invalid refractive index eta=" << eta 
		          << " for layer " << (int)layer.id
		          << ". Must be between 1.0 and 3.0" << std::endl;
		return false;
	}
	if (mua < 0) {
		std::cerr << "Error: Invalid absorption coefficient mua=" << mua 
		          << " for layer " << (int)layer.id
		          << ". Must be >= 0" << std::endl;
		return false;
	}
	if (mus < 0) {
		std::cerr << "Error: Invalid scattering coefficient mus=" << mus 
		          << " for layer " << (int)layer.id
		          << ". Must be >= 0" << std::endl;
		return false;
	}
	
	// Create material with the inline properties
	uint8_t material_id = static_cast<uint8_t>(tissues_.size());
	Material inline_material(material_id, ani, eta, mua, mus);
	tissues_.push_back(std::move(inline_material));
	
	// Update layer to reference the new material
	layer.tissue_id = material_id;
	
	if (verbose_) {
		std::cout << "Created material " << (int)material_id 
		          << " for layer " << (int)layer.id
		          << " (eta=" << eta << ", mua=" << mua << ", mus=" << mus << ", ani=" << ani << ")" << std::endl;
	}
	
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
		else if (equals(line, "material")) {
			extract_config_block(in_config, "material", datamap);
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
