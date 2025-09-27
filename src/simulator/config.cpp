#include "config.hpp"

#include <fstream>
#include <iostream>
#include <sstream>
#include <utility>

#include <toml++/toml.h>

#include "math/triangle.hpp"
#include "debug_logger.hpp"

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
	try {
		// Parse TOML file
		auto config = toml::parse_file(filename);
		
		// Clear existing data
		sources_.clear();
		tissues_.clear();
		layers_.clear();

		// Parse sections
		if (!parse_general_config(config)) {
			std::cerr << "Failed to parse general configuration" << std::endl;
			return false;
		}

		if (!parse_source_config(config)) {
			std::cerr << "Failed to parse source configuration" << std::endl;
			return false;
		}

		if (!parse_layer_configs(config)) {
			std::cerr << "Failed to parse layer configurations" << std::endl;
			return false;
		}

		// Update calculated values
		set_num_layers(static_cast<uint64_t>(layers_.size()));
		set_num_sources(static_cast<uint64_t>(sources_.size()));

		return true;
	}
	catch (const toml::parse_error& err) {
		std::cerr << "TOML parsing error in " << filename << ": " << err.what() << std::endl;
		return false;
	}
	catch (const std::exception& e) {
		std::cerr << "Error parsing config file " << filename << ": " << e.what() << std::endl;
		return false;
	}
}

bool Config::parse_general_config(const toml::table& config) {
	auto general = config["general"];
	if (!general) {
		std::cerr << "No [general] section found in config" << std::endl;
		return false;
	}

	// Parse general parameters
	if (auto photons = general["photons"].value<uint64_t>()) {
		num_photons_ = *photons;
	}
	
	if (auto voxel_size = general["voxel_size"].value<double>()) {
		vox_size_ = *voxel_size;
	}
	
	if (auto log = general["log"].value<bool>()) {
		log_ = *log;
	}
	
	if (auto deterministic = general["deterministic"].value<bool>()) {
		deterministic_ = *deterministic;
	}

	return true;
}

bool Config::parse_source_config(const toml::table& config) {
	auto source_section = config["source"];
	if (!source_section) {
		std::cerr << "No [source] section found in config" << std::endl;
		return false;
	}

	Source src;
	src.id = static_cast<uint64_t>(sources_.size());

	// Parse position
	if (auto position_arr = source_section["position"].as_array()) {
		if (position_arr->size() != 3) {
			std::cerr << "Source position must be an array of 3 numbers" << std::endl;
			return false;
		}
		src.origin = glm::dvec3(
			position_arr->at(0).value_or<double>(0.0),
			position_arr->at(1).value_or<double>(0.0),
			position_arr->at(2).value_or<double>(0.0)
		);
	} else {
		std::cerr << "Source position must be an array of 3 numbers" << std::endl;
		return false;
	}

	// Parse direction
	if (auto direction_arr = source_section["direction"].as_array()) {
		if (direction_arr->size() != 3) {
			std::cerr << "Source direction must be an array of 3 numbers" << std::endl;
			return false;
		}
		src.direction = glm::dvec3(
			direction_arr->at(0).value_or<double>(0.0),
			direction_arr->at(1).value_or<double>(0.0),
			direction_arr->at(2).value_or<double>(0.0)
		);
	} else {
		std::cerr << "Source direction must be an array of 3 numbers" << std::endl;
		return false;
	}

	// Validate direction is not zero
	if (src.direction == glm::dvec3(0)) {
		std::cerr << "Source direction cannot be zero vector" << std::endl;
		return false;
	}

	// normalize direction
	src.direction = glm::normalize(src.direction);

	sources_.push_back(src);
	return true;
}

bool Config::parse_layer_configs(const toml::table& config) {
	// TOML array of tables: [[layer]]
	auto layers_array = config["layer"].as_array();
	if (!layers_array) {
		std::cerr << "No [[layer]] sections found in config" << std::endl;
		return false;
	}

	for (const auto& layer_node : *layers_array) {
		auto layer_table = layer_node.as_table();
		if (!layer_table) {
			std::cerr << "Invalid layer configuration" << std::endl;
			continue;
		}

		Layer layer;
		layer.id = static_cast<uint8_t>(layers_.size());

		// Parse material properties
		bool has_eta = false, has_mua = false, has_mus = false, has_ani = false;
		double eta = 1.37, mua = 1.0, mus = 10.0, ani = 0.0;

		if (auto val = (*layer_table)["eta"].value<double>()) {
			eta = *val;
			has_eta = true;
		}
		if (auto val = (*layer_table)["mua"].value<double>()) {
			mua = *val;
			has_mua = true;
		}
		if (auto val = (*layer_table)["mus"].value<double>()) {
			mus = *val;
			has_mus = true;
		}
		if (auto val = (*layer_table)["ani"].value<double>()) {
			ani = *val;
			has_ani = true;
		}

		// Validate required material properties
		if (!has_eta || !has_mua || !has_mus || !has_ani) {
			std::cerr << "Error: Layer " << (int)layer.id << " is missing required material properties. ";
			std::cerr << "All layers must define: eta, mua, mus, ani" << std::endl;
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
		if (mua < 0 || mus < 0) {
			std::cerr << "Error: Invalid absorption/scattering coefficients (mua=" << mua 
			          << ", mus=" << mus << ") for layer " << (int)layer.id
			          << ". Must be >= 0" << std::endl;
			return false;
		}

		// Parse geometry
		std::vector<glm::dvec3> vertices;
		std::vector<glm::uvec3> faces;
		std::vector<Triangle> triangle_mesh;

		// Parse vertices
		if (auto vertices_arr = (*layer_table)["vertices"].as_array()) {
			vertices = parse_vertices(*vertices_arr);
		} else {
			std::cerr << "Layer " << (int)layer.id << " missing vertices array" << std::endl;
			return false;
		}

		// Parse faces
		if (auto faces_arr = (*layer_table)["faces"].as_array()) {
			faces = parse_faces(*faces_arr);
		} else {
			std::cerr << "Layer " << (int)layer.id << " missing faces array" << std::endl;
			return false;
		}

		// Validate geometry
		if (vertices.size() < 4) {
			std::cerr << "Error: Layer " << (int)layer.id 
			          << " must have at least 4 vertices for a valid polyhedron" << std::endl;
			return false;
		}
		if (faces.size() < 4) {
			std::cerr << "Error: Layer " << (int)layer.id 
			          << " must have at least 4 faces for a valid polyhedron" << std::endl;
			return false;
		}

		// Build triangle mesh
		for (const auto& face : faces) {
			if (face.x >= vertices.size() || face.y >= vertices.size() || face.z >= vertices.size()) {
				std::cerr << "Error: Layer " << (int)layer.id << " has face with invalid vertex index" << std::endl;
				return false;
			}

			Triangle triangle(vertices[face.x], vertices[face.y], vertices[face.z]);
			triangle_mesh.push_back(triangle);
		}

		// Set the layer mesh
		layer.mesh = triangle_mesh;
		layer.set_triangles(triangle_mesh);
		layer.validate_and_fix_normals();

		// Create material
		Material material(ani, eta, mua, mus);
		tissues_.push_back(std::move(material));

		// Assign material to layer (use index as tissue_id)
		layer.tissue_id = static_cast<uint8_t>(tissues_.size() - 1);

		if (log_) {
			std::ostringstream debug_msg;
			debug_msg << "Created material " << (int)layer.tissue_id 
			          << " for layer " << (int)layer.id
			          << " (eta=" << eta << ", mua=" << mua << ", mus=" << mus << ", ani=" << ani << ")";
			DebugLogger::instance().log_info(debug_msg.str());
		}

		layers_.push_back(std::move(layer));
	}

	return true;
}

// Helper methods for parsing TOML arrays
glm::dvec3 Config::parse_vec3(const toml::array& arr) const {
	if (arr.size() != 3) {
		throw std::runtime_error("Vec3 array must have exactly 3 elements");
	}
	return glm::dvec3(
		arr[0].value_or<double>(0.0),
		arr[1].value_or<double>(0.0),
		arr[2].value_or<double>(0.0)
	);
}

glm::uvec3 Config::parse_uvec3(const toml::array& arr) const {
	if (arr.size() != 3) {
		throw std::runtime_error("UVec3 array must have exactly 3 elements");
	}
	return glm::uvec3(
		arr[0].value_or<uint32_t>(0),
		arr[1].value_or<uint32_t>(0),
		arr[2].value_or<uint32_t>(0)
	);
}

std::vector<glm::dvec3> Config::parse_vertices(const toml::array& vertices_array) const {
	std::vector<glm::dvec3> vertices;
	vertices.reserve(vertices_array.size());

	for (const auto& vertex_node : vertices_array) {
		if (auto vertex_arr = vertex_node.as_array()) {
			vertices.push_back(parse_vec3(*vertex_arr));
		} else {
			throw std::runtime_error("Each vertex must be an array of 3 numbers");
		}
	}

	return vertices;
}

std::vector<glm::uvec3> Config::parse_faces(const toml::array& faces_array) const {
	std::vector<glm::uvec3> faces;
	faces.reserve(faces_array.size());

	for (const auto& face_node : faces_array) {
		if (auto face_arr = face_node.as_array()) {
			faces.push_back(parse_uvec3(*face_arr));
		} else {
			throw std::runtime_error("Each face must be an array of 3 vertex indices");
		}
	}

	return faces;
}

std::vector<glm::dvec3> Config::parse_normals(const toml::array& normals_array) const {
	std::vector<glm::dvec3> normals;
	normals.reserve(normals_array.size());

	for (const auto& normal_node : normals_array) {
		if (auto normal_arr = normal_node.as_array()) {
			normals.push_back(glm::normalize(parse_vec3(*normal_arr)));
		} else {
			throw std::runtime_error("Each normal must be an array of 3 numbers");
		}
	}

	return normals;
}
