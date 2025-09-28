#include "config.hpp"

#include <fstream>
#include <iostream>
#include <limits>
#include <sstream>
#include <utility>
#include <filesystem>

#include <toml++/toml.h>
#include "common/error_types.hpp"

#include "logger.hpp"
#include "common/result.hpp"
#include "common/error_handler.hpp"
#include "common/error_types.hpp"
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

bool Config::initialize(const std::string& config_file) {
	if (!initialized_) {
		instance_ = std::unique_ptr<Config>(new Config());
		auto result = instance_->parse_config_file(config_file);
		if (!result.is_ok()) {
			ErrorHandler::instance().report_error(ErrorMessage::format(result.error(), "Config initialization failed"));
			// Reset instance since parsing failed
			instance_.reset();
			return false;
		}
		initialized_ = true;
	}
	return true;
}

bool Config::is_initialized() {
	return initialized_;
}

Config& Config::get() {
	if (!initialized_ || !instance_) {
		// Instead of throwing, create a temporary fallback instance
		// This prevents crashes during error recovery
		
		// Create a minimal fallback config to prevent crashes
		static std::unique_ptr<Config> fallback_instance = nullptr;
		if (!fallback_instance) {
			fallback_instance = std::unique_ptr<Config>(new Config());
			// Initialize with minimal default values to prevent further crashes
			// Set some safe defaults to prevent issues
			fallback_instance->config_filename_ = "No Configuration Loaded";
		}
		return *fallback_instance;
	}
	return *instance_;
}

void Config::shutdown() {
	instance_.reset();
	initialized_ = false;
}

// Structured error handling version
Result<void, ConfigError> Config::parse_config_file(const std::string& filename) {
	// Check if file exists
	if (!std::filesystem::exists(filename)) {
		return Result<void, ConfigError>::error(ConfigError::FileNotFound);
	}

	try {
		// Store the config filename
		config_filename_ = filename;
		
		// Parse TOML file
		toml::table config = toml::parse_file(filename);
		
		// Parse each configuration section
		if (!parse_general_config(config)) {
			return Result<void, ConfigError>::error(ConfigError::ValidationError);
		}
		if (!parse_source_config(config)) {
			return Result<void, ConfigError>::error(ConfigError::ValidationError);
		}
		if (!parse_layer_configs(config)) {
			return Result<void, ConfigError>::error(ConfigError::GeometryError);
		}
		
		return Result<void, ConfigError>::ok();
		
	}
	catch (const toml::parse_error& err) {
		// Store error details for higher-level handling
		return Result<void, ConfigError>::error(ConfigError::ParseError);
	}
	catch (const std::exception& err) {
		// Store error details for higher-level handling
		return Result<void, ConfigError>::error(ConfigError::ParseError);
	}
}

bool Config::parse_general_config(const toml::table& config) {
	auto general = config["general"];
	if (!general) {
		REPORT_CONFIG_ERROR("No [general] section found in config");
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
		REPORT_CONFIG_ERROR("No [source] section found in config");
		return false;
	}

	Source src;
	src.id = static_cast<uint64_t>(sources_.size());

	// Parse position
	if (auto position_arr = source_section["position"].as_array()) {
		if (!validate_array_size(*position_arr, 3, "Source position")) {
			return false;
		}
		src.origin = glm::dvec3(
			position_arr->at(0).value_or<double>(0.0),
			position_arr->at(1).value_or<double>(0.0),
			position_arr->at(2).value_or<double>(0.0)
		);
	} else {
		ErrorHandler::instance().report_error("Source position must be an array of 3 numbers");
		return false;
	}

	// Parse direction
	if (auto direction_arr = source_section["direction"].as_array()) {
		if (!validate_array_size(*direction_arr, 3, "Source direction")) {
			return false;
		}
		src.direction = glm::dvec3(
			direction_arr->at(0).value_or<double>(0.0),
			direction_arr->at(1).value_or<double>(0.0),
			direction_arr->at(2).value_or<double>(0.0)
		);
	} else {
		ErrorHandler::instance().report_error("Source direction must be an array of 3 numbers");
		return false;
	}

	// Validate direction is not zero
	if (src.direction == glm::dvec3(0)) {
		ErrorHandler::instance().report_error("Source direction cannot be zero vector");
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
		ErrorHandler::instance().report_error("No [[layer]] sections found in config");
		return false;
	}

	for (const auto& layer_node : *layers_array) {
		auto layer_table = layer_node.as_table();
		if (!layer_table) {
			ErrorHandler::instance().report_error("Invalid layer configuration");
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
			ErrorHandler::instance().report_error("Layer " + std::to_string((int)layer.id) + " is missing required material properties. All layers must define: eta, mua, mus, ani");
			return false;
		}

		// Validate material property ranges
		if (!validate_range(ani, -1.0, 1.0, "anisotropy factor (ani)", layer.id)) {
			return false;
		}
		if (!validate_range(eta, 1.0, 3.0, "refractive index (eta)", layer.id)) {
			return false;
		}
		if (!validate_range(mua, 0.0, std::numeric_limits<double>::max(), "absorption coefficient (mua)", layer.id) ||
		    !validate_range(mus, 0.0, std::numeric_limits<double>::max(), "scattering coefficient (mus)", layer.id)) {
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
			ErrorHandler::instance().report_error("Layer " + std::to_string((int)layer.id) + " missing vertices array");
			return false;
		}

		// Parse faces
		if (auto faces_arr = (*layer_table)["faces"].as_array()) {
			faces = parse_faces(*faces_arr);
		} else {
			ErrorHandler::instance().report_error("Layer " + std::to_string((int)layer.id) + " missing faces array");
			return false;
		}

		// Validate geometry
		if (!validate_geometry_size(vertices.size(), 4, "vertices", layer.id) ||
		    !validate_geometry_size(faces.size(), 4, "faces", layer.id)) {
			return false;
		}

		// Build triangle mesh
		for (const auto& face : faces) {
			if (face.x >= vertices.size() || face.y >= vertices.size() || face.z >= vertices.size()) {
				ErrorHandler::instance().report_error("Layer " + std::to_string((int)layer.id) + " has face with invalid vertex index");
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
		materials_.push_back(std::move(material));

		// Assign material to layer (use index as tissue_id)
		layer.tissue_id = static_cast<uint8_t>(materials_.size() - 1);

		if (log_) {
			std::ostringstream debug_msg;
			debug_msg << "Created material " << (int)layer.tissue_id 
			          << " for layer " << (int)layer.id
			          << " (eta=" << eta << ", mua=" << mua << ", mus=" << mus << ", ani=" << ani << ")";
			Logger::instance().log_info(debug_msg.str());
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

// Validation helper method implementations
template<typename T>
bool Config::validate_range(T value, T min, T max, const std::string& param_name, int layer_id) const {
	if (value < min || value > max) {
		std::ostringstream error_msg;
		error_msg << "Invalid " << param_name << "=" << value;
		if (layer_id >= 0) {
			error_msg << " for layer " << layer_id;
		}
		if (max == std::numeric_limits<T>::max()) {
			error_msg << ". Must be >= " << min;
		} else {
			error_msg << ". Must be between " << min << " and " << max;
		}
		ErrorHandler::instance().report_error(ErrorMessage::format(ConfigError::ValidationError, error_msg.str()));
		return false;
	}
	return true;
}

bool Config::validate_array_size(const toml::array& arr, size_t expected_size, const std::string& param_name) const {
	if (arr.size() != expected_size) {
		ErrorHandler::instance().report_error(param_name + " must be an array of " + std::to_string(expected_size) + " numbers");
		return false;
	}
	return true;
}

bool Config::validate_geometry_size(size_t actual_size, size_t min_size, const std::string& type, int layer_id) const {
	if (actual_size < min_size) {
		std::ostringstream error_msg;
		error_msg << "Layer " << layer_id << " must have at least " << min_size << " " << type << " for a valid polyhedron";
		ErrorHandler::instance().report_error(ErrorMessage::format(ConfigError::GeometryError, error_msg.str()));
		return false;
	}
	return true;
}

// Structured error validation method implementations
template<typename T>
Result<void, ConfigError> Config::validate_range_structured(T value, T min, T max, const std::string& param_name, int layer_id) const {
	if (value < min || value > max) {
		std::ostringstream error_msg;
		error_msg << "Validation failed for parameter '" << param_name << "'";
		if (layer_id >= 0) error_msg << " in layer " << layer_id;
		error_msg << ": value " << value << " not in range [" << min << ", " << max << "]";
		ErrorHandler::instance().report_error(ErrorMessage::format(ConfigError::InvalidValue, error_msg.str()));
		return Result<void, ConfigError>::error(ConfigError::InvalidValue);
	}
	return Result<void, ConfigError>::ok();
}

Result<void, ConfigError> Config::validate_array_size_structured(const toml::array& arr, size_t expected_size, const std::string& param_name) const {
	if (arr.size() != expected_size) {
		std::cerr << "Array validation failed for parameter '" << param_name 
			  << "': array has " << arr.size() << " elements, expected " << expected_size << std::endl;
		return Result<void, ConfigError>::error(ConfigError::ValidationError);
	}
	return Result<void, ConfigError>::ok();
}

Result<void, ConfigError> Config::validate_geometry_size_structured(size_t actual_size, size_t min_size, const std::string& type, int layer_id) const {
	if (actual_size < min_size) {
		std::ostringstream error_msg;
		error_msg << "Geometry validation failed for " << type << " geometry";
		if (layer_id >= 0) error_msg << " in layer " << layer_id;
		error_msg << ": has " << actual_size << " elements, minimum required is " << min_size;
		ErrorHandler::instance().report_error(ErrorMessage::format(ConfigError::GeometryError, error_msg.str()));
		return Result<void, ConfigError>::error(ConfigError::GeometryError);
	}
	return Result<void, ConfigError>::ok();
}

// Explicit template instantiations
template bool Config::validate_range<double>(double value, double min, double max, const std::string& param_name, int layer_id) const;

// Structured error template instantiations  
template Result<void, ConfigError> Config::validate_range_structured<double>(double, double, double, const std::string&, int) const;
template Result<void, ConfigError> Config::validate_range_structured<int>(int, int, int, const std::string&, int) const;
template Result<void, ConfigError> Config::validate_range_structured<size_t>(size_t, size_t, size_t, const std::string&, int) const;
