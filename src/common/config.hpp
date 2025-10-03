/**
 * @file config.hpp
 * @brief Global configuration management system
 *
 * The Config class provides centralized configuration management for the entire
 * Photron application using TOML format files. It handles simulation parameters,
 * geometry definitions, material properties, and runtime options with validation
 * and default value support.
 */

#pragma once

#include <cstdint>
#include <memory>
#include <span>
#include <stdexcept>
#include <string>
#include <vector>

#include <glm/glm.hpp>
#include <toml++/toml.h>

#include "math/concepts.hpp"

// Forward declarations
class Layer;
struct Source;
class Material;
template<typename T, typename E> class Result;
enum class ConfigError;

/**
 * @class Config
 * @brief Singleton configuration management system with TOML support
 *
 * The Config class serves as the single source of truth for all application
 * configuration data. It provides:
 *
 * **TOML Configuration Loading:**
 * - Hierarchical configuration structure
 * - Type-safe parameter parsing
 * - Comprehensive validation and error reporting
 * - Default value fallback mechanisms
 *
 * **Parameter Categories:**
 * - **Simulation**: Photon count, random seeds, algorithm parameters
 * - **Geometry**: Voxel dimensions, layer definitions, material properties
 * - **Sources**: Light source configurations and emission parameters
 * - **Output**: File paths, logging levels, result formats
 *
 * **Key Features:**
 * - Thread-safe singleton access pattern
 * - Immutable configuration after initialization
 * - Structured error reporting for invalid configurations
 * - Memory-efficient storage with move semantics
 *
 * **Configuration File Structure:**
 * ```toml
 * [general]
 * photons = 1000000
 * voxel_size = 0.01
 *
 * [source]
 * position = [0.0, 0.0, 1.0]
 * direction = [0.0, 0.0, -1.0]
 *
 * [[material]]
 * name = "skin"
 * mua = 0.1
 * mus = 10.0
 * ```
 */
class Config
{
private:
	/// Singleton instance pointer (managed automatically)
	static std::unique_ptr<Config> instance_;
	/// Initialization state flag to prevent double-initialization
	static bool initialized_;

	// Voxel grid parameters
	uint32_t nx_ {0};       ///< Number of voxels in X direction
	uint32_t ny_ {0};       ///< Number of voxels in Y direction
	uint32_t nz_ {0};       ///< Number of voxels in Z direction

	double vox_size_ {0.0}; ///< Uniform voxel edge length (cubic voxels: dx=dy=dz)

	// Simulation parameters
	uint64_t num_layers_ {0};  ///< Number of material layers in simulation
	uint64_t num_voxels_ {0};  ///< Total number of voxels (nx * ny * nz)
	uint64_t num_photons_ {0}; ///< Number of Monte Carlo photons to simulate
	uint64_t num_sources_ {0}; ///< Number of configured light sources

	// Runtime control flags
	bool log_ {false};           ///< Enable logging and progress messages
	bool deterministic_ {false}; ///< Use fixed random seed for reproducible results

	// Configuration metadata
	std::string config_filename_; ///< Path to loaded configuration file (for reference)

	// Parsed configuration objects
	std::vector<Source> sources_;     ///< Light source configurations
	std::vector<Material> materials_; ///< Material property definitions
	std::vector<Layer> layers_;

	// TOML parsing helper methods
	glm::dvec3 parse_vec3(const toml::array& arr) const;
	glm::uvec3 parse_uvec3(const toml::array& arr) const;
	std::vector<glm::dvec3> parse_vertices(const toml::array& vertices_array) const;
	std::vector<glm::uvec3> parse_faces(const toml::array& faces_array) const;
	std::vector<glm::dvec3> parse_normals(const toml::array& normals_array) const;

	// Validation helper methods
	// Enhanced validation with structured errors
	template<typename T>
	Result<void, ConfigError> validate_range_structured(
		T value, T min, T max, const std::string& param_name, int layer_id = -1) const;
	Result<void, ConfigError> validate_array_size_structured(const toml::array& arr,
															 size_t expected_size,
															 const std::string& param_name) const;
	Result<void, ConfigError> validate_geometry_size_structured(size_t actual_size,
																size_t min_size,
																const std::string& type,
																int layer_id) const;

	// Legacy validation methods (for backward compatibility)
	template<typename T>
	bool validate_range(T value, T min, T max, const std::string& param_name, int layer_id = -1) const;
	bool validate_array_size(const toml::array& arr, size_t expected_size, const std::string& param_name) const;
	bool validate_geometry_size(size_t actual_size, size_t min_size, const std::string& type, int layer_id) const;

public:
	// Delete copy constructor and assignment operator for singleton
	Config(const Config&) = delete;
	Config& operator=(const Config&) = delete;

	// Static methods for singleton access
	/**
	 * @brief Initialize the global Config instance
	 */
	static void initialize();

	/**
	 * @brief Initialize the global Config instance with a config file
	 * @param config_file Path to the configuration file
	 * @return true if initialization and parsing succeeded, false otherwise
	 */
	static bool initialize(const std::string& config_file);

	/**
	 * @brief Check if config service has been initialized
	 * @return true if initialized, false otherwise
	 */
	static bool is_initialized();

	/**
	 * @brief Get the current config instance
	 * @return Reference to the current config
	 * @throws std::runtime_error if config has not been initialized
	 */
	static Config& get();

	/**
	 * @brief Reset the config service (for testing or shutdown)
	 */
	static void shutdown();

	// Getters with [[nodiscard]] for safety
	[[nodiscard]] constexpr uint32_t nx() const noexcept { return nx_; }
	[[nodiscard]] constexpr uint32_t ny() const noexcept { return ny_; }
	[[nodiscard]] constexpr uint32_t nz() const noexcept { return nz_; }
	[[nodiscard]] constexpr uint64_t num_layers() const noexcept { return num_layers_; }
	[[nodiscard]] constexpr uint64_t num_voxels() const noexcept { return num_voxels_; }
	[[nodiscard]] constexpr uint64_t num_photons() const noexcept { return num_photons_; }
	[[nodiscard]] constexpr uint64_t num_sources() const noexcept { return num_sources_; }
	[[nodiscard]] constexpr double vox_size() const noexcept { return vox_size_; }
	[[nodiscard]] constexpr double ambient_eta() const noexcept { return 1.0; } // Always air
	[[nodiscard]] constexpr bool partial() const noexcept { return true; }      // Always enabled
	[[nodiscard]] constexpr bool log() const noexcept { return log_; }
	[[nodiscard]] constexpr bool deterministic() const noexcept { return deterministic_; }
	[[nodiscard]] const std::string& config_filename() const noexcept { return config_filename_; }

	// Compatibility accessors returning references
	[[nodiscard]] const std::vector<Source>& sources() const noexcept { return sources_; }
	[[nodiscard]] const std::vector<Material>& materials() const noexcept { return materials_; }
	[[nodiscard]] const std::vector<Layer>& layers() const noexcept { return layers_; }

	// Move versions for transferring ownership (needed for Layer which is move-only)
	[[nodiscard]] std::vector<Layer>&& move_layers() { return std::move(layers_); }
	[[nodiscard]] std::vector<Material>&& move_materials() { return std::move(materials_); }

	// Direct access to vectors is preferred, but move semantics needed for Layer transfer

	// Setters
	void set_nx(uint32_t nx) noexcept { nx_ = nx; }
	void set_ny(uint32_t ny) noexcept { ny_ = ny; }
	void set_nz(uint32_t nz) noexcept { nz_ = nz; }
	void set_num_layers(uint64_t num_layers) noexcept { num_layers_ = num_layers; }
	void set_num_voxels(uint64_t num_voxels) noexcept { num_voxels_ = num_voxels; }
	void set_num_photons(uint64_t num_photons) noexcept { num_photons_ = num_photons; }
	void set_num_sources(uint64_t num_sources) noexcept { num_sources_ = num_sources; }
	void set_vox_size(double vox_size) noexcept { vox_size_ = vox_size; }
	// void set_ambient_eta(double ambient_eta) noexcept { /* always 1.0 */ } // No longer configurable
	// void set_partial(bool partial) noexcept { partial_ = partial; } // No longer configurable
	void set_log(bool log) noexcept { log_ = log; }
	void set_deterministic(bool deterministic) noexcept { deterministic_ = deterministic; }

	// Main parsing interface
	// Configuration parsing with structured error handling
	Result<void, ConfigError> parse_config_file(const std::string& filename);

	// TOML parsing methods
	bool parse_general_config(const toml::table& config);
	bool parse_source_config(const toml::table& config);
	bool parse_layer_configs(const toml::table& config);

private:
	// Private constructor for singleton
	Config() = default;
};
