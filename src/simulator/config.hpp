#pragma once

#include <cstdint>
#include <list>
#include <map>
#include <memory>
#include <sstream>
#include <stdexcept>
#include <string>
#include <string_view>
#include <vector>
#include <span>

#include "math/concepts.hpp"
#include "simulator/layer.hpp"
#include "simulator/photon.hpp"
#include "simulator/material.hpp"

class Config
{
private:
	static std::unique_ptr<Config> instance_;
	static bool initialized_;

	uint32_t nx_ {0};          		// number of voxels in the x direction
	uint32_t ny_ {0};          		// number of voxels in the y direction
	uint32_t nz_ {0};          		// number of voxels in the z direction

	double vox_size_ {0.0};      	// uniform size of each voxel (dx=dy=dz)

	uint64_t num_layers_ {0};  		// number of layers
	uint64_t num_voxels_ {0};  		// number of voxels
	uint64_t num_photons_ {0}; 		// number of photons
	uint64_t num_sources_ {0}; 		// number of light sources

	bool log_ {false};         		// display logging debug and progress messages
	bool deterministic_ {false};   	// use deterministic random seed for reproducible results
	
	// Parsing configuration data
	std::vector<Source> sources_;
	std::vector<Material> tissues_;
	std::vector<Layer> layers_;

	// Utility functions for parsing
	template<Numeric T>
	[[nodiscard]] constexpr T str2num(std::string_view str) noexcept {
		T num{};
		std::stringstream ss{std::string(str)};
		ss >> num;
		// Check for parsing failure or NaN
		if (ss.fail() || (std::floating_point<T> && num != num)) {
			return T{};
		}
		return num;
	}

	bool equals(std::string_view s, std::string_view t) { return s == t; }

	void trim_comment(std::string& str) {
		size_t position = str.find_first_of('#');
		if (position != str.npos) {
			str.erase(position);
		}
	}

	void trim_spaces(std::string& str) {
		static constexpr char delims[] = " \t\r\n"; // space, tab, return, newline
		size_t position = str.find_first_not_of(delims);
		if (position == std::string::npos) {
			str = "";                   // String is all whitespace
			return;
		}
		str.erase(0, position);
		position = str.find_last_not_of(delims);
		if (position != std::string::npos) {
			str.erase(position + 1);
		}
	}

	[[nodiscard]] std::vector<double> split(std::string_view str, char split_char) {
		std::vector<double> dblarray;
		dblarray.clear();
		dblarray.reserve(8); // Reserve space for common case

		auto position = str.find(split_char);
		size_t beg = 0;

		// decompose statement using ranges when possible
		while (position != std::string_view::npos) {
			const auto piece_view = str.substr(beg, position - beg);
			std::string piece{piece_view}; // Convert to string for trim operations
			trim_spaces(piece);
			if (!piece.empty()) {
				const double val = str2num<double>(piece);
				dblarray.push_back(val);
			}
			beg = position + 1;
			position = str.find(split_char, beg);
		}

		// add the last piece if there's anything left
		if (beg < str.size()) {
			const auto piece_view = str.substr(beg);
			std::string piece{piece_view};
			trim_spaces(piece);
			if (!piece.empty()) {
				const double val = str2num<double>(piece);
				dblarray.push_back(val);
			}
		}

		return dblarray;
	}

	template<ConfigContainer Container>
	[[nodiscard]] std::vector<std::pair<std::string, std::string>> parameter_values(const Container& lines) 
		requires std::convertible_to<std::ranges::range_value_t<Container>, std::string> {
		std::vector<std::pair<std::string, std::string>> param_values;
		param_values.reserve(std::ranges::size(lines)); // Pre-allocate for efficiency

		for (const auto& line : lines) {
			const auto delim = line.find_first_of("=");

			// parameter name
			std::string param = line.substr(0, delim);
			trim_spaces(param);

			// parameter value
			std::string value = line.substr(delim + 1);
			trim_comment(value);
			trim_spaces(value);

			param_values.emplace_back(std::move(param), std::move(value));
		}

		return param_values;
	}

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
	 */
	static void initialize(const std::string& config_file);
	
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
	[[nodiscard]] constexpr bool partial() const noexcept { return true; } // Always enabled
	[[nodiscard]] constexpr bool log() const noexcept { return log_; }
	[[nodiscard]] constexpr bool deterministic() const noexcept { return deterministic_; }

	// Access parsed data with span for zero-copy access
	[[nodiscard]] std::span<const Source> sources() const noexcept { 
		return std::span<const Source>(sources_); 
	}
	[[nodiscard]] std::span<const Material> tissues() const noexcept { 
		return std::span<const Material>(tissues_); 
	}
	[[nodiscard]] std::span<const Layer> layers() const noexcept { 
		return std::span<const Layer>(layers_); 
	}
	
	// Compatibility accessors returning references
	[[nodiscard]] const std::vector<Source>& sources_vector() const noexcept { return sources_; }
	[[nodiscard]] const std::vector<Material>& tissues_vector() const noexcept { return tissues_; }
	[[nodiscard]] const std::vector<Layer>& layers_vector() const noexcept { return layers_; }
	
	// Move versions for when transferring ownership
	[[nodiscard]] std::vector<Layer>&& move_layers() { return std::move(layers_); }
	[[nodiscard]] std::vector<Material>&& move_tissues() { return std::move(tissues_); }

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
	bool parse_config_file(const std::string& filename);

	// Section-specific parsing methods
	bool parse_general_config(std::list<std::string>& data);
	bool parse_source_config(std::list<std::string>& data);
	bool parse_tissue_config(std::list<std::string>& data);
	bool parse_layer_config(std::list<std::string>& data);

private:
	// Internal parsing utilities
	bool extract_config_data(const std::string& filename, std::multimap<std::string, std::list<std::string>>& datamap);
	void extract_config_block(std::ifstream& in_config, const std::string& section,
							  std::multimap<std::string, std::list<std::string>>& datamap);

	// Private constructor for singleton
	Config() = default;
};
