#pragma once

#include <cstdint>
#include <list>
#include <map>
#include <sstream>
#include <string>
#include <string_view>
#include <vector>

#include "simulator/layer.hpp"
#include "simulator/source.hpp"
#include "simulator/tissue.hpp"

class Config
{
private:
	uint32_t nx_ = 0;          // number of voxels in the x direction
	uint32_t ny_ = 0;          // number of voxels in the y direction
	uint32_t nz_ = 0;          // number of voxels in the z direction

	uint64_t num_layers_ = 0;  // number of layers
	uint64_t num_voxels_ = 0;  // number of voxels
	uint64_t num_photons_ = 0; // number of photons
	uint64_t num_sources_ = 0; // number of light sources

	double vox_size_ = 0.0;      // uniform size of each voxel (dx=dy=dz)
	double ambient_eta_ = 0.0;   // refractive index of ambient medium

	bool partial_ = true;         // partial reflection
	bool progress_ = false;        // display progress messages

	// Parsing configuration data
	std::vector<Source> sources_;
	std::vector<Tissue> tissues_;
	std::vector<Layer> layers_;

	// Utility functions for parsing
	template<typename T>
	T str2num(const std::string& str) {
		T num;
		std::stringstream ss(str);
		ss >> num;
		// isnan check
		if (num != num) {
			return T();
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

	std::vector<double> split(const std::string& str, char split_char) {
		std::vector<double> dblarray;
		dblarray.clear();

		size_t position = str.find(split_char);
		size_t beg = 0;

		// decompose statement
		while (position != std::string::npos) {
			std::string piece = str.substr(beg, position - beg);
			trim_spaces(piece);
			if (!piece.empty()) {
				double val = str2num<double>(piece);
				dblarray.push_back(val);
			}
			beg = position + 1;
			position = str.find(split_char, beg);
		}

		// add the last piece if there's anything left
		if (beg < str.size()) {
			std::string piece = str.substr(beg);
			trim_spaces(piece);
			if (!piece.empty()) {
				double val = str2num<double>(piece);
				dblarray.push_back(val);
			}
		}

		return dblarray;
	}

	std::vector<std::pair<std::string, std::string>> parameter_values(std::list<std::string>& lines) {
		std::vector<std::pair<std::string, std::string>> param_values;

		for (const auto& line : lines) {
			size_t delim = line.find_first_of("=");

			// parameter name
			std::string param = line.substr(0, delim);
			trim_spaces(param);

			// parameter value
			std::string value = line.substr(delim + 1);
			trim_comment(value);
			trim_spaces(value);

			param_values.push_back(std::make_pair(param, value));
		}

		return param_values;
	}

public:
	Config() = default;

	// Getters
	uint32_t nx() const noexcept { return nx_; }
	uint32_t ny() const noexcept { return ny_; }
	uint32_t nz() const noexcept { return nz_; }
	uint64_t num_layers() const noexcept { return num_layers_; }
	uint64_t num_voxels() const noexcept { return num_voxels_; }
	uint64_t num_photons() const noexcept { return num_photons_; }
	uint64_t num_sources() const noexcept { return num_sources_; }
	double vox_size() const noexcept { return vox_size_; }
	double ambient_eta() const noexcept { return ambient_eta_; }
	bool partial() const noexcept { return partial_; }
	bool progress() const noexcept { return progress_; }

	// Access parsed data
	const std::vector<Source>& sources() const noexcept { return sources_; }
	const std::vector<Tissue>& tissues() const noexcept { return tissues_; }
	const std::vector<Layer>& layers() const noexcept { return layers_; }
	
	// Move versions for when transferring ownership
	std::vector<Layer>&& move_layers() { return std::move(layers_); }

	// Setters
	void set_nx(uint32_t nx) noexcept { nx_ = nx; }
	void set_ny(uint32_t ny) noexcept { ny_ = ny; }
	void set_nz(uint32_t nz) noexcept { nz_ = nz; }
	void set_num_layers(uint64_t num_layers) noexcept { num_layers_ = num_layers; }
	void set_num_voxels(uint64_t num_voxels) noexcept { num_voxels_ = num_voxels; }
	void set_num_photons(uint64_t num_photons) noexcept { num_photons_ = num_photons; }
	void set_num_sources(uint64_t num_sources) noexcept { num_sources_ = num_sources; }
	void set_vox_size(double vox_size) noexcept { vox_size_ = vox_size; }
	void set_ambient_eta(double ambient_eta) noexcept { ambient_eta_ = ambient_eta; }
	void set_partial(bool partial) noexcept { partial_ = partial; }
	void set_progress(bool progress) noexcept { progress_ = progress; }

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
};
