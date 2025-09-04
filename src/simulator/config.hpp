#pragma once

#include <cstdint>
#include <list>
#include <map>
#include <sstream>
#include <string>
#include <vector>

#include "simulator/layer.hpp"
#include "simulator/source.hpp"
#include "simulator/tissue.hpp"

class Config
{
private:
	uint32_t nx_;          // number of voxels in the x direction
	uint32_t ny_;          // number of voxels in the y direction
	uint32_t nz_;          // number of voxels in the z direction

	uint64_t num_layers_;  // number of layers
	uint64_t num_voxels_;  // number of voxels
	uint64_t num_photons_; // number of photons
	uint64_t num_sources_; // number of light sources

	double vox_size_;      // uniform size of each voxel (dx=dy=dz)
	double ambient_eta_;   // refractive index of ambient medium

	bool partial_;         // partial reflection
	bool progress_;        // display progress messages

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

	bool equals(const std::string& s, const std::string& t) { return (s.compare(t) == 0); }

	void trim_comment(std::string& str) {
		size_t position = str.find_first_of('#');
		if (position != str.npos) {
			str.erase(position);
		}
	}

	void trim_spaces(std::string& str) {
		char const* delims = " \t\r\n"; // space, tab, return, newline
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
	Config();

	// Getters
	uint32_t nx() const { return nx_; }
	uint32_t ny() const { return ny_; }
	uint32_t nz() const { return nz_; }
	uint64_t num_layers() const { return num_layers_; }
	uint64_t num_voxels() const { return num_voxels_; }
	uint64_t num_photons() const { return num_photons_; }
	uint64_t num_sources() const { return num_sources_; }
	double vox_size() const { return vox_size_; }
	double ambient_eta() const { return ambient_eta_; }
	bool partial() const { return partial_; }
	bool progress() const { return progress_; }

	// Access parsed data
	const std::vector<Source>& sources() const { return sources_; }
	const std::vector<Tissue>& tissues() const { return tissues_; }
	const std::vector<Layer>& layers() const { return layers_; }
	
	// Move versions for when transferring ownership
	std::vector<Layer>&& move_layers() { return std::move(layers_); }

	// Setters
	void set_nx(uint32_t nx) { nx_ = nx; }
	void set_ny(uint32_t ny) { ny_ = ny; }
	void set_nz(uint32_t nz) { nz_ = nz; }
	void set_num_layers(uint64_t num_layers) { num_layers_ = num_layers; }
	void set_num_voxels(uint64_t num_voxels) { num_voxels_ = num_voxels; }
	void set_num_photons(uint64_t num_photons) { num_photons_ = num_photons; }
	void set_num_sources(uint64_t num_sources) { num_sources_ = num_sources; }
	void set_vox_size(double vox_size) { vox_size_ = vox_size; }
	void set_ambient_eta(double ambient_eta) { ambient_eta_ = ambient_eta; }
	void set_partial(bool partial) { partial_ = partial; }
	void set_progress(bool progress) { progress_ = progress; }

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
