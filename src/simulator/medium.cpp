#include "medium.hpp"

#include <algorithm>
#include <cfloat>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <map>

#include "math/ray.hpp"
#include "voxel.hpp"

Medium::Medium(Config& config) : config_(config) {
	// Reserve space for common use cases
	layers_.reserve(10);
	tissues_.reserve(5);

	// Transfer layer and tissue data from config to medium
	layers_ = config_.move_layers();
	tissues_ = config_.move_tissues();

	// Update geometry for each layer after transfer
	for (auto& layer : layers_) {
		layer.update_geometry();
	}

	// Initialize bounds to extreme values for proper min/max calculation
	bounds_.min_bounds = glm::dvec3(std::numeric_limits<double>::max());
	bounds_.max_bounds = glm::dvec3(-std::numeric_limits<double>::max());
}

bool Medium::initialize() {
	// Initialize volume
	if (!initialize_volume()) {
		std::cerr << "An error occurred while initializing the voxel grid." << std::endl;
		return false;
	}
	std::cout << "Voxel grid initialized successfully." << std::endl;

	// Intialize layers
	if (!initialize_layers()) {
		std::cerr << "An error occurred while initializing the layers." << std::endl;
		return false;
	}
	std::cout << "Layers initialized successfully." << std::endl;

	// Reset simulation data after voxels are created
	reset_simulation_data();
    std::cout << "Simulation data reset completed." << std::endl;

	// voxelize the geometry
	if (!voxelize_layers()) {
		std::cerr << "An error occurred during geometry voxelization." << std::endl;
		return false;
	}
	std::cout << "Geometry voxelization completed successfully." << std::endl;

	// Identify surface voxels
	identify_surface_voxels();
	std::cout << "Surface voxel identification completed." << std::endl;

	return true;
}

bool Medium::initialize_volume() {
	// Initialize bounds to extreme values for proper min/max calculation
	bounds_.min_bounds = glm::dvec3(std::numeric_limits<double>::max());
	bounds_.max_bounds = glm::dvec3(-std::numeric_limits<double>::max());

	// compute grid boundary extent with padding to ensure full coverage
	for (const auto& layer : layers_) {
		// see if a vertex denotes a new boundary
		for (const auto& triangle : layer.mesh) {
			glm::dvec3 v0 = triangle.v0();
			glm::dvec3 v1 = triangle.v1();
			glm::dvec3 v2 = triangle.v2();

			// get the maximum value among the previous maximum or new vertices
			bounds_.min_bounds.x = std::min({bounds_.min_bounds.x, v0.x, v1.x, v2.x}); // left (-x)
			bounds_.max_bounds.x = std::max({bounds_.max_bounds.x, v0.x, v1.x, v2.x}); // right (+x)
			bounds_.min_bounds.y = std::min({bounds_.min_bounds.y, v0.y, v1.y, v2.y}); // top (+y)
			bounds_.max_bounds.y = std::max({bounds_.max_bounds.y, v0.y, v1.y, v2.y}); // bottom (-y)
			bounds_.min_bounds.z = std::min({bounds_.min_bounds.z, v0.z, v1.z, v2.z}); // front (+z)
			bounds_.max_bounds.z = std::max({bounds_.max_bounds.z, v0.z, v1.z, v2.z}); // rear (-z)
		}
	}

	// Add padding to ensure voxels fully cover the mesh geometry
	double padding = config_.vox_size() * 2.0; // Two voxel widths of padding
	bounds_.min_bounds -= glm::dvec3(padding);
	bounds_.max_bounds += glm::dvec3(padding);

	// check for inconsistency and zero width/height/depth
	if ((bounds_.min_bounds.x >= bounds_.max_bounds.x) || (bounds_.min_bounds.y >= bounds_.max_bounds.y)
		|| (bounds_.min_bounds.z >= bounds_.max_bounds.z)) {
		std::cerr << "Invalid bounds detected. Zero or negative dimensions." << std::endl;
		return false;
	}

	// Get the size vector and extract individual dimensions
	auto size_vec = bounds_.size();

	// Number of voxels in each dimension (as float)
	float nx = static_cast<float>(size_vec.x / config_.vox_size());
	float ny = static_cast<float>(size_vec.y / config_.vox_size());
	float nz = static_cast<float>(size_vec.z / config_.vox_size());

	// Initialize the voxel grid with proper dimensions
	uint32_t grid_nx = static_cast<uint32_t>(nx);
	uint32_t grid_ny = static_cast<uint32_t>(ny);
	uint32_t grid_nz = static_cast<uint32_t>(nz);

	// TODO: Avoid setting Config values from here

	// number of voxels in each dimension
	config_.set_nx(grid_nx);
	config_.set_ny(grid_ny);
	config_.set_nz(grid_nz);

	// total number of voxels
	config_.set_num_voxels(grid_nx * grid_ny * grid_nz);

	// Check for voxel sizes that are too large
	if (config_.vox_size() > size_vec.x || config_.vox_size() > size_vec.y || config_.vox_size() > size_vec.z) {
		std::cerr << "Error: Voxel size is too large compared to geometry dimensions" << std::endl;
		return false;
	}

	// check grid dimensions
	if (config_.num_voxels() < 1 || config_.nx() < 1 || config_.ny() < 1 || config_.nz() < 1) {
		return false;
	}

	volume_ = Volume(config_.vox_size(), grid_nx, grid_ny, grid_nz);

	return true;
}

bool Medium::initialize_layers() {
	config_.set_num_layers(static_cast<uint64_t>(layers_.size()));

	// error checking
	if (config_.num_layers() < 1 || config_.num_voxels() < 1) {
		return false;
	}

	// check for duplicate layers
	for (size_t i = 1; i < layers_.size(); ++i) {
		if (layers_[i - 1] == layers_[i]) {
			return false;
		}
	}

	// check for duplicate tissues
	for (uint32_t i = 1; i < tissues_.size(); ++i) {
		if (tissues_[i - 1] == tissues_[i]) {
			return false;
		}
	}

	// check if layer's tissue id is out of range
	for (const auto& layer : layers_) {
		if (layer.tissue_id >= tissues_.size()) {
			return false;
		}
	}

	return true;
}

bool Medium::voxelize_layers() {
	uint32_t nx = config_.nx();
	uint32_t ny = config_.ny();
	uint32_t nz = config_.nz();
	double vox_size = config_.vox_size();

	std::cout << "Voxelizing geometry using point-containment method..." << std::endl;
	std::cout << "Grid dimensions: " 
              << nx << "x" 
              << ny << "x" 
              << nz << " (total: " 
              << (nx * ny * nz) << " voxels)"
			  << std::endl;
	std::cout << "Voxel size: " << vox_size << std::endl;
	std::cout << "Grid bounds: min(" 
              << bounds_.min_bounds.x << "," 
              << bounds_.min_bounds.y << "," 
              << bounds_.min_bounds.z << ") max(" 
              << bounds_.max_bounds.x << "," 
              << bounds_.max_bounds.y << "," 
              << bounds_.max_bounds.z << ")"
			  << std::endl;

	int total_voxels_assigned = 0;
	int partial_voxels = 0;

	for (uint32_t iz = 0; iz < nz; iz++) {
		for (uint32_t iy = 0; iy < ny; iy++) {
			for (uint32_t ix = 0; ix < nx; ix++) {
				// Calculate voxel center in world coordinates
				glm::dvec3 voxel_center = glm::dvec3(bounds_.min_bounds.x + (ix + 0.5) * vox_size,
                                                     bounds_.min_bounds.y + (iy + 0.5) * vox_size,
                                                     bounds_.min_bounds.z + (iz + 0.5) * vox_size);

				// Calculate voxel corners to check for partial intersection
				glm::dvec3 voxel_min = glm::dvec3(bounds_.min_bounds.x + ix * vox_size, 
                                                  bounds_.min_bounds.y + iy * vox_size,
                                                  bounds_.min_bounds.z + iz * vox_size);
				
                glm::dvec3 voxel_max = voxel_min + glm::dvec3(vox_size);

				// Check if voxel intersects with geometry
				bool center_inside = false;
				bool any_corner_inside = false;

				// Test center and corners
				for (const auto& layer : layers_) {
					if (layer.contains_point(voxel_center)) {
						center_inside = true;
						break;
					}
				}

				// Test corners to detect partial intersection
				std::vector<glm::dvec3> corners = {
                    voxel_min,
					{voxel_max.x, voxel_min.y, voxel_min.z},
					{voxel_min.x, voxel_max.y, voxel_min.z},
					{voxel_max.x, voxel_max.y, voxel_min.z},
					{voxel_min.x, voxel_min.y, voxel_max.z},
					{voxel_max.x, voxel_min.y, voxel_max.z},
					{voxel_min.x, voxel_max.y, voxel_max.z},
					voxel_max
                };

				for (const auto& corner : corners) {
					for (const auto& layer : layers_) {
						if (layer.contains_point(corner)) {
							any_corner_inside = true;
							break;
						}
					}
					if (any_corner_inside) {
						break;
                    }
				}

				// Assign tissue if center is inside OR any corner is inside (partial voxel)
				if (center_inside || any_corner_inside) {
					Voxel* voxel = volume_(ix, iy, iz);

					// Find which specific layer contains this voxel (check center point)
					// IMPORTANT: Process layers in reverse order so inner layers override outer layers
					for (const auto& layer : layers_) {
						if (layer.contains_point(voxel_center)) {
							voxel->tissue = &tissues_[layer.tissue_id];
							break; // Use the innermost layer that contains the center
						}
					}

					// If center point didn't hit any layer but corners did,
					// use the innermost layer that contains any corner
					if (voxel->tissue == nullptr) {
						for (const auto& corner : corners) {
							for (const auto& layer : layers_) {
								if (layer.contains_point(corner)) {
									voxel->tissue = &tissues_[layer.tissue_id];
									break;
								}
							}
                            if (voxel->tissue != nullptr) {
								break;
                            }
						}
					}

					// Compute volume fractions for proper boundary physics
					if (center_inside && any_corner_inside) {
						// Check if this is actually a boundary voxel
						bool all_corners_inside = true;
						for (const auto& corner : corners) {
							bool corner_inside = false;
							for (const auto& layer : layers_) {
								if (layer.contains_point(corner)) {
									corner_inside = true;
									break;
								}
							}
							if (!corner_inside) {
								all_corners_inside = false;
								break;
							}
						}

						if (all_corners_inside) {
							// Fully inside
							voxel->volume_fraction_inside = 1.0;
							voxel->volume_fraction_outside = 0.0;
							voxel->is_boundary_voxel = false;
						}
						else {
							// Partial voxel - compute accurate volume fractions
							voxel->volume_fraction_inside = volume_.fraction_inside_fast(voxel_min, voxel_max, layers_, 2);
							voxel->volume_fraction_outside = 1.0 - voxel->volume_fraction_inside;
							voxel->is_boundary_voxel = true;
						}
					}
					else if (center_inside) {
						// Center inside but some corners outside - partial
						voxel->volume_fraction_inside = volume_.fraction_inside_fast(voxel_min, voxel_max, layers_, 2);
						voxel->volume_fraction_outside = 1.0 - voxel->volume_fraction_inside;
						voxel->is_boundary_voxel = true;
					}
					else {
						// Only corners inside - partial
						voxel->volume_fraction_inside = volume_.fraction_inside_fast(voxel_min, voxel_max, layers_, 2);
						voxel->volume_fraction_outside = 1.0 - voxel->volume_fraction_inside;
						voxel->is_boundary_voxel = true;
					}

					total_voxels_assigned++;

					if (voxel->is_boundary_voxel) {
						partial_voxels++;
					}
				}
			}
		}
	}

	std::cout << "Point-containment voxelization completed. Assigned tissue to " << total_voxels_assigned << " voxels ("
			  << partial_voxels << " partial)." << std::endl;
	return true;
}

void Medium::identify_surface_voxels() {
	int surface_voxel_count = 0;
	const glm::uvec3& dimensions = volume_.dimensions();
	
	for (uint32_t z = 0; z < dimensions.z; z++) {
		for (uint32_t y = 0; y < dimensions.y; y++) {
			for (uint32_t x = 0; x < dimensions.x; x++) {
				Voxel* voxel = volume_.at(x, y, z);
				
				// Only consider voxels that have tissue assigned (inside or partially inside geometry)
				if (!voxel->tissue) {
					voxel->is_surface_voxel = false;
					continue;
				}
				
				bool is_surface = false;
				
				// Method 1: Boundary voxels are automatically surface voxels
				if (voxel->is_boundary_voxel) {
					is_surface = true;
				}
				
				// Method 2: Check for empty neighbors (traditional surface detection)
				if (!is_surface) {
					bool has_empty_neighbor = false;
					
					for (int dz = -1; dz <= 1; dz++) {
						for (int dy = -1; dy <= 1; dy++) {
							for (int dx = -1; dx <= 1; dx++) {
								if (dx == 0 && dy == 0 && dz == 0){
                                    continue; // Skip self
                                }
								
								int nx = (int)x + dx;
								int ny = (int)y + dy;
								int nz = (int)z + dz;
								
								// Check bounds
								if (nx < 0 || ny < 0 || nz < 0 
                                || 	nx >= (int)dimensions.x
                                || 	ny >= (int)dimensions.y 
                                || 	nz >= (int)dimensions.z) {
									// Out of bounds = empty neighbor
									has_empty_neighbor = true;
									break;
								}
								
								// Check if neighbor has no tissue (is empty)
								Voxel* neighbor = volume_.at(nx, ny, nz);
								if (!neighbor->tissue) {
									has_empty_neighbor = true;
									break;
								}
							}
							if (has_empty_neighbor) {
                                break;
                            }
						}
						if (has_empty_neighbor) {
                            break;
                        }
					}
					
					if (has_empty_neighbor) {
						is_surface = true;
					}
				}
				
				// Method 3: Detect grazing surface voxels by checking if voxel boundary intersects any mesh surface
				if (!is_surface) {
					// Get voxel bounding box
					Cuboid voxel_box = voxel_corners(voxel);
					glm::dvec3 voxel_min = voxel_box.min_point();
					glm::dvec3 voxel_max = voxel_box.max_point();
					
					// Dynamic tolerance based on voxel size to handle floating-point precision
					double voxel_size = volume_.voxel_size();
					double tolerance = voxel_size * 0.1; // 10% of voxel size for floating-point safety
					
					// Check if any mesh triangle intersects this voxel
					for (const auto& layer : layers_) {
						for (size_t face_idx = 0; face_idx < layer.mesh.size(); face_idx++) {
							const Triangle& triangle = layer.mesh[face_idx];
							
							// Get triangle bounding box with tolerance expansion
							glm::dvec3 tri_min = glm::min(glm::min(triangle.v0(), triangle.v1()), triangle.v2()) - tolerance;
							glm::dvec3 tri_max = glm::max(glm::max(triangle.v0(), triangle.v1()), triangle.v2()) + tolerance;
							
							// Expand voxel bounding box slightly for floating-point safety
							glm::dvec3 safe_voxel_min = voxel_min - tolerance;
							glm::dvec3 safe_voxel_max = voxel_max + tolerance;
							
							// Check if triangle bounding box overlaps with voxel bounding box
							bool overlaps_x = (tri_max.x >= safe_voxel_min.x) && (tri_min.x <= safe_voxel_max.x);
							bool overlaps_y = (tri_max.y >= safe_voxel_min.y) && (tri_min.y <= safe_voxel_max.y);
							bool overlaps_z = (tri_max.z >= safe_voxel_min.z) && (tri_min.z <= safe_voxel_max.z);
							
							if (overlaps_x && overlaps_y && overlaps_z) {
								is_surface = true;
								break;
							}
						}
						if (is_surface) {
                            break;
                        }
					}
				}
				
				voxel->is_surface_voxel = is_surface;
				if (is_surface) {
					surface_voxel_count++;
				}
			}
		}
	}
	
	std::cout << "Identified " << surface_voxel_count << " surface voxels." << std::endl;
}

void Medium::reset_simulation_data() {
	// Reset record
	record_ = Record();

	// Reset voxel data
	for (auto& voxel_ptr : volume_) {
		if (voxel_ptr) {
			voxel_ptr->absorption = 0.0;
			voxel_ptr->emittance = 0.0;
		}
	}

	// Reset metrics
	metrics_ = Metrics();
}

void Medium::normalize() {
	// Normalize globally recorded parameters
	record_.total_absorption /= config_.num_photons() * config_.num_sources();
	record_.diffuse_reflection /= config_.num_photons() * config_.num_sources();
	record_.diffuse_transmission /= config_.num_photons() * config_.num_sources();

	record_.specular_reflection /= config_.num_photons() * config_.num_sources();
	record_.specular_transmission /= config_.num_photons() * config_.num_sources();

	// Normalize voxel data
	for (const auto& voxel_ptr : volume_) {
		if (voxel_ptr && voxel_ptr->tissue) {
			voxel_ptr->absorption /= config_.num_photons() * config_.num_sources();
			voxel_ptr->emittance /= config_.num_photons() * config_.num_sources();
		}
	}
}

void Medium::write_results() const {
	std::string str_abs = "absorption.out";
	std::string str_emi = "emittance.out";

	std::ofstream ofs_abs(str_abs.c_str(), std::ios_base::out);
	std::ofstream ofs_emi(str_emi.c_str(), std::ios_base::out);

	if (!ofs_abs.good() || !ofs_emi.good()) {
		std::cerr << "Error: Could not open output files for writing" << std::endl;
		return;
	}

	// Write voxel absorption and emittance data
	const glm::uvec3& dimensions = volume_.dimensions();

	for (uint32_t z = 0; z < dimensions.z; z++) {
		for (uint32_t y = 0; y < dimensions.y; y++) {
			for (uint32_t x = 0; x < dimensions.x; x++) {
				const Voxel* voxel = volume_(x, y, z);
				if (voxel && voxel->tissue) {
					ofs_abs << voxel->absorption << " ";
					ofs_emi << voxel->emittance << " ";
				}
				else {
					ofs_abs << "0.0 ";
					ofs_emi << "0.0 ";
				}
			}
			ofs_abs << std::endl;
			ofs_emi << std::endl;
		}
		ofs_abs << std::endl;
		ofs_emi << std::endl;
	}
}

Voxel* Medium::voxel_at(glm::dvec3& position) {
	if (!bounds_.includes(position)) {
		return nullptr;
	}

	// Calculate relative position from minimum bounds
	double dx = position.x - bounds_.min_bounds.x;
	double dy = position.y - bounds_.min_bounds.y;
	double dz = position.z - bounds_.min_bounds.z;

	// indices start at minimum boundaries of voxels
	uint32_t ix = static_cast<uint32_t>(std::floor(dx / config_.vox_size()));
	uint32_t iy = static_cast<uint32_t>(std::floor(dy / config_.vox_size()));
	uint32_t iz = static_cast<uint32_t>(std::floor(dz / config_.vox_size()));

	// avoid index overflow
	if (ix >= config_.nx()) {
		ix = config_.nx() - 1;
	}
	if (iy >= config_.ny()) {
		iy = config_.ny() - 1;
	}
	if (iz >= config_.nz()) {
		iz = config_.nz() - 1;
	}

	// retrieve the voxel at the position using the Volume
	if (!volume_.is_valid_coordinate(ix, iy, iz)) {
		return nullptr;
	}

	Voxel* voxel = volume_(ix, iy, iz);

	// if voxel does not have a tissue, it is outside the medium
	return voxel->tissue ? voxel : nullptr;
}

Cuboid Medium::voxel_corners(Voxel* voxel) const {
	uint32_t ix_min = voxel->ix();
	uint32_t iy_min = voxel->iy();
	uint32_t iz_min = voxel->iz();

	uint32_t ix_max = voxel->ix() + 1;
	uint32_t iy_max = voxel->iy() + 1;
	uint32_t iz_max = voxel->iz() + 1;

	// minimum voxel position
	float x_min = static_cast<float>(bounds_.min_bounds.x + (config_.vox_size() * ix_min));
	float y_min = static_cast<float>(bounds_.min_bounds.y + (config_.vox_size() * iy_min));
	float z_min = static_cast<float>(bounds_.min_bounds.z + (config_.vox_size() * iz_min));

	// maximum voxel position
	float x_max = static_cast<float>(bounds_.min_bounds.x + (config_.vox_size() * ix_max));
	float y_max = static_cast<float>(bounds_.min_bounds.y + (config_.vox_size() * iy_max));
	float z_max = static_cast<float>(bounds_.min_bounds.z + (config_.vox_size() * iz_max));

	// round off coordinate values around the origin
	x_min = (std::fabs(x_min) < 1E-10) ? 0 : x_min;
	y_min = (std::fabs(y_min) < 1E-10) ? 0 : y_min;
	z_min = (std::fabs(z_min) < 1E-10) ? 0 : z_min;

	x_max = (std::fabs(x_max) < 1E-10) ? 0 : x_max;
	y_max = (std::fabs(y_max) < 1E-10) ? 0 : y_max;
	z_max = (std::fabs(z_max) < 1E-10) ? 0 : z_max;

	return Cuboid(x_min, y_min, z_min, x_max, y_max, z_max);
}

double Medium::intersection(Source& source) const {
	Ray ray = Ray(source.origin, source.direction);
	Triangle intersected_triangle{};
	glm::dvec3 intersection_point{};

	double closest_distance = std::numeric_limits<double>::max();

	for (const auto& layer : layers_) {
		double distance = ray.intersect_layer(layer, intersection_point);
		if (distance < closest_distance) {
			closest_distance = distance;
		}
	}

    source.intersect = intersection_point;
	source.triangle = intersected_triangle;

	return closest_distance;
}

bool Medium::contains_point(const glm::dvec3& point) const {
	// Test if the point is inside any layer's geometry
	for (const auto& layer : layers_) {
		if (layer.contains_point(point)) {
			return true;
		}
	}

	return false;
}
