#include "medium.hpp"

#include <algorithm>
#include <cfloat>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <map>

#include "math/ray.hpp"
#include "voxel.hpp"
#include "../app.hpp" // For output path utilities

Medium::Medium(Config& config) : config_(config) {
	// Reserve space for common use cases
	layers_.reserve(10);
	tissues_.reserve(5);

	// Transfer layer and material data from config to medium
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
	if (config_.verbose()) {
		std::cout << "Voxel grid initialized successfully." << std::endl;
	}

	// Intialize layers
	if (!initialize_layers()) {
		std::cerr << "An error occurred while initializing the layers." << std::endl;
		return false;
	}
	if (config_.verbose()) {
		std::cout << "Layers initialized successfully." << std::endl;
	}

	// Reset simulation data after voxels are created
	reset_simulation_data();
	if (config_.verbose()) {
		std::cout << "Simulation data reset completed." << std::endl;
	}

	// voxelize the geometry
	if (!voxelize_layers()) {
		std::cerr << "An error occurred during geometry voxelization." << std::endl;
		return false;
	}
	if (config_.verbose()) {
		std::cout << "Geometry voxelization completed successfully." << std::endl;
	}

	// PHASE 2: Skip old surface detection - Distance Field Based voxelization already set correct surface flags
	// disambiguate_surface_voxels();  // DISABLED: SDF voxelization is already accurate
	if (config_.verbose()) {
		std::cout << "Surface voxel disambiguation skipped (using SDF results)." << std::endl;
	}

	// Identify surface voxels - DISABLED: SDF already set is_surface_voxel correctly
	// identify_surface_voxels();  // DISABLED: Trust SDF voxelization results
	if (config_.verbose()) {
		std::cout << "Surface voxel identification skipped (using SDF results)." << std::endl;
	}

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

	// check if layer's material id is out of range
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

	if (config_.verbose()) {
		std::cout << "Voxelizing geometry..." << std::endl;
		std::cout << "Grid dimensions: " 
		          << nx << "x" << ny << "x" << nz 
				  << " (total " << (nx * ny * nz) 
				  << " voxels, size: " << vox_size 
				  << ")" << std::endl;
		std::cout << "Grid bounds: min(" 
		          << bounds_.min_bounds.x << "," 
		          << bounds_.min_bounds.y << "," 
		          << bounds_.min_bounds.z << ") max(" 
		          << bounds_.max_bounds.x << "," 
		          << bounds_.max_bounds.y << "," 
		          << bounds_.max_bounds.z << ")"
				  << std::endl;
	}

	int total_voxels_assigned = 0;
	int partial_voxels = 0;

	// Use Distance Field Based voxelization for robust surface detection
	for (uint32_t iz = 0; iz < nz; iz++) {
		for (uint32_t iy = 0; iy < ny; iy++) {
			for (uint32_t ix = 0; ix < nx; ix++) {
				// Calculate voxel position in world coordinates
				glm::dvec3 voxel_min = glm::dvec3(bounds_.min_bounds.x + ix * vox_size, 
                                                  bounds_.min_bounds.y + iy * vox_size,
                                                  bounds_.min_bounds.z + iz * vox_size);
				glm::dvec3 voxel_max = voxel_min + glm::dvec3(vox_size);

				// Use Distance Field Based voxelization for robust classification
				VoxelClassification classification = volume_.distance_field_voxelization(voxel_min, voxel_max, layers_);

				if (classification.is_inside_geometry) {
					Voxel* voxel = volume_(ix, iy, iz);
					
					// Assign material from dominant layer
					voxel->material = &tissues_[classification.dominant_tissue_id];
					voxel->volume_fraction_inside = classification.volume_fraction;
					voxel->volume_fraction_outside = 1.0 - classification.volume_fraction;
					voxel->is_boundary_voxel = classification.is_boundary_voxel;
					voxel->is_surface_voxel = classification.is_surface_voxel;
					
					total_voxels_assigned++;
					if (voxel->is_boundary_voxel) {
						partial_voxels++;
					}
				}
			}
		}
	}

	if (config_.verbose()) {
		std::cout << "Point-containment voxelization completed. Assigned material to " << total_voxels_assigned << " voxels ("
				  << partial_voxels << " partial)." << std::endl;
	}

	// PHASE 2: Detect external surface voxels (voxels adjacent to ambient)
	detect_external_surface_voxels();
	
	return true;
}

void Medium::disambiguate_surface_voxels() {
	const glm::uvec3& dimensions = volume_.dimensions();
	int false_positives_removed = 0;
	
	std::cout << "Starting surface disambiguation for boundary voxels..." << std::endl;
	
	// Debug: Check specific problematic voxels first
	int problematic_coords[][3] = {{20, 7, 2}, {19, 9, 5}, {13, 8, 12}};
	for (int i = 0; i < 3; i++) {
		int x = problematic_coords[i][0], y = problematic_coords[i][1], z = problematic_coords[i][2];
		if (x < (int)dimensions.x && y < (int)dimensions.y && z < (int)dimensions.z) {
			Voxel* voxel = volume_.at(x, y, z);
			std::cout << "DEBUG: Problematic voxel (" << x << "," << y << "," << z 
					  << ") - material=" << (voxel->material ? "YES" : "NO")
					  << " boundary=" << (voxel->is_boundary_voxel ? "YES" : "NO")
					  << " surface=" << (voxel->is_surface_voxel ? "YES" : "NO")
					  << " volume_fraction=" << voxel->volume_fraction_inside << std::endl;
		}
	}
	
	// Check all voxels with material for potential false positives (not just boundary voxels)
	for (uint32_t z = 0; z < dimensions.z; z++) {
		for (uint32_t y = 0; y < dimensions.y; y++) {
			for (uint32_t x = 0; x < dimensions.x; x++) {
				Voxel* voxel = volume_.at(x, y, z);
				
				// Only process voxels that have material
				if (!voxel->material) {
					continue;
				}
				
				// Debug output for specific problematic voxels
				bool is_problematic = (x == 20 && y == 7 && z == 2) || (x == 19 && y == 9 && z == 5) || (x == 13 && y == 8 && z == 12);
				if (is_problematic) {
					std::cout << "DEBUG: Processing boundary voxel (" << x << "," << y << "," << z 
							  << ") with volume_inside=" << voxel->volume_fraction_inside << std::endl;
				}
				
				// Check 6-connected neighbors (faces only, not edges/corners)
				// These are the potential "opposing" voxels across a surface
				struct FaceOffset { int dx, dy, dz; };
				std::vector<FaceOffset> face_neighbors = {
					{-1, 0, 0}, {1, 0, 0},   // left, right
					{0, -1, 0}, {0, 1, 0},   // down, up  
					{0, 0, -1}, {0, 0, 1}    // back, front
				};
				
				bool should_remove_tissue = false;
				
				for (const auto& offset : face_neighbors) {
					int dx = offset.dx, dy = offset.dy, dz = offset.dz;
					int nx = (int)x + dx;
					int ny = (int)y + dy; 
					int nz = (int)z + dz;
					
					// Check bounds
					if (nx < 0 || ny < 0 || nz < 0 
					|| 	nx >= (int)dimensions.x
					|| 	ny >= (int)dimensions.y 
					|| 	nz >= (int)dimensions.z) {
						continue;
					}
					
					Voxel* neighbor = volume_.at(nx, ny, nz);
					
					// If neighbor is also a boundary voxel with material,
					// check if it has significantly more volume inside
					if (neighbor->material && neighbor->is_boundary_voxel) {
						// Surface disambiguation: if neighbor has much more volume inside,
						// then current voxel is likely a false positive
						double volume_difference = neighbor->volume_fraction_inside - voxel->volume_fraction_inside;
						
						if (volume_difference > 0.3) { // Neighbor has 30%+ more volume inside
							should_remove_tissue = true;
							break;
						}
					}
				}
				
				// Remove material from false positive voxels
				if (should_remove_tissue) {
					voxel->material = nullptr;
					voxel->is_boundary_voxel = false; // No longer a boundary since no material
					false_positives_removed++;
				}
			}
		}
	}
	
	if (config_.verbose()) {
		std::cout << "Surface disambiguation removed " << false_positives_removed << " false positive voxels." << std::endl;
	}
}

void Medium::identify_surface_voxels() {
	int surface_voxel_count = 0;
	const glm::uvec3& dimensions = volume_.dimensions();
	
	for (uint32_t z = 0; z < dimensions.z; z++) {
		for (uint32_t y = 0; y < dimensions.y; y++) {
			for (uint32_t x = 0; x < dimensions.x; x++) {
				Voxel* voxel = volume_.at(x, y, z);
				
				// SIMPLE SURFACE DETECTION: A voxel is surface if:
				// 1. It has material (part of medium)
				// 2. It's a boundary voxel (intersects geometry boundary)
				// 3. OR it has an empty neighbor (adjacent to outside)
				
				bool is_surface = false;
				
				if (voxel->material) {
					// Method 1: Boundary voxel (intersects geometry surface)
					if (voxel->is_boundary_voxel) {
						is_surface = true;
					}
					
					// Method 2: Has empty neighbors (traditional surface detection)
					if (!is_surface) {
						for (int dz = -1; dz <= 1 && !is_surface; dz++) {
							for (int dy = -1; dy <= 1 && !is_surface; dy++) {
								for (int dx = -1; dx <= 1 && !is_surface; dx++) {
									if (dx == 0 && dy == 0 && dz == 0) continue; // Skip self
									
									int nx = (int)x + dx;
									int ny = (int)y + dy;
									int nz = (int)z + dz;
									
									// Check bounds - if neighbor is outside grid, this is surface
									if (nx < 0 || ny < 0 || nz < 0 
									|| 	nx >= (int)dimensions.x
									|| 	ny >= (int)dimensions.y 
									|| 	nz >= (int)dimensions.z) {
										is_surface = true;
										break;
									}
									
									// Check if neighbor has material
									Voxel* neighbor = volume_.at(nx, ny, nz);
									if (!neighbor->material) {
										is_surface = true;
										break;
									}
								}
							}
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
	
	if (config_.verbose()) {
		std::cout << "Identified " << surface_voxel_count << " surface voxels." << std::endl;
	}
}

void Medium::reset_simulation_data() {
	// Reset voxel data
	for (auto& voxel_ptr : volume_) {
		if (voxel_ptr) {
			voxel_ptr->absorption = 0.0;
			voxel_ptr->emittance = 0.0;
			// Reset directional emittance fields as well
			voxel_ptr->emittance_reflected = 0.0;
			voxel_ptr->emittance_transmitted = 0.0;
			voxel_ptr->emittance_diffuse = 0.0;
		}
	}

	// Reset metrics
	metrics_.reset();
}

void Medium::normalize() {
	// Normalize raw accumulator parameters
	double divisor = static_cast<double>(config_.num_photons()) * static_cast<double>(config_.num_sources());
	metrics_.normalize_raw_values(divisor);

	// Normalize voxel data
	for (const auto& voxel_ptr : volume_) {
		if (voxel_ptr && voxel_ptr->material) {
			voxel_ptr->absorption /= config_.num_photons() * config_.num_sources();
			voxel_ptr->emittance /= config_.num_photons() * config_.num_sources();
			
			// CRITICAL FIX: Also normalize the directional emittance fields used by energy conservation
			voxel_ptr->emittance_reflected /= config_.num_photons() * config_.num_sources();
			voxel_ptr->emittance_transmitted /= config_.num_photons() * config_.num_sources();
			voxel_ptr->emittance_diffuse /= config_.num_photons() * config_.num_sources();
		}
	}
}

void Medium::write_results() const {
	std::string str_abs = App::get_output_path("absorption.out");
	std::string str_emi = App::get_output_path("emittance.out");

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
				if (voxel && voxel->material) {
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
	// ROBUST BOUNDARY HANDLING: Use larger epsilon tolerance for boundary inclusion
	static const double BOUNDARY_EPSILON = 1e-6;
	
	// CRITICAL: Check if photon is actually exiting the medium completely
	// Don't try to find voxels for photons that have truly exited
	if (!contains_point(position)) {
		// Position is completely outside medium - this is a true exit
		return nullptr;
	}
	
	// Position is within medium boundaries, but may need voxel coordinate adjustment
	if (!bounds_.includes(position)) {
		// ROBUST BOUNDARY HANDLING: Try epsilon-nudged position
		glm::dvec3 nudged_pos = position;
		bool position_fixed = false;
		
		// If position is very close to bounds, nudge it inside
		if (position.x <= bounds_.min_bounds.x + BOUNDARY_EPSILON) {
			nudged_pos.x = bounds_.min_bounds.x + BOUNDARY_EPSILON;
			position_fixed = true;
		} else if (position.x >= bounds_.max_bounds.x - BOUNDARY_EPSILON) {
			nudged_pos.x = bounds_.max_bounds.x - BOUNDARY_EPSILON;
			position_fixed = true;
		}
		
		if (position.y <= bounds_.min_bounds.y + BOUNDARY_EPSILON) {
			nudged_pos.y = bounds_.min_bounds.y + BOUNDARY_EPSILON;
			position_fixed = true;
		} else if (position.y >= bounds_.max_bounds.y - BOUNDARY_EPSILON) {
			nudged_pos.y = bounds_.max_bounds.y - BOUNDARY_EPSILON;
			position_fixed = true;
		}
		
		if (position.z <= bounds_.min_bounds.z + BOUNDARY_EPSILON) {
			nudged_pos.z = bounds_.min_bounds.z + BOUNDARY_EPSILON;
			position_fixed = true;
		} else if (position.z >= bounds_.max_bounds.z - BOUNDARY_EPSILON) {
			nudged_pos.z = bounds_.max_bounds.z - BOUNDARY_EPSILON;
			position_fixed = true;
		}
		
		if (position_fixed && bounds_.includes(nudged_pos)) {
			position = nudged_pos; // Update the reference position
		} else {
			return nullptr;
		}
	}

	// Calculate relative position from minimum bounds
	double dx = position.x - bounds_.min_bounds.x;
	double dy = position.y - bounds_.min_bounds.y;
	double dz = position.z - bounds_.min_bounds.z;

	// ROBUST COORDINATE CALCULATION: Use consistent floor() with boundary epsilon
	uint32_t ix = static_cast<uint32_t>(std::floor(dx / config_.vox_size() + BOUNDARY_EPSILON));
	uint32_t iy = static_cast<uint32_t>(std::floor(dy / config_.vox_size() + BOUNDARY_EPSILON));
	uint32_t iz = static_cast<uint32_t>(std::floor(dz / config_.vox_size() + BOUNDARY_EPSILON));



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

	// if voxel does not have a material, it is outside the medium
	return voxel->material ? voxel : nullptr;
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

	// Find closest intersection across all layers
	for (const auto& layer : layers_) {
		Triangle temp_triangle{};
		glm::dvec3 temp_intersection{};
		
		// Use the layer's triangles span to find the first intersection
		std::span<const Triangle> triangles = layer.triangles();
		double distance = ray.intersect_first_triangle(
			std::span<Triangle>(const_cast<Triangle*>(triangles.data()), triangles.size()),
			temp_intersection, temp_triangle
		);
		
		if (distance > 0 && distance < closest_distance) {
			closest_distance = distance;
			intersection_point = temp_intersection;
			intersected_triangle = temp_triangle;
		}
	}

    source.intersect = intersection_point;
	source.triangle = intersected_triangle;

	return closest_distance;
}

bool Medium::contains_point(const glm::dvec3& point) const {
	// ROBUST BOUNDARY HANDLING: Use epsilon tolerance for boundary inclusion
	static const double BOUNDARY_EPSILON = 1e-8;
	
	// Test if the point is inside any layer's geometry
	for (const auto& layer : layers_) {
		// First check if point is strictly inside
		if (layer.contains_point(point)) {
			return true;
		}
		
		// If not strictly inside, try epsilon-nudged positions for boundary tolerance
		// This handles floating point precision issues at exact boundaries
		for (double dx : {-BOUNDARY_EPSILON, 0.0, BOUNDARY_EPSILON}) {
			for (double dy : {-BOUNDARY_EPSILON, 0.0, BOUNDARY_EPSILON}) {
				for (double dz : {-BOUNDARY_EPSILON, 0.0, BOUNDARY_EPSILON}) {
					if (dx == 0.0 && dy == 0.0 && dz == 0.0) continue; // Skip original point
					
					glm::dvec3 nudged_point = point + glm::dvec3(dx, dy, dz);
					if (layer.contains_point(nudged_point)) {
						return true;
					}
				}
			}
		}
	}

	return false;
}

void Medium::detect_external_surface_voxels() {
	const glm::uvec3& dimensions = volume_.dimensions();
	int surface_voxel_count = 0;
	
	if (config_.verbose()) {
		std::cout << "Detecting external surface voxels..." << std::endl;
	}
	
	// For each voxel with material, check if it has neighbors without material
	for (uint32_t z = 0; z < dimensions.z; z++) {
		for (uint32_t y = 0; y < dimensions.y; y++) {
			for (uint32_t x = 0; x < dimensions.x; x++) {
				Voxel* voxel = volume_.at(x, y, z);
				
				// Only check voxels that have material
				if (!voxel->material) {
					continue;
				}
				
				bool is_external_surface = false;
				
				// Check all 26 neighbors (including face, edge, and corner neighbors)
				for (int dz = -1; dz <= 1 && !is_external_surface; dz++) {
					for (int dy = -1; dy <= 1 && !is_external_surface; dy++) {
						for (int dx = -1; dx <= 1 && !is_external_surface; dx++) {
							if (dx == 0 && dy == 0 && dz == 0) continue; // Skip self
							
							int nx = (int)x + dx;
							int ny = (int)y + dy;
							int nz = (int)z + dz;
							
							// Check bounds - if neighbor is outside grid, this is external surface
							if (nx < 0 || ny < 0 || nz < 0 
							|| 	nx >= (int)dimensions.x
							|| 	ny >= (int)dimensions.y 
							|| 	nz >= (int)dimensions.z) {
								is_external_surface = true;
								break;
							}
							
							// Check if neighbor has material - if not, this is external surface
							Voxel* neighbor = volume_.at(nx, ny, nz);
							if (!neighbor->material) {
								is_external_surface = true;
								break;
							}
						}
					}
				}
				
				// Set surface flag
				voxel->is_surface_voxel = is_external_surface;
				if (is_external_surface) {
					surface_voxel_count++;
				}
			}
		}
	}
	
	if (config_.verbose()) {
		std::cout << "External surface detection completed. Found " << surface_voxel_count << " surface voxels." << std::endl;
	}
}
