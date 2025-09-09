#include "simulator.hpp"

#include <algorithm>
#include <cfloat>
#include <cstdlib>
#include <ctime>
#include <fstream>
#include <iostream>
#include <limits>
#include <memory>
#include <numbers>
#include <sstream>

#include "math/math.hpp"
#include "math/random.hpp"
#include "math/ray.hpp"

/***********************************************************
 * Simulator constructor.
 ***********************************************************/
Simulator::Simulator() : mcml_random(std::make_shared<Random>()), mcml_weight_threshold(1e-4) {
	// Modern C++20: Use default initialization instead of explicit construction
	paths.clear();
	layers.clear();
	photons.clear();
	tissues.clear();
	sources.clear();

	// Reserve space for common use cases
	paths.reserve(1000);
	photons.reserve(10000);
	layers.reserve(10);
	tissues.reserve(5);
	sources.reserve(5);

	// Volume manages its own memory automatically

	// Initialize MCML random number generator
	mcml_random->seed(static_cast<int>(std::time(nullptr)));
}

/***********************************************************
 * Simulator destructor.
 ***********************************************************/
Simulator::~Simulator() {
	// Volume automatically cleans up its own memory
	// Smart pointers will automatically clean up vertex memory
}

/***********************************************************
 * INITIALIZATION
 ***********************************************************/

/***********************************************************
 * Parse the input file and initializes the data structures.
 ***********************************************************/
bool Simulator::initialize(std::string file) {
	std::cout << "Initializing Photron" << std::endl;

	// read and parse input configuration file
	std::cout << "Parsing configuration file: " << file << std::endl;
	if (!parse(file)) {
		std::cerr << "An error occurred while parsing the input file." << std::endl;
		return false;
	}
	std::cout << "Configuration parsed successfully." << std::endl;

	// initialize voxel grid
	if (!initialize_grid()) {
		std::cerr << "An error occurred while initializing the voxel grid." << std::endl;
		return false;
	}
	std::cout << "Voxel grid initialized successfully." << std::endl;

	// initialize other data structures
	if (!initialize_data()) {
		std::cerr << "An error occurred while initializing the data structures." << std::endl;
		return false;
	}
	std::cout << "Data structures initialized successfully." << std::endl;

	// reset all simulation data from previous runs (after voxels are created)
	reset_simulation_data();

	// voxelize the geometry
	if (!voxelize_layers()) {
		std::cerr << "An error occurred during geometry voxelization." << std::endl;
		return false;
	}
	std::cout << "Geometry voxelization completed successfully." << std::endl;
	
	// Identify surface voxels (voxels at the outer boundary of the voxel cluster)
	identify_surface_voxels();
	std::cout << "Surface voxel identification completed." << std::endl;
	
	return true;
}

/***********************************************************
 *	INITIALIZATION SUBROUTINES
 ***********************************************************/

/***********************************************************
 * Read and parse the given configuration file.
 ***********************************************************/
bool Simulator::parse(const std::string& fconfig) {
	// Clear existing data structures for reinitialization
	layers.clear();
	tissues.clear();
	sources.clear();
	photons.clear();
	paths.clear();

	// Use Config class to parse the entire configuration file
	if (!config.parse_config_file(fconfig)) {
		std::cerr << "Failed to parse configuration file." << std::endl;
		return false;
	}

	// Extract parsed data from Config
	sources = config.sources_vector();
	tissues = config.tissues_vector();
	layers = config.move_layers();

	return true;
}

/***********************************************************
 * Initialize the voxel grid.
 ***********************************************************/
bool Simulator::initialize_grid() {
	// Volume will be initialized after we calculate the bounds and dimensions

	// Initialize bounds to extreme values for proper min/max calculation
	bounds.min_bounds = glm::dvec3(DBL_MAX);
	bounds.max_bounds = glm::dvec3(-DBL_MAX);

	// compute grid boundary extent with padding to ensure full coverage
	for (const auto& layer : layers) {
		// see if a vertex denotes a new boundary
		for (const auto& triangle : layer.mesh) {
			glm::dvec3 v0 = triangle.v0();
			glm::dvec3 v1 = triangle.v1();
			glm::dvec3 v2 = triangle.v2();

			// get the maximum value among the previous maximum or new vertices
			bounds.min_bounds.x = std::min({bounds.min_bounds.x, v0.x, v1.x, v2.x}); // left (-x)
			bounds.max_bounds.x = std::max({bounds.max_bounds.x, v0.x, v1.x, v2.x}); // right (+x)

			bounds.min_bounds.y = std::min({bounds.min_bounds.y, v0.y, v1.y, v2.y}); // top (+y)
			bounds.max_bounds.y = std::max({bounds.max_bounds.y, v0.y, v1.y, v2.y}); // bottom (-y)

			bounds.min_bounds.z = std::min({bounds.min_bounds.z, v0.z, v1.z, v2.z}); // front (+z)
			bounds.max_bounds.z = std::max({bounds.max_bounds.z, v0.z, v1.z, v2.z}); // rear (-z)
		}
	}

	// Add padding to ensure voxels fully cover the mesh geometry
	double padding = config.vox_size() * 2.0; // Two voxel widths of padding
	bounds.min_bounds -= glm::dvec3(padding);
	bounds.max_bounds += glm::dvec3(padding);

	// check for inconsistency and zero width/height/depth
	if (bounds.min_bounds.x >= bounds.max_bounds.x || bounds.min_bounds.y >= bounds.max_bounds.y
		|| bounds.min_bounds.z >= bounds.max_bounds.z) {
		std::cerr << "Invalid bounds detected. Possibly no geometry data or parsing error." << std::endl;
		return false;
	}

	// Note: width, height, depth are now computed properties - no need to set them explicitly

	// Get the size vector and extract individual dimensions
	auto size_vec = bounds.size();

	// number of voxels in each dimension (as float)
	float nx = static_cast<float>(size_vec.x / config.vox_size());
	float ny = static_cast<float>(size_vec.y / config.vox_size());
	float nz = static_cast<float>(size_vec.z / config.vox_size());

	// number of voxels in each dimension
	config.set_nx(static_cast<uint32_t>(nx));
	config.set_ny(static_cast<uint32_t>(ny));
	config.set_nz(static_cast<uint32_t>(nz));

	// total number of voxels
	config.set_num_voxels(static_cast<uint64_t>(config.nx()) * static_cast<uint64_t>(config.ny())
						  * static_cast<uint64_t>(config.nz()));

	// check for voxel sizes that are too large
	if (config.vox_size() > size_vec.x || config.vox_size() > size_vec.y || config.vox_size() > size_vec.z) {
		return false;
	}

	// check grid dimensions
	if (config.num_voxels() < 1 || config.nx() < 1 || config.ny() < 1 || config.nz() < 1) {
		return false;
	}

	// Initialize the voxel grid with proper dimensions
	voxel_grid = Volume(config.vox_size(), config.nx(), config.ny(), config.nz());

	return true;
}

/***********************************************************
 * Initialize configuration properties.
 ***********************************************************/
bool Simulator::initialize_data() {
	// initialize config properties
	config.set_num_layers(static_cast<uint64_t>(layers.size()));
	config.set_num_sources(static_cast<uint64_t>(sources.size()));

	// error checking
	if (config.num_layers() < 1 || config.num_voxels() < 1) {
		return false;
	}

	// check for duplicate layers
	for (uint32_t i = 1; i < layers.size(); ++i) {
		if (layers[i - 1] == layers[i]) {
			return false;
		}
	}

	// check for duplicate tissues
	for (uint32_t i = 1; i < tissues.size(); ++i) {
		if (tissues[i - 1] == tissues[i]) {
			return false;
		}
	}

	// check if layer's tissue id is out of range
	for (const auto& layer : layers) {
		if (layer.tissue_id >= tissues.size()) {
			return false;
		}
	}

	// associate light sources with their geometric intersections
	for (auto& source : sources) {
		// intersection point and triangle that is hit
		glm::dvec3 intersect;
		Triangle triangle;

		// find intersection of ray from this source with geometry (point, triangle, normal)
		Ray ray = Ray(source.origin, source.direction);
		double distance = ray.intersect_first_triangle_from_layers(layers, intersect, triangle);
		if (distance == std::numeric_limits<double>::max()) {
			// Calculate intersection with the top plane of the geometry
			// The top face is at Y=0, source is at Y=0.2 going down
			double t = (0.0 - source.origin.y) / source.direction.y; // intersect Y=0 plane
			if (t > 0) {
				intersect.x = source.origin.x + t * source.direction.x;
				intersect.y = 0.0;
				intersect.z = source.origin.z + t * source.direction.z;

				// Use the first triangle from the first layer as the intersected triangle
				if (!layers.empty() && !layers[0].mesh.empty()) {
					triangle = layers[0].mesh[0];
				}
			}
			else {
				std::cerr << "Error: Cannot calculate intersection for source " << source.id << std::endl;
				return false;
			}
		}

		source.intersect = intersect;
		source.triangle = triangle;
	}

	// initialize photons
	for (uint64_t i = 0; i < config.num_photons(); ++i) {
		photons.emplace_back(i);
	}

	return true;
}

/***********************************************************
 * Voxelize the geometry.
 ***********************************************************/
bool Simulator::voxelize_layers() {
	// Use point-in-mesh testing for each voxel center
	// This is more accurate than ray-casting for complex geometries

	uint32_t nx = config.nx();
	uint32_t ny = config.ny();
	uint32_t nz = config.nz();
	double vox_size = config.vox_size();

	std::cout << "Voxelizing geometry using point-containment method..." << std::endl;
	std::cout << "Grid dimensions: " << nx << "x" << ny << "x" << nz << " (total: " << (nx * ny * nz) << " voxels)"
			  << std::endl;
	std::cout << "Voxel size: " << vox_size << std::endl;
	std::cout << "Grid bounds: min(" << bounds.min_bounds.x << "," << bounds.min_bounds.y << "," << bounds.min_bounds.z
			  << ") max(" << bounds.max_bounds.x << "," << bounds.max_bounds.y << "," << bounds.max_bounds.z << ")"
			  << std::endl;

	int total_voxels_assigned = 0;
	int partial_voxels = 0;

	for (uint32_t iz = 0; iz < nz; iz++) {
		for (uint32_t iy = 0; iy < ny; iy++) {
			for (uint32_t ix = 0; ix < nx; ix++) {
				// Calculate voxel center in world coordinates
				glm::dvec3 voxel_center =
					glm::dvec3(bounds.min_bounds.x + (ix + 0.5) * vox_size, bounds.min_bounds.y + (iy + 0.5) * vox_size,
							   bounds.min_bounds.z + (iz + 0.5) * vox_size);

				// Calculate voxel corners to check for partial intersection
				glm::dvec3 voxel_min =
					glm::dvec3(bounds.min_bounds.x + ix * vox_size, bounds.min_bounds.y + iy * vox_size,
							   bounds.min_bounds.z + iz * vox_size);
				glm::dvec3 voxel_max = voxel_min + glm::dvec3(vox_size);

				// Check if voxel intersects with geometry
				bool center_inside = false;
				bool any_corner_inside = false;

				// Test center and corners
				for (const auto& layer : layers) {
					if (is_point_inside_layer_mesh(voxel_center, layer)) {
						center_inside = true;
						break;
					}
				}

				// Test corners to detect partial intersection
				std::vector<glm::dvec3> corners = {voxel_min,
												   {voxel_max.x, voxel_min.y, voxel_min.z},
												   {voxel_min.x, voxel_max.y, voxel_min.z},
												   {voxel_max.x, voxel_max.y, voxel_min.z},
												   {voxel_min.x, voxel_min.y, voxel_max.z},
												   {voxel_max.x, voxel_min.y, voxel_max.z},
												   {voxel_min.x, voxel_max.y, voxel_max.z},
												   voxel_max};

				for (const auto& corner : corners) {
					for (const auto& layer : layers) {
						if (is_point_inside_layer_mesh(corner, layer)) {
							any_corner_inside = true;
							break;
						}
					}
					if (any_corner_inside)
						break;
				}

				// Assign tissue if center is inside OR any corner is inside (partial voxel)
				if (center_inside || any_corner_inside) {
					Voxel* voxel = voxel_grid(ix, iy, iz);

					// Find which specific layer contains this voxel (check center point)
					// IMPORTANT: Process layers in reverse order so inner layers override outer layers
					for (auto it = layers.rbegin(); it != layers.rend(); ++it) {
						const auto& layer = *it;
						if (is_point_inside_layer_mesh(voxel_center, layer)) {
							voxel->tissue = &tissues[layer.tissue_id];
							break; // Use the innermost layer that contains the center
						}
					}

					// If center point didn't hit any layer but corners did,
					// use the innermost layer that contains any corner
					if (voxel->tissue == nullptr) {
						for (const auto& corner : corners) {
							for (auto it = layers.rbegin(); it != layers.rend(); ++it) {
								const auto& layer = *it;
								if (is_point_inside_layer_mesh(corner, layer)) {
									voxel->tissue = &tissues[layer.tissue_id];
									break;
								}
							}
							if (voxel->tissue != nullptr)
								break;
						}
					}

					// Compute volume fractions for proper boundary physics
					if (center_inside && any_corner_inside) {
						// Check if this is actually a boundary voxel
						bool all_corners_inside = true;
						for (const auto& corner : corners) {
							bool corner_inside = false;
							for (const auto& layer : layers) {
								if (is_point_inside_layer_mesh(corner, layer)) {
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
							voxel->volume_fraction_inside =
								voxel_grid.compute_volume_fraction_inside_fast(voxel_min, voxel_max, *this, 2);
							voxel->volume_fraction_outside = 1.0 - voxel->volume_fraction_inside;
							voxel->is_boundary_voxel = true;
						}
					}
					else if (center_inside) {
						// Center inside but some corners outside - partial
						voxel->volume_fraction_inside =
							voxel_grid.compute_volume_fraction_inside_fast(voxel_min, voxel_max, *this, 2);
						voxel->volume_fraction_outside = 1.0 - voxel->volume_fraction_inside;
						voxel->is_boundary_voxel = true;
					}
					else {
						// Only corners inside - partial
						voxel->volume_fraction_inside =
							voxel_grid.compute_volume_fraction_inside_fast(voxel_min, voxel_max, *this, 2);
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

/***********************************************************
 * Identify voxels that are at the outer surface of the voxel cluster.
 * This includes both voxels with empty neighbors AND voxels that graze the surface.
 ***********************************************************/
void Simulator::identify_surface_voxels() {
	int surface_voxel_count = 0;
	const glm::uvec3& dimensions = voxel_grid.dimensions();
	
	for (uint32_t z = 0; z < dimensions.z; z++) {
		for (uint32_t y = 0; y < dimensions.y; y++) {
			for (uint32_t x = 0; x < dimensions.x; x++) {
				Voxel* voxel = voxel_grid.at(x, y, z);
				
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
								if (dx == 0 && dy == 0 && dz == 0) continue; // Skip self
								
								int nx = (int)x + dx;
								int ny = (int)y + dy;
								int nz = (int)z + dz;
								
								// Check bounds
								if (nx < 0 || ny < 0 || nz < 0 || 
									nx >= (int)dimensions.x || 
									ny >= (int)dimensions.y || 
									nz >= (int)dimensions.z) {
									// Out of bounds = empty neighbor
									has_empty_neighbor = true;
									break;
								}
								
								// Check if neighbor has no tissue (is empty)
								Voxel* neighbor = voxel_grid.at(nx, ny, nz);
								if (!neighbor->tissue) {
									has_empty_neighbor = true;
									break;
								}
							}
							if (has_empty_neighbor) break;
						}
						if (has_empty_neighbor) break;
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
					double voxel_size = voxel_grid.voxel_size();
					double tolerance = voxel_size * 0.1; // 10% of voxel size for floating-point safety
					
					// Check if any mesh triangle intersects this voxel
					for (const auto& layer : layers) {
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
						if (is_surface) break;
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

/***********************************************************
 * Run the Monte Carlo photon transport simulation.
 ***********************************************************/
void Simulator::simulate() {
	std::cout << "Running Monte Carlo simulation" << std::endl;

	metrics.start_clock();

	// For each light source
	for (auto& source : sources) {
		// Compute specular reflection
		specular_reflection(source);

		// for each photon
		uint32_t p = 0;
		for (auto& photon : photons) {
			// progress report
			if (config.progress() && ((p + 1) % 1000) == 0) {
				std::cout << "Photon " << p + 1 << "/" << config.num_photons() << std::endl;
			}

			// launch the photon (create a new path)
			launch(photons[p], source);

			// Safety mechanism to prevent infinite photon loops
			int photon_iteration_counter = 0;
			const int max_photon_iterations = 1000000;

			while (photons[p].alive) {
				photon_iteration_counter++;
				if (photon_iteration_counter > max_photon_iterations) {
					std::cerr << "Warning: Photon " << photons[p].id << " exceeded maximum iterations, terminating."
							  << std::endl;
					photons[p].alive = false;
					break;
				}
				step_size(photon); // Set new step size
				transfer(photon);  // Propagate photon through the medium in substeps
				roulette(photon);  // Determine photon termination
				
				// Only scatter if photon hasn't just crossed a boundary
				// Scattering should only occur in bulk medium, not at interfaces
				if (!photon.cross) {
					scatter(photon);   // Scatter photon into a new direction
				}
			}

			p++;
		}
	}

	// Normalize physical quantities
	normalize();

	metrics.stop_clock();
	metrics.collect_data(record.total_absorption, record.specular_reflection, record.diffuse_reflection,
						 record.specular_transmission, record.diffuse_transmission);
	metrics.write_to_file();
	metrics.print_report();
}

/***********************************************************
 * Simulate a single additional photon for interactive use.
 ***********************************************************/
void Simulator::simulate_single_photon() {
	// Use the first source (there should be at least one)
	if (sources.empty()) {
		std::cerr << "Error: No light sources available for single photon simulation" << std::endl;
		return;
	}

	Source& source = sources[0];

	// Create a new photon with unique ID
	uint64_t new_photon_id = photons.size();
	Photon new_photon(new_photon_id);

	// Launch the photon
	launch(new_photon, source);

	// Safety mechanism to prevent infinite photon loops
	int photon_iteration_counter = 0;
	const int max_photon_iterations = 1000000;

	while (new_photon.alive) {
		photon_iteration_counter++;
		if (photon_iteration_counter > max_photon_iterations) {
			std::cerr << "Warning: Single photon exceeded maximum iterations, terminating." << std::endl;
			new_photon.alive = false;
			break;
		}

		step_size(new_photon); // Set new step size
		transfer(new_photon);  // Propagate photon through the medium in substeps
		roulette(new_photon);  // Determine photon termination
		scatter(new_photon);   // Scatter photon into a new direction
	}

	// Add the completed photon to the photons vector for rendering
	photons.push_back(new_photon);
}

/***********************************************************
 * Set up photon properties for tracing.
 ***********************************************************/
void Simulator::launch(Photon& photon, Source& source) {
	photon.alive = true;
	photon.weight = 1.0 - record.specular_reflection;
	photon.source = source;
	photon.direction = source.direction;
	photon.position = source.intersect;
	photon.voxel = voxel_at(photon.position);

	// create vertices for new light path
	auto light = std::make_shared<PhotonNode>(source.origin, photon.weight);
	auto intersection = std::make_shared<PhotonNode>(source.intersect, photon.weight);
	auto reflection = std::make_shared<PhotonNode>(move(source.intersect, source.specular_direction, 0.1),
												   record.specular_reflection);

	light->next = intersection;      // intersection vertex/node
	intersection->prev = light;      // light source origin
	intersection->emit = reflection; // specular reflection

	// use intersection point as head
	PhotonPath path = PhotonPath(static_cast<long>(photon.id), intersection);
	paths.push_back(path);

	metrics.add_vertex(photon.position.x, photon.position.y, photon.position.z);
}

/***********************************************************
 * Set dimensionless step size for next segment of the
 * random walk using MCML 3.0.0 algorithm.
 ***********************************************************/
void Simulator::step_size(Photon& photon) {
	generate_step_size(photon);
	metrics.add_step_size(photon.step / photon.mu_s());
}

/***********************************************************
 * Transfer a photon through a voxelized medium with
 * individual substeps.
 ***********************************************************/
void Simulator::transfer(Photon& photon) {
	/*
	 * set substep (max = distance to boundary)
	 * deposit weight in current voxel
	 * if (photon crosses boundary)
	 *     if (photon goes outside medium)
	 *         record partial transmission
	 *     else if (photon moves to differing refractive indexed media)
	 *         reflect from or transmit across the boundary
	 *     else if (photon moves to equal refractive indexed media)
	 *         continue normal propagation
	 * decrease step size by traveled distance
	 */

	// Safety mechanism to prevent infinite loops
	int substep_counter = 0;
	const int max_substeps = 100000;

	while (photon.step >= 1E-10 && photon.alive) {
		substep_counter++;
		if (substep_counter > max_substeps) {
			std::cerr << "Warning: Photon exceeded maximum substeps, terminating." << std::endl;
			photon.alive = false;
			break;
		}

		// set substep - this will handle mesh boundary detection
		sub_step(photon);

		// deposit weight
		deposit(photon);

		// possibly cross boundary
		if (photon.cross) {
			cross(photon);
		}
		else {
			photon.position = move(photon.position, photon.direction, photon.sub_step);
		}

		// prevent errors due to crossing to ambient medium
		if (!photon.voxel) {
			photon.alive = false;
			
			// Find the appropriate boundary/surface voxel to record emittance
			Voxel* emittance_voxel = nullptr;
			
			// Calculate tolerance based on photon transport step size to handle numerical precision
			double step_tolerance = photon.sub_step * 0.1; // 10% of last substep size
			double min_tolerance = voxel_grid.voxel_size() * 0.01; // Minimum 1% of voxel size
			double tolerance = std::max(step_tolerance, min_tolerance);
			
			// First try: voxel containing the intersection point
			Voxel* intersect_voxel = voxel_at(photon.intersect);
			if (intersect_voxel && (intersect_voxel->is_boundary_voxel || intersect_voxel->is_surface_voxel)) {
				emittance_voxel = intersect_voxel;
			} 
			
			// Second try: check voxels within tolerance of both intersection point and photon position
			if (!emittance_voxel) {
				// Search in a small region around both intersection and position points
				std::vector<glm::dvec3> search_points = {
					photon.intersect,
					photon.position,
					(photon.intersect + photon.position) * 0.5 // midpoint
				};
				
				for (const auto& point : search_points) {
					// Check the point itself and nearby voxels
					for (double dx = -tolerance; dx <= tolerance; dx += voxel_grid.voxel_size() * 0.5) {
						for (double dy = -tolerance; dy <= tolerance; dy += voxel_grid.voxel_size() * 0.5) {
							for (double dz = -tolerance; dz <= tolerance; dz += voxel_grid.voxel_size() * 0.5) {
								glm::dvec3 test_point = point + glm::dvec3(dx, dy, dz);
								Voxel* candidate = voxel_at(test_point);
								if (candidate && (candidate->is_boundary_voxel || candidate->is_surface_voxel)) {
									emittance_voxel = candidate;
									break;
								}
							}
							if (emittance_voxel) break;
						}
						if (emittance_voxel) break;
					}
					if (emittance_voxel) break;
				}
			}
			
			// Third try: the previous voxel if it's boundary/surface
			if (!emittance_voxel && photon.prev_voxel && (photon.prev_voxel->is_boundary_voxel || photon.prev_voxel->is_surface_voxel)) {
				emittance_voxel = photon.prev_voxel;
			}
			
			// Last resort: use previous voxel even if not boundary/surface (with warning)
			if (!emittance_voxel && photon.prev_voxel) {
				emittance_voxel = photon.prev_voxel;
				std::cerr << "Warning: Using non-boundary voxel for emittance (numerical precision issue)" << std::endl;
			}
			
			// Ensure we have a voxel to record emittance
			if (emittance_voxel) {
				photon.voxel = emittance_voxel;
				radiate(photon, photon.direction, photon.weight);
			} else {
				std::cerr << "ERROR: Could not find any voxel for emittance recording!" << std::endl;
			}
		}

		// update step size
		photon.step -= (photon.sub_step * photon.mu_s());

		metrics.add_vertex(photon.position.x, photon.position.y, photon.position.z);
	}
}

/***********************************************************
 * Set the photon's next substep and initialize the given
 * intersection point and voxel normal if it crosses the
 * voxel boundary.
 ***********************************************************/
void Simulator::sub_step(Photon& photon) {
	// Get voxel boundaries
	Cuboid box = voxel_corners(photon.voxel);

	// Check if photon is exactly on a voxel boundary and nudge if necessary
	bool on_boundary = false;
	glm::dvec3 adjusted_position = photon.position;

	// Check each axis for boundary conditions
	double x_diff_min = std::abs(photon.position.x - box.min_point().x);
	double x_diff_max = std::abs(photon.position.x - box.max_point().x);
	double y_diff_min = std::abs(photon.position.y - box.min_point().y);
	double y_diff_max = std::abs(photon.position.y - box.max_point().y);
	double z_diff_min = std::abs(photon.position.z - box.min_point().z);
	double z_diff_max = std::abs(photon.position.z - box.max_point().z);

	if (x_diff_min < MathConstants::BOUNDARY_EPSILON) {
		adjusted_position.x = box.min_point().x + MathConstants::BOUNDARY_EPSILON;
		on_boundary = true;
	}
	else if (x_diff_max < MathConstants::BOUNDARY_EPSILON) {
		adjusted_position.x = box.max_point().x - MathConstants::BOUNDARY_EPSILON;
		on_boundary = true;
	}

	if (y_diff_min < MathConstants::BOUNDARY_EPSILON) {
		adjusted_position.y = box.min_point().y + MathConstants::BOUNDARY_EPSILON;
		on_boundary = true;
	}
	else if (y_diff_max < MathConstants::BOUNDARY_EPSILON) {
		adjusted_position.y = box.max_point().y - MathConstants::BOUNDARY_EPSILON;
		on_boundary = true;
	}

	if (z_diff_min < MathConstants::BOUNDARY_EPSILON) {
		adjusted_position.z = box.min_point().z + MathConstants::BOUNDARY_EPSILON;
		on_boundary = true;
	}
	else if (z_diff_max < MathConstants::BOUNDARY_EPSILON) {
		adjusted_position.z = box.max_point().z - MathConstants::BOUNDARY_EPSILON;
		on_boundary = true;
	}

	// Create ray from (possibly adjusted) photon position and direction
	Ray ray = Ray(adjusted_position, photon.direction);

	double voxdist = ray.intersect_cuboid_internal(box, photon.intersect, photon.voxel_normal);

	// Check if no voxel intersections were found
	if (voxdist == std::numeric_limits<double>::max()) {
		std::cerr << "Critical error: ray-voxel intersection test failed." << std::endl;
		std::cerr << "  Ray origin: " << ray.origin().x << ", " << ray.origin().y << ", " << ray.origin().z
				  << std::endl;
		std::cerr << "  Ray direction: " << ray.direction().x << ", " << ray.direction().y << ", " << ray.direction().z
				  << std::endl;
		std::cerr << "  Ray direction length: " << glm::length(ray.direction()) << std::endl;
		std::cerr << "  Voxel bounds: min(" << box.min_point().x << ", " << box.min_point().y << ", "
				  << box.min_point().z << ") max(" << box.max_point().x << ", " << box.max_point().y << ", "
				  << box.max_point().z << ")" << std::endl;
		std::cerr << "  Photon position: (" << photon.position.x << ", " << photon.position.y << ", "
				  << photon.position.z << ")" << std::endl;
		if (on_boundary) {
			std::cerr << "  Boundary condition detected and position adjusted." << std::endl;
			std::cerr << "  Original position: (" << photon.position.x << ", " << photon.position.y << ", "
					  << photon.position.z << ")" << std::endl;
			std::cerr << "  Adjusted position: (" << adjusted_position.x << ", " << adjusted_position.y << ", "
					  << adjusted_position.z << ")" << std::endl;
		}
		std::cerr << "  Photon step: " << photon.step << std::endl;
		std::cerr << "  Photon weight: " << photon.weight << std::endl;

		// Instead of crashing, try to handle this gracefully
		std::cerr << "  Attempting graceful recovery..." << std::endl;

		// Terminate this photon instead of crashing the entire program
		photon.alive = false;
		photon.sub_step = 0.0;
		photon.cross = false;
		return;
	}

	// Compute free path for a substep (distance to next scattering event)
	double freepath = photon.step / photon.mu_s();

	// Check if photon is currently inside the mesh
	bool photon_inside_mesh = is_point_inside_geometry(photon.position);

	if (!photon_inside_mesh) {
		// Photon is outside mesh - it should have exited already
		photon.alive = false;
		return;
	}

	// Find the closest mesh boundary intersection
	double mesh_dist = std::numeric_limits<double>::max();
	glm::dvec3 mesh_intersection {0.0, 0.0, 0.0};
	glm::dvec3 mesh_normal {0.0, 0.0, 0.0};
	bool found_mesh_intersection = false;

	// Look for exit points from the mesh
	for (const auto& layer : layers) {
		for (const auto& triangle : layer.mesh) {
			Triangle triangle_copy = triangle;
			glm::dvec3 intersection;

			if (ray.intersect_triangle(triangle_copy, intersection)) {
				double dist = glm::distance(photon.position, intersection);

				// Only consider intersections that are forward and reasonably close
				if (dist > 1e-10 && dist < mesh_dist) {
					// Check if this intersection would take us outside the mesh
					glm::dvec3 test_point = intersection + photon.direction * 1e-6;
					if (!is_point_inside_geometry(test_point)) {
						mesh_dist = dist;
						mesh_intersection = intersection;
						mesh_normal = triangle.normal();
						found_mesh_intersection = true;
					}
				}
			}
		}
	}

	// Determine the substep based on the shortest distance
	if (found_mesh_intersection && mesh_dist <= freepath) {
		// Mesh boundary is closest - photon will exit the medium at geometry boundary
		photon.intersect = mesh_intersection;
		
		// Check if intersection is very close to a vertex (special case)
		const double vertex_threshold = 1e-6;
		bool is_vertex_intersection = false;
		glm::dvec3 averaged_normal = mesh_normal;
		
		// Check all triangles to see if intersection is near any vertices
		std::vector<glm::dvec3> vertex_normals;
		for (const auto& layer : layers) {
			for (const auto& triangle : layer.mesh) {
				// Check distance to each vertex
				double dist_v0 = glm::distance(mesh_intersection, triangle.v0());
				double dist_v1 = glm::distance(mesh_intersection, triangle.v1());
				double dist_v2 = glm::distance(mesh_intersection, triangle.v2());
				
				if (dist_v0 < vertex_threshold || dist_v1 < vertex_threshold || dist_v2 < vertex_threshold) {
					// This triangle shares the vertex - include its normal
					vertex_normals.push_back(triangle.normal());
					is_vertex_intersection = true;
				}
			}
		}
		
		// If vertex intersection, average the normals of adjacent faces
		if (is_vertex_intersection && !vertex_normals.empty()) {
			glm::dvec3 sum_normal(0.0);
			for (const auto& normal : vertex_normals) {
				sum_normal += normal;
			}
			averaged_normal = glm::normalize(sum_normal);
		}
		
		photon.voxel_normal = averaged_normal;
		photon.sub_step = mesh_dist;
		photon.cross = true;
	}
	else if (voxdist <= freepath) {
		// Voxel boundary is closer than scattering event but no mesh exit
		// This handles internal voxel transitions within the geometry
		photon.sub_step = voxdist;
		photon.cross = true;
	}
	else {
		// Scattering event occurs before any boundary
		photon.sub_step = freepath;
		photon.cross = false;
	}
}

/***********************************************************
 * Deposit some of the photon's weight into the geometry.
 ***********************************************************/
void Simulator::deposit(Photon& photon) {
	// Cancel if photon is outside of medium or doesn't have tissue
	if (!photon.voxel || !photon.voxel->tissue) {
		return;
	}

	// For boundary voxels, only deposit in the portion that's inside the geometry
	double effective_volume_fraction = 1.0;
	if (photon.voxel->is_boundary_voxel) {
		// Check if photon is actually inside the geometry at this position
		if (!is_point_inside_geometry(photon.position)) {
			// Photon is in the outside portion of a boundary voxel - no absorption
			return;
		}
		// Scale absorption by the volume fraction inside
		effective_volume_fraction = photon.voxel->volume_fraction_inside;
	}
	else {
		// Also check if photon is actually inside the geometry for non-boundary voxels
		if (!is_point_inside_geometry(photon.position)) {
			return;
		}
	}

	// deposited weight (scaled by effective volume fraction)
	double deltaw = photon.weight * (1 - std::exp(-photon.mu_a() * photon.sub_step)) * effective_volume_fraction;

	// update photon weight
	photon.weight -= deltaw;

	// assign deposited weight to voxel
	photon.voxel->absorption += deltaw;

	// update total absorption
	record.total_absorption += deltaw;
}

/***********************************************************
 * Determine the action to take when a photon is about to
 * traverse a voxel face. It can either cross to the ambient
 * medium, to another medium with a differing refractive
 * index, to another medium with the same refractive index,
 * or within the same medium.
 *
 * Reflection and transmission can be handled partially at
 * external-internal medium boundaries, but is always handled
 * as an all-or-none event at internal medium boundaries.
 *
 * If appropriate, the new photon direction, position and
 * voxel are computed.
 ***********************************************************/
void Simulator::cross(Photon& photon) {
	// Safety check - ensure photon has valid voxel and tissue
	if (!photon.voxel || !photon.voxel->tissue) {
		// Photon is outside medium - kill it
		photon.alive = false;
		return;
	}

	// Additional check: if photon is outside geometry, it should exit
	glm::dvec3 next_position = move(photon.position, photon.direction, photon.sub_step);
	bool next_in_geometry = is_point_inside_geometry(next_position);

	if (!next_in_geometry) {
		// Photon is exiting the medium - ensure intersect point is correct
		// If photon.intersect wasn't set by mesh detection, compute the actual exit point
		bool current_in_geometry = is_point_inside_geometry(photon.position);
		if (current_in_geometry) {
			// Find the exact exit point by intersecting with geometry
			Ray exit_ray(photon.position, photon.direction);
			glm::dvec3 true_exit_point {0.0, 0.0, 0.0};
			bool found_exit = false;
			double min_exit_dist = std::numeric_limits<double>::max();

			// Find closest geometry exit
			for (const auto& layer : layers) {
				for (const auto& triangle : layer.mesh) {
					Triangle triangle_copy = triangle;
					glm::dvec3 intersection;

					if (exit_ray.intersect_triangle(triangle_copy, intersection)) {
						double dist = glm::distance(photon.position, intersection);

						if (dist > 1e-10 && dist < min_exit_dist) {
							// Verify this intersection takes us outside
							glm::dvec3 test_point = intersection + photon.direction * 1e-6;
							if (!is_point_inside_geometry(test_point)) {
								min_exit_dist = dist;
								true_exit_point = intersection;
								found_exit = true;
							}
						}
					}
				}
			}

			// Use the correct exit point
			if (found_exit) {
				photon.intersect = true_exit_point;
			}
			else {
				// Fallback: interpolate to geometry boundary
				photon.intersect = photon.position + photon.direction * (photon.sub_step * 0.999);
			}
		}
		else {
			// Photon is already outside - this shouldn't happen
			photon.intersect = photon.position;
		}

		// Photon is exiting the medium - handle as ambient medium crossing
		double eta = config.ambient_eta();
		glm::dvec3 transmittance, reflectance;

		double reflection = internal_reflection(photon, eta, transmittance, reflectance);

		if (reflection == 0.0) {
			// Total transmission - photon exits
			photon.direction = transmittance;
			photon.alive = false;
			radiate(photon, transmittance, photon.weight);
		}
		else if (reflection == 1.0) {
			// Total internal reflection
			photon.direction = reflectance;
			photon.position = move_delta(photon.intersect, photon.direction);
			photon.voxel = voxel_at(photon.position);
		}
		else {
			// Partial reflection/transmission
			if (config.partial()) {
				radiate(photon, transmittance, photon.weight * (1.0 - reflection));
				photon.weight *= reflection;
				photon.direction = reflectance;
				photon.position = move_delta(photon.intersect, photon.direction);
				photon.voxel = voxel_at(photon.position);
			}
			else {
				// All-or-none
				if (mcml_random->next() > reflection) {
					photon.direction = transmittance;
					photon.alive = false;
					radiate(photon, transmittance, photon.weight);
				}
				else {
					photon.direction = reflectance;
					photon.position = move_delta(photon.intersect, photon.direction);
					photon.voxel = voxel_at(photon.position);
				}
			}
		}
		return;
	}

	// directions of transmission and reflection
	glm::dvec3 transmittance, reflectance;

	// First compute the transmission and reflection directions using current photon state
	double eta_ambient = config.ambient_eta();
	double temp_reflection = internal_reflection(photon, eta_ambient, transmittance, reflectance);

	// retrieve the refractive index of the medium that would be struck by the transmitted ray
	glm::dvec3 newpos = move_delta(photon.intersect, transmittance);
	Voxel* newvox = voxel_at(newpos);
	photon.prev_voxel = photon.voxel;

	// determine refractive index of the medium being entered
	double eta = (newvox == nullptr) ? config.ambient_eta() : newvox->tissue->eta;

	// Recalculate with correct refractive index
	temp_reflection = internal_reflection(photon, eta, transmittance, reflectance);

	// 1. crossing to ambient medium
	if (newvox == nullptr) {
		// Debug: Check if photon is exiting from an interior voxel (shouldn't happen!)
		if (photon.voxel && !photon.voxel->is_surface_voxel) {
			std::cerr << "ERROR: Photon attempting to exit from interior voxel at (" 
					  << photon.voxel->ix() << ", " << photon.voxel->iy() << ", " << photon.voxel->iz() 
					  << ") to ambient medium!" << std::endl;
			std::cerr << "  Photon position: (" << photon.position.x << ", " << photon.position.y << ", " << photon.position.z << ")" << std::endl;
			std::cerr << "  Intersection point: (" << photon.intersect.x << ", " << photon.intersect.y << ", " << photon.intersect.z << ")" << std::endl;
			std::cerr << "  New position: (" << newpos.x << ", " << newpos.y << ", " << newpos.z << ")" << std::endl;
			std::cerr << "  Direction: (" << transmittance.x << ", " << transmittance.y << ", " << transmittance.z << ")" << std::endl;
		}
		
		// At exterior boundary: normal should point INTO the medium (opposite to typical mesh normal)
		// For reflection at exterior boundary, we want normal pointing inward
		glm::dvec3 normal = glm::normalize(photon.voxel_normal);
		glm::dvec3 incident = glm::normalize(photon.direction);
		
		// At exterior boundary, flip normal to point inward if it's pointing outward
		// Check if normal points away from medium center or in same direction as ray
		if (glm::dot(incident, normal) > 0.0) {
			normal = -normal; // Make normal point inward
		}
		
		// Recalculate reflection and transmission with correct normal
		double cos_i = std::abs(glm::dot(incident, normal));
		
		// Reflection: incident ray bounces back into medium
		glm::dvec3 corrected_reflectance = incident - 2.0 * glm::dot(incident, normal) * normal;
		corrected_reflectance = glm::normalize(corrected_reflectance);
		
		// Transmission: use the pre-calculated transmittance direction
		double reflection = temp_reflection;
		double transmission = 1.0 - reflection;

		// total transmission
		if (reflection == 0.0) {
			// photon dies
			photon.direction = transmittance;
			photon.alive = false;

			radiate(photon, transmittance, photon.weight);
		}
		// total internal reflection
		else if (reflection == 1.0) {
			// photon reflects off surface
			photon.direction = corrected_reflectance;
			photon.position = move_delta(photon.intersect, photon.direction);
			photon.voxel = voxel_at(photon.position);
		}
		else {
			// partial reflection
			if (config.partial()) {
				// emit partial reflectance
				radiate(photon, transmittance, photon.weight * transmission);

				// adjust photon weight
				photon.weight *= reflection;

				// update direction/position/voxel
				photon.direction = corrected_reflectance;
				photon.position = move_delta(photon.intersect, photon.direction);
				photon.voxel = voxel_at(photon.position);
			}
			else { // all-or-none transmission/reflection
				// total transmission
				if (mcml_random->next() > reflection) {
					photon.direction = transmittance;
					photon.alive = false;

					radiate(photon, transmittance, photon.weight);
				}
				// total reflection
				else {
					photon.direction = corrected_reflectance;
					photon.position = move_delta(photon.intersect, photon.direction);
					photon.voxel = voxel_at(photon.position);
				}
			}
		}

		metrics.increment_scatters();
	}
	// 2. crossing to another medium
	else if (newvox != nullptr && newvox->tissue != nullptr && photon.voxel->tissue != nullptr
			 && photon.voxel->tissue->eta != newvox->tissue->eta) {
		// Use already computed reflection/transmission
		double reflection = temp_reflection;

		// total transmission
		if (reflection == 0.0) {
			photon.direction = transmittance;
			photon.position = move_delta(photon.intersect, photon.direction);
			photon.voxel = voxel_at(photon.position);
		}
		// total internal reflection
		else if (reflection == 1.0) {
			photon.direction = reflectance;
			photon.position = move_delta(photon.intersect, photon.direction);
			photon.voxel = voxel_at(photon.position);
		}
		else { // all-or-none transmission/reflection
			// total transmission
			if (mcml_random->next() > reflection) {
				photon.direction = transmittance;
				photon.position = move_delta(photon.intersect, photon.direction);
				photon.voxel = voxel_at(photon.position);
			}
			// total reflection
			else {
				photon.direction = reflectance;
				photon.position = move_delta(photon.intersect, photon.direction);
				photon.voxel = voxel_at(photon.position);
			}
		}

		metrics.increment_scatters();
	}
	// 3. crossing within the same medium (total transmission)
	else {
		// direction is unchanged
		photon.position = move_delta(photon.intersect, photon.direction);
		photon.voxel = voxel_at(photon.position);
	}
}

/***********************************************************
 * Record the emittance from a photon (partially) leaving
 * the material.
 ***********************************************************/
void Simulator::radiate(Photon& photon, glm::dvec3& direction, double weight) {
	if (!photon.voxel)
		return;

	// Record emittance at the selected voxel
	Voxel* exit_voxel = photon.voxel;
	
	// Calculate appropriate weight based on voxel type
	double effective_weight = weight;
	if (exit_voxel->is_boundary_voxel) {
		// Scale by the volume fraction outside the geometry for boundary voxels
		effective_weight = weight * exit_voxel->volume_fraction_outside;
	}
	// For surface voxels or other cases, use full weight
	
	// Always record emittance - every photon exit must be recorded somewhere
	exit_voxel->emittance += effective_weight;
	
	// Optional: Log when recording at non-ideal voxels for debugging
	if (!exit_voxel->is_boundary_voxel && !exit_voxel->is_surface_voxel) {
		static int warning_count = 0;
		if (warning_count < 10) { // Limit spam
			std::cerr << "Info: Emittance recorded at interior voxel (" 
					  << exit_voxel->ix() << ", " << exit_voxel->iy() << ", " << exit_voxel->iz() 
					  << ") - may indicate transport precision issues." << std::endl;
			warning_count++;
		}
	}
}

/***********************************************************
 * Determine the survivability of a given photon using MCML 3.0.0 algorithm.
 ***********************************************************/
void Simulator::roulette(Photon& photon) {
	// Use MCML 3.0.0 Russian roulette
	roulette_photon(photon);
}

/***********************************************************
 * Scatter the photon into a new direction based on the
 * Henyey-Greenstein phase function using MCML 3.0.0 algorithm.
 ***********************************************************/
void Simulator::scatter(Photon& photon) {
	if (!photon.alive) {
		return;
	}

	// Get tissue properties for scattering
	Tissue* tissue = photon.voxel->tissue;
	if (!tissue) {
		photon.alive = false;
		return;
	}

	// Use MCML 3.0.0 scattering algorithm
	scatter_photon(photon, *tissue);

	// normalize direction vector (safety check)
	photon.direction = glm::normalize(photon.direction);

	// prevent scattering into ambient medium when close to boundaries
	glm::dvec3 newpos = move_delta(photon.position, photon.direction);
	if (!voxel_at(newpos)) {
		photon.alive = false;
		radiate(photon, photon.direction, photon.weight);
	}

	metrics.increment_scatters();

	// add new internal position to path
	paths.back().add_internal_vertex(std::make_shared<PhotonNode>(photon.position, photon.weight));
}

/***********************************************************
 * Normalize the recorded values based on the number of
 * photons traced.
 ***********************************************************/
void Simulator::normalize() {
	// normalize globally recorded parameters
	record.total_absorption /= config.num_photons() * config.num_sources();
	record.diffuse_reflection /= config.num_photons() * config.num_sources();
	record.diffuse_transmission /= config.num_photons() * config.num_sources();

	record.specular_reflection /= config.num_sources();                          // rs is only computed once per source
	record.specular_transmission /= config.num_photons() * config.num_sources(); // ts is computed per photon

	// normalize voxel data
	for (const auto& voxel_ptr : voxel_grid) {
		auto* voxel = voxel_ptr.get();
		// skip computation for voxels outside the medium
		if (!voxel->tissue) {
			continue;
		}

		voxel->absorption /= (config.num_photons() * config.num_sources());
		voxel->emittance /= (config.num_photons() * config.num_sources());
	}
}

/***********************************************************
 * Compute the specular reflectance from a light source at
 * the surface.
 ***********************************************************/
void Simulator::specular_reflection(Source& source) {
	Voxel* voxel = voxel_at(source.intersect);

	// voxel should never be nullptr at this point
	if (!voxel) {
		std::cerr << "Critical error: specular reflection could not be computed." << std::endl;
		exit(EXIT_FAILURE);
	}

	// refractive indices of ambient medium and medium that is hit
	double n1 = config.ambient_eta();
	double n2 = voxel->tissue->eta;

	// set the specular reflection
	double temp_ratio = (n1 - n2) / (n1 + n2);
	record.specular_reflection = (n2 != n1) ? temp_ratio * temp_ratio : 0;

	// reflection direction: R = V - 2(V . N)N
	// Ensure normal points away from incident ray (outward from surface)
	glm::dvec3 normal = source.triangle.normal();
	double projection_scalar = glm::dot(source.direction, normal);
	
	// If normal points toward incident ray, flip it to point outward
	if (projection_scalar > 0.0) {
		normal = -normal;
		projection_scalar = -projection_scalar;
	}
	
	glm::dvec3 twice_projection = normal * (projection_scalar * 2.0);
	glm::dvec3 rsdir = glm::dvec3(source.direction.x - twice_projection.x, source.direction.y - twice_projection.y,
								  source.direction.z - twice_projection.z);

	source.specular_direction = rsdir;
}

/***********************************************************
 * Compute the fraction of incoming light that is reflected
 * back at an interface between two media. Also compute the
 * directions of transmission and reflection.
 ***********************************************************/
double Simulator::internal_reflection(Photon& photon, double& eta_t, glm::dvec3& transmittance,
									  glm::dvec3& reflectance) {
	// Modern Fresnel equations with improved numerical stability
	double eta_i = photon.voxel->tissue->eta;
	double eta_ratio = eta_t / eta_i;
	double eta_ratio_sq = eta_ratio * eta_ratio;

	// Get normalized vectors
	glm::dvec3 normal = glm::normalize(photon.voxel_normal);
	glm::dvec3 incident = glm::normalize(photon.direction);

	// Ensure normal points toward the incident medium (opposite to ray direction)
	if (glm::dot(incident, normal) > 0.0) {
		normal = -normal;
	}

	// Use dot product with properly oriented normal
	double cos_i = std::abs(glm::dot(incident, normal));
	cos_i = std::clamp(cos_i, 0.0, 1.0); // Ensure valid range

	double sin_i_sq = 1.0 - cos_i * cos_i;
	double sin_t_sq = eta_ratio_sq * sin_i_sq;

	// Total internal reflection check with numerical tolerance
	if (sin_t_sq >= 1.0 - 1e-12) {
		// Total internal reflection
		reflectance = incident - 2.0 * glm::dot(incident, normal) * normal;
		reflectance = glm::normalize(reflectance);
		transmittance = glm::dvec3(0.0); // No transmission
		return 1.0;
	}

	// Calculate transmission angle with numerical stability
	double cos_t = std::sqrt(std::max(0.0, 1.0 - sin_t_sq));

	// Fresnel equations for s and p polarizations
	double r_s_num = eta_i * cos_i - eta_t * cos_t;
	double r_s_den = eta_i * cos_i + eta_t * cos_t;
	double r_p_num = eta_t * cos_i - eta_i * cos_t;
	double r_p_den = eta_t * cos_i + eta_i * cos_t;

	// Avoid division by zero
	double r_s = (std::abs(r_s_den) > 1e-12) ? (r_s_num * r_s_num) / (r_s_den * r_s_den) : 1.0;
	double r_p = (std::abs(r_p_den) > 1e-12) ? (r_p_num * r_p_num) / (r_p_den * r_p_den) : 1.0;

	double reflection = 0.5 * (r_s + r_p);
	reflection = std::clamp(reflection, 0.0, 1.0);

	// Transmission direction (Snell's law in vector form)
	transmittance = eta_ratio * incident + (eta_ratio * cos_i - cos_t) * normal;
	transmittance = glm::normalize(transmittance);

	// Reflection direction
	reflectance = incident - 2.0 * glm::dot(incident, normal) * normal;
	reflectance = glm::normalize(reflectance);

	return reflection;
}

/***********************************************************
 * Return a pointer to a voxel that encapsulates the given
 * position, or nullptr if the position is outside the medium.
 ***********************************************************/
Voxel* Simulator::voxel_at(glm::dvec3& position) {
	if (!bounds.includes(position)) {
		return nullptr;
	}

	// distances from boundaries
	double dx = std::fabs(bounds.min_bounds.x - position.x);
	double dy = std::fabs(bounds.min_bounds.y - position.y);
	double dz = std::fabs(bounds.min_bounds.z - position.z);

	// indices start at minimum boundaries of voxels
	uint32_t ix = static_cast<uint32_t>(std::floor(dx / config.vox_size()));
	uint32_t iy = static_cast<uint32_t>(std::floor(dy / config.vox_size()));
	uint32_t iz = static_cast<uint32_t>(std::floor(dz / config.vox_size()));

	// avoid index overflow
	if (ix >= config.nx()) {
		ix = config.nx() - 1;
	}
	if (iy >= config.ny()) {
		iy = config.ny() - 1;
	}
	if (iz >= config.nz()) {
		iz = config.nz() - 1;
	}

	// retrieve the voxel at the position using the Volume
	if (!voxel_grid.is_valid_coordinate(ix, iy, iz)) {
		return nullptr;
	}

	Voxel* voxel = voxel_grid(ix, iy, iz);

	// if voxel does not have a tissue, it is outside the medium
	return voxel->tissue ? voxel : nullptr;
}

/***********************************************************
 * Return the minimum and maximum positions of the given
 * voxel as a cuboid structure.
 ***********************************************************/
Cuboid Simulator::voxel_corners(Voxel* voxel) {
	uint32_t ix_min = voxel->ix();
	uint32_t iy_min = voxel->iy();
	uint32_t iz_min = voxel->iz();

	uint32_t ix_max = voxel->ix() + 1;
	uint32_t iy_max = voxel->iy() + 1;
	uint32_t iz_max = voxel->iz() + 1;

	// minimum voxel position
	float x_min = static_cast<float>(bounds.min_bounds.x + (config.vox_size() * ix_min));
	float y_min = static_cast<float>(bounds.min_bounds.y + (config.vox_size() * iy_min));
	float z_min = static_cast<float>(bounds.min_bounds.z + (config.vox_size() * iz_min));

	// maximum voxel position
	float x_max = static_cast<float>(bounds.min_bounds.x + (config.vox_size() * ix_max));
	float y_max = static_cast<float>(bounds.min_bounds.y + (config.vox_size() * iy_max));
	float z_max = static_cast<float>(bounds.min_bounds.z + (config.vox_size() * iz_max));

	// round off coordinate values around the origin
	x_min = (std::fabs(x_min) < 1E-10) ? 0 : x_min;
	y_min = (std::fabs(y_min) < 1E-10) ? 0 : y_min;
	z_min = (std::fabs(z_min) < 1E-10) ? 0 : z_min;

	x_max = (std::fabs(x_max) < 1E-10) ? 0 : x_max;
	y_max = (std::fabs(y_max) < 1E-10) ? 0 : y_max;
	z_max = (std::fabs(z_max) < 1E-10) ? 0 : z_max;

	return Cuboid(x_min, y_min, z_min, x_max, y_max, z_max);
}

/***********************************************************
 * Return the destination for a given origin, direction and
 * distance.
 ***********************************************************/
glm::dvec3 Simulator::move(glm::dvec3& position, glm::dvec3& direction, double d) {
	glm::dvec3 point = position;

	// return the end point of a sub-step
	point.x = position.x + direction.x * d;
	point.y = position.y + direction.y * d;
	point.z = position.z + direction.z * d;

	return point;
}

/***********************************************************
 * Return the destination for a given origin and direction
 * after making a small hop.
 ***********************************************************/
glm::dvec3 Simulator::move_delta(glm::dvec3& position, glm::dvec3& direction) {
	glm::dvec3 point = position;

	// delta distance (based on voxel size)
	double d = config.vox_size() * 0.00001;

	point.x = position.x + direction.x * d;
	point.y = position.y + direction.y * d;
	point.z = position.z + direction.z * d;

	return point;
}

/***********************************************************
 *	Write the resulting physical quantities to a file.
 ***********************************************************/
void Simulator::report() {
	std::string str_sim = "simulation.out";
	std::string str_abs = "absorption.out";
	std::string str_emi = "emittance.out";
	std::string str_ptd = "photons.out";

	std::ofstream ofs_rep(str_sim.c_str(), std::ios_base::out); // simulation report
	std::ofstream ofs_abs(str_abs.c_str(), std::ios_base::out); // absorption report
	std::ofstream ofs_emi(str_emi.c_str(), std::ios_base::out); // emittance report
	std::ofstream ofs_ptn(str_ptd.c_str(), std::ios_base::out); // photon exitance report

	if (!ofs_rep.good() || !ofs_abs.good() || !ofs_emi.good() || !ofs_ptn.good()) {
		std::cerr << "Error: an output file could not be opened." << std::endl;
		return;
	}

	if (ofs_rep.good()) {
		ofs_rep.precision(8);
		ofs_rep << "################################################################" << std::endl;
		ofs_rep << "# SIMULATION REPORT" << std::endl;
		ofs_rep << "################################################################" << std::endl;
		ofs_rep << std::endl;

		// write input configuration
		ofs_rep << "Configuration" << std::endl;
		ofs_rep << "################################################################" << std::endl;
		ofs_rep << std::endl;
		ofs_rep << "Number of photons: " << '\t' << config.num_photons() << std::endl;
		ofs_rep << "Number of layers:  " << '\t' << config.num_layers() << std::endl;
		ofs_rep << "Number of voxels:  " << '\t' << config.num_voxels() << std::endl;
		ofs_rep << "Grid dimensions:   " << '\t' << config.nx() << " x " << config.ny() << " x " << config.nz()
				<< std::endl;
		ofs_rep << "Voxel dimensions:  " << '\t' << config.vox_size() << " x " << config.vox_size() << " x "
				<< config.vox_size() << std::endl;
		ofs_rep << std::endl << std::endl;

		// write recorded parameters: a, rs, rd, (ts, td)
		ofs_rep << "Recorded parameters" << std::endl;
		ofs_rep << "################################################################" << std::endl;
		ofs_rep << "Total absorption:      " << '\t' << std::fixed << record.total_absorption << std::endl;
		ofs_rep << "Diffuse reflection:    " << '\t' << std::fixed << record.diffuse_reflection << std::endl;
		ofs_rep << "Specular reflection:   " << '\t' << std::fixed << record.specular_reflection << std::endl;
		ofs_rep << "Diffuse transmission:  " << '\t' << std::fixed << record.diffuse_transmission << std::endl;
		ofs_rep << "Specular transmission: " << '\t' << std::fixed << record.specular_transmission << std::endl;
		ofs_rep << std::endl;
		ofs_rep.close();
	}
	else {
		std::cerr << "Error: file " << str_sim << " could not be opened." << std::endl;
		ofs_rep.close();
	}

	// write voxel absorption (each 'block' is a slice)
	if (ofs_abs.good()) {
		ofs_abs.precision(5);
		ofs_abs << "################################################################" << std::endl;
		ofs_abs << "# ABSORPTION REPORT" << std::endl;
		ofs_abs << "################################################################" << std::endl;
		ofs_abs << std::endl;

		// from rear to front
		for (uint32_t iz = 0; iz < config.nz(); ++iz) {
			ofs_abs << "Slice " << iz + 1 << "/" << config.nz();
			if (iz == 0) {
				ofs_abs << " (rear)";
			}
			if (iz == config.nz() - 1) {
				ofs_abs << " (front)";
			}
			ofs_abs << std::endl;

			for (uint32_t iy = config.ny() - 1; iy > 0; --iy) { // top to bottom
				for (uint32_t ix = 0; ix < config.nx(); ++ix) { // left to right
					Voxel* voxel = voxel_grid(ix, iy, iz);
					ofs_abs << std::fixed << voxel->absorption << '\t';
				}
				ofs_abs << std::endl;
			}
			ofs_abs << std::endl;
		}
		ofs_abs.close();
	}
	else {
		std::cerr << "Error: file " << str_abs << " could not be opened." << std::endl;
		ofs_abs.close();
	}

	// write voxel emittance
	if (ofs_emi.good()) {
		ofs_emi.precision(5);
		ofs_emi << "################################################################" << std::endl;
		ofs_emi << "# EMITTANCE REPORT" << std::endl;
		ofs_emi << "################################################################" << std::endl;
		ofs_emi << std::endl;

		// rear to front
		for (uint32_t iz = 0; iz < config.nz(); ++iz) {
			ofs_emi << "Slice " << iz + 1 << "/" << config.nz();
			if (iz == 0) {
				ofs_emi << " (rear)";
			}
			if (iz == config.nz() - 1) {
				ofs_emi << " (front)";
			}
			ofs_emi << std::endl;

			for (uint32_t iy = config.ny() - 1; iy > 0; --iy) { // top to bottom
				for (uint32_t ix = 0; ix < config.nx(); ++ix) { // left to right
					Voxel* voxel = voxel_grid(ix, iy, iz);
					ofs_emi << std::fixed << voxel->emittance << '\t';
				}
				ofs_emi << std::endl;
			}
			ofs_emi << std::endl;
		}
		ofs_emi.close();
	}
	else {
		std::cerr << "Error: file " << str_emi << " could not be opened." << std::endl;
		ofs_emi.close();
	}

	// exiting photons (position, direction, weight)
	if (ofs_ptn.good()) {
		ofs_ptn.precision(8);
		ofs_ptn << "################################################################" << std::endl;
		ofs_ptn << "# PHOTON REPORT" << std::endl;
		ofs_ptn << "################################################################" << std::endl;
		ofs_ptn << std::endl;

		if (emitters.empty()) {
			ofs_ptn << "No photons exited the medium." << std::endl;
		}

		for (const auto& emitter : emitters) {
			ofs_ptn << "Photon" << std::endl;
			ofs_ptn << "{" << std::endl;
			ofs_ptn << std::fixed << '\t' << "id  = " << emitter.id << std::endl;
			ofs_ptn << std::fixed << '\t' << "position = " << emitter.position.x << ", " << emitter.position.y << ", "
					<< emitter.position.z << std::endl;
			ofs_ptn << std::fixed << '\t' << "direction = " << emitter.direction.x << ", " << emitter.direction.y
					<< ", " << emitter.direction.z << std::endl;
			ofs_ptn << std::fixed << '\t' << "val = " << emitter.weight << std::endl;
			ofs_ptn << "}" << std::endl;
			ofs_ptn << std::endl;
		}
		ofs_ptn.close();
	}
	else {
		std::cerr << "Error: file " << str_ptd << " could not be opened." << std::endl;
		ofs_ptn.close();
	}
}

/***********************************************************
 * MCML 3.0.0 Monte Carlo Methods - Integrated Implementation
 ***********************************************************/

void Simulator::generate_step_size(Photon& photon) {
	// Modern Beer-Lambert law implementation with improved numerical stability
	if (photon.step < 1e-12) { // Higher precision threshold
		double rnd;
		// More robust zero-avoidance for random number generation
		// Use std::numeric_limits for better precision handling
		constexpr double min_random = std::numeric_limits<double>::epsilon();
		do {
			rnd = mcml_random->next();
		}
		while (rnd <= min_random);

		// Use high-precision logarithm for step size calculation
		photon.step = -std::log(rnd);
	}
}

void Simulator::scatter_photon(Photon& photon, const Tissue& tissue) {
	// Modern numerically stable Henyey-Greenstein phase function implementation
	double cos_theta, sin_theta, cos_phi, sin_phi;
	double g = tissue.g; // anisotropy factor
	double rnd = mcml_random->next();

	// More numerically stable Henyey-Greenstein sampling
	if (std::abs(g) > 1e-12) {
		// Use improved numerical stability for extreme anisotropy values
		double g2 = g * g;
		double one_minus_g2 = 1.0 - g2;
		double temp = one_minus_g2 / (1.0 - g + 2.0 * g * rnd);
		cos_theta = (1.0 + g2 - temp * temp) / (2.0 * g);

		// Robust clamping with higher precision threshold
		cos_theta = std::clamp(cos_theta, -1.0, 1.0);
	}
	else {
		// Isotropic scattering for negligible anisotropy
		cos_theta = 2.0 * rnd - 1.0;
	}

	// More robust trigonometric calculation with numerical safety
	sin_theta = std::sqrt(std::max(0.0, 1.0 - cos_theta * cos_theta));

	// Uniform azimuthal angle sampling with better numerical properties
	rnd = mcml_random->next();
	double phi = MathConstants::TWO_PI * rnd;
	cos_phi = std::cos(phi);
	sin_phi = std::sin(phi);

	// Modernized direction update using GLM vector operations and rotation matrices
	// This replaces the manual coordinate transformation with more robust methods
	glm::dvec3 old_direction = photon.direction;

	// Use Rodrigues' rotation formula for more stable direction update
	if (std::abs(old_direction.z) > 0.99999) {
		// Special case: nearly aligned with z-axis - use direct coordinate assignment
		photon.direction =
			glm::dvec3(sin_theta * cos_phi, sin_theta * sin_phi, cos_theta * (old_direction.z > 0 ? 1.0 : -1.0));
	}
	else {
		// General case: Use stable rotation about axis perpendicular to incident direction
		// Create orthonormal basis vectors
		glm::dvec3 w = old_direction; // incident direction
		glm::dvec3 u = glm::normalize(glm::cross(std::abs(w.z) < 0.9 ? glm::dvec3(0, 0, 1) : glm::dvec3(1, 0, 0), w));
		glm::dvec3 v = glm::cross(w, u);

		// Construct scattered direction in local coordinate system
		photon.direction = sin_theta * cos_phi * u + sin_theta * sin_phi * v + cos_theta * w;
	}

	// Ensure direction normalization for numerical stability
	photon.direction = glm::normalize(photon.direction);

	// Mark that photon has scattered
	photon.scatters = true;
}

void Simulator::roulette_photon(Photon& photon) {
	// Modernized Russian roulette with adaptive threshold and survival probability
	if (photon.weight < mcml_weight_threshold) {
		// Use weight-dependent survival probability for better variance reduction
		// This is more efficient than fixed 10% survival rate
		double survival_probability = std::max(0.1, photon.weight / mcml_weight_threshold);
		survival_probability = std::min(survival_probability, 0.5); // Cap at 50% for stability

		if (mcml_random->next() <= survival_probability) {
			// Survive with proper weight normalization
			photon.weight /= survival_probability;
		}
		else {
			// Terminate photon
			photon.alive = false;
		}
	}
}

void Simulator::set_rng_seed(int seed) {
	mcml_random->seed(seed);
}

void Simulator::reset_simulation_data() {
	// Reset global record data
	record.total_absorption = 0.0;      // total absorption
	record.diffuse_reflection = 0.0;    // diffuse reflection
	record.specular_reflection = 0.0;   // specular reflection
	record.diffuse_transmission = 0.0;  // diffuse transmission
	record.specular_transmission = 0.0; // specular transmission

	// Reset all voxel data
	for (const auto& voxel_ptr : voxel_grid) {
		auto* voxel = voxel_ptr.get();
		if (voxel) {
			voxel->absorption = 0.0;
			voxel->emittance = 0.0;
		}
	}

	// Reset experimenter accumulated data
	metrics.reset();

	std::cout << "Simulation data reset completed." << std::endl;
}

bool Simulator::is_point_inside_geometry(const glm::dvec3& point) const {
	// Test if the point is inside any layer's geometry
	for (const auto& layer : layers) {
		if (layer.contains_point(point)) {
			return true;
		}
	}

	return false;
}

bool Simulator::is_point_inside_layer_mesh(const glm::dvec3& point, const Layer& layer) const {
	// Use the Layer's containment test
	return layer.contains_point(point);
}
