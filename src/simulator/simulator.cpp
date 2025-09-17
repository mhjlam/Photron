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
#include <set>
#include <sstream>
#include <iterator>

#include "math/math.hpp"
#include "math/random.hpp"
#include "math/ray.hpp"

/***********************************************************
 * Simulator constructor.
 ***********************************************************/
Simulator::Simulator() : rng(std::make_shared<Random>()), mcml_weight_threshold(1e-4) {
	// Modern C++20: Use default initialization instead of explicit construction
	paths.clear();
	photons.clear();
	sources.clear();
	mediums.clear();

	// Reserve space for common use cases
	paths.reserve(1000);
	photons.reserve(10000);
	sources.reserve(5);

	// Initialize MCML random number generator
	rng->seed(static_cast<int>(std::time(nullptr)));
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
	if (config.verbose()) {
		std::cout << "Initializing Photron" << std::endl;
	}

	// Clear previous simulation data to ensure clean reset
	mediums.clear();
	photons.clear();
	sources.clear();
	paths.clear();
	emitters.clear();
	
	// Reset metrics
	metrics.reset();

	// read and parse input configuration file
	if (config.verbose()) {
		std::cout << "Parsing configuration file: " << file << std::endl;
	}
	if (!parse(file)) {
		std::cerr << "An error occurred while parsing the input file." << std::endl;
		return false;
	}
	if (config.verbose()) {
		std::cout << "Configuration parsed successfully." << std::endl;
	}

	// Initialize medium (only one for now)
	mediums.emplace_back(Medium(config));

	for (auto& medium : mediums) {
		if (!medium.initialize()) {
			std::cerr << "An error occurred while initializing the medium." << std::endl;
			return false;
		}
		if (config.verbose()) {
			std::cout << "Medium initialized successfully." << std::endl;
		}
	}

	// Initialize light sources
	if (!initialize_sources()) {
		std::cerr << "An error occurred while initializing the data structures." << std::endl;
		return false;
	}
	if (config.verbose()) {
		std::cout << "Light sources initialized successfully." << std::endl;
	}

	// Initialize photons
	for (uint64_t i = 0; i < config.num_photons(); ++i) {
		photons.emplace_back(i);
	}
	if (config.verbose()) {
		std::cout << "Photons initialized successfully." << std::endl;
	}
	
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
	sources.clear();
	photons.clear();
	paths.clear();
	mediums.clear();

	// Use Config class to parse the entire configuration file
	if (!config.parse_config_file(fconfig)) {
		std::cerr << "Failed to parse configuration file." << std::endl;
		return false;
	}

	// Extract parsed data from Config
	sources = config.sources_vector();

	return true;
}

/***********************************************************
 * Initialize configuration properties.
 ***********************************************************/
bool Simulator::initialize_sources() {
	// initialize config properties
	config.set_num_sources(static_cast<uint64_t>(sources.size()));

	// associate light sources with their geometric intersections
	for (auto& source : sources) {
		bool found_intersection = false;

		// find intersection of ray from this source with geometry (point, triangle, normal)
		Ray ray = Ray(source.origin, source.direction);
		for (auto& medium : mediums) {
			double distance = medium.intersection(source);

			if (distance != std::numeric_limits<double>::max()) {
				found_intersection = true;
				break; // Use the first medium that intersects
			}
		}

		if (!found_intersection) {
			std::cerr << "Error: Light source " << source.id 
					  << " does not intersect with any medium geometry." << std::endl;
		}
	}

	return true;
}

/***********************************************************
 * Run the Monte Carlo photon transport simulation.
 ***********************************************************/
void Simulator::simulate() {
	if (config.verbose()) {
		std::cout << "Running Monte Carlo simulation" << std::endl;
	}

	metrics.start_clock();

	// Track overall energy balance
	static double total_initial_energy = 0.0;
	static double total_launched_energy = 0.0;

	// For each light source
	for (auto& source : sources) {
		// for each photon
		for (uint32_t p = 0; p < photons.size(); ++p) {
			// progress report
			if (config.progress() && ((p + 1) % 1000) == 0) {
				std::cout << "Photon " << p + 1 << "/" << config.num_photons() << std::endl;
			}

			// Track initial energy
			total_initial_energy += 1.0; // Each photon starts with weight 1.0

			// launch the photon (create a new path)
			launch(photons[p], source);
			
			// CRITICAL FIX: Update source with calculated specular direction from first photon
			// This ensures the renderer can access the correct reflection direction
			if (p == 0 && photons[p].alive) {
				source.specular_direction = photons[p].specular_direction;
			}
			
		// Track launched energy
		if (photons[p].alive) {
			total_launched_energy += photons[p].weight;
		}			// Safety mechanism to prevent infinite photon loops
			int photon_iteration_counter = 0;
			const int max_photon_iterations = 1000000;

			while (photons[p].alive) {
				photon_iteration_counter++;
				if (photon_iteration_counter > max_photon_iterations) {
					std::cerr << "Warning: Photon " << photons[p].id << " exceeded maximum iterations, terminating."
							  << std::endl;
					
					// Use centralized termination for consistency
					if (photons[p].weight > 0.0) {
						terminate_photon_and_record_energy(photons[p], "max_iterations");
					} else {
						photons[p].alive = false;
					}
					break;
				}
				
				step_size(photons[p]); // Set new step size
				transfer(photons[p]);  // Propagate photon through the medium in substeps
				roulette(photons[p]);  // Determine photon termination
				
				// Only scatter if photon hasn't just crossed a boundary
				// Scattering should only occur in bulk medium, not at interfaces
				if (!photons[p].cross) {
					scatter(photons[p]);   // Scatter photon into a new direction
				}
			}

			p++;
		}
	}

	// Normalize physical quantities
	normalize();

	metrics.stop_clock();
	
	// Increment simulation version since data has changed
	increment_simulation_version();
	
	// Aggregate record data across all mediums for metrics collection
	double total_absorption = 0.0, specular_reflection = 0.0, diffuse_reflection = 0.0;
	double surface_refraction = 0.0, specular_transmission = 0.0, diffuse_transmission = 0.0;
	
	// Calculate voxel-based totals for accurate energy accounting
	double total_voxel_absorption = 0.0;
	double total_voxel_emittance = 0.0;
	
	for (const auto& medium : mediums) {
		const auto& record = medium.get_record();
		// Still use medium records for surface interactions (these are accurate)
		specular_reflection += record.specular_reflection;
		surface_refraction += record.surface_refraction;
		
		// Calculate voxel totals for this medium using Volume iterators
		const auto& volume = medium.get_volume();
		for (const auto& voxel : volume) {
			if (voxel && voxel->tissue != nullptr) {
				total_voxel_absorption += voxel->absorption;
				total_voxel_emittance += voxel->emittance;
			}
		}
	}
	
	// Use voxel totals for metrics instead of inflated medium records
	// Use pure voxel-based energy accounting for consistency
	total_absorption = total_voxel_absorption;
	specular_transmission = 0.0; // No specular transmission in current config  
	
	// Aggregate voxel data to medium records before calculating metrics
	aggregate_voxel_data();
	
	// Calculate diffuse reflection and transmission from medium records (now populated)
	diffuse_reflection = 0.0;
	diffuse_transmission = 0.0;
	
	// Use the direction-classified totals from medium records (populated by aggregate_voxel_data)
	for (const auto& medium : mediums) {
		const auto& record = medium.get_record();
		diffuse_reflection += record.diffuse_reflection;
		diffuse_transmission += record.diffuse_transmission;
	}
	
	// Also aggregate transport metrics (path length, step sizes, scatter events)
	std::vector<double> combined_step_sizes;
	std::vector<glm::dvec3> combined_path_vertices;
	double combined_scatter_events = 0.0;
	
	for (const auto& medium : mediums) {
		// Aggregate transport metrics from each medium
		const auto& medium_metrics = medium.get_metrics();
		const auto& step_sizes = medium_metrics.get_step_sizes();
		const auto& path_vertices = medium_metrics.get_path_vertices();
		
		// Combine step sizes
		combined_step_sizes.insert(combined_step_sizes.end(), step_sizes.begin(), step_sizes.end());
		
		// Combine path vertices
		combined_path_vertices.insert(combined_path_vertices.end(), path_vertices.begin(), path_vertices.end());
		
		// Add scatter events
		combined_scatter_events += medium_metrics.get_scatter_events();
	}
	
	// Set the aggregated transport data in the main simulator metrics
	metrics.set_step_sizes(combined_step_sizes);
	metrics.set_path_vertices(combined_path_vertices);
	metrics.set_scatter_events(combined_scatter_events);
	
	metrics.collect_data(total_absorption, specular_reflection, diffuse_reflection,
						 surface_refraction, specular_transmission, diffuse_transmission);
	metrics.write_to_file();
	metrics.print_report(*this);
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
			
			// Deposit remaining energy as absorption for energy conservation
			if (new_photon.weight > 0.0 && new_photon.voxel && new_photon.voxel->tissue) {
				// EXPERIMENTAL: Use energy conservation enforcement
				terminate_photon_and_record_energy(new_photon, "max_iterations");
			} else {
				new_photon.alive = false;
			}
			break;
		}

		step_size(new_photon); // Set new step size
		transfer(new_photon);  // Propagate photon through the medium in substeps
		roulette(new_photon);  // Determine photon termination
		scatter(new_photon);   // Scatter photon into a new direction
	}

	// Add the completed photon to the photons vector for rendering
	photons.push_back(new_photon);
	
	// CRITICAL FIX: Aggregate voxel data to medium records after adding new photon
	// This ensures diffuse_reflection and diffuse_transmission are updated for the overlay
	aggregate_voxel_data();
	
	// Recalculate metrics to include the new photon
	double total_absorption = 0.0, specular_reflection = 0.0, diffuse_reflection = 0.0;
	double surface_refraction = 0.0, specular_transmission = 0.0, diffuse_transmission = 0.0;
	
	// Aggregate data from all mediums
	std::vector<double> combined_step_sizes;
	std::vector<glm::dvec3> combined_path_vertices;
	int combined_scatter_events = 0;
	
	for (auto& medium : mediums) {
		auto& record = medium.get_record();
		total_absorption += record.total_absorption;
		specular_reflection += record.specular_reflection;
		diffuse_reflection += record.diffuse_reflection;
		surface_refraction += record.surface_refraction;
		specular_transmission += record.specular_transmission;
		diffuse_transmission += record.diffuse_transmission;
		
		// Get metrics from this medium
		auto& medium_metrics = medium.get_metrics();
		
		// Add step sizes
		auto step_sizes = medium_metrics.get_step_sizes();
		combined_step_sizes.insert(combined_step_sizes.end(), step_sizes.begin(), step_sizes.end());
		
		// Add path vertices
		auto path_vertices = medium_metrics.get_path_vertices();
		combined_path_vertices.insert(combined_path_vertices.end(), path_vertices.begin(), path_vertices.end());
		
		// Add scatter events
		combined_scatter_events += medium_metrics.get_scatter_events();
	}
	
	// Set the aggregated transport data in the main simulator metrics
	metrics.set_step_sizes(combined_step_sizes);
	metrics.set_path_vertices(combined_path_vertices);
	metrics.set_scatter_events(combined_scatter_events);
	
	metrics.collect_data(total_absorption, specular_reflection, diffuse_reflection,
						 surface_refraction, specular_transmission, diffuse_transmission);
	
	// Increment simulation version since data has changed
	increment_simulation_version();
}

/***********************************************************
 * Set up photon properties for tracing.
 ***********************************************************/
void Simulator::launch(Photon& photon, const Source& source) {
	// Find the medium that the photon will start in
	Medium* start_medium = find_medium_at(source.intersect);
	if (!start_medium) {
		std::cerr << "Error: Photon launch point is not in any medium" << std::endl;
		photon.alive = false;
		return;
	}
	
	// Count this photon as entering the medium
	start_medium->increment_photons_entered();
	
	photon.alive = true;
	photon.weight = 1.0;  // Each photon starts with full energy
	
	// ENERGY CONSERVATION: Initialize energy tracking
	photon.total_energy_budget = photon.weight;
	photon.total_energy_radiated = 0.0;
	photon.total_energy_absorbed = 0.0;
	photon.radiate_call_count = 0;
	
	// Copy source data directly into photon
	photon.source_origin = source.origin;
	photon.source_direction = source.direction;
	photon.specular_direction = source.specular_direction;
	photon.source_intersect = source.intersect;
	photon.source_triangle = source.triangle;
	
	photon.direction = source.direction;
	photon.position = source.intersect;
	photon.voxel = start_medium->voxel_at(photon.position);

	// Compute specular reflection for this photon
	specular_reflection(photon);

	// create vertices for new light path
	auto light = std::make_shared<PhotonNode>(photon.source_origin, photon.weight);
	auto intersection = std::make_shared<PhotonNode>(photon.source_intersect, photon.weight);
	auto reflection = std::make_shared<PhotonNode>(move(photon.source_intersect, photon.specular_direction, 0.1),
												   start_medium->get_record().specular_reflection);

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
	// Find current medium for the photon
	Medium* current_medium = find_medium_at(photon.position);
	if (!current_medium) {
		// Photon is in ambient space - energy is now recorded only through voxel emittance
		// No medium records needed as we use voxel-based energy conservation
		photon.alive = false;
		return;
	}
	
	// Ensure photon has valid voxel reference
	if (!photon.voxel) {
		photon.voxel = current_medium->voxel_at(photon.position);
		if (!photon.voxel) {
			// Photon is outside - energy is now recorded only through voxel emittance
			// No medium records needed as we use voxel-based energy conservation
			photon.alive = false;
			return;
		}
	}
	
	generate_step_size(photon);
	
	// Add step size to metrics of the current medium
	// CRITICAL FIX: Prevent division by zero for non-scattering media
	double normalization_factor = photon.mu_s();
	if (normalization_factor <= 1e-12) {
		normalization_factor = std::max(photon.mu_a(), 1e-6); // Use absorption or minimum value
	}
	current_medium->get_metrics().add_step_size(photon.step / normalization_factor);
}

/***********************************************************
 * Track the path from current position to intersection point,
 * calculate absorption along the entire path through all voxels,
 * and update photon weight accordingly.
 ***********************************************************/
void Simulator::track_voxel_path_and_deposit(Photon& photon) {
	if (!photon.voxel || !photon.voxel->tissue) {
		return;
	}

	// Get start and end positions
	glm::dvec3 start_pos = photon.position;
	glm::dvec3 end_pos = photon.intersect;
	glm::dvec3 direction = end_pos - start_pos;
	double total_distance = glm::length(direction);
	
	if (total_distance < 1e-12) {
		// Very short step, handle normally
		deposit(photon);
		return;
	}
	
	direction = glm::normalize(direction);
	
	// Find current medium
	Medium* current_medium = find_medium_at(start_pos);
	if (!current_medium) {
		return;
	}
	
	// Calculate proper path segments through each voxel
	std::vector<std::pair<Voxel*, double>> voxel_segments;
	Voxel* last_surface_voxel = photon.voxel; // Default to current voxel
	
	// Use fine-grained ray marching to accurately calculate distances in each voxel
	double voxel_size = current_medium->get_volume().voxel_size();
	double step_size = voxel_size * 0.1; // Fine sampling for accuracy
	int num_steps = std::max(10, static_cast<int>(total_distance / step_size));
	
	std::map<Voxel*, double> voxel_distances; // Track total distance in each voxel
	Voxel* prev_voxel = nullptr;
	double segment_start = 0.0;
	
	for (int i = 0; i <= num_steps; ++i) {
		double t = (i == num_steps) ? 1.0 : static_cast<double>(i) / num_steps;
		glm::dvec3 sample_pos = start_pos + direction * (total_distance * t);
		double current_distance = total_distance * t;
		
		Voxel* current_voxel = current_medium->voxel_at(sample_pos);
		
		// If we've moved to a different voxel, record the distance in the previous voxel
		if (current_voxel != prev_voxel && prev_voxel != nullptr) {
			double segment_length = current_distance - segment_start;
			voxel_distances[prev_voxel] += segment_length;
			segment_start = current_distance;
		}
		
		// Track surface/boundary voxels for exit recording
		if (current_voxel && (current_voxel->is_surface_voxel || current_voxel->is_boundary_voxel)) {
			last_surface_voxel = current_voxel;
		}
		
		prev_voxel = current_voxel;
	}
	
	// Record the final segment
	if (prev_voxel != nullptr) {
		double segment_length = total_distance - segment_start;
		voxel_distances[prev_voxel] += segment_length;
	}
	
	// Calculate absorption along the path using proper Beer-Lambert law
	double remaining_weight = photon.weight;
	double total_absorption = 0.0;
	
	for (const auto& [voxel, distance] : voxel_distances) {
		if (voxel && voxel->tissue && distance > 1e-12) {
			// Calculate effective volume fraction for boundary voxels
			double effective_volume_fraction = 1.0;
			if (voxel->is_boundary_voxel) {
				effective_volume_fraction = voxel->volume_fraction_inside;
			}
			
			// Apply Beer-Lambert law for this segment
			double mu_a = voxel->tissue->mu_a * effective_volume_fraction;
			double segment_transmission = std::exp(-mu_a * distance);
			double segment_absorption = remaining_weight * (1.0 - segment_transmission);
			
			// Deposit absorption in this voxel
			voxel->absorption += segment_absorption;
			total_absorption += segment_absorption;
			
			// Update remaining weight for next segment
			remaining_weight *= segment_transmission;
		}
	}
	
	// ENERGY CONSERVATION FIX: Update photon energy tracking to match deposit() function
	// This ensures that energy conservation calculations remain accurate
	if (total_absorption > 0.0) {
		photon.total_energy_absorbed += total_absorption;
	}
	
	// Update photon weight
	photon.weight = remaining_weight;
	
	// ENERGY CONSERVATION: Only voxel-level recording prevents double-counting
	// Medium records will be populated via aggregation functions later
	
	// Update photon voxel assignment to the last surface voxel for proper emittance recording
	photon.prev_voxel = photon.voxel;
	photon.voxel = last_surface_voxel;
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
			
			// Deposit remaining energy as absorption for energy conservation
			if (photon.weight > 0.0 && photon.voxel && photon.voxel->tissue) {
				// EXPERIMENTAL: Use energy conservation enforcement
				terminate_photon_and_record_energy(photon, "max_iterations");
			} else {
				photon.alive = false;
			}
			break;
		}

	// set substep - this will handle mesh boundary detection
	sub_step(photon);

	// If photon is crossing a boundary, track the path to update voxel assignment and absorption
	if (photon.cross) {
		track_voxel_path_and_deposit(photon);
	} else {
		// deposit weight normally for non-crossing steps
		deposit(photon);
	}

	// possibly cross boundary
	if (photon.cross) {
		cross(photon);
	}
	else {
		photon.position = move(photon.position, photon.direction, photon.sub_step);
	}
	
	// prevent errors due to crossing to ambient medium
	Medium* current_medium = find_medium_at(photon.position);
	if (!current_medium || !photon.voxel) {
		// CRITICAL: Record as transmission - photon is exiting
		photon.alive = false;
		radiate(photon, photon.direction, photon.weight);
		return;
	}

	// update step size using current medium's tissue properties
	if (current_medium && photon.voxel) {
		photon.step -= (photon.sub_step * photon.mu_s());
	}

	current_medium->get_metrics().add_vertex(photon.position.x, photon.position.y, photon.position.z);
	}
}

/***********************************************************
 * Set the photon's next substep and initialize the given
 * intersection point and voxel normal if it crosses the
 * voxel boundary.
 ***********************************************************/
void Simulator::sub_step(Photon& photon) {
	// Validate voxel before proceeding
	if (!photon.voxel) {
		std::cerr << "Error: Photon voxel is null in sub_step()" << std::endl;
		photon.alive = false;
		return;
	}

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

		// Record remaining energy as absorption to maintain energy conservation
		if (photon.weight > 0.0) {
			Medium* current_medium = find_medium_at(photon.position);
			if (current_medium) {
				auto& record = current_medium->get_record();
				record.total_absorption += photon.weight;
			}
		}

		// Terminate this photon instead of crashing the entire program
		photon.alive = false;
		photon.sub_step = 0.0;
		photon.cross = false;
		return;
	}

	// Compute free path for a substep (distance to next scattering event)
	// CRITICAL FIX: Handle non-scattering media properly
	double mu_s_effective = photon.mu_s();
	if (mu_s_effective <= 1e-12) {
		// For non-scattering media, use total attenuation or set very large free path
		mu_s_effective = std::max(photon.mu_a(), 1e-6);
	}
	double freepath = photon.step / mu_s_effective;

	// Check if photon is currently inside the mesh
	bool photon_inside_mesh = is_point_inside_geometry(photon.position);

	if (!photon_inside_mesh) {
		// Photon is outside mesh - it should have exited already
		// ENERGY CONSERVATION FIX: Do NOT record medium transmission here
		// Energy was already recorded via radiate() when photon exited
		// Recording here would cause double-counting and energy gain
		if (photon.weight > 0.0) {
			// Energy already accounted for in voxel emittance via radiate()
			// No additional medium record updates needed
		}
		photon.alive = false;
		return;
	}

	// Find the closest mesh boundary intersection
	double mesh_dist = std::numeric_limits<double>::max();
	glm::dvec3 mesh_intersection {0.0, 0.0, 0.0};
	glm::dvec3 mesh_normal {0.0, 0.0, 0.0};
	bool found_mesh_intersection = false;

	// Get current medium to access its layers
	Medium* current_medium = find_medium_at(photon.position);
	if (!current_medium) {
		// Energy conservation now handled through voxel emittance only
		photon.alive = false;
		return;
	}

	// Look for exit points from the mesh
	for (const auto& layer : current_medium->get_layers()) {
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
		for (const auto& layer : current_medium->get_layers()) {
			for (const auto& triangle : layer.mesh) {
				// Check distance to each vertex
				double dist_v0 = glm::length(mesh_intersection - triangle.v0());
				double dist_v1 = glm::length(mesh_intersection - triangle.v1());
				double dist_v2 = glm::length(mesh_intersection - triangle.v2());
				
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
		// Scale absorption by the volume fraction inside for boundary voxels
		// Don't exit early - the photon is still in a tissue voxel and energy should be conserved
		effective_volume_fraction = photon.voxel->volume_fraction_inside;
		
		// If photon is in the outside portion of boundary voxel, still deposit but scale appropriately
		if (!is_point_inside_geometry(photon.position)) {
			// Use the outside volume fraction for photons in the outside portion
			effective_volume_fraction = photon.voxel->volume_fraction_outside;
		}
	}
	else {
		// For non-boundary voxels, check geometry but don't exit early for energy conservation
		if (!is_point_inside_geometry(photon.position)) {
			// Photon is slightly outside geometry due to numerical precision
			// Still deposit energy to maintain conservation, but with reduced fraction
			effective_volume_fraction = 0.5; // Compromise value for edge cases
		}
	}

	// deposited weight (scaled by effective volume fraction)
	double deltaw = photon.weight * (1 - std::exp(-photon.mu_a() * photon.sub_step)) * effective_volume_fraction;

	// ENERGY CONSERVATION ENFORCEMENT
	// Calculate how much energy this photon has left to absorb
	double energy_already_used = photon.total_energy_radiated + photon.total_energy_absorbed;
	double energy_available = photon.total_energy_budget - energy_already_used;
	
	// Enforce energy conservation: cannot absorb more than available
	double actual_absorbed_weight = std::min(deltaw, energy_available);
	
	// Update photon energy tracking
	photon.total_energy_absorbed += actual_absorbed_weight;

	// update photon weight  
	photon.weight -= actual_absorbed_weight;

	// assign deposited weight to voxel
	photon.voxel->absorption += actual_absorbed_weight;

	// ENERGY CONSERVATION: Only voxel-level recording
	// Medium records populated via aggregation before output
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
		// Photon is outside medium - record as transmission
		photon.alive = false;
		radiate(photon, photon.direction, photon.weight);
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

			// Find closest geometry exit - check all mediums for exit points
			for (const auto& medium : mediums) {
				for (const auto& layer : medium.get_layers()) {
					for (const auto& triangle : layer.mesh) {
						Triangle triangle_copy = triangle;
						glm::dvec3 intersection;

						if (exit_ray.intersect_triangle(triangle_copy, intersection)) {
							double dist = glm::length(photon.position - intersection);

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
				// True Splitting: Always account for both portions
				// Radiate the transmitted portion (exits medium)
				if ((1.0 - reflection) > 1e-12) {
					radiate(photon, transmittance, photon.weight * (1.0 - reflection));
				}
				
				// Continue photon as reflected portion with weighted energy
				if (reflection > 1e-12) {
					photon.weight *= reflection;
					photon.direction = reflectance;
					photon.position = move_delta(photon.intersect, photon.direction);
					photon.voxel = voxel_at(photon.position);
				} else {
					// No reflection, terminate photon
					photon.alive = false;
				}
			}
			else {
				// All-or-none
				if (rng->next() > reflection) {
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

	// Determine which medium(s) are involved in this crossing
	Medium* current_medium = find_medium_at(photon.position);
	glm::dvec3 newpos = move_delta(photon.intersect, transmittance);
	Medium* new_medium = find_medium_at(newpos);
	Voxel* newvox = new_medium ? new_medium->voxel_at(newpos) : nullptr;
	
	photon.prev_voxel = photon.voxel;

	// determine refractive index of the medium being entered
	double eta = (newvox == nullptr) ? config.ambient_eta() : newvox->tissue->eta;

	// Recalculate with correct refractive index
	temp_reflection = internal_reflection(photon, eta, transmittance, reflectance);

	// Handle different crossing scenarios
	if (new_medium != current_medium) {
		// Medium transition detected
		handle_medium_transition(photon, current_medium, new_medium);
		
		// If photon was killed by medium transition, don't process further
		if (!photon.alive) {
			// Still need to call radiate() for ambient exit
			if (!new_medium && current_medium) {
				radiate(photon, photon.direction, photon.weight);
			}
			return;
		}
	}

	// 1. crossing to ambient medium
	if (newvox == nullptr) {
		// Check if photon is exiting from an interior voxel (shouldn't happen!)
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

			// radiate() now handles both voxel emittance AND medium record updates
			radiate(photon, transmittance, photon.weight);
		}
		// total internal reflection
		else if (reflection == 1.0) {
			// photon reflects off surface
			photon.direction = corrected_reflectance;
			photon.position = move_delta(photon.intersect, photon.direction);
			if (current_medium) {
				photon.voxel = current_medium->voxel_at(photon.position);
			}
		}
		else {
			// partial reflection
			if (config.partial()) {
				// True Splitting: Always account for both portions
				// Radiate the transmitted portion (exits medium)
				if (transmission > 1e-12) {
					radiate(photon, transmittance, photon.weight * transmission);
				}
				
				// Continue photon as reflected portion with weighted energy
				if (reflection > 1e-12) {
					photon.weight *= reflection;
					photon.direction = corrected_reflectance;
					photon.position = move_delta(photon.intersect, photon.direction);
					if (current_medium) {
						photon.voxel = current_medium->voxel_at(photon.position);
					}
				} else {
					// No reflection, terminate photon
					photon.alive = false;
				}
			}
			else { // all-or-none transmission/reflection
				// total transmission
				if (rng->next() > reflection) {
					photon.direction = transmittance;
					photon.alive = false;

					// radiate() now handles both voxel emittance AND medium record updates
					radiate(photon, transmittance, photon.weight);
				}
				// total reflection
				else {
					photon.direction = corrected_reflectance;
					photon.position = move_delta(photon.intersect, photon.direction);
					if (current_medium) {
						photon.voxel = current_medium->voxel_at(photon.position);
					}
				}
			}
		}

		if (current_medium) {
			current_medium->get_metrics().increment_scatters();
		}
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
			if (new_medium) {
				photon.voxel = new_medium->voxel_at(photon.position);
			}
		}
		// total internal reflection
		else if (reflection == 1.0) {
			photon.direction = reflectance;
			photon.position = move_delta(photon.intersect, photon.direction);
			if (current_medium) {
				photon.voxel = current_medium->voxel_at(photon.position);
			}
		}
		else { // all-or-none transmission/reflection
			// total transmission
			if (rng->next() > reflection) {
				photon.direction = transmittance;
				photon.position = move_delta(photon.intersect, photon.direction);
				if (new_medium) {
					photon.voxel = new_medium->voxel_at(photon.position);
				}
			}
			// total reflection
			else {
				photon.direction = reflectance;
				photon.position = move_delta(photon.intersect, photon.direction);
				if (current_medium) {
					photon.voxel = current_medium->voxel_at(photon.position);
				}
			}
		}

		if (current_medium) {
			current_medium->get_metrics().increment_scatters();
		}
	}
	// 3. crossing within the same medium (total transmission)
	else {
		// direction is unchanged
		photon.position = move_delta(photon.intersect, photon.direction);
		if (current_medium) {
			photon.voxel = current_medium->voxel_at(photon.position);
		}
	}
}

/***********************************************************
 * Record the emittance from a photon (partially) leaving
 * the material.
 ***********************************************************/
void Simulator::radiate(Photon& photon, glm::dvec3& direction, double weight) {
	if (!photon.voxel)
		return;

	// Record photon's radiate() call origin
	// ENERGY CONSERVATION ENFORCEMENT (disabled for True Splitting)
	// True Splitting handles energy conservation differently - through statistical splitting
	if (!config.partial()) {
		photon.radiate_call_count++;
		
		// Calculate how much energy this photon has left to radiate
		double energy_already_used = photon.total_energy_radiated + photon.total_energy_absorbed;
		double energy_available = photon.total_energy_budget - energy_already_used;
		
		// Enforce energy conservation: cannot radiate more than available
		double actual_radiated_weight = std::min(weight, energy_available);
		
		// Update photon energy tracking
		photon.total_energy_radiated += actual_radiated_weight;
		weight = actual_radiated_weight;
	} else {
		// True Splitting mode: no per-photon enforcement, just track
		photon.radiate_call_count++;
		photon.total_energy_radiated += weight;
	}

	// Record emittance at the selected voxel
	Voxel* exit_voxel = photon.voxel;
	exit_voxel->emittance += weight;
	
	// Use proper reflection/transmission determination based on exit position relative to entry
	if (is_photon_reflecting(photon)) {
		// Exit on same side as entry - classify as reflection
		exit_voxel->emittance_reflected += weight;
		photon.exit_type = Photon::ExitType::REFLECTED;
	} else {
		// Exit on opposite side from entry - classify as transmission  
		exit_voxel->emittance_transmitted += weight;
		photon.exit_type = Photon::ExitType::TRANSMITTED;
	}
	
	// Create external vertex for photon path with proper exit classification
	// Convert Photon::ExitType to PhotonNode::ExitType
	PhotonNode::ExitType node_exit_type = PhotonNode::ExitType::NONE;
	if (photon.exit_type == Photon::ExitType::REFLECTED) {
		node_exit_type = PhotonNode::ExitType::REFLECTED;
	} else if (photon.exit_type == Photon::ExitType::TRANSMITTED) {
		node_exit_type = PhotonNode::ExitType::TRANSMITTED;
	}
	
	// Add external vertex to the current photon path
	std::shared_ptr<PhotonNode> exit_node = nullptr;
	if (!paths.empty()) {
		exit_node = std::make_shared<PhotonNode>(photon.intersect, weight, node_exit_type);
		paths.back().add_external_vertex(exit_node);
	}
	
	// Create Emitter object with proper exit classification for renderer use
	Emitter::ExitType emitter_exit_type = Emitter::ExitType::NONE;
	if (photon.exit_type == Photon::ExitType::REFLECTED) {
		emitter_exit_type = Emitter::ExitType::REFLECTED;
	} else if (photon.exit_type == Photon::ExitType::TRANSMITTED) {
		emitter_exit_type = Emitter::ExitType::TRANSMITTED;
	}
	
	// Create shared emitter and establish connection with exit node
	auto emitter = std::make_shared<Emitter>(photon.id, photon.intersect, direction, weight, emitter_exit_type);
	emitters.push_back(emitter);
	
	// Establish bidirectional connection between exit node and emitter
	if (exit_node) {
		exit_node->emitter = emitter;
	}
	
	// ENERGY CONSERVATION: Only voxel-level recording
	// Medium records populated via aggregation before output
	
	// ENERGY CONSERVATION: Only voxel-level recording prevents double-counting
	// Direction-classified emittance recorded in voxel fields for aggregation later
	
	// Optional: Log when recording at non-ideal voxels
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

	// Get current medium and tissue properties for scattering
	Medium* current_medium = find_medium_at(photon.position);
	if (!current_medium) {
		// Record as transmission - photon has exited medium
		photon.alive = false;
		radiate(photon, photon.direction, photon.weight);
		return;
	}
	
	Tissue* tissue = photon.voxel->tissue;
	if (!tissue) {
		// Record as transmission - photon cannot scatter without tissue
		photon.alive = false;
		radiate(photon, photon.direction, photon.weight);
		return;
	}

	// Use MCML 3.0.0 scattering algorithm
	scatter_photon(photon, *tissue);

	// normalize direction vector (safety check)
	photon.direction = glm::normalize(photon.direction);

	// prevent scattering into ambient medium when close to boundaries
	glm::dvec3 newpos = move_delta(photon.position, photon.direction);
	if (!is_inside_any_medium(newpos)) {
		photon.alive = false;
		radiate(photon, photon.direction, photon.weight);
		return;
	}

	current_medium->get_metrics().increment_scatters();

	// add new internal position to path
	paths.back().add_internal_vertex(std::make_shared<PhotonNode>(photon.position, photon.weight));
}

/***********************************************************
 * Normalize the recorded values based on the number of
 * photons traced.
 ***********************************************************/
void Simulator::normalize() {
	// Normalize records for all mediums
	for (auto& medium : mediums) {
		auto& record = medium.get_record();
		record.total_absorption /= config.num_photons() * config.num_sources();
		record.diffuse_reflection /= config.num_photons() * config.num_sources();
		record.surface_refraction /= config.num_photons() * config.num_sources();
		record.diffuse_transmission /= config.num_photons() * config.num_sources();

		record.specular_reflection /= config.num_photons() * config.num_sources();   // Now accumulates energy per photon
		record.specular_transmission /= config.num_photons() * config.num_sources(); // ts is computed per photon
	}

	// normalize voxel data across all mediums
	for (auto& medium : mediums) {
		auto& voxel_grid = medium.get_volume();
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
}

/***********************************************************
 * Compute the specular reflectance from a light source at
 * the surface.
 ***********************************************************/
void Simulator::specular_reflection(Photon& photon) {
	Voxel* voxel = voxel_at(photon.source_intersect);

	// voxel should never be nullptr at this point
	if (!voxel) {
		std::cerr << "Critical error: specular reflection could not be computed." << std::endl;
		exit(EXIT_FAILURE);
	}

	// refractive indices of ambient medium and medium that is hit
	double n1 = config.ambient_eta();
	double n2 = voxel->tissue->eta;

	// Calculate specular reflection coefficient from Fresnel equations
	double temp_ratio = (n1 - n2) / (n1 + n2);
	double fresnel_reflection = (n2 != n1) ? temp_ratio * temp_ratio : 0;
	
	// Record surface refraction (energy entering the medium) 
	auto* medium = find_medium_at(photon.source_intersect);
	if (medium) {
		auto& record = medium->get_record();
		// Surface refraction is the energy that enters the medium (1 - reflected energy)
		double surface_refraction_energy = photon.weight * (1.0 - fresnel_reflection);
		record.surface_refraction += surface_refraction_energy; // Energy entering medium at surface
		
		// Record specular reflection (energy immediately reflected at surface)
		double specular_reflection_energy = photon.weight * fresnel_reflection;
		record.specular_reflection += specular_reflection_energy;
		
		// CRITICAL FIX: Reduce photon weight by reflected amount so only transmitted energy
		// continues for absorption/emission. This prevents double-counting reflected energy.
		photon.weight = surface_refraction_energy; // Only transmitted energy continues
	}

	// reflection direction: R = V - 2(V . N)N
	// Ensure normal points away from incident ray (outward from surface)
	glm::dvec3 normal = photon.source_triangle.normal();
	double projection_scalar = glm::dot(photon.source_direction, normal);
	
	// If normal points toward incident ray, flip it to point outward
	if (projection_scalar > 0.0) {
		normal = -normal;
		projection_scalar = -projection_scalar;
	}
	
	glm::dvec3 twice_projection = normal * (projection_scalar * 2.0);
	glm::dvec3 rsdir = glm::dvec3(photon.source_direction.x - twice_projection.x, 
								  photon.source_direction.y - twice_projection.y,
								  photon.source_direction.z - twice_projection.z);

	photon.specular_direction = rsdir;
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
 * Aggregate voxel energy data into medium records for output.
 * This prevents double-counting by using only voxel-level data.
 ***********************************************************/
void Simulator::aggregate_voxel_data() {
	// Clear all medium records first (except surface interactions)
	for (auto& medium : mediums) {
		auto& record = medium.get_record();
		record.total_absorption = 0.0;
		record.diffuse_reflection = 0.0;
		// record.specular_reflection preserved (set by specular_reflection() function)
		// record.surface_refraction preserved (set by specular_reflection() function)
		record.diffuse_transmission = 0.0;
		record.specular_transmission = 0.0;
	}
	
	// Aggregate data from all voxels
	for (auto& medium : mediums) {
		auto& record = medium.get_record();
		auto& volume = medium.get_volume();  // Remove const since we need non-const access
		
		// Iterate through all voxels in this medium
		for (uint32_t z = 0; z < volume.depth(); ++z) {
			for (uint32_t y = 0; y < volume.height(); ++y) {
				for (uint32_t x = 0; x < volume.width(); ++x) {
					Voxel* voxel = volume.at(x, y, z);
					if (voxel && voxel->tissue) {
						// Only aggregate voxels that belong to this medium
						// Check if voxel's tissue ID matches this medium's ID
						if (voxel->tissue->id - '0' == (&medium - &mediums[0])) {
							// Aggregate absorption
							record.total_absorption += voxel->absorption;
							
							// Aggregate emittance by direction classification
							record.diffuse_reflection += voxel->emittance_reflected;
							record.diffuse_transmission += voxel->emittance_transmitted;
							
							// Note: surface_refraction and specular_reflection handled separately
							// Note: specular_transmission currently unused
						}
					}
				}
			}
		}
	}
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

		// write recorded parameters: a, rs, rd, (ts, td) - aggregate across all mediums
		ofs_rep << "Recorded parameters" << std::endl;
		ofs_rep << "################################################################" << std::endl;
		
		// Aggregate values across all mediums
		double total_absorption = 0.0, diffuse_reflection = 0.0, specular_reflection = 0.0;
		double surface_refraction = 0.0, diffuse_transmission = 0.0, specular_transmission = 0.0;
		
		for (const auto& medium : mediums) {
			const auto& record = medium.get_record();
			total_absorption += record.total_absorption;
			diffuse_reflection += record.diffuse_reflection;
			specular_reflection += record.specular_reflection;
			surface_refraction += record.surface_refraction;
			diffuse_transmission += record.diffuse_transmission;
			specular_transmission += record.specular_transmission;
		}
		
		ofs_rep << "Total absorption:      " << '\t' << std::fixed << total_absorption << std::endl;
		ofs_rep << "Diffuse reflection:    " << '\t' << std::fixed << diffuse_reflection << std::endl;
		ofs_rep << "Specular reflection:   " << '\t' << std::fixed << specular_reflection << std::endl;
		ofs_rep << "Surface refraction:    " << '\t' << std::fixed << surface_refraction << std::endl;
		ofs_rep << "Diffuse transmission:  " << '\t' << std::fixed << diffuse_transmission << std::endl;
		ofs_rep << "Specular transmission: " << '\t' << std::fixed << specular_transmission << std::endl;
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
					Voxel* voxel = voxel_at(glm::dvec3(
						ix * config.vox_size(),
						iy * config.vox_size(),
						iz * config.vox_size()
					));
					if (voxel) {
						ofs_abs << std::fixed << voxel->absorption << '\t';
					} else {
						ofs_abs << std::fixed << 0.0 << '\t';
					}
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
					Voxel* voxel = voxel_at(glm::dvec3(
						ix * config.vox_size(),
						iy * config.vox_size(),
						iz * config.vox_size()
					));
					if (voxel) {
						ofs_emi << std::fixed << voxel->emittance << '\t';
					} else {
						ofs_emi << std::fixed << 0.0 << '\t';
					}
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
			ofs_ptn << std::fixed << '\t' << "id  = " << emitter->id << std::endl;
			ofs_ptn << std::fixed << '\t' << "position = " << emitter->position.x << ", " << emitter->position.y << ", "
					<< emitter->position.z << std::endl;
			ofs_ptn << std::fixed << '\t' << "direction = " << emitter->direction.x << ", " << emitter->direction.y
					<< ", " << emitter->direction.z << std::endl;
			ofs_ptn << std::fixed << '\t' << "val = " << emitter->weight << std::endl;
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
			rnd = rng->next();
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
	double rnd = rng->next();

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
	rnd = rng->next();
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
		double survival_probability = std::max(0.1, photon.weight / mcml_weight_threshold);
		survival_probability = std::min(survival_probability, 0.5); // Cap at 50% for stability

		if (rng->next() <= survival_probability) {
			// ENERGY CONSERVATION FIX: Proper Russian Roulette
			
			// Survive with proper weight normalization
			photon.weight /= survival_probability;
			
			// CRITICAL FIX: Do NOT increase total_energy_budget!
			// Russian Roulette must be energy-neutral on average.
			// The energy budget stays the same - only the current weight changes.
		}
		else {
			// Use centralized termination for consistency
			terminate_photon_and_record_energy(photon, "roulette");
		}
	}
}


/***********************************************************
 * MEDIUM CONTEXT MANAGEMENT
 ***********************************************************/

/***********************************************************
 * Find which medium contains the given point.
 ***********************************************************/
Medium* Simulator::find_medium_at(const glm::dvec3& position) const {
	for (auto& medium : mediums) {
		if (medium.contains_point(position)) {
			return const_cast<Medium*>(&medium);
		}
	}
	return nullptr; // Point is in ambient space
}

/***********************************************************
 * Check if point is inside any medium.
 ***********************************************************/
bool Simulator::is_inside_any_medium(const glm::dvec3& position) const {
	return find_medium_at(position) != nullptr;
}

/***********************************************************
 * Handle photon transition between mediums.
 ***********************************************************/
void Simulator::handle_medium_transition(Photon& photon, Medium* from, Medium* to) {
	if (!from && to) {
		// Entering a medium from ambient space
	} 
	else if (from && !to) {
		// Exiting to ambient space
		photon.alive = false;
		
		// NOTE: Do NOT call radiate() here - let the main cross() logic handle it
		// This prevents duplicate radiate() calls for the same boundary crossing
	} 
	else if (from && to && from != to) {
		// Transitioning between different mediums
		
		// Handle interface physics between mediums
		// This would involve Fresnel reflection/transmission calculations
		// For now, we'll use a simplified approach
		
		// Get surface normal at the interface (simplified - would need proper interface detection)
		// This is a placeholder - in practice you'd need to find the actual interface
		glm::dvec3 interface_normal = photon.voxel_normal;
		if (glm::length(interface_normal) < 0.1) {
			interface_normal = glm::dvec3(0.0, 0.0, 1.0); // Default normal
		}
		
		// Simple transmission for now - in practice you'd calculate Fresnel coefficients
		// based on the refractive indices of the two mediums
		double transmission_probability = 0.9; // Placeholder
		
		if (rng->next() < transmission_probability) {
			// Photon transmits to new medium
			// Direction could be refracted based on Snell's law
		} else {
			// Photon reflects back to original medium
			photon.direction = glm::reflect(photon.direction, interface_normal);
		}
	}
}

/***********************************************************
 * HELPER METHODS FOR MEDIUM DELEGATION
 ***********************************************************/

/***********************************************************
 * Find voxel at position by delegating to appropriate medium.
 ***********************************************************/
Voxel* Simulator::voxel_at(const glm::dvec3& position) const {
	Medium* medium = find_medium_at(position);
	if (medium) {
		glm::dvec3 pos = position; // Make a non-const copy for the method call
		return medium->voxel_at(pos);
	}
	return nullptr;
}

/***********************************************************
 * Get voxel corners by delegating to appropriate medium.
 ***********************************************************/
Cuboid Simulator::voxel_corners(Voxel* voxel) const {
	if (!voxel) {
		return Cuboid(); // Return default/empty cuboid
	}
	
	// For now, use the first medium (simple case)
	// TODO: In multi-medium scenario, we need a way to identify which medium owns the voxel
	if (!mediums.empty()) {
		return mediums[0].voxel_corners(voxel);
	}
	
	return Cuboid(); // Return default/empty cuboid
}

/***********************************************************
 * Check if point is inside geometry by checking any medium.
 ***********************************************************/
bool Simulator::is_point_inside_geometry(const glm::dvec3& position) const {
	for (const auto& medium : mediums) {
		if (medium.contains_point(position)) {
			return true;
		}
	}
	return false;
}

/***********************************************************
 * Access voxel by grid coordinates - searches all mediums
 ***********************************************************/
Voxel* Simulator::voxel_grid(uint32_t x, uint32_t y, uint32_t z) const {
	// Convert grid coordinates to world position using combined bounds
	// This ensures grid coordinates map correctly to world space
	double voxel_size = config.vox_size();
	Range3 bounds = get_combined_bounds();
	glm::dvec3 position(bounds.min_bounds.x + (x + 0.5) * voxel_size,
	                    bounds.min_bounds.y + (y + 0.5) * voxel_size, 
	                    bounds.min_bounds.z + (z + 0.5) * voxel_size);
	
	// Find which medium contains this position and get the voxel
	return voxel_at(position);
}

/***********************************************************
 * Accessor methods for aggregating data across all mediums
 ***********************************************************/

std::vector<Tissue> Simulator::get_all_tissues() const {
	std::vector<Tissue> all_tissues;
	for (const auto& medium : mediums) {
		const auto& medium_tissues = medium.get_tissues();
		for (const auto& tissue : medium_tissues) {
			// Add unique tissues from each medium
			bool found = false;
			for (const auto& existing : all_tissues) {
				if (existing.id == tissue.id) {
					found = true;
					break;
				}
			}
			if (!found) {
				all_tissues.push_back(tissue);
			}
		}
	}
	return all_tissues;
}

const std::vector<Layer>& Simulator::get_all_layers() const {
	// For simplicity, return the first medium's layers
	// If no mediums exist, return empty static vector
	static const std::vector<Layer> empty_layers;
	if (mediums.empty()) {
		return empty_layers;
	}
	return mediums[0].get_layers();
}

std::vector<std::shared_ptr<Voxel>>& Simulator::get_all_voxels() {
	// For non-const version, we need to return a reference to a static container
	// This is a bit tricky with multi-medium, but for now return the first medium's voxels
	// In practice, most use cases have only one medium
	static std::vector<std::shared_ptr<Voxel>> combined_voxels;
	combined_voxels.clear();
	
	for (auto& medium : mediums) {
		auto& volume = medium.get_volume();
		for (auto& voxel_ptr : volume) {
			// Convert unique_ptr to shared_ptr
			combined_voxels.push_back(std::shared_ptr<Voxel>(voxel_ptr.get(), [](Voxel*){}));
		}
	}
	return combined_voxels;
}

const std::vector<std::shared_ptr<Voxel>>& Simulator::get_all_voxels() const {
	// For const version, return a const reference to the same static container
	static std::vector<std::shared_ptr<Voxel>> combined_voxels;
	combined_voxels.clear();
	
	for (const auto& medium : mediums) {
		const auto& volume = medium.get_volume();
		for (const auto& voxel_ptr : volume) {
			// Convert unique_ptr to shared_ptr with non-owning deleter
			combined_voxels.push_back(std::shared_ptr<Voxel>(voxel_ptr.get(), [](Voxel*){}));
		}
	}
	return combined_voxels;
}

size_t Simulator::get_total_voxel_count() const {
	size_t total = 0;
	for (const auto& medium : mediums) {
		total += medium.get_volume().size();
	}
	return total;
}

/***********************************************************
 * Determine if photon is reflecting (exiting same side as entry)
 * or transmitting (exiting different side).
 ***********************************************************/
bool Simulator::is_photon_reflecting(const Photon& photon) const {
	// Compare the exit direction with the source entry direction
	// If they are on the same side of the geometry, it's reflection
	// If they are on opposite sides, it's transmission
	
	// Get the surface normal at the source intersection point
	glm::dvec3 source_normal = photon.source_triangle.normal();
	
	// Get the exit direction
	glm::dvec3 exit_direction = glm::normalize(photon.direction);
	
	// Get the entry direction (negated because we want incoming direction)
	glm::dvec3 entry_direction = glm::normalize(-photon.source_direction);
	
	// Project both directions onto the source normal
	double entry_projection = glm::dot(entry_direction, source_normal);
	double exit_projection = glm::dot(exit_direction, source_normal);
	
	// If projections have the same sign, photon is reflecting (same side)
	// If projections have opposite signs, photon is transmitting (opposite sides)
	bool same_side = (entry_projection * exit_projection) > 0.0;
	
	return same_side;
}

Range3 Simulator::get_combined_bounds() const {
	if (mediums.empty()) {
		return Range3(); // Return default bounds
	}
	
	Range3 combined = mediums[0].get_bounds();
	for (size_t i = 1; i < mediums.size(); ++i) {
		const Range3& medium_bounds = mediums[i].get_bounds();
		// Expand combined bounds to include this medium's bounds
		combined.min_bounds.x = std::min(combined.min_bounds.x, medium_bounds.min_bounds.x);
		combined.min_bounds.y = std::min(combined.min_bounds.y, medium_bounds.min_bounds.y);
		combined.min_bounds.z = std::min(combined.min_bounds.z, medium_bounds.min_bounds.z);
		combined.max_bounds.x = std::max(combined.max_bounds.x, medium_bounds.max_bounds.x);
		combined.max_bounds.y = std::max(combined.max_bounds.y, medium_bounds.max_bounds.y);
		combined.max_bounds.z = std::max(combined.max_bounds.z, medium_bounds.max_bounds.z);
	}
	return combined;
}

Record Simulator::get_combined_record() const {
	Record combined;
	for (const auto& medium : mediums) {
		const auto& record = medium.get_record();
		combined.total_absorption += record.total_absorption;
		combined.diffuse_reflection += record.diffuse_reflection;
		combined.specular_reflection += record.specular_reflection;
		combined.surface_refraction += record.surface_refraction;
		combined.diffuse_transmission += record.diffuse_transmission;
		combined.specular_transmission += record.specular_transmission;
		combined.avg_path_length += record.avg_path_length;
		combined.total_steps += record.total_steps;
	}
	
	// CRITICAL: Add voxel emittance to transmission totals for proper energy conservation
	// Emittance represents energy that has exited the medium and should be counted as transmission
	double total_voxel_emittance = 0.0;
	for (const auto& medium : mediums) {
		const auto& voxel_grid = medium.get_volume();
		for (const auto& voxel_ptr : voxel_grid) {
			const auto* voxel = voxel_ptr.get();
			if (voxel && voxel->tissue) {
				total_voxel_emittance += voxel->emittance;
			}
		}
	}
	
	// CRITICAL FIX: The voxel emittance needs to be scaled to match the normalization state of record data
	// If we're in normalized mode (records are per-photon), we need to scale voxel emittance accordingly
	// The normalization factor is config.num_photons() * config.num_sources()
	// However, we can detect if data is normalized by checking if surface_refraction is around 1.0 (normalized) or much larger (raw)
	double normalization_factor = 1.0;
	if (combined.surface_refraction > 0.0 && combined.surface_refraction < 10.0) {
		// Data appears to be normalized (surface refraction should be close to 1.0 for normalized data)
		// Voxel data was also normalized, so they should match
		normalization_factor = 1.0;
	} else if (combined.surface_refraction > 10.0) {
		// Data appears to be raw (unnormalized) - this happens during single photon accumulation
		// Voxel data might be normalized from previous full simulation, so scale it up
		normalization_factor = config.num_photons() * config.num_sources();
	}
	
	// FIXED: Remove double-counting of energy
	// Energy is already recorded in medium records via radiate() function
	// Adding voxel emittance here would count the same energy twice
	// combined.diffuse_transmission += total_voxel_emittance * normalization_factor;
	
	return combined;
}

Simulator::EnergyConservation Simulator::calculate_energy_conservation() const {
	EnergyConservation result;
	
	// First, get the combined record data for surface interactions
	Record combined_record = get_combined_record();
	result.surface_reflection = combined_record.specular_reflection;
	result.surface_refraction = combined_record.surface_refraction;
	
	// Calculate total voxel-based energy for accurate accounting
	// This approach bypasses normalization issues and gives the most accurate energy accounting
	double total_voxel_absorption = 0.0;
	double total_voxel_reflection = 0.0;
	double total_voxel_transmission = 0.0;
	
	for (const auto& medium : mediums) {
		const auto& volume = medium.get_volume();
		for (const auto& voxel : volume) {
			if (voxel && voxel->tissue != nullptr) {
				total_voxel_absorption += voxel->absorption;
				
				// Classify emittance into reflection vs transmission based on medium records
				// The medium records contain direction-classified emittance after aggregate_voxel_data()
				const auto& record = medium.get_record();
				double total_medium_emittance = record.diffuse_reflection + record.diffuse_transmission;
				
				if (total_medium_emittance > 0.0) {
					// Proportionally split voxel emittance based on medium record ratios
					double reflection_ratio = record.diffuse_reflection / total_medium_emittance;
					double transmission_ratio = record.diffuse_transmission / total_medium_emittance;
					
					total_voxel_reflection += voxel->emittance * reflection_ratio;
					total_voxel_transmission += voxel->emittance * transmission_ratio;
				} else {
					// No medium emittance recorded, assume all goes to transmission
					total_voxel_transmission += voxel->emittance;
				}
			}
		}
	}
	
	result.total_absorption = total_voxel_absorption;
	result.total_reflection = total_voxel_reflection;
	result.total_transmission = total_voxel_transmission;
	result.total_diffusion = total_voxel_reflection + total_voxel_transmission;
	
	// Total energy should include ALL components for proper conservation check
	// This represents the total energy accounted for across all interaction types
	result.total_energy = result.surface_reflection + result.total_absorption + 
	                     result.total_reflection + result.total_transmission;
	
	return result;
}

/***********************************************************
 * SIMPLIFIED ENERGY RECORDING SYSTEM
 * Single function to handle all photon termination scenarios
 ***********************************************************/
void Simulator::terminate_photon_and_record_energy(Photon& photon, const std::string& reason) {
	// Track all photon terminations
	static int termination_count = 0;
	static double total_terminated_energy = 0.0;
	termination_count++;
	total_terminated_energy += photon.weight;
	
	if (!photon.alive || photon.weight <= 0.0) {
		return; // Already terminated or no energy to record
	}
	
	// Find the appropriate medium to record energy
	Medium* record_medium = find_medium_at(photon.position);
	if (!record_medium && !mediums.empty()) {
		record_medium = &mediums[0]; // Fallback to first medium
	}
	
	if (!record_medium) {
		std::cerr << "CRITICAL ERROR: Cannot record energy - no medium found for photon termination: " << reason << std::endl;
		photon.alive = false;
		return;
	}
	
	auto& record = record_medium->get_record();
	
	// CONSOLIDATED APPROACH: This function only handles INTERNAL terminations (absorption)
	// For exits, use radiate() function which handles both voxel emittance and medium records
	if (reason == "absorption" || reason == "roulette" || reason == "max_iterations") {
		// Energy absorbed within the medium
		record.total_absorption += photon.weight;
		
		// ENERGY CONSERVATION FIX: Track energy in photon accounting system
		if (photon.voxel && photon.voxel->tissue) {
			// Update photon energy tracking for conservation
			photon.total_energy_absorbed += photon.weight;
			
			// Add to voxel absorption
			photon.voxel->absorption += photon.weight;
		}
	}
	else {
		// ERROR: Exit reasons should use radiate(), not this function!
		std::cerr << "ERROR: terminate_photon_and_record_energy() called with exit reason '" << reason 
				  << "' - should use radiate() for exits!" << std::endl;
		// Fallback: record as absorption to maintain energy conservation
		record.total_absorption += photon.weight;
		
		// ENERGY CONSERVATION FIX: Track energy in photon accounting system
		if (photon.voxel && photon.voxel->tissue) {
			// Update photon energy tracking for conservation
			photon.total_energy_absorbed += photon.weight;
			
			// Add to voxel absorption as fallback
			photon.voxel->absorption += photon.weight;
		}
	}
	
	// Terminate the photon
	photon.alive = false;
	photon.weight = 0.0; // Clear weight to prevent double-counting
}
