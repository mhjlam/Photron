#include "simulator.hpp"

#include "common/result.hpp"
#include "common/error_types.hpp"
#include "common/optimized_logger.hpp"
#include "common/error_handler.hpp"

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

#include "app.hpp" // For output path utilities
#include "math/math.hpp"
#include "math/random.hpp"
#include "math/ray.hpp"
#include "math/voxel_dda3d.hpp"

/***********************************************************
 * Simulator constructor.
 ***********************************************************/
Simulator::Simulator() : rng(std::make_shared<Random>()), mcml_weight_threshold(1e-4) {
	// Modern C++20: Use default initialization instead of explicit construction
	photons.clear();
	sources.clear();
	mediums.clear();

	// Reserve space for common use cases
	photons.reserve(10000);
	sources.reserve(5);

	// Initialize MCML random number generator 
	// Note: Will be re-seeded in initialize() based on config settings
	rng->seed(static_cast<int>(std::time(nullptr)));
}

/***********************************************************
 * Initialize shared metrics for energy statistics
 ***********************************************************/
void Simulator::initialize_shared_metrics() {
	// Create shared metrics instance that can be used for energy statistics
	shared_metrics_ = std::make_shared<Metrics>();
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
Result<void, SimulationError> Simulator::initialize(std::string file) {
	// Clear previous simulation data to ensure clean reset
	mediums.clear();
	photons.clear();
	sources.clear();
	emitters.clear();
	
	// Reset metrics
	metrics.reset();

	// Re-initialize the Config singleton with the file to ensure fresh data
	// This is critical for rerun functionality since data gets moved out during first run
	Config::shutdown(); // Clear the singleton
	if (!Config::initialize(file)) {
		return Result<void, SimulationError>::error(SimulationError::InvalidConfiguration);
	}

	// read and parse input configuration file
	auto parse_result = parse(file);
	if (!parse_result.is_ok()) {
		return parse_result; // Forward the error
	}
	
	// Only show loading message if logging is enabled
	if (Config::get().log()) {
		std::cout << "Loading configuration: " << file << std::endl << std::endl;
	}

	// Configure random number generator based on deterministic setting
	if (Config::get().deterministic()) {
		// Use fixed seed for reproducible results
		const int deterministic_seed = 12345;
		rng->seed(deterministic_seed);
		if (Config::get().log()) {
			std::cout << "Deterministic mode enabled: Using fixed seed " << deterministic_seed << std::endl;
		}
		FAST_LOG_INFO("Deterministic mode enabled: Using fixed seed " + std::to_string(deterministic_seed));
	} else {
		// Use time-based seed for stochastic behavior
		int time_seed = static_cast<int>(std::time(nullptr));
		rng->seed(time_seed);
		if (Config::get().log()) {
			std::cout << "Stochastic mode: Using time-based seed " << time_seed << std::endl;
		}
		FAST_LOG_INFO("Stochastic mode: Using time-based seed " + std::to_string(time_seed));
	}

	// Initialize medium (only one for now)
	mediums.emplace_back(Medium(Config::get()));

	for (auto& medium : mediums) {
		if (!medium.initialize()) {
			std::cerr << "An error occurred while initializing the medium." << std::endl;
			FAST_LOG_ERROR("Failed to initialize medium");
			return Result<void, SimulationError>::error(SimulationError::MediumInitializationFailure);
		}
		if (Config::get().log()) {
			// Log initialization success to file only
		FAST_LOG_INFO("Medium initialized successfully.");
		}
	}

	// Initialize 3D DDA instances for robust voxel traversal
	initialize_dda_instances();
	if (Config::get().log()) {
		// Log DDA initialization to file only  
		FAST_LOG_INFO("3D DDA voxel traversal initialized successfully.");
	}

	// Initialize debug logger (only in log mode)
	if (Config::get().log()) {
		Logger::instance().initialize(
			App::get_output_path("trace.csv"), 
			App::get_output_path("debug.log"), 
			true
		);
		Logger::instance().log_info("=== Photron Simulation Started ===");
		Logger::instance().log_info("Initializing Photron");
		Logger::instance().log_info("Configuration parsed successfully.");
		std::cout << "Debug logger initialized for photon tracing." << std::endl;
		Logger::instance().log_info("Debug logger initialized for photon tracing.");
	} else {
		Logger::instance().initialize("", "", false);  // Disable logging
	}

	// Initialize light sources
	auto init_sources_result = initialize_sources();
	if (!init_sources_result.is_ok()) {
		return init_sources_result; // Forward the error
	}
	if (Config::get().log()) {
		// Log light source initialization to file only
		FAST_LOG_INFO("Light sources initialized successfully.");
	}

	// Initialize photons
	for (uint64_t i = 0; i < Config::get().num_photons(); ++i) {
		photons.emplace_back(i);
	}
	
	// Shared metrics should be set externally via set_shared_metrics() before initialization
	if (!shared_metrics_) {
		std::cerr << "Warning: shared_metrics_ not set. Please call set_shared_metrics() before initialize()." << std::endl;
		// Fallback: create local metrics instance
		shared_metrics_ = std::make_shared<Metrics>();
	}
	
	if (Config::get().log()) {
		// Simple confirmation message for console
		std::cout << "Photron simulation initialized successfully." << std::endl;
		FAST_LOG_INFO("Photron simulation initialized successfully.");
	}
	
	return Result<void, SimulationError>::ok();
}

/***********************************************************
 *	INITIALIZATION SUBROUTINES
 ***********************************************************/

/***********************************************************
 * Read and parse the given configuration file.
 ***********************************************************/
Result<void, SimulationError> Simulator::parse(const std::string& /* fconfig */) {
	// Clear existing data structures for reinitialization
	sources.clear();
	photons.clear();
	mediums.clear();

	// Config should already be initialized by this point with the same file
	// Just extract parsed data from Config
	if (!Config::is_initialized()) {
		return Result<void, SimulationError>::error(SimulationError::ConfigNotInitialized);
	}
	
	sources = Config::get().sources();
	
	// Only show debug info if logging is enabled
	if (Config::get().log()) {
		std::cout << "Parsed " << sources.size() << " sources from config" << std::endl;
	}
	
	if (sources.empty()) {
		return Result<void, SimulationError>::error(SimulationError::NoSources);
	}

	return Result<void, SimulationError>::ok();
}

/***********************************************************
 * Initialize configuration properties.
 ***********************************************************/
Result<void, SimulationError> Simulator::initialize_sources() {
	// initialize config properties
	Config::get().set_num_sources(static_cast<uint64_t>(sources.size()));

	// associate light sources with their geometric intersections
	for (auto& source : sources) {
		bool found_intersection = false;

		// find intersection of ray from this source with geometry (point, triangle, normal)
		Ray ray = Ray(source.origin, source.direction);
		for (auto& medium : mediums) {
			double distance = medium.intersection(source);

			if (distance != std::numeric_limits<double>::max()) {
				// Verify that the intersection point is in a valid voxel using DDA
				Medium* validation_medium = find_medium_at_with_dda(source.intersect);
				if (validation_medium != nullptr) {
					found_intersection = true;
					break; // Use the first medium that intersects with valid voxel
				} else {
					// Try slightly adjusting the intersection point inward along the direction
					glm::dvec3 adjusted_intersect = source.intersect + source.direction * 1e-6;
					validation_medium = find_medium_at_with_dda(adjusted_intersect);
					if (validation_medium != nullptr) {
						source.intersect = adjusted_intersect; // Use the adjusted point
						found_intersection = true;
						break;
					}
				}
			}
		}

		if (!found_intersection) {
			return Result<void, SimulationError>::error(SimulationError::SourceIntersectionFailure);
		}
	}

	return Result<void, SimulationError>::ok();
}

// Legacy compatibility methods (deprecated)
bool Simulator::parse_legacy(const std::string& fconfig) {
	auto result = parse(fconfig);
	return result.is_ok();
}

bool Simulator::initialize_sources_legacy() {
	auto result = initialize_sources();
	return result.is_ok();
}

/***********************************************************
 * Run the Monte Carlo photon transport simulation.
 ***********************************************************/
Result<void, SimulationError> Simulator::simulate() {
	// Validate simulation state before starting
	if (sources.empty()) {
		return Result<void, SimulationError>::error(SimulationError::InvalidConfiguration);
	}
	
	if (photons.empty()) {
		return Result<void, SimulationError>::error(SimulationError::NoPhotons);
	}

	metrics.start_clock();

	// Track overall energy balance
	static double total_initial_energy = 0.0;
	static double total_launched_energy = 0.0;

	// For each light source
	for (auto& source : sources) {
		// for each photon
		for (uint32_t p = 0; p < photons.size(); ++p) {
			// Track initial energy
			total_initial_energy += 1.0; // Each photon starts with weight 1.0

			// launch the photon (create a new path)
			launch(photons[p], source);
			
			// CRITICAL FIX: Update source with calculated specular direction from first photon
			// This ensures the renderer can access the correct reflection direction
			if (p == 0 && photons[p].alive) {
				source.specular_direction = photons[p].specular_direction();
			}
			
			// Track launched energy
			if (photons[p].alive) {
				total_launched_energy += photons[p].weight;
			}
			
			// Safety mechanism to prevent infinite photon loops
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
			
			// Progress report AFTER photon simulation is complete - show every 2% for better visibility
			uint32_t progress_interval = std::max(1u, static_cast<uint32_t>(photons.size() / 50)); // Show every 2%
			
			if ((p + 1) % progress_interval == 0 || p == 0 || (p + 1) == photons.size()) {
				// Always show console progress - this is essential user feedback
				double progress_percent = ((double)(p + 1) / photons.size()) * 100.0;
				std::cout << "\rProgress: " << (p + 1) << "/" << Config::get().num_photons() 
						  << " (" << std::fixed << std::setprecision(1) << progress_percent << "%)" << std::flush;
				
				// Always call GUI progress callback
				if (progress_callback_) {
					progress_callback_(p + 1, Config::get().num_photons());
				}
			}
		}
	}

	// Always complete the progress line
	std::cout << std::endl;

	// Normalize physical quantities
	normalize();

	// NOTE: Voxel data output now handled by report() function

	metrics.stop_clock();
	
	// Set simulation completion data for accurate reporting
	metrics.set_simulation_completion(photons.size());
	
	// Increment simulation version since data has changed
	increment_simulation_version();
	
	// Aggregate voxel data to medium records before calculating final metrics
	aggregate_voxel_data();
	
    // Use consolidated energy aggregation method
	auto energy_data = aggregate_medium_energy_data();
	
	// Set only the scatter events in the main simulator metrics (path data is now handled by medium metrics)
	metrics.set_scatter_events(energy_data.scatter_events);
	
	metrics.collect_data(energy_data.total_absorption, energy_data.specular_reflection, energy_data.diffuse_reflection,
						 energy_data.surface_refraction, energy_data.specular_transmission, energy_data.diffuse_transmission);
	
	// Update shared metrics for GUI components
	if (shared_metrics_) {
		shared_metrics_->collect_data(energy_data.total_absorption, energy_data.specular_reflection, energy_data.diffuse_reflection,
									  energy_data.surface_refraction, energy_data.specular_transmission, energy_data.diffuse_transmission);
		shared_metrics_->set_scatter_events(energy_data.scatter_events);
		shared_metrics_->stop_clock(); // Update timing for GUI
		shared_metrics_->set_simulation_completion(photons.size()); // Set completion data for GUI
	}
	
	metrics.print_report(*this);
	
	return Result<void, SimulationError>::ok();
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
			if (new_photon.weight > 0.0 && new_photon.voxel && new_photon.voxel->material) {
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
	
	// Use consolidated energy aggregation method
	auto energy_data = aggregate_medium_energy_data();
	
	// Set only the scatter events in the main simulator metrics (path data is now handled by medium metrics)
	metrics.set_scatter_events(energy_data.scatter_events);
	
	metrics.collect_data(energy_data.total_absorption, energy_data.specular_reflection, energy_data.diffuse_reflection,
						 energy_data.surface_refraction, energy_data.specular_transmission, energy_data.diffuse_transmission);
	
	// Update shared metrics if available (for GUI components)
	if (shared_metrics_) {
		shared_metrics_->collect_data(energy_data.total_absorption, energy_data.specular_reflection, energy_data.diffuse_reflection,
									  energy_data.surface_refraction, energy_data.specular_transmission, energy_data.diffuse_transmission);
		shared_metrics_->set_scatter_events(energy_data.scatter_events);
	}
	
	// Increment simulation version since data has changed
	increment_simulation_version();
}

/***********************************************************
 * Set up photon properties for tracing.
 ***********************************************************/
void Simulator::launch(Photon& photon, const Source& source) {
	// Find the medium that the photon will start in
	Medium* start_medium = find_medium_at_with_dda(source.intersect);
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
	
	// Initialize scatter count for detailed tracking
	photon.scatter_count = 0;
	photon.scatters = false;
	
	// Create detailed tracking data for this photon
	DetailedPhotonData detailed_data;
	detailed_data.id = photon.id;
	detailed_data.entrance_position = source.origin;
	detailed_data.entrance_direction = source.direction;
	detailed_data.initial_weight = photon.weight;
	detailed_data.remaining_weight = photon.weight;
	detailed_data.has_exit = false;
	detailed_data.exited_medium = false;
	detailed_data.termination_reason = "initialized";
	
	// Ensure detailed_photon_data has enough space
	if (detailed_photon_data.size() <= photon.id) {
		detailed_photon_data.resize(photon.id + 1);
	}
	detailed_photon_data[photon.id] = detailed_data;
	
	// Copy source data directly into photon
	photon.source.origin = source.origin;
	photon.source.direction = source.direction;
	photon.source.specular_direction = source.specular_direction;
	photon.source.intersect = source.intersect;
	photon.source.triangle = source.triangle;
	
	photon.direction = source.direction;
	photon.position = source.intersect;
	photon.voxel = start_medium->voxel_at(photon.position);

	// Log photon launch
	glm::ivec3 voxel_coords = photon.voxel ? 
		glm::ivec3(photon.voxel->ix(), photon.voxel->iy(), photon.voxel->iz()) : glm::ivec3(-1);
	Logger::instance().log_photon_event(
		static_cast<int>(photon.id), "LAUNCH", photon.position, photon.direction, 
		photon.weight, voxel_coords, 0, 0.0,
		"Photon launched into medium"
	);

	// Compute specular reflection for this photon
	specular_reflection(photon);

	// CRITICAL FIX: If initial voxel_at failed but specular_reflection succeeded,
	// try to get the voxel again using the same nudging approach
	if (!photon.voxel && photon.alive) {
		// Use the same nudging approach as in specular_reflection
		const double epsilon = 1e-6;
		glm::dvec3 nudged_position = photon.position + epsilon * photon.direction;
		
		auto* nudged_medium = find_medium_at(nudged_position);
		if (nudged_medium) {
			photon.voxel = nudged_medium->voxel_at(nudged_position);
			if (photon.voxel) {
				// Update position to the nudged position for consistent transport
				photon.position = nudged_position;
			}
		}
	}

	// CRITICAL FIX: After specular reflection processing, move photon slightly
	// inside the medium to avoid boundary issues in subsequent transport steps
	if (photon.alive && photon.voxel) {
		// Move photon a small distance along the original direction to get off the exact surface
		const double surface_epsilon = 1e-6;
		glm::dvec3 nudged_position = photon.position + surface_epsilon * photon.direction;
		
		// Verify the nudged position is still in the medium
		auto* nudged_medium = find_medium_at(nudged_position);
		if (nudged_medium) {
			photon.position = nudged_position;
			photon.voxel = nudged_medium->voxel_at(photon.position);
		}
	}

	// create vertices for new light path
	auto light = std::make_shared<PhotonNode>(photon.source.origin, photon.weight);
	auto intersection = std::make_shared<PhotonNode>(photon.source.intersect, photon.weight);
	auto reflection = std::make_shared<PhotonNode>(move(photon.source.intersect, photon.source.specular_direction, 0.1),
												   start_medium->get_metrics().get_surface_reflection());

	light->next = intersection;      // intersection vertex/node
	intersection->prev = light;      // light source origin
	intersection->emit = reflection; // specular reflection

	// Initialize photon's internal path tracking
	photon.path_head = intersection;
	photon.path_last = intersection;
	photon.num_seg_int = 1;
	photon.num_seg_ext = 1;

	metrics.add_vertex(photon.position.x, photon.position.y, photon.position.z);
}

/***********************************************************
 * Set dimensionless step size for next segment of the
 * random walk using MCML 3.0.0 algorithm.
 ***********************************************************/
void Simulator::step_size(Photon& photon) {
	// Find current medium for the photon
	Medium* current_medium = find_medium_at_with_dda(photon.position);
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
	if (!photon.voxel || !photon.voxel->material) {
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
	Medium* current_medium = find_medium_at_with_dda(start_pos);
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
		
		// Track surface voxels for exit recording (ONLY external surfaces, not internal boundaries)
		if (current_voxel && current_voxel->is_surface_voxel) {
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
		if (voxel && voxel->material && distance > 1e-12) {
			// Calculate effective volume fraction for boundary voxels
			double effective_volume_fraction = 1.0;
			if (voxel->is_boundary_voxel) {
				effective_volume_fraction = voxel->volume_fraction_inside;
			}
			
			// Apply Beer-Lambert law for this segment
			double mu_a = voxel->material->mu_a() * effective_volume_fraction;
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
		
		// DEBUG: Track photon state at each step (optimized logging)
		if (substep_counter <= 50) { // Limit debug output
			std::ostringstream oss;
			oss << "Step " << substep_counter << ": Photon " << photon.id 
			    << " weight=" << photon.weight << " step=" << photon.step 
			    << " pos=(" << photon.position.x << "," << photon.position.y << "," << photon.position.z << ")"
			    << " alive=" << (photon.alive ? "yes" : "no");
			FAST_LOG_DEBUG(oss.str());
		}
		
		if (substep_counter > max_substeps) {
			FAST_LOG_WARNING("Photon exceeded maximum substeps, terminating.");
			
			// Deposit remaining energy as absorption for energy conservation
			if (photon.weight > 0.0 && photon.voxel && photon.voxel->material) {
				// EXPERIMENTAL: Use energy conservation enforcement
				terminate_photon_and_record_energy(photon, "max_iterations");
			} else {
				photon.alive = false;
			}
			break;
		}

	// set substep - this will handle mesh boundary detection
	sub_step(photon);

	// TESTING: Disable double absorption to test interface energy splitting fix
	// UNIFIED ABSORPTION TRACKING: Use actual photon path segments for maximum accuracy
	// This ensures absorption follows the same path that gets rendered
	// track_photon_path_segments_for_absorption(photon);
	
	// TRADITIONAL ABSORPTION: Per-step absorption
	deposit(photon);

	// possibly cross boundary
	if (photon.cross) {
		cross(photon);
	}
	else {
		photon.position = move(photon.position, photon.direction, photon.sub_step);
	}
	
	// prevent errors due to crossing to ambient medium
	Medium* current_medium = find_medium_at_with_dda(photon.position);
	if (!current_medium) {
		// CRITICAL: Record as transmission - photon is exiting
		glm::ivec3 last_voxel_coords = photon.voxel ? 
			glm::ivec3(photon.voxel->ix(), photon.voxel->iy(), photon.voxel->iz()) : glm::ivec3(-1);
		
		Logger::instance().log_photon_event(
			static_cast<int>(photon.id), "EXIT", photon.position, photon.direction, 
			photon.weight, last_voxel_coords, -1, photon.weight,
			"Photon exiting medium - calling radiate"
		);
		
		photon.alive = false;
		radiate(photon, photon.direction, photon.weight);
		return;
	}
	
	// ROBUST BOUNDARY HANDLING: Update photon's voxel to match new position
	// BUT preserve material voxel assignment when photon is at exit boundaries
	Voxel* new_voxel = current_medium->voxel_at(photon.position);
	
	if (new_voxel) {
		// Normal case: photon is in a valid material voxel
		photon.voxel = new_voxel;
	} else {
		// CRITICAL EXIT BOUNDARY FIX: Photon position maps to null voxel (air)
		// Keep photon assigned to its last material voxel until actual medium exit
		// This prevents premature assignment to air voxels during exit process
		
		if (photon.voxel && photon.voxel->material) {
			// Photon still has a valid material voxel from previous step
			// Check if photon is actually exiting the medium geometry
			if (!current_medium->contains_point(photon.position)) {
				// Photon has truly exited the medium - proceed with exit logic
				
				glm::ivec3 last_voxel_coords(-1);
				Logger::instance().log_photon_event(
					static_cast<int>(photon.id), "EXIT", photon.position, photon.direction, 
					photon.weight, last_voxel_coords, -1, photon.weight,
					"Photon moved to invalid voxel - calling radiate"
				);
				
				photon.alive = false;
				radiate(photon, photon.direction, photon.weight);
				return;
			} else {
				// Photon is still within medium bounds but in transition zone
				// Keep using last material voxel - this prevents air voxel assignment
				// The photon will properly exit on the next iteration
			}
		} else {
			// Photon has no previous material voxel - this is an error state
			
			glm::ivec3 last_voxel_coords(-1);
			Logger::instance().log_photon_event(
				static_cast<int>(photon.id), "EXIT", photon.position, photon.direction, 
				photon.weight, last_voxel_coords, -1, photon.weight,
				"Photon has no material voxel - calling radiate"
			);
			
			photon.alive = false;
			radiate(photon, photon.direction, photon.weight);
			return;
		}
	}

	// update step size using current medium's material properties
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
	// ROBUST VOXEL SELECTION: Always find the correct voxel at current position
	Medium* photon_medium = find_medium_at(photon.position);
	if (!photon_medium) {
		FAST_LOG_ERROR("No medium found at photon position in sub_step()");
		if (photon.weight > 0.0) {
			terminate_photon_and_record_energy(photon, "no_material_properties");
		} else {
			photon.alive = false;
		}
		return;
	}
	
	// Get the correct voxel at the current position
	Voxel* current_voxel = photon_medium->voxel_at(photon.position);
	if (!current_voxel) {
		FAST_LOG_ERROR("No voxel found at photon position in sub_step()");
		// Energy conservation: Deposit remaining energy as absorption before terminating
		if (photon.weight > 0.0) {
			terminate_photon_and_record_energy(photon, "no_voxel_found");
		} else {
			photon.alive = false;
		}
		return;
	}
	
	// Update photon's voxel reference to the correct current voxel
	photon.voxel = current_voxel;
	
	// Validate voxel has material properties
	if (!photon.voxel->material) {
		FAST_LOG_ERROR("Photon voxel has no material properties in sub_step()");
		// Energy conservation: Deposit remaining energy as absorption before terminating
		if (photon.weight > 0.0) {
			terminate_photon_and_record_energy(photon, "no_material_properties");
		} else {
			photon.alive = false;
		}
		return;
	}

	// Get voxel boundaries for the CORRECT voxel
	Cuboid box = voxel_corners(photon.voxel);

	// ROBUST BOUNDARY HANDLING: Check if photon is exactly on a voxel boundary
	bool on_boundary = false;
	glm::dvec3 adjusted_position = photon.position;

	// Check each axis for boundary conditions with proper epsilon tolerance
	const double EPSILON = MathConstants::BOUNDARY_EPSILON;
	double x_diff_min = std::abs(photon.position.x - box.min_point().x);
	double x_diff_max = std::abs(photon.position.x - box.max_point().x);
	double y_diff_min = std::abs(photon.position.y - box.min_point().y);
	double y_diff_max = std::abs(photon.position.y - box.max_point().y);
	double z_diff_min = std::abs(photon.position.z - box.min_point().z);
	double z_diff_max = std::abs(photon.position.z - box.max_point().z);

	if (x_diff_min < EPSILON) {
		adjusted_position.x = box.min_point().x + EPSILON;
		on_boundary = true;
	}
	else if (x_diff_max < EPSILON) {
		adjusted_position.x = box.max_point().x - EPSILON;
		on_boundary = true;
	}

	if (y_diff_min < EPSILON) {
		adjusted_position.y = box.min_point().y + EPSILON;
		on_boundary = true;
	}
	else if (y_diff_max < EPSILON) {
		adjusted_position.y = box.max_point().y - EPSILON;
		on_boundary = true;
	}

	if (z_diff_min < EPSILON) {
		adjusted_position.z = box.min_point().z + EPSILON;
		on_boundary = true;
	}
	else if (z_diff_max < EPSILON) {
		adjusted_position.z = box.max_point().z - EPSILON;
		on_boundary = true;
	}

	// VALIDATION: Ensure adjusted position is actually inside the voxel
	if (adjusted_position.x < box.min_point().x || adjusted_position.x > box.max_point().x ||
		adjusted_position.y < box.min_point().y || adjusted_position.y > box.max_point().y ||
		adjusted_position.z < box.min_point().z || adjusted_position.z > box.max_point().z) {
		if (Config::get().log()) {
			Logger::instance().log_warning("Adjusted position outside voxel bounds. Using fallback.");
		}
		// Fallback: place photon at voxel center
		adjusted_position = glm::dvec3(
			(box.min_point().x + box.max_point().x) * 0.5,
			(box.min_point().y + box.max_point().y) * 0.5,
			(box.min_point().z + box.max_point().z) * 0.5
		);
		on_boundary = true;
	}

	// Create ray from (robustly adjusted) photon position and direction
	Ray ray = Ray(adjusted_position, photon.direction);

	double voxdist = ray.intersect_cuboid_internal(box, photon.intersect, photon.voxel_normal);

	// ROBUST ERROR HANDLING: Multiple fallback strategies if intersection fails
	if (voxdist == std::numeric_limits<double>::max()) {
		if (Config::get().log()) {
			std::cerr << "WARNING: Primary ray-voxel intersection failed. Attempting fallbacks..." << std::endl;
			std::cerr << "  Ray origin: " << ray.origin().x << ", " << ray.origin().y << ", " << ray.origin().z << std::endl;
			std::cerr << "  Ray direction: " << ray.direction().x << ", " << ray.direction().y << ", " << ray.direction().z << std::endl;
			std::cerr << "  Voxel bounds: min(" << box.min_point().x << ", " << box.min_point().y << ", " << box.min_point().z 
					  << ") max(" << box.max_point().x << ", " << box.max_point().y << ", " << box.max_point().z << ")" << std::endl;
		}
		
		// FALLBACK 1: Try from exact voxel center
		glm::dvec3 voxel_center = glm::dvec3(
			(box.min_point().x + box.max_point().x) * 0.5,
			(box.min_point().y + box.max_point().y) * 0.5,
			(box.min_point().z + box.max_point().z) * 0.5
		);
		Ray fallback_ray1 = Ray(voxel_center, photon.direction);
		voxdist = fallback_ray1.intersect_cuboid_internal(box, photon.intersect, photon.voxel_normal);
		
		if (voxdist != std::numeric_limits<double>::max()) {
			if (Config::get().log()) {
				std::cerr << "  Fallback 1 (voxel center) succeeded." << std::endl;
			}
		} else {
			// FALLBACK 2: Use manual boundary calculation
			if (Config::get().log()) {
				std::cerr << "  Fallback 1 failed. Using manual boundary calculation." << std::endl;
			}
			
			// Find which boundary the ray will hit first
			double t_min = std::numeric_limits<double>::max();
			glm::dvec3 hit_point{0.0, 0.0, 0.0};
			glm::dvec3 hit_normal{0.0, 0.0, 0.0};
			
			// Check each face of the voxel cuboid
			std::vector<std::pair<glm::dvec3, glm::dvec3>> faces = {
				{{box.min_point().x, 0, 0}, {-1, 0, 0}}, // Left face
				{{box.max_point().x, 0, 0}, {1, 0, 0}},  // Right face
				{{0, box.min_point().y, 0}, {0, -1, 0}}, // Bottom face
				{{0, box.max_point().y, 0}, {0, 1, 0}},  // Top face
				{{0, 0, box.min_point().z}, {0, 0, -1}}, // Back face
				{{0, 0, box.max_point().z}, {0, 0, 1}}   // Front face
			};
			
			for (const auto& face : faces) {
				glm::dvec3 face_point = face.first;
				glm::dvec3 face_normal = face.second;
				
				double denom = glm::dot(photon.direction, face_normal);
				if (std::abs(denom) > 1e-12) { // Ray not parallel to face
					double t = glm::dot(face_point - adjusted_position, face_normal) / denom;
					if (t > 1e-12 && t < t_min) { // Valid forward intersection
						glm::dvec3 test_point = adjusted_position + t * photon.direction;
						
						// Check if intersection point is within face bounds
						bool within_bounds = true;
						if (face_normal.x != 0) { // YZ face
							within_bounds = (test_point.y >= box.min_point().y - EPSILON && test_point.y <= box.max_point().y + EPSILON &&
											test_point.z >= box.min_point().z - EPSILON && test_point.z <= box.max_point().z + EPSILON);
						} else if (face_normal.y != 0) { // XZ face
							within_bounds = (test_point.x >= box.min_point().x - EPSILON && test_point.x <= box.max_point().x + EPSILON &&
											test_point.z >= box.min_point().z - EPSILON && test_point.z <= box.max_point().z + EPSILON);
						} else { // XY face
							within_bounds = (test_point.x >= box.min_point().x - EPSILON && test_point.x <= box.max_point().x + EPSILON &&
											test_point.y >= box.min_point().y - EPSILON && test_point.y <= box.max_point().y + EPSILON);
						}
						
						if (within_bounds) {
							t_min = t;
							hit_point = test_point;
							hit_normal = face_normal;
						}
					}
				}
			}
			
			if (t_min < std::numeric_limits<double>::max()) {
				voxdist = t_min;
				photon.intersect = hit_point;
				photon.voxel_normal = hit_normal;
				if (Config::get().log()) {
					std::cerr << "  Fallback 2 (manual calculation) succeeded." << std::endl;
				}
			} else {
				// FALLBACK 3: Emergency exit - force photon to exit current voxel
				if (Config::get().log()) {
					std::cerr << "  All fallbacks failed. Forcing photon exit." << std::endl;
				}
				photon.intersect = adjusted_position + photon.direction * EPSILON;
				photon.voxel_normal = -photon.direction; // Opposite to ray direction
				voxdist = EPSILON;
			}
		}
	}

	// Now we have a valid intersection - continue with sub_step logic
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
		// Energy conservation: Deposit remaining energy as absorption before terminating
		if (photon.weight > 0.0) {
			terminate_photon_and_record_energy(photon, "no_medium_found");
		} else {
			photon.alive = false;
		}
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
	// Cancel if photon is outside of medium or doesn't have material
	if (!photon.voxel || !photon.voxel->material) {
		return;
	}

	// For boundary voxels, only deposit in the portion that's inside the geometry
	double effective_volume_fraction = 1.0;
	if (photon.voxel->is_boundary_voxel) {
		// Scale absorption by the volume fraction inside for boundary voxels
		// Don't exit early - the photon is still in a material voxel and energy should be conserved
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

	// DEBUG: Track absorption details for single photon debug (log mode only, limited output)
	static int deposit_debug_count = 0;
	if (Config::get().log() && deposit_debug_count < 10 && deltaw > 0.001) { // Limit to first 10 steps with significant absorption
		std::ostringstream oss;
		oss << "DEPOSIT: Photon " << photon.id << " weight=" << photon.weight 
		    << " mu_a=" << photon.mu_a() << " sub_step=" << photon.sub_step 
		    << " deltaw=" << deltaw << " effective_volume=" << effective_volume_fraction;
		Logger::instance().log_debug(oss.str());
		deposit_debug_count++;
	}

	// ENERGY CONSERVATION ENFORCEMENT
	// Calculate how much energy this photon has left to absorb
	double energy_already_used = photon.total_energy_radiated + photon.total_energy_absorbed;
	double energy_available = photon.total_energy_budget - energy_already_used;
	
	// CRITICAL FIX: Prevent negative energy calculations
	if (energy_available < 0.0) {
		if (Config::get().log()) {
			std::ostringstream oss;
			oss << "Energy available became negative (" << energy_available 
			    << "), budget=" << photon.total_energy_budget 
			    << ", used=" << energy_already_used 
			    << " (radiated=" << photon.total_energy_radiated 
			    << ", absorbed=" << photon.total_energy_absorbed << ")";
			Logger::instance().log_warning(oss.str());
		}
		energy_available = 0.0;
	}
	
	// Enforce energy conservation: cannot absorb more than available
	double actual_absorbed_weight = std::min(deltaw, energy_available);
	
	// DEBUG: Track absorption calculation
	if (Config::get().log() && (actual_absorbed_weight != deltaw || photon.weight < actual_absorbed_weight)) {
		std::ostringstream oss;
		oss << "ABSORPTION: deltaw=" << deltaw 
		    << ", actual=" << actual_absorbed_weight 
		    << ", current_weight=" << photon.weight 
		    << ", after_weight=" << (photon.weight - actual_absorbed_weight);
		Logger::instance().log_debug(oss.str());
	}
	
	// Update photon energy tracking
	photon.total_energy_absorbed += actual_absorbed_weight;

	// update photon weight  
	photon.weight -= actual_absorbed_weight;
	
	// CRITICAL FIX: Prevent negative weights
	if (photon.weight < 0.0) {
		if (Config::get().log()) {
			std::ostringstream debug_msg;
			debug_msg << "WARNING: Photon weight became negative (" << photon.weight << "), setting to 0";
			Logger::instance().log_warning(debug_msg.str());
		}
		photon.weight = 0.0;
	}
	
	// CRITICAL FIX: Terminate zero-weight photons immediately
	if (photon.weight <= 0.0) {
		photon.alive = false;
		return;
	}

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
	// DEBUG: Track boundary crossings (log mode only, limited output)
	static int crossing_count = 0;
	if (Config::get().log() && crossing_count < 20) { // Limit to first 20 crossings
		std::ostringstream debug_msg;
		debug_msg << "=== BOUNDARY CROSSING === Photon " << photon.id << " at pos=(" 
		          << photon.position.x << "," << photon.position.y << "," << photon.position.z 
		          << ") weight=" << photon.weight;
		Logger::instance().log_debug(debug_msg.str());
		crossing_count++;
	}
	
	// Safety check - ensure photon has valid voxel and material
	if (!photon.voxel || !photon.voxel->material) {
		if (Config::get().log()) {
			Logger::instance().log_debug("  -> Outside medium, recording transmission");
		}
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
		double eta = Config::get().ambient_eta();
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
		return;
	}

	// directions of transmission and reflection
	glm::dvec3 transmittance, reflectance;

	// First compute the transmission and reflection directions using current photon state
	double eta_ambient = Config::get().ambient_eta();
	double temp_reflection = internal_reflection(photon, eta_ambient, transmittance, reflectance);

	// Determine which medium(s) are involved in this crossing
	Medium* current_medium = find_medium_at_with_dda(photon.position);
	glm::dvec3 newpos = move_delta(photon.intersect, transmittance);
	Medium* new_medium = find_medium_at_with_dda(newpos);
	Voxel* newvox = new_medium ? new_medium->voxel_at(newpos) : nullptr;
	
	// ROBUST VOXEL LOOKUP: Handle edge cases where move_delta places photon at voxel boundaries
	if (!newvox && new_medium) {
		// Try the intersect point directly first
		Voxel* intersect_voxel = new_medium->voxel_at(photon.intersect);
		if (intersect_voxel) {
			newvox = intersect_voxel;
		} else {
			// Try a slightly larger delta movement
			glm::dvec3 larger_newpos = photon.intersect + transmittance * (Config::get().vox_size() * 0.001);
			if (new_medium->contains_point(larger_newpos)) {
				Voxel* larger_voxel = new_medium->voxel_at(larger_newpos);
				if (larger_voxel) {
					newvox = larger_voxel;
				}
			}
		}
	}
	
	photon.prev_voxel = photon.voxel;

	// determine refractive index of the medium being entered
	double eta = (newvox == nullptr) ? Config::get().ambient_eta() : newvox->material->eta();

	// Recalculate with correct refractive index
	temp_reflection = internal_reflection(photon, eta, transmittance, reflectance);

	// Handle different crossing scenarios
	if (new_medium != current_medium) {
		// Check if this is a true medium exit or just a layer boundary
		if (!new_medium && current_medium) {
			// TRUE MEDIUM EXIT: Photon exiting to ambient space
			handle_medium_transition(photon, current_medium, new_medium);
			
			// If photon was killed by medium transition, call radiate() for actual exit
			if (!photon.alive) {
				radiate(photon, photon.direction, photon.weight);
			}
			return;
		} else if (current_medium && new_medium && current_medium != new_medium) {
			// DIFFERENT MEDIUM TRANSITION: Should not happen with current config
			handle_medium_transition(photon, current_medium, new_medium);
			
			if (!photon.alive) {
				radiate(photon, photon.direction, photon.weight);
			}
			return;
		}
	}
	
	// Check for LAYER BOUNDARY within same medium (Fresnel reflection)
	if (Config::get().log()) {
		std::ostringstream oss;
		oss << "Checking layer boundary - current_medium=" << (current_medium ? "yes" : "no") 
		    << " new_medium=" << (new_medium ? "yes" : "no") 
		    << " same=" << (current_medium == new_medium ? "yes" : "no");
		Logger::instance().log_debug(oss.str());
	}
	
	if (current_medium && new_medium && current_medium == new_medium) {
		Voxel* current_voxel = photon.voxel;
		
		if (Config::get().log()) {
			Logger::instance().log_debug("Same medium detected - checking voxels");
			std::ostringstream oss;
			oss << "  current_voxel=" << (current_voxel ? "yes" : "no")
			    << " newvox=" << (newvox ? "yes" : "no");
			Logger::instance().log_debug(oss.str());
			
			if (current_voxel && newvox) {
				std::ostringstream oss0;
				oss0 << "  current_tissue=" << (current_voxel->material ? "yes" : "no")
				     << " new_tissue=" << (newvox->material ? "yes" : "no");
				Logger::instance().log_debug(oss0.str());
				
			if (current_voxel->material && newvox->material) {
				bool same_optical = current_voxel->material->has_same_optical_properties(*newvox->material);
				
				std::ostringstream oss1;
				oss1 << "  current_material_hash=" << current_voxel->material->get_optical_properties_hash()
				     << " new_material_hash=" << newvox->material->get_optical_properties_hash();
				Logger::instance().log_debug(oss1.str());					std::ostringstream oss2;
					oss2 << "  current_properties: eta=" << current_voxel->material->eta()
					     << " mua=" << current_voxel->material->mu_a()
					     << " mus=" << current_voxel->material->mu_s()  
					     << " g=" << current_voxel->material->g() 
					     << " hash=" << current_voxel->material->get_optical_properties_hash();
					Logger::instance().log_debug(oss2.str());
					
					std::ostringstream oss3;
					oss3 << "  new_properties: eta=" << newvox->material->eta()
					     << " mua=" << newvox->material->mu_a()
					     << " mus=" << newvox->material->mu_s()
					     << " g=" << newvox->material->g()
					     << " hash=" << newvox->material->get_optical_properties_hash();
					Logger::instance().log_debug(oss3.str());
					
					std::ostringstream oss4;
					oss4 << "  optical_properties_same=" << (same_optical ? "yes" : "no");
					Logger::instance().log_debug(oss4.str());
				}
			}
		}
		
		if (current_voxel && newvox && 
			current_voxel->material && newvox->material && 
			!current_voxel->material->has_same_optical_properties(*newvox->material)) {
			
			// INTERFACE ENERGY SPLITTING - Simple implementation
			double n1 = current_voxel->material->eta();  // From medium
			double n2 = newvox->material->eta();          // To medium
			
			// Calculate angle of incidence
			double cos_theta_i = -glm::dot(photon.direction, photon.voxel_normal);
			cos_theta_i = glm::clamp(cos_theta_i, 0.0, 1.0); // Ensure valid range
			
			// Check for total internal reflection
			double n_ratio = n1 / n2;
			double sin_theta_t_sq = n_ratio * n_ratio * (1.0 - cos_theta_i * cos_theta_i);
			
			if (sin_theta_t_sq > 1.0) {
				// TOTAL INTERNAL REFLECTION - all energy stays in current medium
				// Reflect photon direction
				photon.direction = photon.direction - 2.0 * glm::dot(photon.direction, photon.voxel_normal) * photon.voxel_normal;
				photon.direction = glm::normalize(photon.direction);
				
				if (Config::get().log()) {
					std::ostringstream debug_msg;
					debug_msg << "  -> TOTAL INTERNAL REFLECTION: n1=" << n1 << ", n2=" << n2;
					Logger::instance().log_debug(debug_msg.str());
				}
				return; // Photon reflects back, no interface crossing
			}
			
			// Calculate Fresnel reflection coefficient
			double cos_theta_t = std::sqrt(1.0 - sin_theta_t_sq);
			double R_fresnel;
			
			if (cos_theta_i < 1e-6) {
				// Normal incidence (simplified)
				R_fresnel = std::pow((n1 - n2) / (n1 + n2), 2.0);
			} else {
				// General case - Fresnel equations for s and p polarizations
				double Rs = std::pow((n1 * cos_theta_i - n2 * cos_theta_t) / (n1 * cos_theta_i + n2 * cos_theta_t), 2.0);
				double Rp = std::pow((n2 * cos_theta_i - n1 * cos_theta_t) / (n2 * cos_theta_i + n1 * cos_theta_t), 2.0);
				R_fresnel = 0.5 * (Rs + Rp); // Average for unpolarized light
			}
			
			// Ensure valid reflection coefficient (critical safety check)
			R_fresnel = glm::clamp(R_fresnel, 0.0, 1.0);
			double T_fresnel = 1.0 - R_fresnel;
			
			// ENERGY SPLITTING
			double initial_weight = photon.weight;
			double reflected_energy = initial_weight * R_fresnel;   // Deposited as absorption
			double transmitted_energy = initial_weight * T_fresnel; // Photon continues
			
			// Deposit reflected energy as absorption in current voxel (last voxel before interface)
			if (current_voxel && reflected_energy > 0.0) {
				current_voxel->absorption += reflected_energy;
			}
			
			// Continue photon with transmitted energy
			photon.weight = transmitted_energy;
			
			// Calculate refracted direction using Snell's law
			glm::dvec3 incident = photon.direction;
			glm::dvec3 normal = photon.voxel_normal;
			
			if (cos_theta_i > 0.9999) {
				// Near-normal incidence - no significant refraction
				photon.direction = incident;
			} else {
				// Apply Snell's law for refraction
				glm::dvec3 refracted_tangent = n_ratio * (incident - cos_theta_i * normal);
				glm::dvec3 refracted_direction = refracted_tangent + cos_theta_t * normal;
				photon.direction = glm::normalize(refracted_direction);
			}
			
			if (Config::get().log()) {
				std::ostringstream debug_msg;
				debug_msg << "  -> INTERFACE ENERGY SPLITTING: n1=" << n1 << ", n2=" << n2 
				          << ", R=" << R_fresnel << ", T=" << T_fresnel
				          << ", reflected=" << reflected_energy << ", transmitted=" << transmitted_energy;
				Logger::instance().log_debug(debug_msg.str());
			}
		}
	}

	// 1. crossing to ambient medium
	if (newvox == nullptr) {
		// CRITICAL FIX: Only treat as ambient exit if photon is actually leaving the geometry
		// Check if the new position is truly outside all medium geometries
		bool truly_exiting_geometry = true;
		for (const auto& medium : mediums) {
			if (medium.contains_point(newpos)) {
				truly_exiting_geometry = false;
				break;
			}
		}
		
		if (!truly_exiting_geometry) {
			// FALSE AMBIENT EXIT: Photon is still inside geometry but newvox is null
			// This is a DDA traversal issue, not a true exit - continue transport normally
			photon.position = newpos;
			photon.voxel = voxel_at(photon.position);  // Find correct voxel at new position
			return;
		}
		
		// TRUE AMBIENT EXIT: Photon is genuinely leaving the geometry
		// CRITICAL: For exit detection, photon MUST remain assigned to surface voxel
		// The voxel assignment that happens during transport can assign photons to air voxels
		// But for emittance recording, we need the photon to be associated with the surface voxel it's exiting from
		
		// If photon.voxel is already a surface voxel, keep it - this is correct
		// If photon.voxel is not a surface voxel, try to find the correct surface voxel
		if (photon.voxel && !photon.voxel->is_surface_voxel) {
			// Photon was assigned to wrong voxel during transport - need to find correct exit voxel
			// Use the intersection point to find the surface voxel we're actually exiting from
			Medium* exit_medium = find_medium_at_with_dda(photon.position);
			if (exit_medium) {
				// Look for a surface voxel near the intersection point
				Voxel* intersection_voxel = exit_medium->voxel_at(photon.intersect);
				if (intersection_voxel && intersection_voxel->is_surface_voxel) {
					photon.voxel = intersection_voxel;  // Correct the assignment
				} else {
					// Fallback: search for a nearby surface voxel
					Voxel* surface_voxel = find_last_surface_voxel_with_dda(photon, transmittance);
					if (surface_voxel) {
						photon.voxel = surface_voxel;
					}
				}
			}
		}
		
		// Now check if we have a proper surface voxel for exit
		if (photon.voxel && !photon.voxel->is_surface_voxel) {
			// VALIDATION: Check if photon position is actually outside the medium geometry
			// This helps us understand if the problem is voxelization or transport
			bool position_outside_geometry = true;
			
			// Check if current photon position is inside any medium
			for (const auto& medium : mediums) {
				if (medium.contains_point(photon.position)) {
					position_outside_geometry = false;
					break;
				}
			}
			
			// Only show detailed warnings in log mode
			if (Config::get().log()) {
				if (position_outside_geometry) {
					// Photon position is legitimately outside - voxelization error
					Logger::instance().log_error("VOXELIZATION ERROR: Photon at position outside geometry but voxel marked as interior!");
					Logger::instance().log_error("Voxel (" 
						+ std::to_string(photon.voxel->ix()) + ", " 
						+ std::to_string(photon.voxel->iy()) + ", " 
						+ std::to_string(photon.voxel->iz()) + ") should be surface but isn't.");
				}
				else {
					// Photon position is inside geometry - transport/exit detection issue
					Logger::instance().log_warning("TRANSPORT ISSUE: Photon trying to exit from position still inside geometry!");
					Logger::instance().log_warning("Position (" 
						+ std::to_string(photon.position.x) + ", " 
						+ std::to_string(photon.position.y) + ", " 
						+ std::to_string(photon.position.z) + ") is inside medium but exit attempted.");
				}
				
				Logger::instance().log_warning("Warning: Photon attempting to exit from interior voxel at (" 
					+ std::to_string(photon.voxel->ix()) + ", " 
					+ std::to_string(photon.voxel->iy()) + ", " 
					+ std::to_string(photon.voxel->iz()) + ") to ambient medium!");
				
				Logger::instance().log_warning("New position: (" 
					+ std::to_string(newpos.x) + ", " 
					+ std::to_string(newpos.y) + ", " 
					+ std::to_string(newpos.z) + ")");
				
				Logger::instance().log_warning("Direction: (" 
					+ std::to_string(transmittance.x) + ", " 
					+ std::to_string(transmittance.y) + ", " 
					+ std::to_string(transmittance.z) + ")");
			}
		}
		
		// Handle outward-pointing surface normals correctly
		// For Monte Carlo transport, we need to determine if photon is entering or exiting
		// Note: voxel_normal is already computed and available for use
		
		// For exiting photon, voxel_normal should point outward (as defined in mesh)
		// For entering photon, we may need to consider the geometry context
		// Since we're at an exterior boundary, photon is likely exiting
		// No normal flipping needed - outward normals are correct for exit calculations
		
		// Recalculate reflection and transmission with correct normal
		glm::dvec3 surface_normal = glm::normalize(photon.voxel_normal);
		glm::dvec3 incident_dir = glm::normalize(photon.direction);
		
		// Reflection: incident ray bounces back into medium
		glm::dvec3 corrected_reflectance = incident_dir - 2.0 * glm::dot(incident_dir, surface_normal) * surface_normal;
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
				// For reflection, photon stays in same medium - safe to update voxel
				photon.voxel = current_medium->voxel_at(photon.position);
			}
		}
		else {
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
					// For reflection, photon stays in same medium - safe to update voxel
					photon.voxel = current_medium->voxel_at(photon.position);
				}
			} else {
				// No reflection, terminate photon
				photon.alive = false;
			}
		}

		if (current_medium) {
			current_medium->get_metrics().increment_scatters();
		}
	}
	// 2. crossing to another medium
	else if (newvox != nullptr 
		&& newvox->material != nullptr 
		&& photon.voxel->material != nullptr
		&& photon.voxel->material->eta() != newvox->material->eta())
	{
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
			// CRITICAL: Smart voxel assignment to preserve surface voxel information
			Voxel* new_voxel = current_medium->voxel_at(photon.position);
			
			// Only update voxel assignment if we're sure it's correct
			// If current voxel is a surface voxel and new voxel is null/air, preserve current
			if (new_voxel && new_voxel->material) {
				photon.voxel = new_voxel;  // Safe to assign - it's a material voxel
			} else if (!photon.voxel || !photon.voxel->is_surface_voxel) {
				// Only assign null/air voxels if current voxel isn't a surface voxel
				photon.voxel = new_voxel;
			}
			// If new_voxel is null/air but photon.voxel is surface voxel, preserve it
		}
	}
}

/***********************************************************
 * Find the last surface voxel that a photon traveled through
 * using DDA traversal from photon's last known position.
 ***********************************************************/
Voxel* Simulator::find_last_surface_voxel_with_dda(const Photon& photon, const glm::dvec3& exit_direction) {
	// Find the medium the photon was in
	Medium* current_medium = find_medium_at(photon.position);
	if (!current_medium) {
		return nullptr;
	}
	
	// Find the DDA instance for this medium
	size_t medium_index = 0;
	for (size_t i = 0; i < mediums.size(); ++i) {
		if (&mediums[i] == current_medium) {
			medium_index = i;
			break;
		}
	}
	
	if (medium_index >= medium_ddas_.size()) {
		// Fallback: try to find voxel with step-back approach
		glm::dvec3 step_back_pos = photon.position - exit_direction * 1e-6;
		Medium* last_medium = find_medium_at_with_dda(step_back_pos);
		if (last_medium) {
			return last_medium->voxel_at(step_back_pos);
		}
		return nullptr;
	}
	
	VoxelDDA3D* dda = medium_ddas_[medium_index].get();
	
	// Trace backward from current position to find the path
	glm::dvec3 start_pos = photon.position - exit_direction * 0.01; // Start slightly inside medium
	glm::dvec3 end_pos = photon.position;
	glm::dvec3 direction = glm::normalize(end_pos - start_pos);
	double total_distance = glm::length(end_pos - start_pos);
	
	if (total_distance < 1e-12) {
		// Very short distance, use current voxel
		return photon.voxel;
	}
	
	// Initialize DDA for this ray
	dda->initialize_ray(start_pos, direction);
	
	// Traverse voxels using DDA
	VoxelDDA3D::TraversalResult result = dda->traverse(total_distance);
	
	// Find the last surface voxel in the traversal
	Voxel* last_surface_voxel = nullptr;
	
	Logger::instance().log_photon_event(
		static_cast<int>(photon.id), "DDA_SEARCH", photon.position, direction, 
		0.0, glm::ivec3(-1), -1, total_distance,
		"Starting DDA traversal for last surface voxel, found " + std::to_string(result.voxels.size()) + " voxels"
	);
	
	for (size_t i = 0; i < result.voxels.size(); ++i) {
		const auto& step = result.voxels[i];
		// Get voxel at this DDA position
		glm::dvec3 mutable_pos = step.world_position; // Create mutable copy for voxel_at
		Voxel* voxel = current_medium->voxel_at(mutable_pos);
		
		if (voxel && voxel->material) {
			glm::ivec3 voxel_coords = glm::ivec3(voxel->ix(), voxel->iy(), voxel->iz());
			bool is_surface = voxel->is_surface_voxel; // ONLY true external surface voxels
			
			Logger::instance().log_photon_event(
				static_cast<int>(photon.id), "DDA_VOXEL", step.world_position, direction, 
				0.0, voxel_coords, static_cast<int>(i), step.distance_traveled,
				"Voxel " + std::to_string(i) + "/" + std::to_string(result.voxels.size()) + 
				(is_surface ? std::string(" SURFACE") : std::string(" INTERIOR"))
			);
			
			if (is_surface) {
				last_surface_voxel = voxel;
			}
		}
	}
	
	// If no surface voxel found in DDA traversal, fall back to photon's current voxel
	Voxel* selected_voxel = last_surface_voxel ? last_surface_voxel : photon.voxel;
	
	if (selected_voxel) {
		glm::ivec3 selected_coords = glm::ivec3(selected_voxel->ix(), selected_voxel->iy(), selected_voxel->iz());
		Logger::instance().log_photon_event(
			static_cast<int>(photon.id), "DDA_RESULT", photon.position, direction, 
			0.0, selected_coords, -1, 0.0,
			"Selected voxel: " + (last_surface_voxel ? std::string("DDA_FOUND") : std::string("FALLBACK_TO_CURRENT"))
		);
	} else {
		Logger::instance().log_photon_event(
			static_cast<int>(photon.id), "DDA_ERROR", photon.position, direction, 
			0.0, glm::ivec3(-1), -1, 0.0,
			"No voxel selected - both DDA and fallback failed"
		);
	}
	
	return selected_voxel;
}

/***********************************************************
 * Record the emittance from a photon (partially) leaving
 * the material.
 ***********************************************************/
void Simulator::radiate(Photon& photon, glm::dvec3& direction, double weight) {
	// CRITICAL FIX: Find the LAST material VOXEL before exit using robust boundary handling
	// The intersection point tells us exactly where the photon crossed the boundary
	// We need to find the material voxel that is most inside the medium near this intersection
	
	Medium* exit_medium = find_medium_at_with_dda(photon.position);
	if (!exit_medium) {
		return; // No medium found
	}
	
	// ROBUST VOXEL SELECTION: Handle numerical instability at voxel boundaries
	// When intersection is at boundary, we need the voxel that's most "inside" the medium
	
	Voxel* exit_voxel = nullptr;
	double voxel_size = exit_medium->get_volume().voxel_size();
	
	// Strategy 1: Sample multiple points around the intersection to find best material voxel
	std::vector<std::pair<Voxel*, double>> candidate_voxels;
	
	// Sample points slightly inside the medium from intersection
	glm::dvec3 reverse_direction = -glm::normalize(direction);
	for (int i = 1; i <= 5; ++i) {
		double epsilon = (voxel_size * 0.1) * i; // Progressive steps back into medium
		glm::dvec3 sample_pos = photon.intersect + reverse_direction * epsilon;
		
		Voxel* candidate = exit_medium->voxel_at(sample_pos);
		if (candidate && candidate->material) {
			// Calculate how "deep" this voxel is inside the medium
			// Voxels closer to intersection but still inside get higher priority
			double depth_score = 1.0 / (epsilon + 1e-9); // Higher score for smaller epsilon
			candidate_voxels.push_back({candidate, depth_score});
		}
	}
	
	// Strategy 2: If no good candidates, check voxel neighbors around intersection
	if (candidate_voxels.empty()) {
		// Get voxel coordinates of intersection point
		glm::dvec3 grid_origin = exit_medium->get_bounds().min_bounds;
		glm::ivec3 intersection_coords = glm::ivec3(
			(photon.intersect.x - grid_origin.x) / voxel_size,
			(photon.intersect.y - grid_origin.y) / voxel_size,
			(photon.intersect.z - grid_origin.z) / voxel_size
		);
		
		// Check neighboring voxels (3x3x3 neighborhood)
		for (int dx = -1; dx <= 1; ++dx) {
			for (int dy = -1; dy <= 1; ++dy) {
				for (int dz = -1; dz <= 1; ++dz) {
					glm::ivec3 neighbor_coords = intersection_coords + glm::ivec3(dx, dy, dz);
					
					// Check bounds
					if (neighbor_coords.x >= 0 && static_cast<uint32_t>(neighbor_coords.x) < exit_medium->get_volume().width() &&
						neighbor_coords.y >= 0 && static_cast<uint32_t>(neighbor_coords.y) < exit_medium->get_volume().height() &&
						neighbor_coords.z >= 0 && static_cast<uint32_t>(neighbor_coords.z) < exit_medium->get_volume().depth()) {
						
						Voxel* neighbor = exit_medium->get_volume().at(neighbor_coords.x, neighbor_coords.y, neighbor_coords.z);
						if (neighbor && neighbor->material) {
							// Calculate distance from intersection to voxel center
							glm::dvec3 voxel_center = grid_origin + glm::dvec3(
								(neighbor_coords.x + 0.5) * voxel_size,
								(neighbor_coords.y + 0.5) * voxel_size,
								(neighbor_coords.z + 0.5) * voxel_size
							);
							double distance = glm::length(photon.intersect - voxel_center);
							double proximity_score = 1.0 / (distance + 1e-9);
							candidate_voxels.push_back({neighbor, proximity_score});
						}
					}
				}
			}
		}
	}
	
	// Strategy 3: Fallback to photon's last known position
	if (candidate_voxels.empty()) {
		exit_voxel = exit_medium->voxel_at(photon.position);
		if (!exit_voxel || !exit_voxel->material) {
			exit_voxel = photon.voxel;
		}
	} else {
		// Select the candidate with highest score (most inside the medium)
		std::sort(candidate_voxels.begin(), candidate_voxels.end(), 
			[](const auto& a, const auto& b) { return a.second > b.second; });
		exit_voxel = candidate_voxels[0].first;
	}
	
	if (!exit_voxel || !exit_voxel->material) {
		std::cerr << "ERROR: Cannot determine last material voxel for emittance recording" << std::endl;
		return;
	}
	
	glm::ivec3 surface_coords = glm::ivec3(exit_voxel->ix(), exit_voxel->iy(), exit_voxel->iz());
	
	// Log the radiate event with the determined exit voxel
	Logger::instance().log_photon_event(
		static_cast<int>(photon.id), "RADIATE", photon.position, direction, 
		weight, surface_coords, -1, weight,
		"Using robust last material voxel selection"
	);
	
	// REMOVED: Old validation that blocked emittance recording
	// We now trust that exit_voxel is the correct last material voxel
	
	// Record photon's radiate() call origin
	// ENERGY CONSERVATION ENFORCEMENT (always True Splitting mode now)
	// True Splitting handles energy conservation through statistical splitting
	photon.radiate_call_count++;
	photon.total_energy_radiated += weight;

	// Record emittance at the LAST material VOXEL (before exit)
	double old_emittance = exit_voxel->emittance;
	(void)old_emittance; // Suppress unused variable warning - kept for debugging
	exit_voxel->emittance += weight;
	
	// Log voxel emittance recording
	Logger::instance().log_voxel_emittance(
		static_cast<int>(photon.id), photon.position, direction, weight, surface_coords, weight, "Surface voxel emittance"
	);
	
	// Use proper reflection/transmission determination based on exit position relative to entry
	bool is_reflecting = is_photon_reflecting(photon);
	if (is_reflecting) {
		// Exit on same side as entry - classify as diffuse reflection (not specular)
		exit_voxel->diffuse_reflection += weight;
		photon.exit_type = Photon::ExitType::REFLECTED;
		Logger::instance().log_photon_event(
			static_cast<int>(photon.id), "REFLECT", photon.position, direction, 
			weight, surface_coords, -1, weight,
			"Photon classified as reflection"
		);
	} else {
		// Exit on opposite side from entry - classify as transmission  
		exit_voxel->diffuse_transmission += weight;
		photon.exit_type = Photon::ExitType::TRANSMITTED;
		Logger::instance().log_photon_event(
			static_cast<int>(photon.id), "TRANSMIT", photon.position, direction, 
			weight, surface_coords, -1, weight,
			"Photon classified as transmission"
		);
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
	if (photon.path_last) {
		exit_node = std::make_shared<PhotonNode>(photon.intersect, weight, 
			static_cast<PhotonNode::ExitType>(node_exit_type));
		photon.add_external_vertex(exit_node);
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
	
	// Update detailed photon tracking data for exit
	if (photon.id < detailed_photon_data.size()) {
		DetailedPhotonData& detailed = detailed_photon_data[photon.id];
		detailed.has_exit = true;
		detailed.exit_position = photon.intersect;
		detailed.exit_direction = direction;
		detailed.exited_medium = true;
		detailed.scatter_count = photon.scatter_count;
		detailed.total_absorption_deposited = photon.total_energy_absorbed;
		detailed.remaining_weight = weight;
		detailed.termination_reason = (photon.exit_type == Photon::ExitType::REFLECTED) ? "reflected" : "transmitted";
	}
	
	// Establish bidirectional connection between exit node and emitter
	if (exit_node) {
		exit_node->emitter = emitter;
	}
}

/***********************************************************
 * Determine the survivability of a given photon using MCML 3.0.0 algorithm.
 ***********************************************************/
void Simulator::roulette(Photon& photon) {
	// CRITICAL FIX: Never apply Russian Roulette to negative weights
	if (photon.weight < 0.0) {
		if (Config::get().log()) {
			std::ostringstream debug_msg;
			debug_msg << "ERROR: Attempted Russian Roulette on negative weight (" << photon.weight << "), terminating photon";
			Logger::instance().log_error(debug_msg.str());
		}
		terminate_photon_and_record_energy(photon, "negative_weight");
		return;
	}
	
	// Modernized Russian roulette with adaptive threshold and survival probability
	if (photon.weight < mcml_weight_threshold) {
		// Use weight-dependent survival probability for better variance reduction
		double survival_probability = std::max(0.1, photon.weight / mcml_weight_threshold);
		survival_probability = std::min(survival_probability, 0.5); // Cap at 50% for stability

		if (rng->next() <= survival_probability) {
			// ENERGY CONSERVATION FIX: Proper Russian Roulette
			
			// Survive with proper weight normalization
			double old_weight = photon.weight;
			photon.weight /= survival_probability;
			
			// CRITICAL FIX: Update energy budget to match new weight
			// Russian Roulette increases the weight, so budget must increase proportionally
			photon.total_energy_budget = photon.weight;
			
			if (Config::get().log()) {
				std::ostringstream oss;
				oss << "ROULETTE: weight " << old_weight << " -> " << photon.weight 
				    << " (survival_prob=" << survival_probability << ", new_budget=" << photon.total_energy_budget << ")";
				Logger::instance().log_debug(oss.str());
			}
		}
		else {
			// Use centralized termination for consistency
			terminate_photon_and_record_energy(photon, "roulette");
		}
	}
}

/***********************************************************
 * Scatter the photon into a new direction based on the
 * Henyey-Greenstein phase function using MCML 3.0.0 algorithm.
 ***********************************************************/
void Simulator::scatter(Photon& photon) {
	if (!photon.alive) {
		return;
	}

	// Get current medium and material properties for scattering
	Medium* current_medium = find_medium_at(photon.position);
	if (!current_medium) {
		// Record as transmission - photon has exited medium
		photon.alive = false;
		radiate(photon, photon.direction, photon.weight);
		return;
	}
	
	Material* material = photon.voxel->material;
	if (!material) {
		// Record as transmission - photon cannot scatter without material
		photon.alive = false;
		radiate(photon, photon.direction, photon.weight);
		return;
	}

	// Modern numerically stable Henyey-Greenstein phase function implementation
	double cos_theta, sin_theta, cos_phi, sin_phi;
	double g = material->g(); // anisotropy factor
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
	
	// Increment scatter count for this photon
	photon.scatter_count++;
	photon.scatters = true;  // Mark that this photon has scattered at least once

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
	photon.add_internal_vertex(std::make_shared<PhotonNode>(photon.position, photon.weight));
}

/***********************************************************
 * Normalize the recorded values based on the number of
 * photons traced.
 ***********************************************************/
void Simulator::normalize() {
	// Normalize records for all mediums
	double normalization_factor = static_cast<double>(Config::get().num_photons()) * static_cast<double>(Config::get().num_sources());
	for (auto& medium : mediums) {
		medium.get_metrics().normalize_raw_values(normalization_factor);
	}

	// normalize voxel data across all mediums
	for (auto& medium : mediums) {
		auto& voxel_grid = medium.get_volume();
		for (const auto& voxel_ptr : voxel_grid) {
			auto* voxel = voxel_ptr.get();
			// skip computation for voxels outside the medium
			if (!voxel->material) {
				continue;
			}

			voxel->absorption /= (Config::get().num_photons() * Config::get().num_sources());
			voxel->emittance /= (Config::get().num_photons() * Config::get().num_sources());
			
			// CRITICAL FIX: Also normalize the directional emittance fields used by energy conservation
			voxel->specular_reflection /= (Config::get().num_photons() * Config::get().num_sources());
			voxel->diffuse_transmission /= (Config::get().num_photons() * Config::get().num_sources());
			voxel->diffuse_reflection /= (Config::get().num_photons() * Config::get().num_sources());
		}
	}
}

/***********************************************************
 * Compute the specular reflectance from a light source at
 * the surface.
 ***********************************************************/
void Simulator::specular_reflection(Photon& photon) {
	Voxel* voxel = nullptr;
	
	// Try the more robust DDA-based voxel lookup first
	Medium* dda_medium = find_medium_at_with_dda(photon.source.intersect);
	if (dda_medium) {
		glm::dvec3 pos = photon.source.intersect; // Make a non-const copy
		voxel = dda_medium->voxel_at(pos);
	}
	
	// Fallback to regular voxel_at if DDA didn't work
	if (!voxel) {
		voxel = voxel_at(photon.source.intersect);
	}
	
	// Final fallback: nudge the intersection point slightly into the medium
	if (!voxel) {
		glm::dvec3 nudged_pos = photon.source.intersect + photon.source.direction * 1e-6;
		Medium* nudged_medium = find_medium_at_with_dda(nudged_pos);
		if (nudged_medium) {
			voxel = nudged_medium->voxel_at(nudged_pos);
		}
	}

	// voxel should never be nullptr at this point, but handle gracefully if it occurs
	if (!voxel) {
		// Log the error but continue simulation - assume no reflection occurs
		if (Config::is_initialized() && Config::get().log()) {
			std::cerr << "Warning: specular reflection could not be computed - photon will continue unreflected." << std::endl;
		}
		// Graceful fallback: no reflection, photon continues with full weight
		return;
	}

	// refractive indices of ambient medium and medium that is hit
	double n1 = Config::get().ambient_eta();
	double n2 = voxel->material->eta();

	// Calculate specular reflection coefficient from Fresnel equations
	double temp_ratio = (n1 - n2) / (n1 + n2);
	double fresnel_reflection = (n2 != n1) ? temp_ratio * temp_ratio : 0;
	
	// Record surface refraction (energy entering the medium) 
	auto* medium = find_medium_at(photon.source.intersect);
	if (medium) {
		// Surface refraction is the energy that enters the medium (1 - reflected energy)
		double surface_refraction_energy = photon.weight * (1.0 - fresnel_reflection);
		medium->get_metrics().add_surface_refraction(surface_refraction_energy); // Energy entering medium at surface
		
		// Record specular reflection (energy immediately reflected at surface)
		double specular_reflection_energy = photon.weight * fresnel_reflection;
		medium->get_metrics().add_surface_reflection(specular_reflection_energy);
		
		// ADD SPECULAR REFLECTION TO ENTRY VOXEL EMITTANCE FOR RENDERING
		// This makes specular reflection visible in the renderer without breaking energy conservation
		// in the simulator metrics (which remain separate for physics accuracy)
		voxel->specular_reflection += specular_reflection_energy;
		
		// CRITICAL FIX: Reduce photon weight by reflected amount so only transmitted energy
		// continues for absorption/emission. This prevents double-counting reflected energy.
		photon.weight = surface_refraction_energy; // Only transmitted energy continues
	}

	// reflection direction: R = V - 2(V . N)N  
	// With outward-pointing normals, calculate reflection properly
	glm::dvec3 normal = photon.source.triangle.normal();
	glm::dvec3 incident = photon.source.direction;
	double projection_scalar = glm::dot(incident, normal);
	
	// For outward-pointing normal and incident ray:
	// If projection_scalar < 0: ray approaching from outside (normal case)
	// If projection_scalar > 0: ray approaching from inside (less common)
	// Standard reflection formula: R = I - 2(I·N)N works with outward normals
	// when the incident ray is pointing toward the surface
	
	glm::dvec3 reflection_direction = incident - 2.0 * projection_scalar * normal;
	photon.source.specular_direction = glm::normalize(reflection_direction);
}

/***********************************************************
 * Compute the fraction of incoming light that is reflected
 * back at an interface between two media. Also compute the
 * directions of transmission and reflection.
 ***********************************************************/
double Simulator::internal_reflection(Photon& photon, double& eta_t, glm::dvec3& transmittance,
									  glm::dvec3& reflectance) {
	// Modern Fresnel equations with improved numerical stability
	double eta_i = photon.voxel->material->eta();
	double eta_ratio = eta_t / eta_i;
	double eta_ratio_sq = eta_ratio * eta_ratio;

	// Get normalized vectors
	glm::dvec3 normal = glm::normalize(photon.voxel_normal);
	glm::dvec3 incident = glm::normalize(photon.direction);

	// For proper Fresnel calculations with outward-pointing normals:
	// If incident ray and normal point in same direction, photon is exiting (inside -> outside)
	// If incident ray and normal point in opposite directions, photon is entering (outside -> inside)
	double cos_i = glm::dot(incident, normal);
	
	// For exiting photon (cos_i > 0), use normal as-is
	// For entering photon (cos_i < 0), flip normal to point toward incident medium
	if (cos_i < 0.0) {
		normal = -normal;
		cos_i = -cos_i;
	}
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
	double d = Config::get().vox_size() * 0.00001;

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
		medium.reset_record_absorption_and_diffuse();
	}
	
	// Aggregate data from all voxels
	for (auto& medium : mediums) {
		auto& volume = medium.get_volume();  // Remove const since we need non-const access
		
		// Iterate through all voxels in this medium
		for (uint32_t z = 0; z < volume.depth(); ++z) {
			for (uint32_t y = 0; y < volume.height(); ++y) {
				for (uint32_t x = 0; x < volume.width(); ++x) {
					Voxel* voxel = volume.at(x, y, z);
					if (voxel && voxel->material) {
						// All voxels in this medium's volume belong to this medium
						// Aggregate absorption
						medium.get_metrics().add_total_absorption(voxel->absorption);
						
						// Aggregate emittance by direction classification
						medium.get_metrics().add_diffuse_reflection(voxel->specular_reflection);
						medium.get_metrics().add_diffuse_transmission(voxel->diffuse_transmission);
						
						// Note: surface_refraction and specular_reflection handled separately
						// Note: specular_transmission currently unused
					}
				}
			}
		}
	}
}

/***********************************************************
 *	Write the resulting physical quantities to a file.
 ***********************************************************/
/***********************************************************
 *	Write the resulting physical quantities to a file.
 *	Now uses Metrics class for unified statistics and export functionality.
 ***********************************************************/
void Simulator::report(bool generate_csv) {
	// Use the simulator's own metrics for log export to ensure correct timing
	metrics.export_results(*this, generate_csv);
}

/***********************************************************
 * Check if point is inside any medium.
 ***********************************************************/
bool Simulator::is_inside_any_medium(const glm::dvec3& position) const {
	return find_medium_at(position) != nullptr;
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
 * Handle photon transition between mediums.
 ***********************************************************/
void Simulator::handle_medium_transition(Photon& photon, Medium* from, Medium* to) {
	if (!from && to) {
		// Entering a medium from ambient space
		to->increment_photons_entered();
	} 
	else if (from && !to) {
		// Exiting to ambient space
		photon.alive = false;
	} 
	else if (from && to && from != to) {
		// COMPREHENSIVE MULTI-LAYER INTERFACE PHYSICS
		// Transitioning between different mediums with proper Fresnel calculations
		
		// Get material properties for both media
		Material* from_material = nullptr;
		Material* to_material = nullptr;
		
		// Find representative voxels to get material properties
		glm::dvec3 from_pos = photon.position - photon.direction * 1e-6; // Slightly behind
		glm::dvec3 to_pos = photon.position + photon.direction * 1e-6;   // Slightly ahead
		
		Voxel* from_voxel = from->voxel_at(from_pos);
		Voxel* to_voxel = to->voxel_at(to_pos);
		
		if (from_voxel && from_voxel->material) from_material = from_voxel->material;
		if (to_voxel && to_voxel->material) to_material = to_voxel->material;
		
		if (!from_material || !to_material) {
			std::cerr << "Warning: Interface transition without proper material properties" << std::endl;
			// Fallback to simple transmission
			return;
		}
		
		// Get refractive indices
		double n1 = from_material->eta();  // Incident medium
		double n2 = to_material->eta();    // Transmitted medium
		
		// Get interface normal (use voxel normal or calculate from geometry)
		glm::dvec3 interface_normal = photon.voxel_normal;
		
		// For multi-layer case, we know layers are horizontal (Y-axis boundaries)
		// So interface normal should be primarily in Y direction
		if (glm::length(interface_normal) < 0.1) {
			// Fallback: calculate normal from layer geometry
			interface_normal = glm::dvec3(0.0, 1.0, 0.0); // Upward normal for Y-boundaries
		}
		
		// Ensure normal points into the transmitted medium
		if (glm::dot(interface_normal, photon.direction) > 0) {
			interface_normal = -interface_normal;
		}
		
		// Calculate incident angle using Snell's law
		glm::dvec3 incident_dir = -photon.direction; // Direction toward interface
		double cos_theta_i = glm::dot(incident_dir, interface_normal);
		cos_theta_i = glm::clamp(cos_theta_i, -1.0, 1.0);
		double sin_theta_i = std::sqrt(1.0 - cos_theta_i * cos_theta_i);
		(void)sin_theta_i; // Suppress unused variable warning - kept for debugging
		
		// Check for total internal reflection
		double n_ratio = n1 / n2;
		double sin_theta_t_squared = n_ratio * n_ratio * (1.0 - cos_theta_i * cos_theta_i);
		
		if (sin_theta_t_squared > 1.0) {
			// TOTAL INTERNAL REFLECTION
			photon.direction = glm::reflect(photon.direction, interface_normal);
			photon.direction = glm::normalize(photon.direction);
			
			// For medium transitions, this is true boundary physics - no energy splitting
			// Photon continues in original medium with full weight
			
			if (Config::get().log()) {
				std::ostringstream debug_msg;
				debug_msg << "Total internal reflection at medium interface (n1=" << n1 << ", n2=" << n2 << ")";
				Logger::instance().log_debug(debug_msg.str());
			}
			return;
		}
		
		// Calculate transmitted angle
		double cos_theta_t = std::sqrt(1.0 - sin_theta_t_squared);
		
		// Calculate Fresnel reflection coefficient
		double Rs, Rp, R_fresnel;
		
		if (cos_theta_i < 1e-6) {
			// Normal incidence
			R_fresnel = std::pow((n1 - n2) / (n1 + n2), 2.0);
		} else {
			// General case - calculate s and p polarization components
			Rs = std::pow((n1 * cos_theta_i - n2 * cos_theta_t) / (n1 * cos_theta_i + n2 * cos_theta_t), 2.0);
			Rp = std::pow((n1 * cos_theta_t - n2 * cos_theta_i) / (n1 * cos_theta_t + n2 * cos_theta_i), 2.0);
			R_fresnel = 0.5 * (Rs + Rp); // Average for unpolarized light
		}
		
		// Ensure valid reflection coefficient
		R_fresnel = glm::clamp(R_fresnel, 0.0, 1.0);
		double T_fresnel = 1.0 - R_fresnel; // Transmission coefficient
		
		// Apply probabilistic Fresnel reflection/transmission
		// For medium-to-medium transitions, use standard Monte Carlo Fresnel physics
		// No energy splitting - photon either reflects OR transmits probabilistically
		
		if (rng->next() < R_fresnel) {
			// FRESNEL REFLECTION
			photon.direction = glm::reflect(photon.direction, interface_normal);
			photon.direction = glm::normalize(photon.direction);
			
			if (Config::get().log()) {
				std::ostringstream debug_msg;
				debug_msg << "Fresnel reflection at medium interface (R=" << R_fresnel << ")";
				Logger::instance().log_debug(debug_msg.str());
			}
			
		} else {
		
			// FRESNEL TRANSMISSION WITH REFRACTION
			
			// Calculate refracted direction using Snell's law
			glm::dvec3 transmitted_dir;
			
			if (cos_theta_i > 0.9999) {
				// Near-normal incidence - no direction change
				transmitted_dir = photon.direction;
			} else {
				// General refraction using vector form of Snell's law
				glm::dvec3 tangential = photon.direction - cos_theta_i * interface_normal;
				transmitted_dir = n_ratio * tangential + (n_ratio * cos_theta_i - cos_theta_t) * interface_normal;
			}
			
			photon.direction = glm::normalize(transmitted_dir);
			
			// For medium transitions, photon continues with full weight
			
			if (Config::get().log()) {
				std::ostringstream debug_msg;
				debug_msg << "Fresnel transmission at medium interface (T=" << T_fresnel 
						  << ", angle_i=" << std::acos(cos_theta_i) * 180.0 / std::numbers::pi
						  << "°, angle_t=" << std::acos(cos_theta_t) * 180.0 / std::numbers::pi << "°)";
				Logger::instance().log_debug(debug_msg.str());
			}
		}
		
		// VALIDATION: Ensure photon direction is physically reasonable
		if (glm::length(photon.direction) < 0.9 || glm::length(photon.direction) > 1.1) {
			std::cerr << "Warning: Invalid photon direction after interface transition. Normalizing." << std::endl;
			photon.direction = glm::normalize(photon.direction);
		}
		
		// ENERGY CONSERVATION CHECK
		if (photon.weight < 0.0 || photon.weight > 1.0) {
			std::cerr << "Warning: Invalid photon weight after interface transition: " << photon.weight << std::endl;
			photon.weight = glm::clamp(photon.weight, 0.0, 1.0);
		}
		
		// POST-INTERFACE VALIDATION: Ensure photon is in valid state
		validate_photon_state_after_interface_transition(photon, from, to);
	}
}

/***********************************************************
 * INTERFACE TRANSITION VALIDATION
 ***********************************************************/

/***********************************************************
 * Validate photon state after interface transition to prevent
 * invalid states that could cause ray-voxel intersection failures.
 ***********************************************************/
void Simulator::validate_photon_state_after_interface_transition(Photon& photon, Medium* from_medium, Medium* to_medium) {
	// VALIDATION 1: Ensure photon position is within valid bounds
	bool position_valid = true;
	
	// Check if photon is in either the from or to medium
	bool in_from_medium = from_medium && from_medium->contains_point(photon.position);
	bool in_to_medium = to_medium && to_medium->contains_point(photon.position);
	
	if (!in_from_medium && !in_to_medium) {
		std::cerr << "Error: Photon position invalid after interface transition" << std::endl;
		std::cerr << "  Position: (" << photon.position.x << ", " << photon.position.y << ", " << photon.position.z << ")" << std::endl;
		position_valid = false;
	}
	
	// VALIDATION 2: Ensure photon can find a valid voxel
	Medium* current_medium = find_medium_at(photon.position);
	if (!current_medium) {
		std::cerr << "Error: No medium found at photon position after interface transition" << std::endl;
		position_valid = false;
	}
	else {
		Voxel* current_voxel = current_medium->voxel_at(photon.position);
		if (!current_voxel || !current_voxel->material) {
			std::cerr << "Error: No valid voxel/material at photon position after interface transition" << std::endl;
			position_valid = false;
		}
	}
	
	// VALIDATION 3: Direction vector validation
	double dir_length = glm::length(photon.direction);
	if (dir_length < 0.9 || dir_length > 1.1) {
		std::cerr << "Warning: Invalid direction vector length after interface transition: " << dir_length << std::endl;
		photon.direction = glm::normalize(photon.direction);
	}
	
	// Check for NaN or infinite components
	if (std::isnan(photon.direction.x) || std::isnan(photon.direction.y) || std::isnan(photon.direction.z) ||
		std::isinf(photon.direction.x) || std::isinf(photon.direction.y) || std::isinf(photon.direction.z)) {
		std::cerr << "Error: Invalid direction components (NaN/Inf) after interface transition" << std::endl;
		// Emergency fallback direction
		photon.direction = glm::dvec3(0.0, -1.0, 0.0); // Downward
	}
	
	// RECOVERY ACTIONS for invalid position
	if (!position_valid) {
		std::cerr << "Attempting position recovery..." << std::endl;
		
		// RECOVERY 1: Try small perturbations around current position
		const double RECOVERY_EPSILON = 1e-6;
		std::vector<glm::dvec3> recovery_offsets = {
			glm::dvec3(0, RECOVERY_EPSILON, 0),   // Slightly up
			glm::dvec3(0, -RECOVERY_EPSILON, 0),  // Slightly down
			glm::dvec3(RECOVERY_EPSILON, 0, 0),   // Slightly right
			glm::dvec3(-RECOVERY_EPSILON, 0, 0),  // Slightly left
			glm::dvec3(0, 0, RECOVERY_EPSILON),   // Slightly forward
			glm::dvec3(0, 0, -RECOVERY_EPSILON)   // Slightly back
		};
		
		bool recovery_successful = false;
		for (const auto& offset : recovery_offsets) {
			glm::dvec3 test_pos = photon.position + offset;
			Medium* test_medium = find_medium_at(test_pos);
			if (test_medium) {
				Voxel* test_voxel = test_medium->voxel_at(test_pos);
				if (test_voxel && test_voxel->material) {
					photon.position = test_pos;
					photon.voxel = test_voxel;
					recovery_successful = true;
					std::cerr << "Position recovery successful with offset (" 
							  << offset.x << ", " << offset.y << ", " << offset.z << ")" << std::endl;
					break;
				}
			}
		}
		
		if (!recovery_successful) {
			// RECOVERY 2: Emergency termination with energy conservation
			std::cerr << "Position recovery failed. Terminating photon with energy conservation." << std::endl;
			
			// Record remaining energy as absorption in the last valid medium
			if (photon.weight > 0.0) {
				Medium* emergency_medium = from_medium ? from_medium : to_medium;
				if (emergency_medium) {
					emergency_medium->get_metrics().add_total_absorption(photon.weight);
				}
			}
			
			photon.alive = false;
			photon.weight = 0.0;
			return;
		}
	}
	
	// VALIDATION 4: Ensure photon has correct voxel reference
	if (current_medium) {
		Voxel* correct_voxel = current_medium->voxel_at(photon.position);
		if (correct_voxel && correct_voxel->material) {
			photon.voxel = correct_voxel;
		}
	}
	
	// VALIDATION 5: Final state check
	if (photon.alive && (!photon.voxel || !photon.voxel->material)) {
		std::cerr << "Final validation failed: photon has invalid voxel reference" << std::endl;
		photon.alive = false;
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
	double voxel_size = Config::get().vox_size();
	Range3 bounds = get_combined_bounds();
	glm::dvec3 position(bounds.min_bounds.x + (x + 0.5) * voxel_size,
	                    bounds.min_bounds.y + (y + 0.5) * voxel_size, 
	                    bounds.min_bounds.z + (z + 0.5) * voxel_size);
	
	// Find which medium contains this position and get the voxel
	return voxel_at(position);
}

/***********************************************************
 * Get PhotonPaths from photons for backward compatibility
 ***********************************************************/
/***********************************************************
 * Accessor methods for aggregating data across all mediums
 ***********************************************************/

std::vector<Material> Simulator::get_all_tissues() const {
	std::vector<Material> all_tissues;
	for (const auto& medium : mediums) {
		const auto& medium_tissues = medium.get_tissues();
		for (const auto& material : medium_tissues) {
			// Add unique materials from each medium based on optical properties
			bool found = false;
			for (const auto& existing : all_tissues) {
				if (existing.has_same_optical_properties(material)) {
					found = true;
					break;
				}
			}
			if (!found) {
				all_tissues.push_back(material);
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
 * Determine if photon is reflecting (exiting back toward light source)
 * or transmitting (exiting away from light source).
 * 
 * GENERAL MESH APPROACH: Use surface normal to determine if photon exits
 * on the same side as the light source (reflection) or opposite side (transmission).
 * This works for arbitrary mesh geometries, not just axis-aligned boxes.
 ***********************************************************/
bool Simulator::is_photon_reflecting(const Photon& photon) const {
	// Get the surface normal at exit point
	glm::dvec3 surface_normal = photon.voxel_normal;
	
	// Get incident light direction (from source to exit point)
	glm::dvec3 source_pos = photon.source.origin;
	glm::dvec3 exit_pos = photon.intersect;
	glm::dvec3 incident_direction = glm::normalize(exit_pos - source_pos);
	
	// Determine which side of the surface the light is incident from
	// If incident_direction dot surface_normal < 0, light hits from "outside" (normal side)
	// If incident_direction dot surface_normal > 0, light hits from "inside" (opposite side)
	double incident_dot_normal = glm::dot(incident_direction, surface_normal);
	
	// Get photon exit direction 
	glm::dvec3 exit_direction = photon.direction;
	
	// Determine which side of the surface the photon exits toward
	double exit_dot_normal = glm::dot(exit_direction, surface_normal);
	
	// Classification logic:
	// - If incident and exit are on same side of surface normal → REFLECTION
	// - If incident and exit are on opposite sides of surface normal → TRANSMISSION
	
	bool incident_from_normal_side = (incident_dot_normal < 0);
	bool exit_toward_normal_side = (exit_dot_normal > 0);
	
	// Photon is reflecting if it exits back toward the same side light came from
	return (incident_from_normal_side == exit_toward_normal_side);
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

// Combined metrics access methods (delegated to Metrics)
double Simulator::get_combined_total_absorption() const {
	if (!shared_metrics_) return 0.0;
	return shared_metrics_->aggregate_medium_energy_data(*this).total_absorption;
}

double Simulator::get_combined_diffuse_reflection() const {
	if (!shared_metrics_) return 0.0;
	return shared_metrics_->aggregate_medium_energy_data(*this).diffuse_reflection;
}

double Simulator::get_combined_specular_reflection() const {
	if (!shared_metrics_) return 0.0;
	return shared_metrics_->aggregate_medium_energy_data(*this).specular_reflection;
}

double Simulator::get_combined_surface_refraction() const {
	if (!shared_metrics_) return 0.0;
	return shared_metrics_->aggregate_medium_energy_data(*this).surface_refraction;
}

double Simulator::get_combined_diffuse_transmission() const {
	if (!shared_metrics_) return 0.0;
	return shared_metrics_->aggregate_medium_energy_data(*this).diffuse_transmission;
}

double Simulator::get_combined_specular_transmission() const {
	if (!shared_metrics_) return 0.0;
	return shared_metrics_->aggregate_medium_energy_data(*this).specular_transmission;
}

double Simulator::get_combined_avg_path_length() const {
	double total = 0.0;
	for (const auto& medium : mediums) {
		total += medium.get_metrics().get_path_length();
	}
	return total;
}

int Simulator::get_combined_total_steps() const {
	int total = 0;
	for (const auto& medium : mediums) {
		total += medium.get_metrics().get_total_steps();
	}
	return total;
}

int Simulator::get_combined_photons_entered() const {
	int total = 0;
	for (const auto& medium : mediums) {
		total += medium.get_metrics().get_photons_entered();
	}
	return total;
}

/***********************************************************
 * ENERGY STATISTICS DELEGATION METHODS
 * These methods delegate to the shared Metrics instance
 ***********************************************************/
Metrics::MediumEnergyData Simulator::aggregate_medium_energy_data() const {
	if (!shared_metrics_) {
		// Return empty data if metrics not initialized
		return Metrics::MediumEnergyData{};
	}
	return shared_metrics_->aggregate_medium_energy_data(*this);
}

Metrics::EnergyConservation Simulator::calculate_energy_conservation() const {
	if (!shared_metrics_) {
		// Return empty data if metrics not initialized
		return Metrics::EnergyConservation{};
	}
	return shared_metrics_->calculate_energy_conservation(*this);
}

Metrics::EnergyConservationPercentages Simulator::calculate_energy_percentages() const {
	if (!shared_metrics_) {
		// Return empty data if metrics not initialized
		return Metrics::EnergyConservationPercentages{};
	}
	return shared_metrics_->calculate_energy_percentages(*this);
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
	
	// CONSOLIDATED APPROACH: This function only handles INTERNAL terminations (absorption)
	// For exits, use radiate() function which handles both voxel emittance and medium records
	if (reason == "absorption" || reason == "roulette" || reason == "max_iterations" || 
		reason == "no_voxel_found" || reason == "no_material_properties" || reason == "no_medium_found") {
		// Energy absorbed within the medium
		record_medium->get_metrics().add_total_absorption(photon.weight);
		
		// ENERGY CONSERVATION FIX: Track energy in photon accounting system
		if (photon.voxel && photon.voxel->material) {
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
		record_medium->get_metrics().add_total_absorption(photon.weight);
		
		// ENERGY CONSERVATION FIX: Track energy in photon accounting system
		if (photon.voxel && photon.voxel->material) {
			// Update photon energy tracking for conservation
			photon.total_energy_absorbed += photon.weight;
			
			// Add to voxel absorption as fallback
			photon.voxel->absorption += photon.weight;
		}
	}
	
	// Terminate the photon
	photon.alive = false;
	
	// Update detailed photon tracking data
	if (photon.id < detailed_photon_data.size()) {
		DetailedPhotonData& detailed = detailed_photon_data[photon.id];
		detailed.termination_position = photon.position;
		detailed.termination_direction = photon.direction;
		detailed.scatter_count = photon.scatter_count;
		detailed.total_absorption_deposited = photon.total_energy_absorbed;
		detailed.remaining_weight = photon.weight;
		detailed.termination_reason = reason;
	}
	
	photon.weight = 0.0; // Clear weight to prevent double-counting
}

/***********************************************************
 * Initialize 3D DDA instances for robust voxel traversal
 ***********************************************************/
void Simulator::initialize_dda_instances() {
	medium_ddas_.clear();
	medium_ddas_.reserve(mediums.size());
	
	for (size_t i = 0; i < mediums.size(); ++i) {
		const auto& medium = mediums[i];
		
		// Get grid parameters from medium
		glm::ivec3 grid_dimensions(Config::get().nx(), Config::get().ny(), Config::get().nz());
		glm::dvec3 grid_origin(medium.get_bounds().min_bounds);
		double voxel_size = Config::get().vox_size();
		
		// Create DDA instance for this medium
		auto dda = std::make_unique<VoxelDDA3D>(grid_dimensions, grid_origin, voxel_size);
		medium_ddas_.push_back(std::move(dda));
		
		if (Config::get().log()) {
			std::ostringstream debug_msg;
			debug_msg << "  DDA " << i << ": Grid " << grid_dimensions.x << "x" << grid_dimensions.y << "x" << grid_dimensions.z 
					  << ", Origin " << grid_origin.x << "," << grid_origin.y << "," << grid_origin.z 
					  << ", Voxel size " << voxel_size;
			Logger::instance().log_debug(debug_msg.str());
		}
	}
}

/***********************************************************
 * Track absorption along actual photon path segments for maximum accuracy
 ***********************************************************/
void Simulator::track_photon_path_segments_for_absorption(Photon& photon) {
	if (!photon.voxel || !photon.voxel->material) {
		return;
	}

	// Use the same path segment that will be rendered
	glm::dvec3 start_pos = photon.position;
	glm::dvec3 end_pos = photon.intersect;
	double total_distance = glm::length(end_pos - start_pos);
	
	if (total_distance < 1e-12) {
		return;
	}
	
	glm::dvec3 direction = glm::normalize(end_pos - start_pos);
	
	// Find current medium
	Medium* current_medium = find_medium_at(start_pos);
	if (!current_medium) {
		return;
	}
	
	// Calculate absorption using the same DDA traversal but on the actual photon path
	// Find the DDA instance for this medium
	size_t medium_index = 0;
	for (size_t i = 0; i < mediums.size(); ++i) {
		if (&mediums[i] == current_medium) {
			medium_index = i;
			break;
		}
	}
	
	if (medium_index >= medium_ddas_.size()) {
		// DDA not available - this should not happen in production
		std::cerr << "WARNING: DDA not available for medium " << medium_index << std::endl;
		return;
	}

	VoxelDDA3D* dda = medium_ddas_[medium_index].get();
	
	// Initialize DDA for this ray
	dda->initialize_ray(start_pos, direction);
	
	// Traverse voxels using DDA on the actual photon path
	VoxelDDA3D::TraversalResult result = dda->traverse(total_distance);

	// Calculate absorption along the path using DDA results
	double remaining_weight = photon.weight;
	double total_absorption = 0.0;
	Voxel* last_surface_voxel = photon.voxel; // Default to current voxel

	for (const auto& step : result.voxels) {
		// Get voxel at this DDA position - but ensure coordinate consistency
		glm::dvec3 mutable_pos = step.world_position; // Create mutable copy for voxel_at
		Voxel* voxel = current_medium->voxel_at(mutable_pos);
		
		// COORDINATE FIX: If voxel_at returns wrong coordinates, use DDA coordinates directly
		if (voxel) {
			glm::ivec3 dda_coords = step.voxel_coords;
			glm::ivec3 medium_coords = glm::ivec3(voxel->ix(), voxel->iy(), voxel->iz());
			if (dda_coords != medium_coords) {
				// Use DDA coordinates directly when there's a mismatch
				if (dda_coords.x >= 0 && dda_coords.x < static_cast<int>(Config::get().nx()) &&
					dda_coords.y >= 0 && dda_coords.y < static_cast<int>(Config::get().ny()) &&
					dda_coords.z >= 0 && dda_coords.z < static_cast<int>(Config::get().nz())) {
					voxel = voxel_grid(static_cast<uint32_t>(dda_coords.x), 
									   static_cast<uint32_t>(dda_coords.y), 
									   static_cast<uint32_t>(dda_coords.z));
				}
			}
		}
		
		if (voxel && voxel->material) {
			// Calculate distance for this voxel segment
			double segment_distance = 0.0;
			if (!result.voxels.empty()) {
				// Calculate distance between consecutive steps
				auto it = std::find_if(result.voxels.begin(), result.voxels.end(),
					[&step](const VoxelDDA3D::StepResult& s) { return s.voxel_coords == step.voxel_coords; });
				
				if (it != result.voxels.end()) {
					size_t index = std::distance(result.voxels.begin(), it);
					if (index < result.voxels.size() - 1) {
						segment_distance = result.voxels[index + 1].distance_traveled - step.distance_traveled;
					} else {
						segment_distance = total_distance - step.distance_traveled;
					}
				}
			}
			
			if (segment_distance > 1e-12) {
				// Calculate effective volume fraction for boundary voxels
				double effective_volume_fraction = 1.0;
				if (voxel->is_boundary_voxel) {
					effective_volume_fraction = voxel->volume_fraction_inside;
				}
				
				// Apply Beer-Lambert law for this segment
				double mu_a = voxel->material->mu_a() * effective_volume_fraction;
				double segment_transmission = std::exp(-mu_a * segment_distance);
				double segment_absorption = remaining_weight * (1.0 - segment_transmission);
				
				// Deposit absorption in this voxel
				voxel->absorption += segment_absorption;
				total_absorption += segment_absorption;
				
				// Update remaining weight for next segment
				remaining_weight *= segment_transmission;
				
				// Track surface voxels for exit recording (ONLY external surfaces, not internal boundaries)
				if (voxel->is_surface_voxel) {
					last_surface_voxel = voxel;
				}
			}
		}
	}
	
	// Update photon weight and energy tracking
	photon.weight = remaining_weight;
	
	// ENERGY CONSERVATION FIX: Update photon energy tracking
	if (total_absorption > 0.0) {
		photon.total_energy_absorbed += total_absorption;
	}
	
	// Store the last surface voxel for emittance recording
	if (last_surface_voxel) {
		photon.last_surface_voxel = last_surface_voxel;
	}
	
	// Update photon's voxel reference to last surface voxel for proper exit recording
	if (last_surface_voxel) {
		photon.voxel = last_surface_voxel;
	}
}

/***********************************************************
 * Robust voxel traversal using 3D DDA algorithm
 ***********************************************************/
void Simulator::track_voxel_path_with_dda(Photon& photon) {
	if (!photon.voxel || !photon.voxel->material) {
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
	
	// Find current medium and its DDA instance
	Medium* current_medium = find_medium_at(start_pos);
	if (!current_medium) {
		return;
	}
	
	// Find the DDA instance for this medium
	size_t medium_index = 0;
	for (size_t i = 0; i < mediums.size(); ++i) {
		if (&mediums[i] == current_medium) {
			medium_index = i;
			break;
		}
	}
	
	if (medium_index >= medium_ddas_.size()) {
		// DDA not available - this should not happen in production
		std::cerr << "WARNING: DDA not available for medium " << medium_index << std::endl;
		return;
	}
	
	VoxelDDA3D* dda = medium_ddas_[medium_index].get();
	
	// Initialize DDA for this ray
	dda->initialize_ray(start_pos, direction);
	
	// Traverse voxels using DDA
	VoxelDDA3D::TraversalResult result = dda->traverse(total_distance);
	
	// Calculate absorption along the path using DDA results
	double remaining_weight = photon.weight;
	double total_absorption = 0.0;
	Voxel* last_surface_voxel = photon.voxel; // Default to current voxel

	for (const auto& step : result.voxels) {
		// Get voxel at this DDA position
		glm::dvec3 mutable_pos = step.world_position; // Create mutable copy for voxel_at
		Voxel* voxel = current_medium->voxel_at(mutable_pos);
		
		if (voxel && voxel->material) {
			// Calculate distance for this voxel segment
			double segment_distance = 0.0;
			if (!result.voxels.empty()) {
				// Calculate distance between consecutive steps
				auto it = std::find_if(result.voxels.begin(), result.voxels.end(),
					[&step](const VoxelDDA3D::StepResult& s) { return s.voxel_coords == step.voxel_coords; });
				
				if (it != result.voxels.end()) {
					size_t index = std::distance(result.voxels.begin(), it);
					if (index < result.voxels.size() - 1) {
						segment_distance = result.voxels[index + 1].distance_traveled - step.distance_traveled;
					} else {
						segment_distance = total_distance - step.distance_traveled;
					}
				}
			}
			
			if (segment_distance > 1e-12) {
				// Calculate effective volume fraction for boundary voxels
				double effective_volume_fraction = 1.0;
				if (voxel->is_boundary_voxel) {
					effective_volume_fraction = voxel->volume_fraction_inside;
				}
				
				// Apply Beer-Lambert law for this segment
				double mu_a = voxel->material->mu_a() * effective_volume_fraction;
				double segment_transmission = std::exp(-mu_a * segment_distance);
				double segment_absorption = remaining_weight * (1.0 - segment_transmission);
				
				// Deposit absorption in this voxel
				voxel->absorption += segment_absorption;
				total_absorption += segment_absorption;
				
				// Update remaining weight for next segment
				remaining_weight *= segment_transmission;
				
				// Track surface voxels for exit recording (ONLY external surfaces, not internal boundaries)
				if (voxel->is_surface_voxel) {
					last_surface_voxel = voxel;
				}
			}
		}
	}	// Update photon weight and energy tracking
	photon.weight = remaining_weight;
	
	// ENERGY CONSERVATION FIX: Update photon energy tracking
	if (total_absorption > 0.0) {
		photon.total_energy_absorbed += total_absorption;
	}
	
	// Store the last surface voxel for emittance recording
	if (last_surface_voxel) {
		photon.last_surface_voxel = last_surface_voxel;
	}
	
	// Update photon's voxel reference to last surface voxel for proper exit recording
	if (last_surface_voxel) {
		photon.voxel = last_surface_voxel;
	}
}

/***********************************************************
 * Robust medium detection using 3D DDA
 ***********************************************************/
Medium* Simulator::find_medium_at_with_dda(const glm::dvec3& position) const {
	// NUMERICAL STABILITY: Add epsilon tolerance for boundary detection
	static const double EPSILON = 1e-9;
	
	// Try each medium's DDA for precise boundary detection
	for (size_t i = 0; i < mediums.size() && i < medium_ddas_.size(); ++i) {
		const auto& medium = mediums[i];
		VoxelDDA3D* dda = medium_ddas_[i].get();
		Medium* mutable_medium = const_cast<Medium*>(&medium);
		
		// STAGE 1: Standard DDA validation
		glm::ivec3 voxel_coords = dda->world_to_voxel(position);
		if (dda->is_valid_voxel(voxel_coords)) {
			glm::dvec3 mutable_pos = position;
			Voxel* voxel = mutable_medium->voxel_at(mutable_pos);
			if (voxel && voxel->material) {
				return mutable_medium;
			}
		}
		
		// STAGE 2: Epsilon-nudged DDA validation (multiple attempts)
		glm::dvec3 robust_position = position;
		bool dda_success = false;
		
		for (double eps_scale = 1.0; eps_scale <= 1000.0 && !dda_success; eps_scale *= 10.0) {
			// Try nudging in all 8 directions (corners of epsilon cube)
			std::vector<glm::dvec3> nudge_directions = {
				glm::dvec3(-EPSILON * eps_scale, -EPSILON * eps_scale, -EPSILON * eps_scale),
				glm::dvec3(-EPSILON * eps_scale, -EPSILON * eps_scale,  EPSILON * eps_scale),
				glm::dvec3(-EPSILON * eps_scale,  EPSILON * eps_scale, -EPSILON * eps_scale),
				glm::dvec3(-EPSILON * eps_scale,  EPSILON * eps_scale,  EPSILON * eps_scale),
				glm::dvec3( EPSILON * eps_scale, -EPSILON * eps_scale, -EPSILON * eps_scale),
				glm::dvec3( EPSILON * eps_scale, -EPSILON * eps_scale,  EPSILON * eps_scale),
				glm::dvec3( EPSILON * eps_scale,  EPSILON * eps_scale, -EPSILON * eps_scale),
				glm::dvec3( EPSILON * eps_scale,  EPSILON * eps_scale,  EPSILON * eps_scale)
			};
			
			for (const auto& nudge : nudge_directions) {
				glm::dvec3 nudged_pos = position + nudge;
				glm::ivec3 nudged_coords = dda->world_to_voxel(nudged_pos);
				
				if (dda->is_valid_voxel(nudged_coords)) {
					glm::dvec3 test_pos = nudged_pos;
					Voxel* voxel = mutable_medium->voxel_at(test_pos);
					if (voxel && voxel->material) {
						return mutable_medium;
					}
				}
			}
		}
		
		// STAGE 3: BYPASS DDA - Direct voxel_at() check (most robust)
		glm::dvec3 bypass_pos = position;
		Voxel* direct_voxel = mutable_medium->voxel_at(bypass_pos);
		if (direct_voxel && direct_voxel->material) {
			// Direct voxel check succeeded - DDA validation was the problem
			return mutable_medium;
		}
		
		// STAGE 4: Geometric containment check - bypass all coordinate calculations
		if (medium.contains_point(position)) {
			// Position is geometrically inside this medium
			// Try to find ANY valid voxel near this position
			for (double search_radius = EPSILON; search_radius <= 0.01; search_radius *= 10.0) {
				for (int dx = -1; dx <= 1; dx++) {
					for (int dy = -1; dy <= 1; dy++) {
						for (int dz = -1; dz <= 1; dz++) {
							glm::dvec3 search_pos = position + glm::dvec3(dx * search_radius, dy * search_radius, dz * search_radius);
							glm::dvec3 test_pos = search_pos;
							Voxel* search_voxel = mutable_medium->voxel_at(test_pos);
							if (search_voxel && search_voxel->material) {
								return mutable_medium;
							}
						}
					}
				}
			}
		}
	}
	
	return nullptr; // All fallback stages failed - truly in ambient space
}