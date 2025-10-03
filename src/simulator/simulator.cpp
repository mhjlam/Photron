/**
 * @file simulator.cpp
 * @brief Implementation of Monte Carlo photon transport simulation engine
 *
 * Contains the complete implementation of the Simulator class, including photon
 * transport, medium interactions, scattering calculations, and energy conservation
 * tracking for multi-layered tissue modeling.
 */

#include "simulator.hpp"

#include <algorithm>
#include <ctime>
#include <fstream>
#include <iostream>
#include <iterator>
#include <limits>
#include <memory>
#include <numbers>
#include <set>
#include <sstream>

#include "app.hpp"
#include "common/config.hpp"
#include "common/error_handler.hpp"
#include "common/error_types.hpp"
#include "common/logger.hpp"
#include "common/result.hpp"
#include "math/cuboid.hpp"
#include "math/dda.hpp"
#include "math/math.hpp"
#include "math/random.hpp"
#include "math/range.hpp"
#include "math/ray.hpp"
#include "math/triangle.hpp"
#include "simulator/aggregator.hpp"
#include "simulator/geom_lookup.hpp"
#include "simulator/layer.hpp"
#include "simulator/material.hpp"
#include "simulator/medium.hpp"
#include "simulator/photon.hpp"
#include "simulator/volume.hpp"
#include "simulator/voxel.hpp"

/**
 * @brief Construct new Simulator with default Monte Carlo parameters
 *
 * Initializes core data structures and random number generator with
 * sensible defaults for Monte Carlo photon transport simulation.
 */
Simulator::Simulator() : rng(std::make_shared<Random>()), mcml_weight_threshold(MathConstants::MCML_WEIGHT_THRESHOLD) {
	// Initialize core simulation containers
	photons.clear();
	sources.clear();
	mediums.clear();

	// Reserve space for typical simulation sizes to avoid reallocations
	photons.reserve(10000); // Typical photon path storage
	sources.reserve(5);     // Few light sources expected

	// Seed random number generator with current time
	// Will be re-seeded during initialize() based on configuration
	rng->seed(static_cast<int>(std::time(nullptr)));
}

// =============================================================================
// INITIALIZATION
// =============================================================================

/***********************************************************
 * Parse the input file and initializes the data structures.
 ***********************************************************/
Result<void, SimulationError> Simulator::initialize(std::string file) {
	// Reset simulation state for clean initialization
	mediums.clear();
	photons.clear();
	sources.clear();
	emitters.clear();

	// Reset accumulated metrics from previous runs
	metrics.reset();

	// Reinitialize configuration system for rerun support
	Config::shutdown();
	if (!Config::initialize(file)) {
		return Result<void, SimulationError>::error(SimulationError::InvalidConfiguration);
	}

	// Parse and validate configuration file
	auto parse_result = parse(file);
	if (!parse_result.is_ok()) {
		return parse_result;
	}

	// Display initialization progress
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
	}
	else {
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
			ErrorHandler::instance().report_error("An error occurred while initializing the medium.");
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

	// Initialize component instances for modular functionality
	initialize_components();
	if (Config::get().log()) {
		// Log component initialization to file only
		FAST_LOG_INFO("Simulator components (GeomLookup, Aggregator) initialized successfully.");
	}

	// Initialize logger if logging is enabled
	if (Config::get().log()) {
		Logger::instance().initialize(App::get_output_path("trace.csv"), App::get_output_path("debug.log"), true);
		Logger::instance().log_info("=== Photron Simulation Started ===");
		Logger::instance().log_info("Initializing Photron");
		Logger::instance().log_info("Configuration parsed successfully.");
		std::cout << "Debug logger initialized for photon tracing." << std::endl;
		Logger::instance().log_info("Debug logger initialized for photon tracing.");
	}
	else {
		Logger::instance().initialize("", "", false); // Disable logging
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
		ErrorHandler::instance().report_warning(
			"shared_metrics_ not set. Please call set_shared_metrics() before initialize().");
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

// =============================================================================
// INITIALIZATION SUBROUTINES
// =============================================================================

/***********************************************************
 * Read and parse the given configuration file.
 ***********************************************************/
Result<void, SimulationError> Simulator::parse(const std::string& /* fconfig */) {
	// Clear existing data structures for reinitialization
	sources.clear();
	photons.clear();
	mediums.clear();

	// Verify configuration is properly initialized
	if (!Config::is_initialized()) {
		return Result<void, SimulationError>::error(SimulationError::ConfigNotInitialized);
	}

	sources = Config::get().sources();

	// Only show info if logging is enabled
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
				}
				else {
					// Try slightly adjusting the intersection point inward along the direction
					glm::dvec3 adjusted_intersect =
						source.intersect + source.direction * MathConstants::PHOTON_NUDGE_EPSILON;
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

/***********************************************************
 * Run the Monte Carlo photon transport simulation.
 ***********************************************************/
Result<void, SimulationError> Simulator::simulate() {
	// Validate simulation preconditions
	if (sources.empty()) {
		return Result<void, SimulationError>::error(SimulationError::InvalidConfiguration);
	}

	if (photons.empty()) {
		return Result<void, SimulationError>::error(SimulationError::NoPhotons);
	}

	// Initialize simulation timing and energy tracking
	metrics.start_clock();

	// Track overall energy balance for conservation validation
	static double total_initial_energy = 0.0;
	static double total_launched_energy = 0.0;

	// Main Monte Carlo simulation loop structure
	// For each light source (typically one for most simulations)
	for (auto& source : sources) {
		// Process each photon through complete transport history
		for (uint32_t p = 0; p < photons.size(); ++p) {
			// Initialize photon energy accounting
			total_initial_energy += 1.0; // Each photon starts with weight 1.0

			// Launch photon from light source
			photon_transport_->launch(photons[p], source);

			// Set source properties from first photon's optical calculations
			if (p == 0 && photons[p].alive) {
				source.specular_direction = photons[p].specular_direction();
			}

			// Track successfully launched energy
			if (photons[p].alive) {
				total_launched_energy += photons[p].weight;
			}

			// Monte Carlo transport loop with safety limits
			int photon_iteration_counter = 0;
			const int max_photon_iterations = 1000000;

			while (photons[p].alive) {
				// Safety mechanism to prevent infinite loops
				photon_iteration_counter++;
				if (photon_iteration_counter > max_photon_iterations) {
					std::string warning_msg =
						"Photon " + std::to_string(photons[p].id) + " exceeded maximum iterations, terminating.";
					ErrorHandler::instance().report_warning(warning_msg);
					FAST_LOG_WARNING(warning_msg);

					// Use centralized energy-conserving termination
					if (photons[p].weight > 0.0) {
						photon_transport_->terminate_photon_and_record_energy(photons[p], "max_iterations");
					}
					else {
						photons[p].alive = false;
					}
					break;
				}

				// Core Monte Carlo transport steps
				photon_transport_->step_size(photons[p]); // Sample step size from extinction coefficient
				photon_transport_->transfer(photons[p]);  // Propagate through medium with absorption
				photon_transport_->roulette(photons[p]);  // Russian roulette for photon termination

				// Skip scattering after boundary crossing (handled in cross())
				if (!photons[p].cross) {
					photon_transport_->scatter(photons[p]); // Sample new direction from phase function
				}
			}

			// Progress reporting
			uint32_t progress_interval = std::max(1u, static_cast<uint32_t>(photons.size() / 50)); // Every 2%

			if ((p + 1) % progress_interval == 0 || p == 0 || (p + 1) == photons.size()) {
				double progress_percent = ((double)(p + 1) / photons.size()) * 100.0;
				std::cout << "\rProgress: " << (p + 1) << "/" << Config::get().num_photons() << " (" << std::fixed
						  << std::setprecision(1) << progress_percent << "%)" << std::flush;

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
	if (aggregator_) {
		aggregator_->normalize();
	}

	// Voxel data output is handled by report() function

	metrics.stop_clock();

	// Set simulation completion data for accurate reporting
	metrics.set_simulation_completion(photons.size());

	// Increment simulation version since data has changed
	increment_simulation_version();

	// Aggregate voxel data to medium records before calculating final metrics
	if (aggregator_) {
		aggregator_->aggregate_voxel_data();
	}

	// Use consolidated energy aggregation method
	auto energy_data = metrics.aggregate_medium_energy_data(*this);

	// Set only the scatter events in the main simulator metrics (path data is now handled by medium metrics)
	metrics.set_scatter_events(energy_data.scatter_events);

	metrics.collect_data(energy_data.total_absorption,
						 energy_data.specular_reflection,
						 energy_data.diffuse_reflection,
						 energy_data.surface_refraction,
						 energy_data.specular_transmission,
						 energy_data.diffuse_transmission);

	// Set shared metrics for GUI components
	if (shared_metrics_) {
		shared_metrics_->collect_data(energy_data.total_absorption,
									  energy_data.specular_reflection,
									  energy_data.diffuse_reflection,
									  energy_data.surface_refraction,
									  energy_data.specular_transmission,
									  energy_data.diffuse_transmission);
		shared_metrics_->set_scatter_events(energy_data.scatter_events);
		shared_metrics_->stop_clock();                              // Set timing for GUI
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
		ErrorHandler::instance().report_error("No light sources available for single photon simulation");
		return;
	}

	Source& source = sources[0];

	// Create a new photon with unique ID
	uint64_t new_photon_id = photons.size();
	Photon new_photon(new_photon_id);

	// Launch the photon
	photon_transport_->launch(new_photon, source);

	// Safety mechanism to prevent infinite photon loops
	int photon_iteration_counter = 0;
	const int max_photon_iterations = 1000000;

	while (new_photon.alive) {
		photon_iteration_counter++;
		if (photon_iteration_counter > max_photon_iterations) {
			ErrorHandler::instance().report_warning("Single photon exceeded maximum iterations, terminating.");
			FAST_LOG_WARNING("Single photon exceeded maximum iterations, terminating.");

			// Deposit remaining energy as absorption for energy conservation
			if (new_photon.weight > 0.0 && new_photon.voxel && new_photon.voxel->material) {
				// Use energy conservation enforcement
				photon_transport_->terminate_photon_and_record_energy(new_photon, "max_iterations");
			}
			else {
				new_photon.alive = false;
			}
			break;
		}

		photon_transport_->step_size(new_photon); // Set new step size
		photon_transport_->transfer(new_photon);  // Propagate photon through the medium in substeps
		photon_transport_->roulette(new_photon);  // Determine photon termination
		photon_transport_->scatter(new_photon);   // Scatter photon into a new direction
	}

	// Add the completed photon to the photons vector for rendering
	photons.push_back(new_photon);

	// Aggregate voxel data to medium records after adding new photon
	// This ensures diffuse_reflection and diffuse_transmission are current for the overlay
	if (aggregator_) {
		aggregator_->aggregate_voxel_data();
	}

	// Use consolidated energy aggregation method
	auto energy_data = metrics.aggregate_medium_energy_data(*this);

	// Set only the scatter events in the main simulator metrics (path data is now handled by medium metrics)
	metrics.set_scatter_events(energy_data.scatter_events);

	metrics.collect_data(energy_data.total_absorption,
						 energy_data.specular_reflection,
						 energy_data.diffuse_reflection,
						 energy_data.surface_refraction,
						 energy_data.specular_transmission,
						 energy_data.diffuse_transmission);

	// Set shared metrics if available (for GUI components)
	if (shared_metrics_) {
		shared_metrics_->collect_data(energy_data.total_absorption,
									  energy_data.specular_reflection,
									  energy_data.diffuse_reflection,
									  energy_data.surface_refraction,
									  energy_data.specular_transmission,
									  energy_data.diffuse_transmission);
		shared_metrics_->set_scatter_events(energy_data.scatter_events);
	}

	// Increment simulation version since data has changed
	increment_simulation_version();
}

/***********************************************************
 * Compute the fraction of incoming light that is reflected
 * back at an interface between two media. Also compute the
 * directions of transmission and reflection.
 ***********************************************************/
double Simulator::internal_reflection(Photon& photon, double& eta_t, glm::dvec3& transm, glm::dvec3& refl) {
	// Modern Fresnel equations with improved numerical stability
	double eta_i = photon.voxel->material->eta();
	double eta_ratio = eta_t / eta_i;
	double eta_ratio_sq = eta_ratio * eta_ratio;

	// Get normalized vectors
	glm::dvec3 normal = glm::normalize(photon.voxel_normal);
	glm::dvec3 incident = glm::normalize(photon.direction);

	// Fresnel calculations: same direction = exiting, opposite = entering
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
		refl = incident - 2.0 * glm::dot(incident, normal) * normal;
		refl = glm::normalize(refl);
		transm = glm::dvec3(0.0); // No transmission
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
	transm = eta_ratio * incident + (eta_ratio * cos_i - cos_t) * normal;
	transm = glm::normalize(transm);

	// Reflection direction
	refl = incident - 2.0 * glm::dot(incident, normal) * normal;
	refl = glm::normalize(refl);

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
 * Write the resulting physical quantities to a file.
 ***********************************************************/
void Simulator::report(bool generate_csv) {
	// Process simulation results for reporting
	if (aggregator_) {
		aggregator_->aggregate_voxel_data();
		aggregator_->normalize();
	}

	// Always use the simulator's own metrics for export since it needs full Simulator reference
	metrics.export_results(*this, generate_csv);
}

/***********************************************************
 * Check if point is inside any medium.
 ***********************************************************/
bool Simulator::is_inside_any_medium(const glm::dvec3& position) const {
	if (geom_lookup_) {
		return geom_lookup_->is_inside_any_medium(position);
	}
	// Fallback if geom_lookup not initialized
	return find_medium_at(position) != nullptr;
}

// =============================================================================
// MEDIUM CONTEXT MANAGEMENT
// =============================================================================

/***********************************************************
 * Find which medium contains the given point.
 ***********************************************************/
Medium* Simulator::find_medium_at(const glm::dvec3& position) const {
	if (geom_lookup_) {
		return geom_lookup_->find_medium_at(position);
	}
	// Fallback if geom_lookup not initialized
	for (auto& medium : mediums) {
		if (medium.contains_point(position)) {
			return const_cast<Medium*>(&medium);
		}
	}
	return nullptr; // Point is in ambient space
}

// =============================================================================
// INTERFACE TRANSITION VALIDATION
// =============================================================================

/***********************************************************
 * Validate photon state after interface transition to prevent
 * invalid states that could cause ray-voxel intersection failures.
 ***********************************************************/
void Simulator::validate_photon_state_after_interface_transition(Photon& photon,
																 Medium* from_medium,
																 Medium* to_medium) {
	// VALIDATION 1: Ensure photon position is within valid bounds
	bool position_valid = true;

	// Check if photon is in either the from or to medium
	bool in_from_medium = from_medium && from_medium->contains_point(photon.position);
	bool in_to_medium = to_medium && to_medium->contains_point(photon.position);

	if (!in_from_medium && !in_to_medium) {
		std::ostringstream pos_error;
		pos_error << "Photon position invalid after interface transition at (" << photon.position.x << ", "
				  << photon.position.y << ", " << photon.position.z << ")";
		ErrorHandler::instance().report_error(
			ErrorMessage::format(SimulationError::InvalidPhotonState, pos_error.str()));
		FAST_LOG_ERROR(pos_error.str());
		position_valid = false;
	}

	// VALIDATION 2: Ensure photon can find a valid voxel
	Medium* current_medium = find_medium_at(photon.position);
	if (!current_medium) {
		ErrorHandler::instance().report_error(ErrorMessage::format(
			SimulationError::NoMediumFound, "No medium found at photon position after interface transition"));
		FAST_LOG_ERROR("No medium found at photon position after interface transition");
		position_valid = false;
	}
	else {
		Voxel* current_voxel = current_medium->voxel_at(photon.position);
		if (!current_voxel || !current_voxel->material) {
			ErrorHandler::instance().report_error(
				ErrorMessage::format(SimulationError::NoVoxelFound,
									 "No valid voxel/material at photon position after interface transition"));
			FAST_LOG_ERROR("No valid voxel/material at photon position after interface transition");
			position_valid = false;
		}
	}

	// VALIDATION 3: Direction vector validation
	double dir_length = glm::length(photon.direction);
	if (dir_length < 0.9 || dir_length > 1.1) {
		std::string dir_msg =
			"Invalid direction vector length after interface transition: " + std::to_string(dir_length);
		ErrorHandler::instance().report_warning(dir_msg);
		FAST_LOG_WARNING(dir_msg);
		photon.direction = glm::normalize(photon.direction);
	}

	// Check for NaN or infinite components
	if (std::isnan(photon.direction.x) || std::isnan(photon.direction.y) || std::isnan(photon.direction.z)
		|| std::isinf(photon.direction.x) || std::isinf(photon.direction.y) || std::isinf(photon.direction.z)) {
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
			glm::dvec3(0, RECOVERY_EPSILON, 0),  // Slightly up
			glm::dvec3(0, -RECOVERY_EPSILON, 0), // Slightly down
			glm::dvec3(RECOVERY_EPSILON, 0, 0),  // Slightly right
			glm::dvec3(-RECOVERY_EPSILON, 0, 0), // Slightly left
			glm::dvec3(0, 0, RECOVERY_EPSILON),  // Slightly forward
			glm::dvec3(0, 0, -RECOVERY_EPSILON)  // Slightly back
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
					std::cerr << "Position recovery successful with offset (" << offset.x << ", " << offset.y << ", "
							  << offset.z << ")" << std::endl;
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

// =============================================================================
// HELPER METHODS FOR MEDIUM DELEGATION
// =============================================================================

/***********************************************************
 * Find voxel at position by delegating to appropriate medium.
 ***********************************************************/
Voxel* Simulator::voxel_at(const glm::dvec3& position) const {
	if (geom_lookup_) {
		return geom_lookup_->voxel_at(position);
	}
	// Fallback if geom_lookup not initialized
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
	if (geom_lookup_) {
		return geom_lookup_->voxel_corners(voxel);
	}
	// Fallback if geom_lookup not initialized
	if (!voxel) {
		return Cuboid(); // Return default/empty cuboid
	}

	// For now, use the first medium (simple case)
	// In multi-medium scenario, we need a way to identify which medium owns the voxel
	if (!mediums.empty()) {
		return mediums[0].voxel_corners(voxel);
	}

	return Cuboid(); // Return default/empty cuboid
}

/***********************************************************
 * Check if point is inside geometry by checking any medium.
 ***********************************************************/
bool Simulator::is_point_inside_geometry(const glm::dvec3& position) const {
	if (geom_lookup_) {
		return geom_lookup_->is_point_inside_geometry(position);
	}
	// Fallback if geom_lookup not initialized
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
	if (geom_lookup_) {
		return geom_lookup_->voxel_grid(x, y, z);
	}
	// Fallback if geom_lookup not initialized
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
 * Accessor methods for aggregating data across all mediums
 ***********************************************************/

std::vector<Material> Simulator::get_all_tissues() const {
	if (aggregator_) {
		return aggregator_->get_all_tissues();
	}
	// Fallback if aggregator not initialized
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
	if (aggregator_) {
		return aggregator_->get_all_layers();
	}

	// Fallback if aggregator not initialized
	// For simplicity, return the first medium's layers
	// If no mediums exist, return empty static vector
	static const std::vector<Layer> empty_layers;
	if (mediums.empty()) {
		return empty_layers;
	}
	return mediums[0].get_layers();
}

std::vector<std::shared_ptr<Voxel>>& Simulator::get_all_voxels() {
	if (aggregator_) {
		return aggregator_->get_all_voxels();
	}

	// Fallback if aggregator not initialized
	// For non-const version, we need to return a reference to a static container
	// This is a bit tricky with multi-medium, but for now return the first medium's voxels
	// In practice, most use cases have only one medium
	static std::vector<std::shared_ptr<Voxel>> combined_voxels;
	combined_voxels.clear();

	for (auto& medium : mediums) {
		auto& volume = medium.get_volume();
		for (auto& voxel_ptr : volume) {
			// Convert unique_ptr to shared_ptr
			combined_voxels.push_back(std::shared_ptr<Voxel>(voxel_ptr.get(), [](Voxel*) {}));
		}
	}
	return combined_voxels;
}

const std::vector<std::shared_ptr<Voxel>>& Simulator::get_all_voxels() const {
	if (aggregator_) {
		return aggregator_->get_all_voxels();
	}
	// Fallback if aggregator not initialized
	// For const version, return a const reference to the same static container
	static std::vector<std::shared_ptr<Voxel>> combined_voxels;
	combined_voxels.clear();

	for (const auto& medium : mediums) {
		const auto& volume = medium.get_volume();
		for (const auto& voxel_ptr : volume) {
			// Convert unique_ptr to shared_ptr with non-owning deleter
			combined_voxels.push_back(std::shared_ptr<Voxel>(voxel_ptr.get(), [](Voxel*) {}));
		}
	}
	return combined_voxels;
}

size_t Simulator::get_total_voxel_count() const {
	if (aggregator_) {
		return aggregator_->get_total_voxel_count();
	}
	// Fallback if aggregator not initialized
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
	if (aggregator_) {
		return aggregator_->get_combined_bounds();
	}
	// Fallback if aggregator not initialized
	// Return the first medium's bounds for simplicity
	if (!mediums.empty()) {
		return mediums[0].get_bounds();
	}

	// Create default bounds
	return Range3(glm::dvec3(0.0), glm::dvec3(1.0));
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
		auto dda = std::make_unique<DDA>(grid_dimensions, grid_origin, voxel_size);
		medium_ddas_.push_back(std::move(dda));

		if (Config::get().log()) {
			std::ostringstream debug_msg;
			debug_msg << "  DDA " << i << ": Grid " << grid_dimensions.x << "x" << grid_dimensions.y << "x"
					  << grid_dimensions.z << ", Origin " << grid_origin.x << "," << grid_origin.y << ","
					  << grid_origin.z << ", Voxel size " << voxel_size;
			Logger::instance().log_debug(debug_msg.str());
		}
	}
}

/***********************************************************
 * Initialize component instances for modular functionality.
 ***********************************************************/
void Simulator::initialize_components() {
	// Create GeomLookup instance with references to mediums and DDAs
	geom_lookup_ = std::make_unique<GeomLookup>(mediums, medium_ddas_);

	// Create Aggregator instance with references to metrics, mediums, photons, and simulator
	aggregator_ = std::make_unique<Aggregator>(metrics, mediums, photons, *this);

	// Create Transport instance with references to simulation data
	photon_transport_ = std::make_unique<Transport>(mediums, metrics, medium_ddas_, rng, this);

	if (Config::get().log()) {
		FAST_LOG_INFO("Component initialization complete - GeomLookup, Aggregator, and Transport ready");
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
		// Very short step, handle normally via Transport
		if (photon_transport_) {
			photon_transport_->deposit(photon);
		}
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

	DDA* dda = medium_ddas_[medium_index].get();

	// Initialize DDA for this ray
	dda->initialize_ray(start_pos, direction);

	// Traverse voxels using DDA
	DDA::TraversalResult result = dda->traverse(total_distance);

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
				auto it = std::find_if(result.voxels.begin(), result.voxels.end(), [&step](const DDA::StepResult& s) {
					return s.voxel_coords == step.voxel_coords;
				});

				if (it != result.voxels.end()) {
					size_t index = std::distance(result.voxels.begin(), it);
					if (index < result.voxels.size() - 1) {
						segment_distance = result.voxels[index + 1].distance_traveled - step.distance_traveled;
					}
					else {
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
	} // Update photon weight and energy tracking
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
	if (geom_lookup_) {
		return geom_lookup_->find_medium_at_with_dda(position);
	}
	// Fallback if geom_lookup not initialized - use the old complex implementation
	// NUMERICAL STABILITY: Add epsilon tolerance for boundary detection
	static const double EPSILON = 1e-9;

	// Try each medium's DDA for precise boundary detection
	for (size_t i = 0; i < mediums.size() && i < medium_ddas_.size(); ++i) {
		const auto& medium = mediums[i];
		DDA* dda = medium_ddas_[i].get();
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
				glm::dvec3(-EPSILON * eps_scale, -EPSILON * eps_scale, EPSILON * eps_scale),
				glm::dvec3(-EPSILON * eps_scale, EPSILON * eps_scale, -EPSILON * eps_scale),
				glm::dvec3(-EPSILON * eps_scale, EPSILON * eps_scale, EPSILON * eps_scale),
				glm::dvec3(EPSILON * eps_scale, -EPSILON * eps_scale, -EPSILON * eps_scale),
				glm::dvec3(EPSILON * eps_scale, -EPSILON * eps_scale, EPSILON * eps_scale),
				glm::dvec3(EPSILON * eps_scale, EPSILON * eps_scale, -EPSILON * eps_scale),
				glm::dvec3(EPSILON * eps_scale, EPSILON * eps_scale, EPSILON * eps_scale)};

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
							glm::dvec3 search_pos =
								position + glm::dvec3(dx * search_radius, dy * search_radius, dz * search_radius);
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
