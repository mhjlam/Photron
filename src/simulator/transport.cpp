/**
 * @file transport.cpp
 * @brief Implementation of photon transport and interaction handlers
 */

#include "transport.hpp"

#include <algorithm>
#include <cmath>
#include <iostream>
#include <sstream>

#include "common/config.hpp"
#include "common/error_handler.hpp"
#include "common/logger.hpp"
#include "math/dda.hpp"
#include "math/math.hpp"
#include "math/ray.hpp"
#include "simulator/medium.hpp"
#include "simulator/metrics.hpp"
#include "simulator/photon.hpp"
#include "simulator/simulator.hpp"
#include "simulator/voxel.hpp"

Transport::Transport(std::vector<Medium>& mediums,
					 Metrics& metrics,
					 const std::vector<std::unique_ptr<DDA>>& dda_instances,
					 std::shared_ptr<Random> rng,
					 Simulator* simulator) :
	mediums_(&mediums), metrics_(&metrics), dda_instances_(&dda_instances), rng_(rng), simulator_(simulator) {
}

void Transport::cross(Photon& photon) {
	// Handle photon boundary crossing and medium transitions
	// Track boundary crossings (log mode only, limited output)
	static int crossing_count = 0;
	if (Config::get().log() && crossing_count < 20) { // Limit to first 20 crossings
		std::ostringstream debug_msg;
		debug_msg << "=== BOUNDARY CROSSING === Photon " << photon.id << " at pos=(" << photon.position.x << ","
				  << photon.position.y << "," << photon.position.z << ") weight=" << photon.weight;
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
	glm::dvec3 next_position = simulator_->move(photon.position, photon.direction, photon.sub_step);
	bool next_in_geometry = simulator_->is_point_inside_geometry(next_position);

	if (!next_in_geometry) {
		// Photon is exiting the medium - ensure intersect point is correct
		// If photon.intersect wasn't set by mesh detection, compute the actual exit point
		bool current_in_geometry = simulator_->is_point_inside_geometry(photon.position);
		if (current_in_geometry) {
			// Find the exact exit point by intersecting with geometry
			Ray exit_ray(photon.position, photon.direction);
			glm::dvec3 true_exit_point {0.0, 0.0, 0.0};
			bool found_exit = false;
			double min_exit_dist = std::numeric_limits<double>::max();

			// Find closest geometry exit - check all mediums for exit points
			for (const auto& medium : delegate_get_mediums()) {
				for (const auto& layer : medium.get_layers()) {
					for (const auto& triangle : layer.mesh) {
						Triangle triangle_copy = triangle;
						glm::dvec3 intersection;

						if (exit_ray.intersect_triangle(triangle_copy, intersection)) {
							double dist = glm::length(photon.position - intersection);

							if (dist > 1e-10 && dist < min_exit_dist) {
								// Verify this intersection takes us outside
								glm::dvec3 test_point = intersection + photon.direction * 1e-6;
								if (!simulator_->is_point_inside_geometry(test_point)) {
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

		double reflection = simulator_->internal_reflection(photon, eta, transmittance, reflectance);

		if (reflection == 0.0) {
			// Total transmission - photon exits
			photon.direction = transmittance;
			photon.alive = false;
			radiate(photon, transmittance, photon.weight);
		}
		else if (reflection == 1.0) {
			// Total internal reflection
			photon.direction = reflectance;
			photon.position = simulator_->move_delta(photon.intersect, photon.direction);
			photon.voxel = simulator_->voxel_at(photon.position);
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
				photon.position = simulator_->move_delta(photon.intersect, photon.direction);
				photon.voxel = simulator_->voxel_at(photon.position);
			}
			else {
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
	double temp_reflection = simulator_->internal_reflection(photon, eta_ambient, transmittance, reflectance);

	// Determine which medium(s) are involved in this crossing
	Medium* current_medium = delegate_find_medium_at_with_dda(photon.position);
	glm::dvec3 newpos = simulator_->move_delta(photon.intersect, transmittance);
	Medium* new_medium = delegate_find_medium_at_with_dda(newpos);
	Voxel* newvox = new_medium ? new_medium->voxel_at(newpos) : nullptr;

	// ROBUST VOXEL LOOKUP: Handle edge cases where move_delta places photon at voxel boundaries
	if (!newvox && new_medium) {
		// Try the intersect point directly first
		Voxel* intersect_voxel = new_medium->voxel_at(photon.intersect);
		if (intersect_voxel) {
			newvox = intersect_voxel;
		}
		else {
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
	temp_reflection = simulator_->internal_reflection(photon, eta, transmittance, reflectance);

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
		}
		else if (current_medium && new_medium && current_medium != new_medium) {
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
			oss << "  current_voxel=" << (current_voxel ? "yes" : "no") << " newvox=" << (newvox ? "yes" : "no");
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
					Logger::instance().log_debug(oss1.str());
					std::ostringstream oss2;
					oss2 << "  current_properties: eta=" << current_voxel->material->eta()
						 << " mua=" << current_voxel->material->mu_a() << " mus=" << current_voxel->material->mu_s()
						 << " g=" << current_voxel->material->g()
						 << " hash=" << current_voxel->material->get_optical_properties_hash();
					Logger::instance().log_debug(oss2.str());

					std::ostringstream oss3;
					oss3 << "  new_properties: eta=" << newvox->material->eta() << " mua=" << newvox->material->mu_a()
						 << " mus=" << newvox->material->mu_s() << " g=" << newvox->material->g()
						 << " hash=" << newvox->material->get_optical_properties_hash();
					Logger::instance().log_debug(oss3.str());

					std::ostringstream oss4;
					oss4 << "  optical_properties_same=" << (same_optical ? "yes" : "no");
					Logger::instance().log_debug(oss4.str());
				}
			}
		}

		if (current_voxel && newvox && current_voxel->material && newvox->material
			&& !current_voxel->material->has_same_optical_properties(*newvox->material)) {
			// INTERFACE ENERGY SPLITTING - Simple implementation
			double n1 = current_voxel->material->eta(); // From medium
			double n2 = newvox->material->eta();        // To medium

			// Calculate angle of incidence
			double cos_theta_i = -glm::dot(photon.direction, photon.voxel_normal);
			cos_theta_i = glm::clamp(cos_theta_i, 0.0, 1.0); // Ensure valid range

			// Check for total internal reflection
			double n_ratio = n1 / n2;
			double sin_theta_t_sq = n_ratio * n_ratio * (1.0 - cos_theta_i * cos_theta_i);

			if (sin_theta_t_sq > 1.0) {
				// TOTAL INTERNAL REFLECTION - all energy stays in current medium
				// Reflect photon direction
				photon.direction =
					photon.direction - 2.0 * glm::dot(photon.direction, photon.voxel_normal) * photon.voxel_normal;
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
			}
			else {
				// General case - Fresnel equations for s and p polarizations
				double Rs =
					std::pow((n1 * cos_theta_i - n2 * cos_theta_t) / (n1 * cos_theta_i + n2 * cos_theta_t), 2.0);
				double Rp =
					std::pow((n2 * cos_theta_i - n1 * cos_theta_t) / (n2 * cos_theta_i + n1 * cos_theta_t), 2.0);
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
			}
			else {
				// Apply Snell's law for refraction
				glm::dvec3 refracted_tangent = n_ratio * (incident - cos_theta_i * normal);
				glm::dvec3 refracted_direction = refracted_tangent + cos_theta_t * normal;
				photon.direction = glm::normalize(refracted_direction);
			}

			if (Config::get().log()) {
				std::ostringstream debug_msg;
				debug_msg << "  -> INTERFACE ENERGY SPLITTING: n1=" << n1 << ", n2=" << n2 << ", R=" << R_fresnel
						  << ", T=" << T_fresnel << ", reflected=" << reflected_energy
						  << ", transmitted=" << transmitted_energy;
				Logger::instance().log_debug(debug_msg.str());
			}
		}
	}

	// 1. crossing to ambient medium
	if (newvox == nullptr) {
		// Only treat as ambient exit if photon is actually leaving the geometry
		// Check if the new position is truly outside all medium geometries
		bool truly_exiting_geometry = true;
		for (const auto& medium : delegate_get_mediums()) {
			if (medium.contains_point(newpos)) {
				truly_exiting_geometry = false;
				break;
			}
		}

		if (!truly_exiting_geometry) {
			// FALSE AMBIENT EXIT: Photon is still inside geometry but newvox is null
			// This is a DDA traversal issue, not a true exit - continue transport normally
			photon.position = newpos;
			photon.voxel = simulator_->voxel_at(photon.position); // Find correct voxel at new position
			return;
		}

		// True ambient exit: maintain surface voxel assignment for proper emittance recording

		// Find correct surface voxel for emittance recording
		if (photon.voxel && !photon.voxel->is_surface_voxel) {
			// Photon was assigned to wrong voxel during transport - need to find correct exit voxel
			// Use the intersection point to find the surface voxel we're actually exiting from
			Medium* exit_medium = delegate_find_medium_at_with_dda(photon.position);
			if (exit_medium) {
				// Look for a surface voxel near the intersection point
				Voxel* intersection_voxel = exit_medium->voxel_at(photon.intersect);
				if (intersection_voxel && intersection_voxel->is_surface_voxel) {
					photon.voxel = intersection_voxel; // Correct the assignment
				}
				else {
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
			for (const auto& medium : delegate_get_mediums()) {
				if (medium.contains_point(photon.position)) {
					position_outside_geometry = false;
					break;
				}
			}

			// Only show detailed warnings in log mode
			if (Config::get().log()) {
				if (position_outside_geometry) {
					// Photon position is legitimately outside - voxelization error
					Logger::instance().log_error(
						"VOXELIZATION ERROR: Photon at position outside geometry but voxel marked as interior!");
					Logger::instance().log_error(
						"Voxel (" + std::to_string(photon.voxel->ix()) + ", " + std::to_string(photon.voxel->iy())
						+ ", " + std::to_string(photon.voxel->iz()) + ") should be surface but isn't.");
				}
				else {
					// Photon position is inside geometry - transport/exit detection issue
					Logger::instance().log_warning(
						"TRANSPORT ISSUE: Photon trying to exit from position still inside geometry!");
					Logger::instance().log_warning(
						"Position (" + std::to_string(photon.position.x) + ", " + std::to_string(photon.position.y)
						+ ", " + std::to_string(photon.position.z) + ") is inside medium but exit attempted.");
				}

				Logger::instance().log_warning("Warning: Photon attempting to exit from interior voxel at ("
											   + std::to_string(photon.voxel->ix()) + ", "
											   + std::to_string(photon.voxel->iy()) + ", "
											   + std::to_string(photon.voxel->iz()) + ") to ambient medium!");

				Logger::instance().log_warning("New position: (" + std::to_string(newpos.x) + ", "
											   + std::to_string(newpos.y) + ", " + std::to_string(newpos.z) + ")");

				Logger::instance().log_warning("Direction: (" + std::to_string(transmittance.x) + ", "
											   + std::to_string(transmittance.y) + ", "
											   + std::to_string(transmittance.z) + ")");
			}
		}

		// Use outward-pointing normals for exit calculations
		// Recalculate reflection and transmission with outward normal
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
			photon.position = simulator_->move_delta(photon.intersect, photon.direction);
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
				photon.position = simulator_->move_delta(photon.intersect, photon.direction);
				if (current_medium) {
					// For reflection, photon stays in same medium - safe to update voxel
					photon.voxel = current_medium->voxel_at(photon.position);
				}
			}
			else {
				// No reflection, terminate photon
				photon.alive = false;
			}
		}

		if (current_medium) {
			current_medium->get_metrics().increment_scatters();
		}
	}
	// 2. crossing to another medium
	else if (newvox != nullptr && newvox->material != nullptr && photon.voxel->material != nullptr
			 && photon.voxel->material->eta() != newvox->material->eta()) {
		// Use already computed reflection/transmission
		double reflection = temp_reflection;

		// total transmission
		if (reflection == 0.0) {
			photon.direction = transmittance;
			photon.position = simulator_->move_delta(photon.intersect, photon.direction);
			if (new_medium) {
				photon.voxel = new_medium->voxel_at(photon.position);
			}
		}
		// total internal reflection
		else if (reflection == 1.0) {
			photon.direction = reflectance;
			photon.position = simulator_->move_delta(photon.intersect, photon.direction);
			if (current_medium) {
				photon.voxel = current_medium->voxel_at(photon.position);
			}
		}
		else { // all-or-none transmission/reflection
			// total transmission
			if (rng_->next() > reflection) {
				photon.direction = transmittance;
				photon.position = simulator_->move_delta(photon.intersect, photon.direction);
				if (new_medium) {
					photon.voxel = new_medium->voxel_at(photon.position);
				}
			}
			// total reflection
			else {
				photon.direction = reflectance;
				photon.position = simulator_->move_delta(photon.intersect, photon.direction);
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
		photon.position = simulator_->move_delta(photon.intersect, photon.direction);
		if (current_medium) {
			// Smart voxel assignment to preserve surface voxel information
			Voxel* new_voxel =
				current_medium->voxel_at(photon.position); // Only update voxel assignment if we're sure it's correct
			// If current voxel is a surface voxel and new voxel is null/air, preserve current
			if (new_voxel && new_voxel->material) {
				photon.voxel = new_voxel; // Safe to assign - it's a material voxel
			}
			else if (!photon.voxel || !photon.voxel->is_surface_voxel) {
				// Only assign null/air voxels if current voxel isn't a surface voxel
				photon.voxel = new_voxel;
			}
			// If new_voxel is null/air but photon.voxel is surface voxel, preserve it
		}
	}
}

void Transport::sub_step(Photon& photon) {
	// ROBUST VOXEL SELECTION: Always find the correct voxel at current position
	Medium* photon_medium = delegate_find_medium_at(photon.position);
	if (!photon_medium) {
		FAST_LOG_ERROR("No medium found at photon position in sub_step()");
		if (photon.weight > 0.0) {
			terminate_photon_and_record_energy(photon, "no_material_properties");
		}
		else {
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
		}
		else {
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
		}
		else {
			photon.alive = false;
		}
		return;
	}

	// Get voxel boundaries for the CORRECT voxel
	Cuboid box = simulator_->voxel_corners(photon.voxel);

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
	if (adjusted_position.x < box.min_point().x || adjusted_position.x > box.max_point().x
		|| adjusted_position.y < box.min_point().y || adjusted_position.y > box.max_point().y
		|| adjusted_position.z < box.min_point().z || adjusted_position.z > box.max_point().z) {
		if (Config::get().log()) {
			Logger::instance().log_warning("Adjusted position outside voxel bounds. Using fallback.");
		}
		// Fallback: place photon at voxel center
		adjusted_position = glm::dvec3((box.min_point().x + box.max_point().x) * 0.5,
									   (box.min_point().y + box.max_point().y) * 0.5,
									   (box.min_point().z + box.max_point().z) * 0.5);
		on_boundary = true;
	}

	// Create ray from (robustly adjusted) photon position and direction
	Ray ray = Ray(adjusted_position, photon.direction);

	double voxdist = ray.intersect_cuboid_internal(box, photon.intersect, photon.voxel_normal);

	// ROBUST ERROR HANDLING: Multiple fallback strategies if intersection fails
	if (voxdist == std::numeric_limits<double>::max()) {
		if (Config::get().log()) {
			std::ostringstream debug_info;
			debug_info << "Primary ray-voxel intersection failed. Ray origin: (" << ray.origin().x << ", "
					   << ray.origin().y << ", " << ray.origin().z << "), direction: (" << ray.direction().x << ", "
					   << ray.direction().y << ", " << ray.direction().z << "), voxel bounds: [(" << box.min_point().x
					   << ", " << box.min_point().y << ", " << box.min_point().z << ") to (" << box.max_point().x
					   << ", " << box.max_point().y << ", " << box.max_point().z << ")]";
			ErrorHandler::instance().report_warning(debug_info.str());
		}

		// FALLBACK 1: Try from exact voxel center
		glm::dvec3 voxel_center = glm::dvec3((box.min_point().x + box.max_point().x) * 0.5,
											 (box.min_point().y + box.max_point().y) * 0.5,
											 (box.min_point().z + box.max_point().z) * 0.5);
		Ray fallback_ray1 = Ray(voxel_center, photon.direction);
		voxdist = fallback_ray1.intersect_cuboid_internal(box, photon.intersect, photon.voxel_normal);

		if (voxdist != std::numeric_limits<double>::max()) {
			if (Config::get().log()) {
				ErrorHandler::instance().report_info("Fallback 1 (voxel center) succeeded.");
			}
		}
		else {
			// FALLBACK 2: Use manual boundary calculation
			if (Config::get().log()) {
				ErrorHandler::instance().report_info("Fallback 1 failed. Using manual boundary calculation.");
			}

			// Find which boundary the ray will hit first
			double t_min = std::numeric_limits<double>::max();
			glm::dvec3 hit_point {0.0, 0.0, 0.0};
			glm::dvec3 hit_normal {0.0, 0.0, 0.0};

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
				if (std::abs(denom) > 1e-12) {    // Ray not parallel to face
					double t = glm::dot(face_point - adjusted_position, face_normal) / denom;
					if (t > 1e-12 && t < t_min) { // Valid forward intersection
						glm::dvec3 test_point = adjusted_position + t * photon.direction;

						// Check if intersection point is within face bounds
						bool within_bounds = true;
						if (face_normal.x != 0) { // YZ face
							within_bounds = (test_point.y >= box.min_point().y - EPSILON
											 && test_point.y <= box.max_point().y + EPSILON
											 && test_point.z >= box.min_point().z - EPSILON
											 && test_point.z <= box.max_point().z + EPSILON);
						}
						else if (face_normal.y != 0) { // XZ face
							within_bounds = (test_point.x >= box.min_point().x - EPSILON
											 && test_point.x <= box.max_point().x + EPSILON
											 && test_point.z >= box.min_point().z - EPSILON
											 && test_point.z <= box.max_point().z + EPSILON);
						}
						else { // XY face
							within_bounds = (test_point.x >= box.min_point().x - EPSILON
											 && test_point.x <= box.max_point().x + EPSILON
											 && test_point.y >= box.min_point().y - EPSILON
											 && test_point.y <= box.max_point().y + EPSILON);
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
					ErrorHandler::instance().report_info("Fallback 2 (manual calculation) succeeded.");
				}
			}
			else {
				// FALLBACK 3: Emergency exit - force photon to exit current voxel
				if (Config::get().log()) {
					ErrorHandler::instance().report_warning("All fallbacks failed. Forcing photon exit.");
				}
				photon.intersect = adjusted_position + photon.direction * EPSILON;
				photon.voxel_normal = -photon.direction; // Opposite to ray direction
				voxdist = EPSILON;
			}
		}
	}

	// Compute free path substep for scattering event
	double mu_s_effective = photon.mu_s();
	if (mu_s_effective <= 1e-12) {
		// For non-scattering media, use total attenuation or set very large free path
		mu_s_effective = std::max(photon.mu_a(), 1e-6);
	}
	double freepath = photon.step / mu_s_effective;

	// Check if photon is currently inside the mesh
	bool photon_inside_mesh = simulator_->is_point_inside_geometry(photon.position);

	if (!photon_inside_mesh) {
		// Photon already exited - avoid double-counting energy
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
	Medium* current_medium = delegate_find_medium_at(photon.position);
	if (!current_medium) {
		// Energy conservation: Deposit remaining energy as absorption before terminating
		if (photon.weight > 0.0) {
			terminate_photon_and_record_energy(photon, "no_medium_found");
		}
		else {
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
					if (!simulator_->is_point_inside_geometry(test_point)) {
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
		const double vertex_threshold = MathConstants::VERTEX_THRESHOLD;
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

void Transport::radiate(Photon& photon, glm::dvec3& direction, double weight) {
	// Record photon exit energy and find last material voxel before exit

	Medium* exit_medium = delegate_find_medium_at_with_dda(photon.position);
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
		glm::ivec3 intersection_coords = glm::ivec3((photon.intersect.x - grid_origin.x) / voxel_size,
													(photon.intersect.y - grid_origin.y) / voxel_size,
													(photon.intersect.z - grid_origin.z) / voxel_size);

		// Check neighboring voxels (3x3x3 neighborhood)
		for (int dx = -1; dx <= 1; ++dx) {
			for (int dy = -1; dy <= 1; ++dy) {
				for (int dz = -1; dz <= 1; ++dz) {
					glm::ivec3 neighbor_coords = intersection_coords + glm::ivec3(dx, dy, dz);

					// Check bounds
					if (neighbor_coords.x >= 0
						&& static_cast<uint32_t>(neighbor_coords.x) < exit_medium->get_volume().width()
						&& neighbor_coords.y >= 0
						&& static_cast<uint32_t>(neighbor_coords.y) < exit_medium->get_volume().height()
						&& neighbor_coords.z >= 0
						&& static_cast<uint32_t>(neighbor_coords.z) < exit_medium->get_volume().depth()) {
						Voxel* neighbor =
							exit_medium->get_volume().at(neighbor_coords.x, neighbor_coords.y, neighbor_coords.z);
						if (neighbor && neighbor->material) {
							// Calculate distance from intersection to voxel center
							glm::dvec3 voxel_center = grid_origin
													  + glm::dvec3((neighbor_coords.x + 0.5) * voxel_size,
																   (neighbor_coords.y + 0.5) * voxel_size,
																   (neighbor_coords.z + 0.5) * voxel_size);
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
	}
	else {
		// Select the candidate with highest score (most inside the medium)
		std::sort(candidate_voxels.begin(), candidate_voxels.end(), [](const auto& a, const auto& b) {
			return a.second > b.second;
		});
		exit_voxel = candidate_voxels[0].first;
	}

	if (!exit_voxel || !exit_voxel->material) {
		ErrorHandler::instance().report_error(ErrorMessage::format(
			SimulationError::NoVoxelFound, "Cannot determine last material voxel for emittance recording"));
		FAST_LOG_ERROR("Cannot determine last material voxel for emittance recording");
		return;
	}

	glm::ivec3 surface_coords = glm::ivec3(exit_voxel->ix(), exit_voxel->iy(), exit_voxel->iz());

	// Log the radiate event with the determined exit voxel
	Logger::instance().log_photon_event(static_cast<int>(photon.id),
										"RADIATE",
										photon.position,
										direction,
										weight,
										surface_coords,
										-1,
										weight,
										"Using robust last material voxel selection");

	// Record photon exit with energy conservation enforcement
	// True Splitting handles energy conservation through statistical splitting
	photon.radiate_call_count++;
	photon.total_energy_radiated += weight;

	// Record emittance at the LAST material VOXEL (before exit)
	double old_emittance = exit_voxel->emittance;
	(void)old_emittance; // Suppress unused variable warning - kept for debugging
	exit_voxel->emittance += weight;

	// Log voxel emittance recording
	Logger::instance().log_voxel_emittance(static_cast<int>(photon.id),
										   photon.position,
										   direction,
										   weight,
										   surface_coords,
										   weight,
										   "Surface voxel emittance");

	// Use proper reflection/transmission determination based on exit position relative to entry
	bool is_reflecting = simulator_->is_photon_reflecting(photon);
	if (is_reflecting) {
		// Exit on same side as entry - classify as diffuse reflection (not specular)
		exit_voxel->diffuse_reflection += weight;
		photon.exit_type = Photon::ExitType::REFLECTED;
		Logger::instance().log_photon_event(static_cast<int>(photon.id),
											"REFLECT",
											photon.position,
											direction,
											weight,
											surface_coords,
											-1,
											weight,
											"Photon classified as reflection");
	}
	else {
		// Exit on opposite side from entry - classify as transmission
		exit_voxel->diffuse_transmission += weight;
		photon.exit_type = Photon::ExitType::TRANSMITTED;
		Logger::instance().log_photon_event(static_cast<int>(photon.id),
											"TRANSMIT",
											photon.position,
											direction,
											weight,
											surface_coords,
											-1,
											weight,
											"Photon classified as transmission");
	}

	// Create external vertex for photon path with proper exit classification
	// Convert Photon::ExitType to PhotonNode::ExitType
	PhotonNode::ExitType node_exit_type = PhotonNode::ExitType::NONE;
	if (photon.exit_type == Photon::ExitType::REFLECTED) {
		node_exit_type = PhotonNode::ExitType::REFLECTED;
	}
	else if (photon.exit_type == Photon::ExitType::TRANSMITTED) {
		node_exit_type = PhotonNode::ExitType::TRANSMITTED;
	}

	// Add external vertex to the current photon path
	std::shared_ptr<PhotonNode> exit_node = nullptr;
	if (photon.path_last) {
		exit_node =
			std::make_shared<PhotonNode>(photon.intersect, weight, static_cast<PhotonNode::ExitType>(node_exit_type));
		photon.add_external_vertex(exit_node);
	}

	// Create Emitter object with proper exit classification for renderer use
	Emitter::ExitType emitter_exit_type = Emitter::ExitType::NONE;
	if (photon.exit_type == Photon::ExitType::REFLECTED) {
		emitter_exit_type = Emitter::ExitType::REFLECTED;
	}
	else if (photon.exit_type == Photon::ExitType::TRANSMITTED) {
		emitter_exit_type = Emitter::ExitType::TRANSMITTED;
	}

	// Create shared emitter and establish connection with exit node
	auto emitter = std::make_shared<Emitter>(photon.id, photon.intersect, direction, weight, emitter_exit_type);
	simulator_->get_emitters().push_back(emitter);

	// Establish bidirectional connection between exit node and emitter
	if (exit_node) {
		exit_node->emitter = emitter;
	}
}

void Transport::transfer(Photon& photon) {
	// Move photon through medium until next scattering event or boundary
	// This method handles photon transport through the medium including boundary crossings

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

	while (photon.step >= MathConstants::STEP_SIZE_THRESHOLD && photon.alive) {
		substep_counter++;

		// Track photon state at each step
		if (substep_counter <= 50) { // Limit output
			std::ostringstream oss;
			oss << "Step " << substep_counter << ": Photon " << photon.id << " weight=" << photon.weight
				<< " step=" << photon.step << " pos=(" << photon.position.x << "," << photon.position.y << ","
				<< photon.position.z << ")"
				<< " alive=" << (photon.alive ? "yes" : "no");
			FAST_LOG_DEBUG(oss.str());
		}

		if (substep_counter > max_substeps) {
			FAST_LOG_WARNING("Photon exceeded maximum substeps, terminating.");

			// Deposit remaining energy as absorption for energy conservation
			if (photon.weight > 0.0 && photon.voxel && photon.voxel->material) {
				// Use energy conservation enforcement
				terminate_photon_and_record_energy(photon, "max_iterations");
			}
			else {
				photon.alive = false;
			}
			break;
		}

		// set substep - this will handle mesh boundary detection
		sub_step(photon);

		// Per-step absorption
		deposit(photon);

		// possibly cross boundary
		if (photon.cross) {
			cross(photon);
		}
		else {
			photon.position = simulator_->move(photon.position, photon.direction, photon.sub_step);
		}

		// Prevent errors due to crossing to ambient medium
		Medium* current_medium = delegate_find_medium_at_with_dda(photon.position);
		if (!current_medium) {
			// Record as transmission - photon is exiting
			glm::ivec3 last_voxel_coords =
				photon.voxel ? glm::ivec3(photon.voxel->ix(), photon.voxel->iy(), photon.voxel->iz()) : glm::ivec3(-1);

			Logger::instance().log_photon_event(static_cast<int>(photon.id),
												"EXIT",
												photon.position,
												photon.direction,
												photon.weight,
												last_voxel_coords,
												-1,
												photon.weight,
												"Photon exiting medium - calling radiate");

			photon.alive = false;
			radiate(photon, photon.direction, photon.weight);
			return;
		}

		// Set photon's voxel to match new position but preserve material assignment at exit boundaries
		Voxel* new_voxel = current_medium->voxel_at(photon.position);

		if (new_voxel) {
			// Normal case: photon is in a valid material voxel
			photon.voxel = new_voxel;
		}
		else {
			// Keep photon assigned to last material voxel until actual medium exit
			if (photon.voxel && photon.voxel->material) {
				// Photon still has a valid material voxel from previous step
				// Check if photon is actually exiting the medium geometry
				if (!current_medium->contains_point(photon.position)) {
					// Photon has truly exited the medium - proceed with exit logic

					glm::ivec3 last_voxel_coords(-1);
					Logger::instance().log_photon_event(static_cast<int>(photon.id),
														"EXIT",
														photon.position,
														photon.direction,
														photon.weight,
														last_voxel_coords,
														-1,
														photon.weight,
														"Photon moved to invalid voxel - calling radiate");

					photon.alive = false;
					radiate(photon, photon.direction, photon.weight);
					return;
				}
			}
			else {
				// Photon has no previous material voxel (error state)
				glm::ivec3 last_voxel_coords(-1);
				Logger::instance().log_photon_event(static_cast<int>(photon.id),
													"EXIT",
													photon.position,
													photon.direction,
													photon.weight,
													last_voxel_coords,
													-1,
													photon.weight,
													"Photon has no material voxel - calling radiate");

				photon.alive = false;
				radiate(photon, photon.direction, photon.weight);
				return;
			}
		}

		// Update step size using current medium's material properties
		if (current_medium && photon.voxel) {
			photon.step -= (photon.sub_step * photon.mu_s());
		}

		current_medium->get_metrics().add_vertex(photon.position.x, photon.position.y, photon.position.z);
	}
}

// handle_medium_transition() method moved to end of file after extraction

Voxel* Transport::find_last_surface_voxel_with_dda(const Photon& photon, const glm::dvec3& exit_direction) {
	if (simulator_->geom_lookup_) {
		return simulator_->geom_lookup_->find_last_surface_voxel_with_dda(photon, exit_direction);
	}
	// Fallback if geom_lookup not initialized
	// Find the medium the photon was in
	Medium* current_medium = simulator_->find_medium_at(photon.position);
	if (!current_medium) {
		return nullptr;
	}

	// Find the DDA instance for this medium
	size_t medium_index = 0;
	for (size_t i = 0; i < simulator_->mediums.size(); ++i) {
		if (&simulator_->mediums[i] == current_medium) {
			medium_index = i;
			break;
		}
	}

	if (medium_index >= simulator_->medium_ddas_.size()) {
		// Fallback: try to find voxel with step-back approach
		glm::dvec3 step_back_pos = photon.position - exit_direction * 1e-6;
		Medium* last_medium = simulator_->find_medium_at_with_dda(step_back_pos);
		if (last_medium) {
			glm::dvec3 mutable_pos = step_back_pos;
			return last_medium->voxel_at(mutable_pos);
		}
		return nullptr;
	}

	DDA* dda = simulator_->medium_ddas_[medium_index].get();

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
	DDA::TraversalResult result = dda->traverse(total_distance);

	// Find the last surface voxel in the traversal
	Voxel* last_surface_voxel = nullptr;

	Logger::instance().log_photon_event(static_cast<int>(photon.id),
										"DDA_SEARCH",
										photon.position,
										direction,
										0.0,
										glm::ivec3(-1),
										-1,
										total_distance,
										"Starting DDA traversal for last surface voxel, found "
											+ std::to_string(result.voxels.size()) + " voxels");

	for (size_t i = 0; i < result.voxels.size(); ++i) {
		const auto& step = result.voxels[i];
		// Get voxel at this DDA position
		glm::dvec3 mutable_pos = step.world_position;  // Create mutable copy for voxel_at
		Voxel* voxel = current_medium->voxel_at(mutable_pos);

		if (voxel && voxel->material) {
			glm::ivec3 voxel_coords = glm::ivec3(voxel->ix(), voxel->iy(), voxel->iz());
			bool is_surface = voxel->is_surface_voxel; // ONLY true external surface voxels

			Logger::instance().log_photon_event(
				static_cast<int>(photon.id),
				"DDA_VOXEL",
				step.world_position,
				direction,
				0.0,
				voxel_coords,
				static_cast<int>(i),
				step.distance_traveled,
				"Voxel " + std::to_string(i) + "/" + std::to_string(result.voxels.size())
					+ (is_surface ? std::string(" SURFACE") : std::string(" INTERIOR")));

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
			static_cast<int>(photon.id),
			"DDA_RESULT",
			photon.position,
			direction,
			0.0,
			selected_coords,
			-1,
			0.0,
			"Selected voxel: " + (last_surface_voxel ? std::string("DDA_FOUND") : std::string("FALLBACK_TO_CURRENT")));
	}
	else {
		Logger::instance().log_photon_event(static_cast<int>(photon.id),
											"DDA_ERROR",
											photon.position,
											direction,
											0.0,
											glm::ivec3(-1),
											-1,
											0.0,
											"No voxel selected - both DDA and fallback failed");
	}

	return selected_voxel;
}

// Private helper methods

double Transport::calculate_fresnel_reflection(double cos_i, double eta_i, double eta_t) const {
	// Basic Fresnel reflection calculation
	double eta = eta_i / eta_t;
	double sin_t_sq = eta * eta * (1.0 - cos_i * cos_i);

	if (sin_t_sq >= 1.0) {
		return 1.0; // Total internal reflection
	}

	double cos_t = std::sqrt(1.0 - sin_t_sq);

	double r_s = (eta_i * cos_i - eta_t * cos_t) / (eta_i * cos_i + eta_t * cos_t);
	double r_p = (eta_t * cos_i - eta_i * cos_t) / (eta_t * cos_i + eta_i * cos_t);

	return 0.5 * (r_s * r_s + r_p * r_p);
}

glm::dvec3 Transport::calculate_specular_reflection_direction(const glm::dvec3& incident,
															  const glm::dvec3& normal) const {
	return incident - 2.0 * glm::dot(incident, normal) * normal;
}

glm::dvec3 Transport::calculate_refraction_direction(const glm::dvec3& incident,
													 const glm::dvec3& normal,
													 double eta) const {
	double cos_i = -glm::dot(incident, normal);
	double sin_t_sq = eta * eta * (1.0 - cos_i * cos_i);

	if (sin_t_sq >= 1.0) {
		return glm::dvec3(0.0); // Total internal reflection
	}

	double cos_t = std::sqrt(1.0 - sin_t_sq);
	return eta * incident + (eta * cos_i - cos_t) * normal;
}

// =============================================================================
// =============================================================================
// Delegation methods for accessing simulator functionality
// =============================================================================

Medium* Transport::delegate_find_medium_at_with_dda(const glm::dvec3& position) {
	if (simulator_) {
		return simulator_->find_medium_at_with_dda(position);
	}
	return nullptr;
}

Medium* Transport::delegate_find_medium_at(const glm::dvec3& position) {
	if (simulator_) {
		return simulator_->find_medium_at(position);
	}
	return nullptr;
}

const std::vector<Medium>& Transport::delegate_get_mediums() const {
	if (simulator_) {
		return simulator_->mediums;
	}
	// Return empty reference - this is a fallback that should not be used
	static const std::vector<Medium> empty_mediums;
	return empty_mediums;
}

// ============================================================================
// EXTRACTED METHOD IMPLEMENTATIONS
// ============================================================================

void Transport::terminate_photon_and_record_energy(Photon& photon, const std::string& reason) {
	// Track all photon terminations
	static int termination_count = 0;
	static double total_terminated_energy = 0.0;
	termination_count++;
	total_terminated_energy += photon.weight;

	if (!photon.alive || photon.weight <= 0.0) {
		return; // Already terminated or no energy to record
	}

	// Find the appropriate medium to record energy
	Medium* record_medium = delegate_find_medium_at(photon.position);
	const auto& mediums = delegate_get_mediums();
	if (!record_medium && !mediums.empty()) {
		record_medium = const_cast<Medium*>(&mediums[0]); // Fallback to first medium
	}

	if (!record_medium) {
		std::cerr << "CRITICAL ERROR: Cannot record energy - no medium found for photon termination: " << reason
				  << std::endl;
		photon.alive = false;
		return;
	}

	// CONSOLIDATED APPROACH: This function only handles INTERNAL terminations (absorption)
	// For exits, use radiate() function which handles both voxel emittance and medium records
	if (reason == "absorption" || reason == "roulette" || reason == "max_iterations" || reason == "no_voxel_found"
		|| reason == "no_material_properties" || reason == "no_medium_found") {
		// Energy absorbed within the medium
		record_medium->get_metrics().add_total_absorption(photon.weight);

		// Track energy in photon accounting system
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

		// Track energy in photon accounting system
		if (photon.voxel && photon.voxel->material) {
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

void Transport::handle_medium_transition(Photon& photon, Medium* from, Medium* to) {
	if (!from && to) {
		// Entering a medium from ambient space
		to->increment_photons_entered();
	}
	else if (from && !to) {
		// Exiting to ambient space
		photon.alive = false;
	}
	else if (from && to && from != to) {
		// Multi-layer interface physics with Fresnel calculations

		// Get material properties for both media
		Material* from_material = nullptr;
		Material* to_material = nullptr;

		// Find representative voxels to get material properties
		glm::dvec3 from_pos = photon.position - photon.direction * 1e-6; // Slightly behind
		glm::dvec3 to_pos = photon.position + photon.direction * 1e-6;   // Slightly ahead

		Voxel* from_voxel = from->voxel_at(from_pos);
		Voxel* to_voxel = to->voxel_at(to_pos);

		if (from_voxel && from_voxel->material)
			from_material = from_voxel->material;
		if (to_voxel && to_voxel->material)
			to_material = to_voxel->material;

		if (!from_material || !to_material) {
			ErrorHandler::instance().report_warning("Interface transition without proper material properties");
			FAST_LOG_WARNING("Interface transition without proper material properties");
			// Fallback to simple transmission
			return;
		}

		// Get refractive indices
		double n1 = from_material->eta(); // Incident medium
		double n2 = to_material->eta();   // Transmitted medium

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
		(void)sin_theta_i;                           // Suppress unused variable warning - kept for debugging

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
		}
		else {
			// General case - calculate s and p polarization components
			Rs = std::pow((n1 * cos_theta_i - n2 * cos_theta_t) / (n1 * cos_theta_i + n2 * cos_theta_t), 2.0);
			Rp = std::pow((n1 * cos_theta_t - n2 * cos_theta_i) / (n1 * cos_theta_t + n2 * cos_theta_i), 2.0);
			R_fresnel = 0.5 * (Rs + Rp); // Average for unpolarized light
		}

		// Ensure valid reflection coefficient
		R_fresnel = glm::clamp(R_fresnel, 0.0, 1.0);
		double T_fresnel = 1.0 - R_fresnel; // Transmission coefficient

		// Apply probabilistic Fresnel reflection/transmission (no energy splitting)

		if (rng_->next() < R_fresnel) {
			// FRESNEL REFLECTION
			photon.direction = glm::reflect(photon.direction, interface_normal);
			photon.direction = glm::normalize(photon.direction);

			if (Config::get().log()) {
				std::ostringstream debug_msg;
				debug_msg << "Fresnel reflection at medium interface (R=" << R_fresnel << ")";
				Logger::instance().log_debug(debug_msg.str());
			}
		}
		else {
			// FRESNEL TRANSMISSION WITH REFRACTION

			// Calculate refracted direction using Snell's law
			glm::dvec3 transmitted_dir;

			if (cos_theta_i > 0.9999) {
				// Near-normal incidence - no direction change
				transmitted_dir = photon.direction;
			}
			else {
				// General refraction using vector form of Snell's law
				glm::dvec3 tangential = photon.direction - cos_theta_i * interface_normal;
				transmitted_dir = n_ratio * tangential + (n_ratio * cos_theta_i - cos_theta_t) * interface_normal;
			}

			photon.direction = glm::normalize(transmitted_dir);

			// For medium transitions, photon continues with full weight

			if (Config::get().log()) {
				std::ostringstream debug_msg;
				debug_msg << "Fresnel transmission at medium interface (T=" << T_fresnel
						  << ", angle_i=" << std::acos(cos_theta_i) * MathConstants::RAD_TO_DEG
						  << "°, angle_t=" << std::acos(cos_theta_t) * MathConstants::RAD_TO_DEG << "°)";
				Logger::instance().log_debug(debug_msg.str());
			}
		}

		// VALIDATION: Ensure photon direction is physically reasonable
		if (glm::length(photon.direction) < 0.9 || glm::length(photon.direction) > 1.1) {
			ErrorHandler::instance().report_warning(
				"Invalid photon direction after interface transition. Normalizing.");
			FAST_LOG_WARNING("Invalid photon direction after interface transition. Normalizing.");
			photon.direction = glm::normalize(photon.direction);
		}

		// ENERGY CONSERVATION CHECK
		if (photon.weight < 0.0 || photon.weight > 1.0) {
			std::string weight_msg =
				"Invalid photon weight after interface transition: " + std::to_string(photon.weight);
			ErrorHandler::instance().report_warning(weight_msg);
			FAST_LOG_WARNING(weight_msg);
			photon.weight = glm::clamp(photon.weight, 0.0, 1.0);
		}

		// POST-INTERFACE VALIDATION: Ensure photon is in valid state
		simulator_->validate_photon_state_after_interface_transition(photon, from, to);
	}
}

void Transport::scatter(Photon& photon) {
	if (!photon.alive) {
		return;
	}

	// Validate scattering medium and material properties
	Medium* current_medium = delegate_find_medium_at(photon.position);
	if (!current_medium) {
		// Photon has exited medium during transport - record as transmission
		photon.alive = false;
		radiate(photon, photon.direction, photon.weight);
		return;
	}

	Material* material = photon.voxel->material;
	if (!material) {
		// No material properties available - cannot scatter, record as transmission
		photon.alive = false;
		radiate(photon, photon.direction, photon.weight);
		return;
	}

	// Sample scattering angles using Henyey-Greenstein phase function
	double cos_theta, sin_theta, cos_phi, sin_phi;
	double g = material->g(); // anisotropy factor [-1, 1]
	double rnd = rng_->next();

	// Numerically stable Henyey-Greenstein sampling for polar angle
	if (std::abs(g) > 1e-12) {
		// Use improved numerical stability for extreme anisotropy values
		double g2 = g * g;
		double one_minus_g2 = 1.0 - g2;
		double temp = one_minus_g2 / (1.0 - g + 2.0 * g * rnd);
		cos_theta = (1.0 + g2 - temp * temp) / (2.0 * g);

		// Robust clamping to prevent numerical errors
		cos_theta = std::clamp(cos_theta, -1.0, 1.0);
	}
	else {
		// Isotropic scattering for negligible anisotropy
		cos_theta = 2.0 * rnd - 1.0;
	}

	// Calculate sine component with numerical safety
	sin_theta = std::sqrt(std::max(0.0, 1.0 - cos_theta * cos_theta));

	// Sample uniform azimuthal angle
	rnd = rng_->next();
	double phi = MathConstants::TWO_PI * rnd;
	cos_phi = std::cos(phi);
	sin_phi = std::sin(phi);

	// Transform scattered direction to world coordinates
	// Uses robust coordinate transformation avoiding singularities
	glm::dvec3 old_direction = photon.direction;

	// Special handling for directions nearly aligned with z-axis
	if (std::abs(old_direction.z) > 0.99999) {
		// Direct coordinate assignment for near-vertical directions
		photon.direction =
			glm::dvec3(sin_theta * cos_phi, sin_theta * sin_phi, cos_theta * (old_direction.z > 0 ? 1.0 : -1.0));
	}
	else {
		// General case: construct orthonormal basis for stable rotation
		glm::dvec3 w = old_direction; // incident direction
		glm::dvec3 u = glm::normalize(glm::cross(std::abs(w.z) < 0.9 ? glm::dvec3(0, 0, 1) : glm::dvec3(1, 0, 0), w));
		glm::dvec3 v = glm::cross(w, u);

		// Construct scattered direction in local coordinate system
		photon.direction = sin_theta * cos_phi * u + sin_theta * sin_phi * v + cos_theta * w;
	}

	// Increment scatter count for this photon
	photon.scatter_count++;
	photon.scatters = true; // Mark that this photon has scattered at least once

	// normalize direction vector (safety check)
	photon.direction = glm::normalize(photon.direction);

	// prevent scattering into ambient medium when close to boundaries
	glm::dvec3 temp_pos = photon.position;
	glm::dvec3 temp_dir = photon.direction;
	glm::dvec3 newpos = simulator_->move_delta(temp_pos, temp_dir);
	if (!simulator_->is_inside_any_medium(newpos)) {
		photon.alive = false;
		radiate(photon, photon.direction, photon.weight);
		return;
	}

	current_medium->get_metrics().increment_scatters();

	// add new internal position to path
	photon.add_internal_vertex(std::make_shared<PhotonNode>(photon.position, photon.weight));
}

void Transport::roulette(Photon& photon) {
	// Never apply Russian Roulette to negative weights
	if (photon.weight < 0.0) {
		if (Config::get().log()) {
			std::ostringstream debug_msg;
			debug_msg << "ERROR: Attempted Russian Roulette on negative weight (" << photon.weight
					  << "), terminating photon";
			Logger::instance().log_error(debug_msg.str());
		}
		terminate_photon_and_record_energy(photon, "negative_weight");
		return;
	}

	// Modernized Russian roulette with adaptive threshold and survival probability
	if (photon.weight < MCML_WEIGHT_THRESHOLD) {
		// Use weight-dependent survival probability for better variance reduction
		double survival_probability = std::max(0.1, photon.weight / MCML_WEIGHT_THRESHOLD);
		survival_probability = std::min(survival_probability, 0.5); // Cap at 50% for stability

		if (rng_->next() <= survival_probability) {
			// Proper Russian Roulette

			// Survive with proper weight normalization
			double old_weight = photon.weight;
			photon.weight /= survival_probability;

			// Update energy budget to match new weight
			// Russian Roulette increases the weight, so budget must increase proportionally
			photon.total_energy_budget = photon.weight;

			if (Config::get().log()) {
				std::ostringstream oss;
				oss << "ROULETTE: weight " << old_weight << " -> " << photon.weight
					<< " (survival_prob=" << survival_probability << ", new_budget=" << photon.total_energy_budget
					<< ")";
				Logger::instance().log_debug(oss.str());
			}
		}
		else {
			// Use extracted termination for consistency
			terminate_photon_and_record_energy(photon, "roulette");
		}
	}
}

void Transport::step_size(Photon& photon) {
	// Find current medium for the photon
	Medium* current_medium = delegate_find_medium_at_with_dda(photon.position);
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
	if (photon.step < MathConstants::PRECISION_THRESHOLD) { // Higher precision threshold
		double rnd;
		// More robust zero-avoidance for random number generation
		// Use std::numeric_limits for better precision handling
		constexpr double min_random = std::numeric_limits<double>::epsilon();
		do {
			rnd = rng_->next();
		}
		while (rnd <= min_random);

		// Use high-precision logarithm for step size calculation
		photon.step = -std::log(rnd);
	}

	// Add step size to metrics of the current medium
	// Prevent division by zero for non-scattering media
	double normalization_factor = photon.mu_s();
	if (normalization_factor <= MathConstants::PRECISION_THRESHOLD) {
		normalization_factor =
			std::max(photon.mu_a(), MathConstants::PHOTON_NUDGE_EPSILON); // Use absorption or minimum value
	}
	current_medium->get_metrics().add_step_size(photon.step / normalization_factor);
}

void Transport::launch(Photon& photon, const Source& source) {
	// Determine starting medium for photon transport
	Medium* start_medium = delegate_find_medium_at_with_dda(source.intersect);
	if (!start_medium) {
		ErrorHandler::instance().report_error(
			ErrorMessage::format(SimulationError::NoMediumFound, "Photon launch point is not in any medium"));
		FAST_LOG_ERROR("Photon launch point is not in any medium");
		photon.alive = false;
		return;
	}

	// Set medium statistics for photon entry tracking
	start_medium->increment_photons_entered();

	// Initialize basic photon transport properties
	photon.alive = true;
	photon.weight = 1.0; // Each photon starts with full energy

	// Initialize energy conservation tracking system
	photon.total_energy_budget = photon.weight;
	photon.total_energy_radiated = 0.0;
	photon.total_energy_absorbed = 0.0;
	photon.radiate_call_count = 0;

	// Initialize scattering event counters for statistical analysis
	photon.scatter_count = 0;
	photon.scatters = false;

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
	glm::ivec3 voxel_coords =
		photon.voxel ? glm::ivec3(photon.voxel->ix(), photon.voxel->iy(), photon.voxel->iz()) : glm::ivec3(-1);
	Logger::instance().log_photon_event(static_cast<int>(photon.id),
										"LAUNCH",
										photon.position,
										photon.direction,
										photon.weight,
										voxel_coords,
										0,
										0.0,
										"Photon launched into medium");

	// Compute specular reflection for this photon
	specular_reflection(photon);

	// If initial voxel_at failed but specular_reflection succeeded,
	// try to get the voxel again using the same nudging approach
	if (!photon.voxel && photon.alive) {
		// Use the same nudging approach as in specular_reflection
		const double epsilon = MathConstants::PHOTON_NUDGE_EPSILON;
		glm::dvec3 nudged_position = photon.position + epsilon * photon.direction;

		auto* nudged_medium = delegate_find_medium_at(nudged_position);
		if (nudged_medium) {
			photon.voxel = nudged_medium->voxel_at(nudged_position);
			if (photon.voxel) {
				// Set position to the nudged position for consistent transport
				photon.position = nudged_position;
			}
		}
	}

	// After specular reflection processing, move photon slightly
	// inside the medium to avoid boundary issues in subsequent transport steps
	if (photon.alive && photon.voxel) {
		// Move photon a small distance along the original direction to get off the exact surface
		const double surface_epsilon = MathConstants::SURFACE_NUDGE_EPSILON;
		glm::dvec3 nudged_position = photon.position + surface_epsilon * photon.direction;

		// Verify the nudged position is still in the medium
		auto* nudged_medium = delegate_find_medium_at(nudged_position);
		if (nudged_medium) {
			photon.position = nudged_position;
			photon.voxel = nudged_medium->voxel_at(photon.position);
		}
	}

	// create vertices for new light path
	auto light = std::make_shared<PhotonNode>(photon.source.origin, photon.weight);
	auto intersection = std::make_shared<PhotonNode>(photon.source.intersect, photon.weight);
	auto reflection =
		std::make_shared<PhotonNode>(simulator_->move(photon.source.intersect, photon.source.specular_direction, 0.1),
									 start_medium->get_metrics().get_surface_reflection());

	light->next = intersection;      // intersection vertex/node
	intersection->prev = light;      // light source origin
	intersection->emit = reflection; // specular reflection

	// Initialize photon's internal path tracking
	photon.path_head = intersection;
	photon.path_last = intersection;
	photon.num_seg_int = 1;
	photon.num_seg_ext = 1;

	metrics_->add_vertex(photon.position.x, photon.position.y, photon.position.z);
}

void Transport::specular_reflection(Photon& photon) {
	Voxel* voxel = nullptr;

	// Try the more robust DDA-based voxel lookup first
	Medium* dda_medium = simulator_->find_medium_at_with_dda(photon.source.intersect);
	if (dda_medium) {
		glm::dvec3 pos = photon.source.intersect; // Make a non-const copy
		voxel = dda_medium->voxel_at(pos);
	}

	// Fallback to regular voxel_at if DDA didn't work
	if (!voxel) {
		voxel = simulator_->voxel_at(photon.source.intersect);
	}

	// Final fallback: nudge the intersection point slightly into the medium
	if (!voxel) {
		glm::dvec3 nudged_pos = photon.source.intersect + photon.source.direction * 1e-6;
		Medium* nudged_medium = simulator_->find_medium_at_with_dda(nudged_pos);
		if (nudged_medium) {
			voxel = nudged_medium->voxel_at(nudged_pos);
		}
	}

	// voxel should never be nullptr at this point, but handle gracefully if it occurs
	if (!voxel) {
		// Log the error but continue simulation - assume no reflection occurs
		if (Config::is_initialized() && Config::get().log()) {
			ErrorHandler::instance().report_warning(
				"Specular reflection could not be computed - photon will continue unreflected.");
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
	auto* medium = simulator_->find_medium_at(photon.source.intersect);
	if (medium) {
		// Surface refraction is the energy that enters the medium (1 - reflected energy)
		double surface_refraction_energy = photon.weight * (1.0 - fresnel_reflection);
		medium->get_metrics().add_surface_refraction(surface_refraction_energy); // Energy entering medium at surface

		// Record specular reflection (energy immediately reflected at surface)
		double specular_reflection_energy = photon.weight * fresnel_reflection;
		medium->get_metrics().add_surface_reflection(specular_reflection_energy);

		// Add specular reflection to entry voxel for rendering visibility
		voxel->specular_reflection += specular_reflection_energy;

		// Reduce photon weight by reflected amount so only transmitted energy
		// continues for absorption/emission. This prevents double-counting reflected energy.
		photon.weight = surface_refraction_energy; // Only transmitted energy continues
	}

	// reflection direction: R = V - 2(V . N)N
	// With outward-pointing normals, calculate reflection properly
	glm::dvec3 normal = photon.source.triangle.normal();
	glm::dvec3 incident = photon.source.direction;
	double projection_scalar = glm::dot(incident, normal);

	// Standard reflection: R = I - 2(I·N)N with outward normals
	// when the incident ray is pointing toward the surface

	glm::dvec3 reflection_direction = incident - 2.0 * projection_scalar * normal;
	photon.source.specular_direction = glm::normalize(reflection_direction);
}

void Transport::track_photon_path_segments_for_absorption(Photon& photon) {
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
	Medium* current_medium = simulator_->find_medium_at(start_pos);
	if (!current_medium) {
		return;
	}

	// Calculate absorption using the same DDA traversal but on the actual photon path
	// Find the DDA instance for this medium
	size_t medium_index = 0;
	for (size_t i = 0; i < simulator_->get_mediums().size(); ++i) {
		if (&simulator_->get_mediums()[i] == current_medium) {
			medium_index = i;
			break;
		}
	}

	if (medium_index >= simulator_->get_medium_ddas().size()) {
		// DDA not available - this should not happen in production
		std::cerr << "WARNING: DDA not available for medium " << medium_index << std::endl;
		return;
	}

	DDA* dda = simulator_->get_medium_ddas()[medium_index].get();

	// Initialize DDA for this ray
	dda->initialize_ray(start_pos, direction);

	// Traverse voxels using DDA on the actual photon path
	DDA::TraversalResult result = dda->traverse(total_distance);

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
				if (dda_coords.x >= 0 && dda_coords.x < static_cast<int>(Config::get().nx()) && dda_coords.y >= 0
					&& dda_coords.y < static_cast<int>(Config::get().ny()) && dda_coords.z >= 0
					&& dda_coords.z < static_cast<int>(Config::get().nz())) {
					voxel = simulator_->voxel_grid(static_cast<uint32_t>(dda_coords.x),
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

void Transport::track_voxel_path_with_dda(Photon& photon) {
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
	Medium* current_medium = simulator_->find_medium_at(start_pos);
	if (!current_medium) {
		return;
	}

	// Find the DDA instance for this medium
	size_t medium_index = 0;
	for (size_t i = 0; i < simulator_->get_mediums().size(); ++i) {
		if (&simulator_->get_mediums()[i] == current_medium) {
			medium_index = i;
			break;
		}
	}

	if (medium_index >= simulator_->get_medium_ddas().size()) {
		// DDA not available - this should not happen in production
		std::cerr << "WARNING: DDA not available for medium " << medium_index << std::endl;
		return;
	}

	DDA* dda = simulator_->get_medium_ddas()[medium_index].get();

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
}

void Transport::deposit(Photon& photon) {
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
		if (!simulator_->is_point_inside_geometry(photon.position)) {
			// Use the outside volume fraction for photons in the outside portion
			effective_volume_fraction = photon.voxel->volume_fraction_outside;
		}
	}
	else {
		// For non-boundary voxels, check geometry but don't exit early for energy conservation
		if (!simulator_->is_point_inside_geometry(photon.position)) {
			// Photon is slightly outside geometry due to numerical precision
			// Still deposit energy to maintain conservation, but with reduced fraction
			effective_volume_fraction = 0.5; // Compromise value for edge cases
		}
	}

	// deposited weight (scaled by effective volume fraction)
	double deltaw = photon.weight * (1 - std::exp(-photon.mu_a() * photon.sub_step)) * effective_volume_fraction;

	// Track absorption details for debugging (log mode only, limited output)
	static int deposit_debug_count = 0;
	if (Config::get().log() && deposit_debug_count < 10
		&& deltaw > 0.001) { // Limit to first 10 steps with significant absorption
		std::ostringstream oss;
		oss << "DEPOSIT: Photon " << photon.id << " weight=" << photon.weight << " mu_a=" << photon.mu_a()
			<< " sub_step=" << photon.sub_step << " deltaw=" << deltaw
			<< " effective_volume=" << effective_volume_fraction;
		Logger::instance().log_debug(oss.str());
		deposit_debug_count++;
	}

	// ENERGY CONSERVATION ENFORCEMENT
	// Calculate how much energy this photon has left to absorb
	double energy_already_used = photon.total_energy_radiated + photon.total_energy_absorbed;
	double energy_available = photon.total_energy_budget - energy_already_used;

	// Prevent negative energy calculations
	if (energy_available < 0.0) {
		if (Config::get().log()) {
			std::ostringstream oss;
			oss << "Energy available became negative (" << energy_available
				<< "), budget=" << photon.total_energy_budget << ", used=" << energy_already_used
				<< " (radiated=" << photon.total_energy_radiated << ", absorbed=" << photon.total_energy_absorbed
				<< ")";
			Logger::instance().log_warning(oss.str());
		}
		energy_available = 0.0;
	}

	// Enforce energy conservation: cannot absorb more than available
	double actual_absorbed_weight = std::min(deltaw, energy_available);

	// Track absorption calculation
	if (Config::get().log() && (actual_absorbed_weight != deltaw || photon.weight < actual_absorbed_weight)) {
		std::ostringstream oss;
		oss << "ABSORPTION: deltaw=" << deltaw << ", actual=" << actual_absorbed_weight
			<< ", current_weight=" << photon.weight << ", after_weight=" << (photon.weight - actual_absorbed_weight);
		Logger::instance().log_debug(oss.str());
	}

	// Update photon energy tracking
	photon.total_energy_absorbed += actual_absorbed_weight;

	// Update photon weight
	photon.weight -= actual_absorbed_weight;

	// Prevent negative weights
	if (photon.weight < 0.0) {
		if (Config::get().log()) {
			std::ostringstream debug_msg;
			debug_msg << "WARNING: Photon weight became negative (" << photon.weight << "), setting to 0";
			Logger::instance().log_warning(debug_msg.str());
		}
		photon.weight = 0.0;
	}

	// Terminate zero-weight photons immediately
	if (photon.weight <= 0.0) {
		photon.alive = false;
		return;
	}

	// assign deposited weight to voxel
	photon.voxel->absorption += actual_absorbed_weight;
}
