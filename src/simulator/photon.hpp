#pragma once

#include <vector>

#include <glm/glm.hpp>

#include "math/concepts.hpp"
#include "math/triangle.hpp"
#include "simulator/voxel.hpp"

// Emitted photon
struct Emitter
{
	uint64_t id;
	glm::dvec3 position;
	glm::dvec3 direction;
	double weight;
	
	// Exit classification for accurate energy labeling
	enum class ExitType { NONE, REFLECTED, TRANSMITTED };
	ExitType exit_type {ExitType::NONE};

	template<Point3D P, Vector3D V>
	explicit constexpr Emitter(uint64_t i, const P& p, const V& d, double w, ExitType exit_classification = ExitType::NONE) noexcept :
		id(i), position(static_cast<glm::dvec3>(p)), direction(static_cast<glm::dvec3>(d)), weight(w), exit_type(exit_classification) {}
};

struct Photon
{
	uint64_t id {0};               		// identifier
	bool alive {true};             		// true if photon still propagates
	bool cross {false};            		// true if photon crosses voxel boundary in substep
	bool exits {false};		   			// true if photon crosses medium boundary in substep

	// ENERGY CONSERVATION TRACKING
	double total_energy_budget {0.0};	// total energy available to this photon (includes Russian Roulette amplifications)
	double total_energy_radiated {0.0};	// cumulative energy radiated through all radiate() calls
	double total_energy_absorbed {0.0};	// cumulative energy absorbed through all absorption events
	int radiate_call_count {0};		// number of times radiate() was called for this photon

	double step {0.0};             		// step distance between scattering points
	double sub_step {0.0};         		// step distance inside a voxel
	double weight {0.0};           		// remaining weight of the photon packet

	glm::dvec3 position {0.0};     		// position at the start of a substep
	glm::dvec3 direction {0.0};    		// propagation direction
	Voxel* voxel {nullptr};        		// resident voxel at start of substep (non-owning)
	Voxel* prev_voxel {nullptr};   		// voxel before crossing an interface (non-owning)

	glm::dvec3 intersect {0.0};    		// voxel boundary intersection
	glm::dvec3 voxel_normal {0.0}; 		// voxel boundary intersection normal

	bool scatters {false};         		// true if photon path scatters at least once

	// Exit classification (set when photon exits the medium)
	enum class ExitType { NONE, REFLECTED, TRANSMITTED };
	ExitType exit_type {ExitType::NONE};	// how the photon exited the medium

	// Source properties merged into photon
	glm::dvec3 source_origin {0.0};             // origin of the source
	glm::dvec3 source_direction {0.0};          // initial direction of incidence from source
	glm::dvec3 specular_direction {0.0};        // direction of specular reflectance
	glm::dvec3 source_intersect {0.0};          // intersection point with geometry
	Triangle source_triangle;                   // triangle at intersection point

	// Medium context is now queried dynamically via Simulator::find_medium_at()


	// Default constructor uses default member initialization
	Photon() = default;

	// Constructor with ID (other members use default initialization)
	explicit Photon(uint64_t i) noexcept : id(i) {}

	// Constructor with source data merged in
	explicit Photon(uint64_t i, const glm::dvec3& src_origin, const glm::dvec3& src_direction) noexcept :
		id(i), source_origin(src_origin), source_direction(src_direction) {}

	// Constructor with full source data
	explicit Photon(uint64_t i, const glm::dvec3& src_origin, const glm::dvec3& src_direction,
					const glm::dvec3& spec_direction, const glm::dvec3& src_intersect, 
					const Triangle& src_triangle) noexcept :
		id(i), source_origin(src_origin), source_direction(src_direction),
		specular_direction(spec_direction), source_intersect(src_intersect), source_triangle(src_triangle) {}

	// Tissue property accessors with null safety (for backward compatibility)
	[[nodiscard]] double g() const noexcept { return (voxel && voxel->tissue) ? voxel->tissue->g : 0.0; }

	[[nodiscard]] double eta() const noexcept { return (voxel && voxel->tissue) ? voxel->tissue->eta : 0.0; }

	[[nodiscard]] double mu_a() const noexcept { return (voxel && voxel->tissue) ? voxel->tissue->mu_a : 0.0; }

	[[nodiscard]] double mu_s() const noexcept { return (voxel && voxel->tissue) ? voxel->tissue->mu_s : 0.0; }
};
