#pragma once

#include <cstdint>
#include "math/glm_types.hpp"
#include "math/concepts.hpp"
#include "simulator/source.hpp"
#include "simulator/voxel.hpp"

// Emitted photon
struct Emitter
{
	uint64_t id;
	glm::dvec3 position;
	glm::dvec3 direction;
	double weight;

	template<Point3D P, Vector3D V>
	explicit constexpr Emitter(uint64_t i, const P& p, const V& d, double w) noexcept
		: id(i), position(static_cast<glm::dvec3>(p)), direction(static_cast<glm::dvec3>(d)), weight(w) {
	}
};

struct Photon
{
	uint64_t id = 0;                     // identifier
	bool alive = true;                   // true if photon still propagates
	bool cross = false;                  // true if photon crosses into another voxel in substep

	double step = 0.0;                   // step distance between scattering
	double sub_step = 0.0;               // step distance inside a voxel
	double weight = 0.0;                 // remaining weight of the photon packet

	glm::dvec3 position{0.0};            // position at the start of a substep
	glm::dvec3 direction{0.0};           // propagation direction
	Voxel* voxel = nullptr;              // resident voxel at start of substep (non-owning)
	Voxel* prev_voxel = nullptr;         // voxel before crossing an interface (non-owning)

	glm::dvec3 intersect{0.0};           // voxel boundary intersection
	Source source{};                     // source from which the photon was shot
	glm::dvec3 voxel_normal{0.0};        // voxel boundary intersection normal

	bool scatters = false;               // true if photon path scatters at least once

	// Default constructor uses default member initialization
	Photon() = default;

	// Constructor with ID (other members use default initialization)
	explicit Photon(uint64_t i) noexcept : id(i) {}

	// Tissue property accessors with null safety
	[[nodiscard]] double g() const noexcept {
		return (voxel && voxel->tissue) ? voxel->tissue->g : 0.0;
	}

	[[nodiscard]] double eta() const noexcept {
		return (voxel && voxel->tissue) ? voxel->tissue->eta : 0.0;
	}

	[[nodiscard]] double mu_a() const noexcept {
		return (voxel && voxel->tissue) ? voxel->tissue->mu_a : 0.0;
	}

	[[nodiscard]] double mu_s() const noexcept {
		return (voxel && voxel->tissue) ? voxel->tissue->mu_s : 0.0;
	}
};
