#pragma once

#include "glm_types.hpp"
#include "source.hpp"
#include "voxel.hpp"

#include <cstdint>

// Emitted photon
struct Emitter {
	uint64_t id;
	glm::dvec3 position;
	glm::dvec3 direction;
	double weight;

	Emitter(uint64_t i, glm::dvec3& p, glm::dvec3& d, double w) {
		id = i;
		position = p;
		direction = d;
		weight = w;
	}
};

struct Photon {
	uint64_t id;          // identifier
	bool alive;           // true if photon still propagates
	bool cross;           // true if photon crosses into another voxel in substep

	double step;          // step distance between scattering
	double sub_step;      // step distance inside a voxel
	double weight;        // remaining weight of the photon packet

	glm::dvec3 position;      // position at the start of a substep
	glm::dvec3 direction;     // propagation direction
	Voxel* voxel;         // resident voxel at start of substep
	Voxel* prev_voxel;    // voxel before crossing an interface

	glm::dvec3 intersect;     // voxel boundary intersection
	Source source;        // source from which the photon was shot
	glm::dvec3 voxel_normal; // voxel boundary intersection normal

	bool scatters;        // true if photon path scatters at least once

	Photon() {
		id = 0;
		alive = true;
		cross = false;

		step = 0;
		sub_step = 0;
		weight = 0;

		position = glm::dvec3(0);
		direction = glm::dvec3(0);
		voxel = nullptr;
		prev_voxel = nullptr;

		intersect = glm::dvec3(0);
		voxel_normal = glm::dvec3(0);
		source = Source();

		scatters = false;
	}

	Photon(uint64_t i) {
		id = i;
		alive = true;
		cross = false;

		step = 0;
		sub_step = 0;
		weight = 0;

		position = glm::dvec3(0);
		direction = glm::dvec3(0);
		voxel = nullptr;
		prev_voxel = nullptr;

		intersect = glm::dvec3(0);
		voxel_normal = glm::dvec3(0);
		source = Source();

		scatters = false;
	}

	double ani() {
		if (voxel == nullptr || voxel->tissue == nullptr) {
			return 0;
		}
		return voxel->tissue->ani;
	}

	double eta() {
		if (voxel == nullptr || voxel->tissue == nullptr) {
			return 0;
		}
		return voxel->tissue->eta;
	}

	double mua() {
		if (voxel == nullptr || voxel->tissue == nullptr) {
			return 0;
		}
		return voxel->tissue->mua;
	}

	double mus() {
		if (voxel == nullptr || voxel->tissue == nullptr) {
			return 0;
		}
		return voxel->tissue->mus;
	}
};
