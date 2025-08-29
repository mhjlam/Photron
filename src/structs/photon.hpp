#pragma once

#include "point3.hpp"
#include "source.hpp"
#include "vector3.hpp"
#include "voxel.hpp"

#include <cstdint>

// Emitted photon
struct Emitter {
	uint64_t id;
	Point3 position;
	Vector3 direction;
	double weight;

	Emitter(uint64_t i, Point3& p, Vector3& d, double w) {
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

	Point3 position;      // position at the start of a substep
	Vector3 direction;    // propagation direction
	Voxel* voxel;         // resident voxel at start of substep
	Voxel* prev_voxel;    // voxel before crossing an interface

	Point3 intersect;     // voxel boundary intersection
	Source source;        // source from which the photon was shot
	Vector3 voxel_normal; // voxel boundary intersection normal

	bool scatters;        // true if photon path scatters at least once

	Photon() {
		id = 0;
		alive = true;
		cross = false;

		step = 0;
		sub_step = 0;
		weight = 0;

		position = Point3();
		direction = Vector3();
		voxel = nullptr;
		prev_voxel = nullptr;

		intersect = Point3();
		voxel_normal = Vector3();
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

		position = Point3();
		direction = Vector3();
		voxel = nullptr;
		prev_voxel = nullptr;

		intersect = Point3();
		voxel_normal = Vector3();
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
