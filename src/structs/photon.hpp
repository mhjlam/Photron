#pragma once

#include "point3.hpp"
#include "source.hpp"
#include "vector3.hpp"
#include "voxel.hpp"

struct Emitter // emitted photon
{
	unsigned long id;
	Point3 pos;
	Vector3 dir;
	double weight;

	Emitter(unsigned long i, Point3& p, Vector3& d, double w) {
		id = i;
		pos = p;
		dir = d;
		weight = w;
	}
};

struct Photon {
	unsigned long id; // identifier
	bool alive;       // true if photon still propagates
	bool cross;       // true if photon crosses into another voxel in substep

	double step;      // step distance between scattering
	double substep;   // step distance inside a voxel
	double weight;    // remaining weight of the photon packet

	Point3 pos;       // position at the start of a substep
	Vector3 dir;      // propagation direction
	Voxel* vox;       // resident voxel at start of substep
	Voxel* prevvox;   // voxel before crossing an interface

	Point3 inter;     // voxel boundary intersection
	Source source;    // source from which the photon was shot
	Vector3 vnorm;    // voxel boundary intersection normal

	bool scatters;    // true if photon path scatters at least once

	Photon() {
		id = 0;
		alive = true;
		cross = false;

		step = 0;
		substep = 0;
		weight = 0;

		pos = Point3();
		dir = Vector3();
		vox = NULL;
		prevvox = NULL;

		inter = Point3();
		vnorm = Vector3();
		source = Source();

		scatters = false;
	}

	Photon(unsigned long i) {
		id = i;
		alive = true;
		cross = false;

		step = 0;
		substep = 0;
		weight = 0;

		pos = Point3();
		dir = Vector3();
		vox = NULL;
		prevvox = NULL;

		inter = Point3();
		vnorm = Vector3();
		source = Source();

		scatters = false;
	}

	double ani() {
		if (vox == NULL || vox->tissue == NULL)
			return 0;

		return vox->tissue->ani;
	}

	double eta() {
		if (vox == NULL || vox->tissue == NULL)
			return 0;

		return vox->tissue->eta;
	}

	double mua() {
		if (vox == NULL || vox->tissue == NULL)
			return 0;

		return vox->tissue->mua;
	}

	double mus() {
		if (vox == NULL || vox->tissue == NULL)
			return 0;

		return vox->tissue->mus;
	}
};
