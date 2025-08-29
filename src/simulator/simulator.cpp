#include "simulator.hpp"
#include "random.hpp"

#include "../utilities/experimenter.hpp"
#include "../utilities/utilities.hpp"

#include <cmath>
#include <cfloat>
#include <cstdlib>
#include <ctime>
#include <fstream>
#include <iostream>
#include <memory>

typedef unsigned int uint;
typedef unsigned long ulong;
typedef unsigned short ushort;

using namespace std;

/***********************************************************
 * Run the Monte Carlo photon transport simulation.
 ***********************************************************/
void Simulator::Simulate() {
	cout << "Running Monte Carlo simulation" << endl;
	
	Experimenter::StartClock();

	// for each light source
	for (ushort s = 0; s < sources.size(); ++s) {
		// compute specular reflection
		SpecularReflection(sources[s]);

		// for each photon
		for (ulong p = 0; p < photons.size(); ++p) {
			// progress report
			if (config.progress && ((p + 1) % 1000) == 0) {
				cout << "Photon " << p + 1 << "/" << config.numphotons << endl;
			}

			// launch the photon (create a new path)
			Launch(photons[p], sources[s]);

			while (photons[p].alive) {
				// set new step size
				StepSize(photons[p]);

				// propagate photon through the medium in substeps
				Transfer(photons[p]);

				// determine photon termination
				Roulette(photons[p]);

				// scatter photon into a new direction
				Scatter(photons[p]);
			}
		}
	}

	// normalize physical quantities
	Normalize();

	Experimenter::StopClock();
	Experimenter::CollectData(record.at, record.rs, record.rd, record.ts, record.td);
	Experimenter::WriteToFile();
	Experimenter::Print();
}

/***********************************************************
 * Set up photon properties for tracing.
 ***********************************************************/
void Simulator::Launch(Photon& photon, Source& source) {
	photon.alive = true;
	photon.weight = 1.0 - record.rs;

	photon.source = source;

	photon.dir = source.dir;
	photon.pos = source.inter;
	photon.vox = VoxelAt(photon.pos);

	// create vertices for new light path
	Vertex* light = new Vertex(source.orig, photon.weight);
	Vertex* isect = new Vertex(source.inter, photon.weight);
	Vertex* refls = new Vertex(Move(source.inter, source.dirspec, 0.1), record.rs);

	light->next = isect; // intersection vertex/node
	isect->prev = light; // light source origin
	isect->emit = refls; // specular reflection

	// use intersection point as head
	Graph path = Graph(static_cast<long>(photon.id), isect);
	paths.push_back(path);

	Experimenter::AddVertex(photon.pos.x, photon.pos.y, photon.pos.z);
}

/***********************************************************
 * Set dimensionless step size for next segment of the
 * random walk using MCML 3.0.0 algorithm.
 ***********************************************************/
void Simulator::StepSize(Photon& photon) {
	// Use MCML 3.0.0 step size generation
	GenerateStepSize(photon);
	
	Experimenter::AddStepSize(photon.step / photon.mus());
}

/***********************************************************
 * Transfer a photon through a voxelized medium with
 * individual substeps.
 ***********************************************************/
void Simulator::Transfer(Photon& photon) {
	/*
	 * set substep (max = distance to voxel boundary)
	 * deposit weight in current voxel
	 * if (photon crosses voxel boundary)
	 *     if (photon goes outside medium)
	 *         record partial transmission
	 *     else if (photon moves to differing refractive indexed media)
	 *         reflect from or transmit across the voxel boundary
	 *     else if (photon moves to equal refractive indexed media)
	 *         continue normal propagation
	 * decrease step size by traveled distance
	 */

	while (photon.step >= 1E-10 && photon.alive) {
		// set substep
		Substep(photon);

		// deposit weight
		Deposit(photon);

		// possibly cross voxel boundary
		if (photon.cross) {
			Cross(photon);
		}
		else {
			photon.pos = Move(photon.pos, photon.dir, photon.substep);
		}

		// prevent errors due to crossing to ambient medium
		if (!photon.vox) {
			photon.alive = false;
			photon.vox = photon.prevvox; // for Radiate
			Radiate(photon, photon.dir, photon.weight);
		}

		// update step size
		photon.step -= (photon.substep * photon.mus());

		Experimenter::AddVertex(photon.pos.x, photon.pos.y, photon.pos.z);
	}
}

/***********************************************************
 * Set the photon's next substep and initialize the given
 * intersection point and voxel normal if it crosses the
 * voxel boundary.
 ***********************************************************/
void Simulator::Substep(Photon& photon) {
	// create ray and get voxel vertices
	Ray ray = Ray(photon.pos, photon.dir);
	Cuboid box = VoxelCorners(photon.vox);

	// find first intersection with voxel faces and get the distance
	double voxdist = first_ray_cuboid_intersect_internal(ray, box, photon.inter, photon.vnorm);

	// check if no intersections were found
	if (voxdist == DBL_MAX) {
		cerr << "Critical error: ray-voxel intersection test failed." << endl;
		exit(EXIT_FAILURE);
	}

	// compute free path for a substep (distance to next scattering event)
	double freepath = photon.step / photon.mus();

	// see if the free path crosses the voxel
	if (voxdist == 0) {
		// already on voxel face; move just beyond
		photon.substep = config.voxsize * 0.001;
		photon.cross = true;
	}
	// crosses the voxel boundary
	else if (freepath > voxdist) {
		photon.substep = voxdist;
		photon.cross = true;
	}
	// does not cross
	else {
		photon.substep = freepath;
		photon.cross = false;
	}
}

/***********************************************************
 * Deposit some of the photon's weight into the geometry.
 ***********************************************************/
void Simulator::Deposit(Photon& photon) {
	// cancel if photon is outside of medium
	if (!photon.vox->tissue) {
		return;
	}

	// deposited weight
	double deltaw = photon.weight * (1 - std::exp(-photon.mua() * photon.substep));

	// update photon weight
	photon.weight -= deltaw;

	// assign deposited weight to voxel
	photon.vox->absorption += deltaw;

	// update total absorption
	record.at += deltaw;
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
void Simulator::Cross(Photon& photon) {
	// directions of transmission and reflection
	Vector3 tran, refl;

	// retrieve the refractive index of the medium that is struck
	Point3 newpos = MoveDelta(photon.inter, photon.dir);
	Voxel* newvox = VoxelAt(newpos);
	photon.prevvox = photon.vox;

	// determine refractive index
	double eta = (newvox == NULL) ? config.ambienteta : newvox->tissue->eta;

	// 1. crossing to ambient medium
	if (newvox == NULL) {
		// compute reflected/transmitted fraction and reflection/transmission directions
		double reflection = InternalReflection(photon, eta, tran, refl);
		double transmission = 1.0 - reflection;

		// total transmission
		if (reflection == 0.0) {
			// photon dies
			photon.dir = tran;
			photon.alive = false;

			Radiate(photon, tran, photon.weight);
		}
		// total internal reflection
		else if (reflection == 1.0) {
			// photon reflects off surface
			photon.dir = refl;
			photon.pos = MoveDelta(photon.inter, photon.dir);
			photon.vox = VoxelAt(photon.pos);
		}
		else {
			// partial reflection
			if (config.partial) {
				// emit partial reflectance
				Radiate(photon, tran, photon.weight * transmission);

				// adjust photon weight
				photon.weight *= reflection;

				// update direction/position/voxel
				photon.dir = refl;
				photon.pos = MoveDelta(photon.inter, photon.dir);
				photon.vox = VoxelAt(photon.pos);
			}
			else // all-or-none transmission/reflection
			{
				// total transmission
				if (randomnum() > reflection) {
					photon.dir = tran;
					photon.alive = false;

					Radiate(photon, tran, photon.weight);
				}
				// total reflection
				else {
					photon.dir = refl;
					photon.pos = MoveDelta(photon.inter, photon.dir);
					photon.vox = VoxelAt(photon.pos);
				}
			}
		}

		Experimenter::IncrementScatters();
	}
	// 2. crossing to another medium
	else if (newvox != NULL && photon.vox->tissue->eta != newvox->tissue->eta) {
		// compute reflected/transmitted fraction and reflection/transmission directions
		double reflection = InternalReflection(photon, eta, tran, refl);

		// total transmission
		if (reflection == 0.0) {
			photon.dir = tran;
			photon.pos = MoveDelta(photon.inter, photon.dir);
			photon.vox = VoxelAt(photon.pos);
		}
		// total internal reflection
		else if (reflection == 1.0) {
			photon.dir = refl;
			photon.pos = MoveDelta(photon.inter, photon.dir);
			photon.vox = VoxelAt(photon.pos);
		}
		else // all-or-none transmission/reflection
		{
			// total transmission
			if (randomnum() > reflection) {
				photon.dir = tran;
				photon.pos = MoveDelta(photon.inter, photon.dir);
				photon.vox = VoxelAt(photon.pos);
			}
			// total reflection
			else {
				photon.dir = refl;
				photon.pos = MoveDelta(photon.inter, photon.dir);
				photon.vox = VoxelAt(photon.pos);
			}
		}

		Experimenter::IncrementScatters();
	}
	// 3. crossing within the same medium (total transmission)
	else {
		// direction is unchanged
		photon.pos = MoveDelta(photon.inter, photon.dir);
		photon.vox = VoxelAt(photon.pos);
	}
}

/***********************************************************
 * Record the emittance from a photon (partially) leaving
 * the material.
 ***********************************************************/
void Simulator::Radiate(Photon& photon, Vector3& dir, double weight) {
	// set voxel emittance
	photon.vox->emittance += weight;

	// add regular vertex that stops at boundary
	paths.back().AddInternalVertex(new Vertex(photon.inter, weight));
	paths.back().AddExternalVertex(new Vertex(Move(photon.inter, dir, 0.1), weight));
	Experimenter::AddVertex(photon.inter.x, photon.inter.y, photon.inter.z);

	// add emitter at this position
	emitters.push_back(Emitter(photon.id, photon.inter, dir, weight));

	// specular or diffuse transmission
	if (!photon.scatters) {
		record.ts += weight;
	}
	else {
		// cos(theta) = a.b / |a||b|
		double costheta = dot(photon.source.dir, dir) / (photon.source.dir.length() * dir.length());
		double theta = acos(costheta);

		// count as transmission or reflection
		if (theta < HALFPI) {
			record.td += weight;
		}
		else {
			record.rd += weight;
		}
	}
}

/***********************************************************
 * Determine the survivability of a given photon using MCML 3.0.0 algorithm.
 ***********************************************************/
void Simulator::Roulette(Photon& photon) {
	// Use MCML 3.0.0 Russian roulette
	RoulettePhoton(photon);
}

/***********************************************************
 * Scatter the photon into a new direction based on the
 * Henyey-Greenstein phase function using MCML 3.0.0 algorithm.
 ***********************************************************/
void Simulator::Scatter(Photon& photon) {
	if (!photon.alive) {
		return;
	}

	// Get tissue properties for scattering
	Tissue* tissue = photon.vox->tissue;
	if (!tissue) {
		photon.alive = false;
		return;
	}

	// Use MCML 3.0.0 scattering algorithm
	ScatterPhoton(photon, *tissue);

	// normalize direction vector (safety check)
	photon.dir.normalize();

	// prevent scattering into ambient medium when close to boundaries
	Point3 newpos = MoveDelta(photon.pos, photon.dir);
	if (!VoxelAt(newpos)) {
		photon.alive = false;
		Radiate(photon, photon.dir, photon.weight);
	}

	Experimenter::IncrementScatters();

	// add new internal position to path
	paths.back().AddInternalVertex(new Vertex(photon.pos, photon.weight));
}

/***********************************************************
 * Normalize the recorded values based on the number of
 * photons traced.
 ***********************************************************/
void Simulator::Normalize() {
	// normalize globally recorded parameters
	record.at /= config.numphotons * config.numsources;
	record.rd /= config.numphotons * config.numsources;
	record.td /= config.numphotons * config.numsources;

	record.rs /= config.numsources;                     // rs is only computed once per source
	record.ts /= config.numphotons * config.numsources; // ts is computed per photon

	// normalize voxel data
	for (ulong i = 0; i < voxels.size(); ++i) {
		// skip computation for voxels outside the medium
		if (!voxels[i]->tissue) {
			continue;
		}

		voxels[i]->absorption /= (config.numphotons * config.numsources);
		voxels[i]->emittance /= (config.numphotons * config.numsources);
	}
}

/***********************************************************
 * Compute the specular reflectance from a light source at
 * the surface.
 ***********************************************************/
void Simulator::SpecularReflection(Source& source) {
	Voxel* voxel = VoxelAt(source.inter);

	// voxel should never be NULL at this point
	if (!voxel) {
		cerr << "Critical error: specular reflection could not be computed." << endl;
		exit(EXIT_FAILURE);
	}

	// refractive indices of ambient medium and medium that is hit
	double n1 = config.ambienteta;
	double n2 = voxel->tissue->eta;

	// set the specular reflection
	record.rs = (n2 != n1) ? sq((n1 - n2) / (n1 + n2)) : 0;

	// reflection dir: R = -2(V . N)N + V
	Vector3 normal = source.intertri.normal;
	Vector3 projection = normal * dot(source.dir, normal);
	Vector3 rsdir = Vector3((projection * -2.0) + source.dir, true);

	source.dirspec = rsdir;
}

/***********************************************************
 * Compute the fraction of incoming light that is reflected
 * back at an interface between two media. Also compute the
 * directions of transmission and reflection.
 ***********************************************************/
double Simulator::InternalReflection(Photon& photon, double& nt, Vector3& tran, Vector3& refl) {
	// fraction reflected
	double reflection;

	// indices of refraction
	double etai = photon.vox->tissue->eta;
	double etat = nt;

	// critical angle
	double cosc = sqrt(1 - sq(etat) / sq(etai));

	// angles of reflection
	double cost = dot(photon.dir, photon.vnorm);
	double cosp = cost;

	// compute the fraction of reflected light

	// matched boundary
	if (etai == etat) {
		cosp = cost;
		reflection = 0.0;
	}
	// near-perpendicular incidence: cos(0)
	else if (cost >= 1.0 - 1E-10) {
		reflection = sq(etat - etai) / sq(etat + etai);
		cosp = cost;
	}
	// near-parallel incidence: cos(90)
	else if (cost < 1E-10) {
		reflection = 1.0;
		cosp = 0;
	}
	// total internal reflection
	else if ((etat < etai) && cost < cosc) {
		reflection = 1.0;
		cosp = 0;
	}
	// general case
	else {
		// angle phi of reflection
		cosp = sqrt(1 - sq(etai) / sq(etat) * (1 - sq(cost)));

		double temp1 = etai * cosp;
		double temp2 = etat * cost;
		double temp3 = etat * cosp;
		double temp4 = etai * cost;

		// fraction reflected
		reflection = 0.5 * (sq(temp1 - temp2) / sq(temp1 + temp2) + sq(temp3 - temp4) / sq(temp3 + temp4));
	}

	// direction of transmission
	double temp1 = etai / etat;
	double temp2 = cosp - cost * temp1;
	tran.x = photon.dir.x * temp1 + temp2 * photon.vnorm.x;
	tran.y = photon.dir.y * temp1 + temp2 * photon.vnorm.y;
	tran.z = photon.dir.z * temp1 + temp2 * photon.vnorm.z;

	// direction of reflection
	refl.x = photon.dir.x - 2 * cost * photon.vnorm.x;
	refl.y = photon.dir.y - 2 * cost * photon.vnorm.y;
	refl.z = photon.dir.z - 2 * cost * photon.vnorm.z;

	tran.normalize();
	refl.normalize();

	return reflection;
}

/***********************************************************
 * Return a pointer to a voxel that encapsulates the given
 * position, or NULL if the position is outside the medium.
 ***********************************************************/
Voxel* Simulator::VoxelAt(Point3& pos) {
	if (!bounds.includes(pos.x, pos.y, pos.z)) {
		return NULL;
	}

	// distances from boundaries
	double dx = fabs(bounds.xmin - pos.x);
	double dy = fabs(bounds.ymin - pos.y);
	double dz = fabs(bounds.zmin - pos.z);

	// indices start at minimum boundaries of voxels
	ushort ix = static_cast<ushort>(floorf(static_cast<float>(dx / config.voxsize)));
	ushort iy = static_cast<ushort>(floorf(static_cast<float>(dy / config.voxsize)));
	ushort iz = static_cast<ushort>(floorf(static_cast<float>(dz / config.voxsize)));

	// avoid index overflow
	if (ix >= config.nx) {
		ix = config.nx - 1;
	}
	if (iy >= config.ny) {
		iy = config.ny - 1;
	}
	if (iz >= config.nz) {
		iz = config.nz - 1;
	}

	// retrieve the voxel at the position
	size_t voxel_index = static_cast<size_t>(iz) * static_cast<size_t>(config.nx) * static_cast<size_t>(config.ny) + 
	                     static_cast<size_t>(iy) * static_cast<size_t>(config.nx) + static_cast<size_t>(ix);
	Voxel* voxel = voxels.at(voxel_index);

	// if voxel does not have a tissue, it is outside the medium
	return (voxel->tissue) ? voxel : NULL;
}

/***********************************************************
 * Return the center position of the given voxel (unused).
 ***********************************************************/
Point3 Simulator::VoxelCenter(Voxel* vox) {
	double vsz = config.voxsize;
	double hvs = vsz / 2.0;

	// infer from indices
	double x = bounds.xmin + (vsz * vox->ix) + hvs;
	double y = bounds.ymin + (vsz * vox->iy) + hvs;
	double z = bounds.zmin + (vsz * vox->iz) + hvs;

	return Point3(x, y, z);
}

/***********************************************************
 * Return the minimum and maximum positions of the given
 * voxel as a cuboid structure.
 ***********************************************************/
Cuboid Simulator::VoxelCorners(Voxel* vox) {
	uint ixmin = vox->ix;
	uint iymin = vox->iy;
	uint izmin = vox->iz;

	uint ixmax = vox->ix + 1;
	uint iymax = vox->iy + 1;
	uint izmax = vox->iz + 1;

	// minimum voxel position
	float xmin = static_cast<float>(bounds.xmin + (config.voxsize * ixmin));
	float ymin = static_cast<float>(bounds.ymin + (config.voxsize * iymin));
	float zmin = static_cast<float>(bounds.zmin + (config.voxsize * izmin));

	// maximum voxel position
	float xmax = static_cast<float>(bounds.xmin + (config.voxsize * ixmax));
	float ymax = static_cast<float>(bounds.ymin + (config.voxsize * iymax));
	float zmax = static_cast<float>(bounds.zmin + (config.voxsize * izmax));

	// round off coordinate values around the origin
	xmin = (fabsf(xmin) < 1E-10) ? 0 : xmin;
	ymin = (fabsf(ymin) < 1E-10) ? 0 : ymin;
	zmin = (fabsf(zmin) < 1E-10) ? 0 : zmin;

	xmax = (fabsf(xmax) < 1E-10) ? 0 : xmax;
	ymax = (fabsf(ymax) < 1E-10) ? 0 : ymax;
	zmax = (fabsf(zmax) < 1E-10) ? 0 : zmax;

	return Cuboid(xmin, ymin, zmin, xmax, ymax, zmax);
}

/***********************************************************
 * Return the destination for a given origin, direction and
 * distance.
 ***********************************************************/
Point3 Simulator::Move(Point3& pos, Vector3& dir, double d) {
	Point3 point = pos;

	// return the end point of a sub-step
	point.x = pos.x + dir.x * d;
	point.y = pos.y + dir.y * d;
	point.z = pos.z + dir.z * d;

	return point;
}

/***********************************************************
 * Return the destination for a given origin and direction
 * after making a small hop.
 ***********************************************************/
Point3 Simulator::MoveDelta(Point3& pos, Vector3& dir) {
	Point3 point = pos;

	// delta distance (based on voxel size)
	double d = config.voxsize * 0.00001;

	point.x = pos.x + dir.x * d;
	point.y = pos.y + dir.y * d;
	point.z = pos.z + dir.z * d;

	return point;
}

/***********************************************************
 *	Write the resulting physical quantities to a file.
 ***********************************************************/
void Simulator::Report() {
	string strsim = "simulation.out";
	string strabs = "absorption.out";
	string stremi = "emittance.out";
	string strptd = "photons.out";

	ofstream ofsrep(strsim.c_str(), ios_base::out); // simulation report
	ofstream ofsabs(strabs.c_str(), ios_base::out); // absorption report
	ofstream ofsemi(stremi.c_str(), ios_base::out); // emittance report
	ofstream ofsptn(strptd.c_str(), ios_base::out); // photon exitance report

	if (!ofsrep.good() || !ofsabs.good() || !ofsemi.good() || !ofsptn.good()) {
		cerr << "Error: an output file could not be opened." << endl;
		return;
	}

	if (ofsrep.good()) {
		ofsrep.precision(8);
		ofsrep << "################################################################" << endl;
		ofsrep << "# SIMULATION REPORT" << endl;
		ofsrep << "################################################################" << endl;
		ofsrep << endl;

		// write input configuration
		ofsrep << "Configuration" << endl;
		ofsrep << "################################################################" << endl;
		ofsrep << endl;
		ofsrep << "Number of photons: " << tab << config.numphotons << endl;
		ofsrep << "Number of layers:  " << tab << config.numlayers << endl;
		ofsrep << "Number of voxels:  " << tab << config.numvoxels << endl;
		ofsrep << "Grid dimensions:   " << tab << config.nx << " x " << config.ny << " x " << config.nz << endl;
		ofsrep << "Voxel dimensions:  " << tab << config.voxsize << " x " << config.voxsize << " x " << config.voxsize
			   << endl;
		ofsrep << endl << endl;

		// write recorded parameters: a, rs, rd, (ts, td)
		ofsrep << "Recorded parameters" << endl;
		ofsrep << "################################################################" << endl;
		ofsrep << "Total absorption:      " << tab << fixed << record.at << endl;
		ofsrep << "Diffuse reflection:    " << tab << fixed << record.rd << endl;
		ofsrep << "Specular reflection:   " << tab << fixed << record.rs << endl;
		ofsrep << "Diffuse transmission:  " << tab << fixed << record.td << endl;
		ofsrep << "Specular transmission: " << tab << fixed << record.ts << endl;
		ofsrep << endl;
		ofsrep.close();
	}
	else {
		cerr << "Error: file " << strsim << " could not be opened." << endl;
		ofsrep.close();
	}

	// write voxel absorption (each 'block' is a slice)
	if (ofsabs.good()) {
		ofsabs.precision(5);
		ofsabs << "################################################################" << endl;
		ofsabs << "# ABSORPTION REPORT" << endl;
		ofsabs << "################################################################" << endl;
		ofsabs << endl;

		// from rear to front
		for (ushort iz = 0; iz < config.nz; ++iz) {
			ofsabs << "Slice " << iz + 1 << "/" << config.nz;
			if (iz == 0) {
				ofsabs << " (rear)";
			}
			if (iz == config.nz - 1) {
				ofsabs << " (front)";
			}
			ofsabs << endl;

			// top to bottom
			for (ushort iy = config.ny - 1; iy > 0; --iy) {
				// left to right
				for (ushort ix = 0; ix < config.nx; ++ix) {
					size_t voxel_index = static_cast<size_t>(iz) * static_cast<size_t>(config.nx) * static_cast<size_t>(config.ny) + 
					                     static_cast<size_t>(iy) * static_cast<size_t>(config.nx) + static_cast<size_t>(ix);
					Voxel* voxel = voxels.at(voxel_index);
					ofsabs << fixed << voxel->absorption << tab;
				}

				ofsabs << endl;
			}

			ofsabs << endl;
		}

		ofsabs.close();
	}
	else {
		cerr << "Error: file " << strabs << " could not be opened." << endl;
		ofsabs.close();
	}

	// write voxel emittance
	if (ofsemi.good()) {
		ofsemi.precision(5);
		ofsemi << "################################################################" << endl;
		ofsemi << "# EMITTANCE REPORT" << endl;
		ofsemi << "################################################################" << endl;
		ofsemi << endl;

		// rear to front
		for (ushort iz = 0; iz < config.nz; ++iz) {
			ofsemi << "Slice " << iz + 1 << "/" << config.nz;
			if (iz == 0) {
				ofsemi << " (rear)";
			}
			if (iz == config.nz - 1) {
				ofsemi << " (front)";
			}
			ofsemi << endl;

			// top to bottom
			for (ushort iy = config.ny - 1; iy > 0; --iy) {
				// left to right
				for (ushort ix = 0; ix < config.nx; ++ix) {
					size_t voxel_index = static_cast<size_t>(iz) * static_cast<size_t>(config.nx) * static_cast<size_t>(config.ny) + 
					                     static_cast<size_t>(iy) * static_cast<size_t>(config.nx) + static_cast<size_t>(ix);
					Voxel* voxel = voxels.at(voxel_index);
					ofsemi << fixed << voxel->emittance << tab;
				}

				ofsemi << endl;
			}

			ofsemi << endl;
		}

		ofsemi.close();
	}
	else {
		cerr << "Error: file " << stremi << " could not be opened." << endl;
		ofsemi.close();
	}

	// exiting photons (position, dir, weight)
	if (ofsptn.good()) {
		ofsptn.precision(8);
		ofsptn << "################################################################" << endl;
		ofsptn << "# PHOTON REPORT" << endl;
		ofsptn << "################################################################" << endl;
		ofsptn << endl;

		if (emitters.empty()) {
			ofsptn << "No photons exited the medium." << endl;
		}

		for (vector<Emitter>::iterator ep = emitters.begin(); ep != emitters.end(); ++ep) {
			ofsptn << "Photon" << endl;
			ofsptn << "{" << endl;
			ofsptn << fixed << tab << "id  = " << ep->id << endl;
			ofsptn << fixed << tab << "pos = " << ep->pos.x << ", " << ep->pos.y << ", " << ep->pos.z << endl;
			ofsptn << fixed << tab << "dir = " << ep->dir.x << ", " << ep->dir.y << ", " << ep->dir.z << endl;
			ofsptn << fixed << tab << "val = " << ep->weight << endl;
			ofsptn << "}" << endl;
			ofsptn << endl;
		}

		ofsptn.close();
	}
	else {
		cerr << "Error: file " << strptd << " could not be opened." << endl;
		ofsptn.close();
	}
}

/***********************************************************
 * MCML 3.0.0 Monte Carlo Methods - Integrated Implementation
 ***********************************************************/

void Simulator::GenerateStepSize(Photon& photon) {
    // MCML 3.0.0 step size generation using Beer-Lambert law
    if (photon.step < 1e-10) {
        double rnd;
        // Avoid zero random number
        while ((rnd = mcml_random->next()) <= 0.0) {
            // Keep generating until we get a non-zero value
        }
        photon.step = -std::log(rnd);
    }
}

void Simulator::ScatterPhoton(Photon& photon, const Tissue& tissue) {
    // MCML 3.0.0 scattering implementation using Henyey-Greenstein phase function
    double cos_theta, sin_theta, cos_phi, sin_phi;
    double g = tissue.ani; // anisotropy factor
    double rnd;
    
    // Sample scattering angle using Henyey-Greenstein phase function
    rnd = mcml_random->next();
    if (std::abs(g) > 1e-6) {
        double temp = (1.0 - g * g) / (1.0 - g + 2.0 * g * rnd);
        cos_theta = (1.0 + g * g - temp * temp) / (2.0 * g);
        
        // Ensure cos_theta is within valid range
        if (cos_theta > 1.0) cos_theta = 1.0;
        else if (cos_theta < -1.0) cos_theta = -1.0;
    }
    else {
        // Isotropic scattering
        cos_theta = 2.0 * rnd - 1.0;
    }
    
    sin_theta = std::sqrt(1.0 - cos_theta * cos_theta);
    
    // Sample azimuthal angle
    rnd = mcml_random->next();
    cos_phi = std::cos(2.0 * PI * rnd);
    sin_phi = std::sin(2.0 * PI * rnd);
    
    // Update direction using MCML 3.0.0 spinning algorithm
    if (std::abs(photon.dir.z) > 0.99999) {
        // Special case: nearly perpendicular to z-axis
        photon.dir.x = sin_theta * cos_phi;
        photon.dir.y = sin_theta * sin_phi;
        photon.dir.z = cos_theta * (photon.dir.z > 0 ? 1.0 : -1.0);
    }
    else {
        // General case
        double temp = std::sqrt(1.0 - photon.dir.z * photon.dir.z);
        double temp_x = sin_theta * (photon.dir.x * photon.dir.z * cos_phi - photon.dir.y * sin_phi) / temp + photon.dir.x * cos_theta;
        double temp_y = sin_theta * (photon.dir.y * photon.dir.z * cos_phi + photon.dir.x * sin_phi) / temp + photon.dir.y * cos_theta;
        double temp_z = -sin_theta * cos_phi * temp + photon.dir.z * cos_theta;
        
        photon.dir.x = temp_x;
        photon.dir.y = temp_y;
        photon.dir.z = temp_z;
    }
    
    // Mark that photon has scattered
    photon.scatters = true;
}

void Simulator::RoulettePhoton(Photon& photon) {
    // MCML 3.0.0 Russian roulette implementation
    if (photon.weight < mcml_weight_threshold) {
        if (mcml_random->next() <= 0.1) {
            // Survive with weight boost
            photon.weight *= 10.0;
        }
        else {
            // Terminate photon
            photon.alive = false;
        }
    }
}

void Simulator::SetMcmlSeed(int seed) {
    mcml_random->seed(seed);
}
