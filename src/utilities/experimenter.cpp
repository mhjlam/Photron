#include <cmath>
#include <fstream>
#include <iostream>
#include <string>

#ifdef _WIN32
#include <windows.h> // for GetTickCount
#endif

#include "experimenter.hpp"
#include "utilities.hpp"

unsigned long Experimenter::t0;
unsigned long Experimenter::t1;

double Experimenter::totalAbsorption;
double Experimenter::totalReflection;
double Experimenter::totalTransmission;
double Experimenter::totalDiffusion;

double Experimenter::surfaceReflection;
double Experimenter::surfaceRefraction;

double Experimenter::pathLength;
double Experimenter::scatterEvents;
double Experimenter::diffusionDistance;
double Experimenter::averageStepSize;

std::vector<double> Experimenter::stepSizes;
std::vector<Point3> Experimenter::pathVertices;

void Experimenter::IncrementScatters() {
	scatterEvents++;
}

void Experimenter::AddVertex(double x, double y, double z) {
	pathVertices.push_back(Point3(x, y, z));
}

void Experimenter::AddStepSize(double s) {
	stepSizes.push_back(s);
}

void Experimenter::CollectData(double at, double rs, double rd, double ts, double td) {
	totalAbsorption = at;
	totalReflection = rs + rd;
	totalTransmission = ts + td;
	totalDiffusion = rd + ts + td;

	surfaceReflection = rs;
	surfaceRefraction = 1.0 - rs;

	pathLength = ComputePathLength();
	averageStepSize = ComputeAverageStepSize();
	diffusionDistance = ComputeDiffusionDistance();
}

double Experimenter::ComputeDiffusionDistance() {
	if (pathVertices.empty())
		return 0;

	std::vector<Point3>::iterator it = pathVertices.begin();

	double minx = it->x;
	double maxx = it->x;
	double miny = it->y;
	double maxy = it->y;
	double minz = it->z;
	double maxz = it->z;

	for (it = pathVertices.begin() + 1; it != pathVertices.end(); ++it) {
		if (it->x < minx) {
			minx = it->x;
		}
		if (it->x > maxx) {
			maxx = it->x;
		}
		if (it->y < miny) {
			miny = it->y;
		}
		if (it->y > maxy) {
			maxy = it->y;
		}
		if (it->z < minz) {
			minz = it->z;
		}
		if (it->z > maxz) {
			maxz = it->z;
		}
	}

	double dx = maxx - minx;
	double dy = maxy - miny;
	double dz = maxz - minz;

	return sqrt(dx * dx + dy * dy + dz * dz);
}

double Experimenter::ComputeAverageStepSize() {
	if (stepSizes.empty()) {
		return 0;
	}

	double totalstepsize = 0;
	size_t numsteps = stepSizes.size();

	for (size_t i = 0; i < numsteps; ++i) {
		totalstepsize += stepSizes[i];
	}

	return totalstepsize / static_cast<double>(numsteps);
}

double Experimenter::ComputePathLength() {
	if (pathVertices.empty()) {
		return 0;
	}

	double pathlength = 0;

	for (unsigned int i = 1; i < pathVertices.size(); ++i) {
		pathlength += dist(pathVertices[i - 1], pathVertices[i]);
	}

	return pathlength;
}

void Experimenter::StartClock() {
#ifdef _WIN32
	t0 = GetTickCount();
#endif
}

void Experimenter::StopClock() {
#ifdef _WIN32
	t1 = GetTickCount() - t0;
#endif
}

void Experimenter::WriteToFile() {
	std::ofstream out("experiment.out", std::ios_base::app);
	if (!out.good()) {
		return;
	}

	out << "Path length        " << pathLength << std::endl;
	out << "Scatter events     " << scatterEvents << std::endl;
	out << "Average step size  " << averageStepSize << std::endl;
	out << "Diffusion distance " << diffusionDistance << std::endl;

	out << "Total absorption   " << totalAbsorption << std::endl;
	out << "Total reflection   " << totalReflection << std::endl;
	out << "Total transmission " << totalTransmission << std::endl;
	out << "Total diffusion    " << totalDiffusion << std::endl;

	out << "Surface reflection " << surfaceReflection << std::endl;
	out << "Surface refraction " << surfaceRefraction << std::endl;

	out << "Total time taken   " << t1 << " ms" << std::endl;

	out << std::endl;
}

void Experimenter::Print() {
	std::cout << "Path length        " << pathLength << std::endl;
	std::cout << "Scatter events     " << scatterEvents << std::endl;
	std::cout << "Average step size  " << averageStepSize << std::endl;
	std::cout << "Diffusion distance " << diffusionDistance << std::endl;

	std::cout << "Total absorption   " << totalAbsorption << std::endl;
	std::cout << "Total reflection   " << totalReflection << std::endl;
	std::cout << "Total transmission " << totalTransmission << std::endl;
	std::cout << "Total diffusion    " << totalDiffusion << std::endl;

	std::cout << "Surface reflection " << surfaceReflection << std::endl;
	std::cout << "Surface refraction " << surfaceRefraction << std::endl;

	std::cout << "Total time taken   " << t1 << " ms" << std::endl;
}
