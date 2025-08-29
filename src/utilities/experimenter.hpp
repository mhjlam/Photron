#pragma once

#include "../structs/point3.hpp"

#include <vector>

class Experimenter
{
private:
	static unsigned long t0, t1;             // start time and time passed

	static double totalAbsorption;           // total absorption
	static double totalReflection;           // total reflection (diffuse + specular)
	static double totalTransmission;         // total transmission (diffuse + specular)
	static double totalDiffusion;            // total diffusion (diffuse reflection + total transmission)

	static double surfaceReflection;         // amount of weight that reflects at the surface
	static double surfaceRefraction;         // amount of weight that refracts at the surface

	static double pathLength;                // path length of photon migration
	static double scatterEvents;             // number of scattering events
	static double averageStepSize;           // average distance between interaction sites
	static double diffusionDistance;         // maximum distance between two points on the path

	static std::vector<double> stepSizes;    // step sizes
	static std::vector<Point3> pathVertices; // path vertices

public:
	static void AddVertex(double x, double y, double z);
	static void AddStepSize(double s);
	static void IncrementScatters();

	static void CollectData(double at, double rs, double rd, double ts, double td);
	static double ComputeDiffusionDistance();
	static double ComputeAverageStepSize();
	static double ComputePathLength();

	static void StartClock();
	static void StopClock();

	static void WriteToFile();
	static void Print();
};
