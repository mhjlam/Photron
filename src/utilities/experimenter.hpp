#pragma once

#include "../structs/glm_types.hpp"

#include <vector>

class Experimenter
{
public:
	static void add_vertex(double x, double y, double z);
	static void add_step_size(double s);
	static void increment_scatters();

	static void collect_data(double at, double rs, double rd, double ts, double td);
	static double compute_diffusion_distance();
	static double compute_average_step_size();
	static double compute_path_length();

	static void start_clock();
	static void stop_clock();

	static void write_to_file();
	static void print_report();
	static void reset();  // Reset all accumulated data
	
	// Getters for UI display
	static double get_path_length() { return path_length_; }
	static double get_scatter_events() { return scatter_events_; }
	static double get_average_step_size() { return average_step_size_; }
	static double get_diffusion_distance() { return diffusion_distance_; }
	static double get_total_absorption() { return total_absorption_; }
	static double get_total_reflection() { return total_reflection_; }
	static double get_total_transmission() { return total_transmission_; }
	static double get_total_diffusion() { return total_diffusion_; }
	static double get_surface_reflection() { return surface_reflection_; }
	static double get_surface_refraction() { return surface_refraction_; }

private:
	static uint64_t t0_, t1_;                  // start time and time passed

	static double total_absorption_;           // total absorption
	static double total_reflection_;           // total reflection (diffuse + specular)
	static double total_transmission_;         // total transmission (diffuse + specular)
	static double total_diffusion_;            // total diffusion (diffuse reflection + total transmission)

	static double surface_reflection_;         // amount of weight that reflects at the surface
	static double surface_refraction_;         // amount of weight that refracts at the surface

	static double path_length_;                // path length of photon migration
	static double scatter_events_;             // number of scattering events
	static double average_step_size_;          // average distance between interaction sites
	static double diffusion_distance_;         // maximum distance between two points on the path

	static std::vector<double> step_sizes_;    // step sizes
	static std::vector<glm::dvec3> path_vertices_; // path vertices
};
