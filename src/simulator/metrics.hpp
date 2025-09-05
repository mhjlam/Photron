#pragma once

#include <vector>

#include <glm/glm.hpp>

class Metrics
{
public:
	Metrics() = default;
	~Metrics() = default;

	// Data collection methods
	void add_vertex(double x, double y, double z);
	void add_step_size(double s);
	void increment_scatters();

	void collect_data(double at, double rs, double rd, double ts, double td);
	double compute_diffusion_distance();
	double compute_average_step_size();
	double compute_path_length();

	void start_clock();
	void stop_clock();

	void write_to_file();
	void print_report();
	void reset(); // Reset all accumulated data

	// Getters for UI display
	double get_path_length() const { return path_length_; }
	double get_scatter_events() const { return scatter_events_; }
	double get_average_step_size() const { return average_step_size_; }
	double get_diffusion_distance() const { return diffusion_distance_; }
	double get_total_absorption() const { return total_absorption_; }
	double get_total_reflection() const { return total_reflection_; }
	double get_total_transmission() const { return total_transmission_; }
	double get_total_diffusion() const { return total_diffusion_; }
	double get_surface_reflection() const { return surface_reflection_; }
	double get_surface_refraction() const { return surface_refraction_; }

private:
	uint64_t t0_ {0}, t1_ {0};              // start time and time passed

	double total_absorption_ {0.0};         // total absorption
	double total_reflection_ {0.0};         // total reflection (diffuse + specular)
	double total_transmission_ {0.0};       // total transmission (diffuse + specular)
	double total_diffusion_ {0.0};          // total diffusion (diffuse reflection + total transmission)

	double surface_reflection_ {0.0};       // amount of weight that reflects at the surface
	double surface_refraction_ {0.0};       // amount of weight that refracts at the surface

	double path_length_ {0.0};              // path length of photon migration
	double scatter_events_ {0.0};           // number of scattering events
	double average_step_size_ {0.0};        // average distance between interaction sites
	double diffusion_distance_ {0.0};       // maximum distance between two points on the path

	std::vector<double> step_sizes_;        // step sizes
	std::vector<glm::dvec3> path_vertices_; // path vertices
};
