#include "experimenter.hpp"

#include <algorithm>
#include <cmath>
#include <cstdint>
#include <fstream>
#include <iostream>
#include <string>

#ifdef _WIN32
#include <windows.h> // for GetTickCount
#endif

#include "utilities.hpp"

uint64_t Experimenter::t0_;
uint64_t Experimenter::t1_;

double Experimenter::total_absorption_;
double Experimenter::total_reflection_;
double Experimenter::total_transmission_;
double Experimenter::total_diffusion_;

double Experimenter::surface_reflection_;
double Experimenter::surface_refraction_;

double Experimenter::path_length_;
double Experimenter::scatter_events_;
double Experimenter::diffusion_distance_;
double Experimenter::average_step_size_;

std::vector<double> Experimenter::step_sizes_;
std::vector<glm::dvec3> Experimenter::path_vertices_;

void Experimenter::increment_scatters() {
	scatter_events_++;
}

void Experimenter::add_vertex(double x, double y, double z) {
	path_vertices_.push_back(glm::dvec3(x, y, z));
}

void Experimenter::add_step_size(double s) {
	step_sizes_.push_back(s);
}

void Experimenter::collect_data(double at, double rs, double rd, double ts, double td) {
	total_absorption_ = at;
	total_reflection_ = rs + rd;
	total_transmission_ = ts + td;
	total_diffusion_ = rd + ts + td;

	surface_reflection_ = rs;
	surface_refraction_ = 1.0 - rs;

	path_length_ = compute_path_length();
	average_step_size_ = compute_average_step_size();
	diffusion_distance_ = compute_diffusion_distance();
}

double Experimenter::compute_diffusion_distance() {
	if (path_vertices_.empty()) {
		return 0;
	}

	auto [min_x, max_x] = std::minmax_element(path_vertices_.begin(), path_vertices_.end(),
											  [](const glm::dvec3& a, const glm::dvec3& b) { return a.x < b.x; });
	auto [min_y, max_y] = std::minmax_element(path_vertices_.begin(), path_vertices_.end(),
											  [](const glm::dvec3& a, const glm::dvec3& b) { return a.y < b.y; });
	auto [min_z, max_z] = std::minmax_element(path_vertices_.begin(), path_vertices_.end(),
											  [](const glm::dvec3& a, const glm::dvec3& b) { return a.z < b.z; });

	double dx = max_x->x - min_x->x;
	double dy = max_y->y - min_y->y;
	double dz = max_z->z - min_z->z;

	return std::sqrt(dx * dx + dy * dy + dz * dz);
}

double Experimenter::compute_average_step_size() {
	if (step_sizes_.empty()) {
		return 0;
	}

	double total_step_size = 0;
	size_t num_steps = step_sizes_.size();

	for (size_t i = 0; i < num_steps; ++i) {
		total_step_size += step_sizes_[i];
	}
	return total_step_size / static_cast<double>(num_steps);
}

double Experimenter::compute_path_length() {
	if (path_vertices_.empty()) {
		return 0;
	}

	double path_len = 0;
	for (std::size_t i = 1; i < path_vertices_.size(); ++i) {
		path_len += distribution_(path_vertices_[i - 1], path_vertices_[i]);
	}
	return path_len;
}

void Experimenter::start_clock() {
#ifdef _WIN32
	t0_ = GetTickCount();
#endif
}

void Experimenter::stop_clock() {
#ifdef _WIN32
	t1_ = GetTickCount() - t0_;
#endif
}

void Experimenter::write_to_file() {
	std::ofstream out("experiment.out", std::ios_base::app);
	if (!out.good()) {
		return;
	}

	out << "Path length        " << path_length_ << std::endl;
	out << "Scatter events     " << scatter_events_ << std::endl;
	out << "Average step size  " << average_step_size_ << std::endl;
	out << "Diffusion distance " << diffusion_distance_ << std::endl;

	out << "Total absorption   " << total_absorption_ << std::endl;
	out << "Total reflection   " << total_reflection_ << std::endl;
	out << "Total transmission " << total_transmission_ << std::endl;
	out << "Total diffusion    " << total_diffusion_ << std::endl;

	out << "Surface reflection " << surface_reflection_ << std::endl;
	out << "Surface refraction " << surface_refraction_ << std::endl;

	out << "Total time taken   " << t1_ << " ms" << std::endl;

	out << std::endl;
}

void Experimenter::print_report() {
	std::cout << "Path length        " << path_length_ << std::endl;
	std::cout << "Scatter events     " << scatter_events_ << std::endl;
	std::cout << "Average step size  " << average_step_size_ << std::endl;
	std::cout << "Diffusion distance " << diffusion_distance_ << std::endl;

	std::cout << "Total absorption   " << total_absorption_ << std::endl;
	std::cout << "Total reflection   " << total_reflection_ << std::endl;
	std::cout << "Total transmission " << total_transmission_ << std::endl;
	std::cout << "Total diffusion    " << total_diffusion_ << std::endl;

	std::cout << "Surface reflection " << surface_reflection_ << std::endl;
	std::cout << "Surface refraction " << surface_refraction_ << std::endl;

	std::cout << "Total time taken   " << t1_ << " ms" << std::endl;
}

void Experimenter::reset() {
	// Reset timing
	t0_ = 0;
	t1_ = 0;

	// Reset accumulated values
	total_absorption_ = 0.0;
	total_reflection_ = 0.0;
	total_transmission_ = 0.0;
	total_diffusion_ = 0.0;

	surface_reflection_ = 0.0;
	surface_refraction_ = 0.0;

	path_length_ = 0.0;
	scatter_events_ = 0.0;
	diffusion_distance_ = 0.0;
	average_step_size_ = 0.0;

	// Clear accumulated data collections
	step_sizes_.clear();
	path_vertices_.clear();
}
