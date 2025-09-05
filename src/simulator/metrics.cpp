#include "metrics.hpp"

#include <algorithm>
#include <ranges>
#include <cmath>
#include <cstdint>
#include <fstream>
#include <iostream>
#include <string>

#ifdef _WIN32
#include <windows.h> // for GetTickCount
#endif

void Metrics::increment_scatters() {
	scatter_events_++;
}

void Metrics::add_vertex(double x, double y, double z) {
	path_vertices_.push_back(glm::dvec3(x, y, z));
}

void Metrics::add_step_size(double s) {
	step_sizes_.push_back(s);
}

void Metrics::collect_data(double at, double rs, double rd, double ts, double td) {
	total_absorption_ = at;
	total_reflection_ = rs + rd;
	total_transmission_ = ts + td;
	total_diffusion_ = rd + total_transmission_;

	surface_reflection_ = rs;
	surface_refraction_ = ts;

	path_length_ = compute_path_length();
	average_step_size_ = compute_average_step_size();
	diffusion_distance_ = compute_diffusion_distance();
}

double Metrics::compute_diffusion_distance() {
	if (path_vertices_.empty()) {
		return 0.0;
	}

	// Modern C++20: Use ranges::minmax_element for better performance
	const auto [min_x, max_x] = std::ranges::minmax_element(path_vertices_, 
		[](const glm::dvec3& a, const glm::dvec3& b) noexcept { return a.x < b.x; });
	const auto [min_y, max_y] = std::ranges::minmax_element(path_vertices_, 
		[](const glm::dvec3& a, const glm::dvec3& b) noexcept { return a.y < b.y; });
	const auto [min_z, max_z] = std::ranges::minmax_element(path_vertices_, 
		[](const glm::dvec3& a, const glm::dvec3& b) noexcept { return a.z < b.z; });

	const double dx = max_x->x - min_x->x;
	const double dy = max_y->y - min_y->y;
	const double dz = max_z->z - min_z->z;

	return std::sqrt(dx * dx + dy * dy + dz * dz);
}

double Metrics::compute_average_step_size() {
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

double Metrics::compute_path_length() {
	if (path_vertices_.empty()) {
		return 0;
	}

	double path_len = 0;
	for (std::size_t i = 1; i < path_vertices_.size(); ++i) {
		glm::dvec3 diff = path_vertices_[i] - path_vertices_[i - 1];
		path_len += glm::length(diff);
	}

	return path_len;
}

void Metrics::start_clock() {
#ifdef _WIN32
	t0_ = GetTickCount64();
#else
	// Linux implementation could use clock_gettime or similar
	t0_ = 0;
#endif
}

void Metrics::stop_clock() {
#ifdef _WIN32
	t1_ = GetTickCount64() - t0_;
#else
	// Linux implementation
	t1_ = 0;
#endif
}

void Metrics::write_to_file() {
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

void Metrics::print_report() {
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

	std::cout << std::endl;
}

void Metrics::reset() {
	t0_ = 0;
	t1_ = 0;

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

	step_sizes_.clear();
	path_vertices_.clear();
}
