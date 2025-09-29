#include "volume.hpp"

// Add includes for complete type definitions
#include "simulator/voxel.hpp"
#include "../common/result.hpp"
#include "../common/error_types.hpp"

#include <algorithm>
#include <iostream>
#include <memory>
#include <random>
#include <stdexcept>
#include <sstream>

#include "math/random.hpp"
#include "simulator/simulator.hpp"
#include "logger.hpp"

// Default constructor
Volume::Volume() : dimensions_ {0, 0, 0}, total_voxels_ {0}, voxel_size_ {0.0} {
}

// Parameterized constructor
Volume::Volume(double voxel_size, uint32_t num_x, uint32_t num_y, uint32_t num_z) :
	dimensions_(num_x, num_y, num_z), total_voxels_(static_cast<uint64_t>(num_x) * num_y * num_z),
	voxel_size_(voxel_size), voxels_() {
	if (voxel_size <= 0.0) {
		throw std::invalid_argument("Voxel size must be positive");
	}

	if (num_x == 0 || num_y == 0 || num_z == 0) {
		throw std::invalid_argument("Grid dimensions must be positive");
	}

	// Check for potential overflow
	if (total_voxels_ / num_x / num_y != num_z) {
		throw std::overflow_error("Grid dimensions would cause integer overflow");
	}

	initialize_voxels();
}

// Destructor
Volume::~Volume() {
	cleanup_voxels();
}

// Move constructor
Volume::Volume(Volume&& other) noexcept :
	dimensions_(other.dimensions_), total_voxels_(other.total_voxels_), voxel_size_(other.voxel_size_),
	voxels_(std::move(other.voxels_)) {
	// Reset the moved-from object
	other.dimensions_ = glm::uvec3(0, 0, 0);
	other.total_voxels_ = 0;
	other.voxel_size_ = 0.0;
}

// Move assignment operator
Volume& Volume::operator=(Volume&& other) noexcept {
	if (this != &other) {
		// Clean up current resources
		cleanup_voxels();

		// Move data from other
		dimensions_ = other.dimensions_;
		total_voxels_ = other.total_voxels_;
		voxel_size_ = other.voxel_size_;
		voxels_ = std::move(other.voxels_);

		// Reset the moved-from object
		other.dimensions_ = glm::uvec3(0, 0, 0);
		other.total_voxels_ = 0;
		other.voxel_size_ = 0.0;
	}
	return *this;
}

// Grid access operators
Voxel* Volume::operator()(uint32_t x, uint32_t y, uint32_t z) {
	if (!is_valid_coordinate(x, y, z)) {
		throw std::out_of_range("Voxel coordinate out of bounds");
	}
	return voxels_[calculate_index(x, y, z)].get();
}

const Voxel* Volume::operator()(uint32_t x, uint32_t y, uint32_t z) const {
	if (!is_valid_coordinate(x, y, z)) {
		throw std::out_of_range("Voxel coordinate out of bounds");
	}
	return voxels_[calculate_index(x, y, z)].get();
}

// Safe access methods with bounds checking
Voxel* Volume::at(uint32_t x, uint32_t y, uint32_t z) {
	if (!is_valid_coordinate(x, y, z)) {
		throw std::out_of_range("Voxel coordinate out of bounds");
	}
	return voxels_.at(calculate_index(x, y, z)).get();
}

const Voxel* Volume::at(uint32_t x, uint32_t y, uint32_t z) const {
	if (!is_valid_coordinate(x, y, z)) {
		throw std::out_of_range("Voxel coordinate out of bounds");
	}
	return voxels_.at(calculate_index(x, y, z)).get();
}

// Linear indexing
uint32_t Volume::calculate_index(uint32_t x, uint32_t y, uint32_t z) const {
	// Standard 3D to 1D mapping: index = x + y*nx + z*nx*ny
	return x + y * dimensions_.x + z * dimensions_.x * dimensions_.y;
}

Voxel* Volume::at(uint32_t linear_index) {
	if (!is_valid_index(linear_index)) {
		throw std::out_of_range("Linear index out of bounds");
	}
	return voxels_.at(linear_index).get();
}

const Voxel* Volume::at(uint32_t linear_index) const {
	if (!is_valid_index(linear_index)) {
		throw std::out_of_range("Linear index out of bounds");
	}
	return voxels_.at(linear_index).get();
}

// Grid management
void Volume::clear() {
	cleanup_voxels();
	dimensions_ = glm::uvec3(0, 0, 0);
	total_voxels_ = 0;
	voxel_size_ = 0.0;
}

// Private helper methods for grid management
void Volume::initialize_voxels() {
	// Modern C++20: Pre-allocate and use emplace_back for better performance
	voxels_.clear();
	voxels_.reserve(total_voxels_);

	// Initialize voxels in z-y-x order for better cache locality
	// Use emplace_back to construct in-place and avoid unnecessary copies
	for (uint32_t z = 0; z < dimensions_.z; ++z) {
		for (uint32_t y = 0; y < dimensions_.y; ++y) {
			for (uint32_t x = 0; x < dimensions_.x; ++x) {
				voxels_.emplace_back(std::make_unique<Voxel>(x, y, z));
			}
		}
	}

	// Shrink to fit for optimal memory usage
	voxels_.shrink_to_fit();
}

void Volume::cleanup_voxels() {
	// Smart pointers automatically clean up when the vector is cleared
	voxels_.clear();
}

// Volume calculation methods (integrated from VoxelVolumeCalculator)
double Volume::fraction_inside(const glm::dvec3& voxel_min, const glm::dvec3& voxel_max, 
							   const std::vector<Layer>& layers, int samples) const {
	if (samples <= 0) {
		return 0.0;
	}

	// Monte Carlo sampling
	std::random_device rd;
	std::mt19937 gen(rd());
	std::uniform_real_distribution<double> dis(0.0, 1.0);

	int inside_count = 0;

	for (int i = 0; i < samples; i++) {
		// Generate random point within voxel
		glm::dvec3 point = glm::dvec3(voxel_min.x + dis(gen) * (voxel_max.x - voxel_min.x),
									  voxel_min.y + dis(gen) * (voxel_max.y - voxel_min.y),
									  voxel_min.z + dis(gen) * (voxel_max.z - voxel_min.z));

		for (const auto& layer : layers) {
			if (layer.contains_point(point)) {
				inside_count++;
			}
		}
	}

	return static_cast<double>(inside_count) / static_cast<double>(samples);
}

double Volume::fraction_inside_fast(const glm::dvec3& voxel_min, const glm::dvec3& voxel_max,
									const std::vector<Layer>& layers, int max_subdivisions) const {
	// Use robust distance field based voxelization and return volume fraction
	(void)max_subdivisions; // Unused parameter - suppress warning
	VoxelClassification classification = distance_field_voxelization(voxel_min, voxel_max, layers);
	return classification.volume_fraction;
}

double Volume::fraction(uint32_t x, uint32_t y, uint32_t z, const glm::dvec3& grid_min,
						const std::vector<Layer>& layers, int samples) const {
	if (!is_valid_coordinate(x, y, z))
		return 0.0;

	// Calculate voxel bounds
	glm::dvec3 voxel_min = grid_min + glm::dvec3(x * voxel_size_, y * voxel_size_, z * voxel_size_);
	glm::dvec3 voxel_max = voxel_min + glm::dvec3(voxel_size_, voxel_size_, voxel_size_);
	return fraction_inside(voxel_min, voxel_max, layers, samples);
}

double Volume::subdivide_test(const glm::dvec3& min_corner, const glm::dvec3& max_corner,
							  const std::vector<Layer>& layers, int depth, int max_depth) const {
	// Test corner points
	std::vector<glm::dvec3> corners = {
		min_corner,
		{max_corner.x, min_corner.y, min_corner.z},
		{min_corner.x, max_corner.y, min_corner.z},
		{max_corner.x, max_corner.y, min_corner.z},
		{min_corner.x, min_corner.y, max_corner.z},
		{max_corner.x, min_corner.y, max_corner.z},
		{min_corner.x, max_corner.y, max_corner.z},
		max_corner
	};

	int inside_corners = 0;
	for (const auto& corner : corners) {
		for (const auto& layer : layers) {
			if (layer.contains_point(corner)) {
				inside_corners++;
			}
		}
	}

	// If all corners have same state and we're not at max depth, we're done
	if (inside_corners == 0) {
		return 0.0; // Fully outside
	}
	else if (inside_corners == 8) {
		return 1.0; // Fully inside
	}
	else if (depth >= max_depth) {
		// At max depth, estimate based on corner ratio
		return static_cast<double>(inside_corners) / 8.0;
	}

	// Mixed state - subdivide
	glm::dvec3 center = (min_corner + max_corner) * 0.5;

	// 8 sub-cubes
	std::vector<std::pair<glm::dvec3, glm::dvec3>> sub_cubes = {
		{min_corner, center},
		{{center.x, min_corner.y, min_corner.z}, {max_corner.x, center.y, center.z}},
		{{min_corner.x, center.y, min_corner.z}, {center.x, max_corner.y, center.z}},
		{{center.x, center.y, min_corner.z}, {max_corner.x, max_corner.y, center.z}},
		{{min_corner.x, min_corner.y, center.z}, {center.x, center.y, max_corner.z}},
		{{center.x, min_corner.y, center.z}, {max_corner.x, center.y, max_corner.z}},
		{{min_corner.x, center.y, center.z}, {center.x, max_corner.y, max_corner.z}},
		{center, max_corner}};

	double total_fraction = 0.0;
	for (const auto& sub_cube : sub_cubes) {
		total_fraction += subdivide_test(sub_cube.first, sub_cube.second, layers, depth + 1, max_depth);
	}

	return total_fraction / 8.0; // Average of 8 sub-cubes
}

// Distance Field Based Voxelization Algorithm
// This is the most robust approach for complex mesh voxelization
VoxelClassification Volume::distance_field_voxelization(const glm::dvec3& voxel_min, const glm::dvec3& voxel_max,
													   const std::vector<Layer>& layers) const {
	
	VoxelClassification result;
	glm::dvec3 voxel_size = voxel_max - voxel_min;
	glm::dvec3 voxel_center = (voxel_min + voxel_max) * 0.5;
	
	// Key parameters for distance field voxelization
	const double SURFACE_THRESHOLD = voxel_size.x * 0.5;  // Half voxel size
	
	// DEBUG: Check for voxels near the problematic area
	bool is_debug_voxel = Config::get().log() && 
		voxel_center.x > -0.1 && voxel_center.x < 0.0 && 
		voxel_center.y > -0.15 && voxel_center.y < -0.1 && 
		voxel_center.z > -0.35 && voxel_center.z < -0.3;
	
	if (is_debug_voxel) {
		std::ostringstream debug_msg;
		debug_msg << "DEBUG VOXEL - Center: (" 
				  << voxel_center.x << ", " 
				  << voxel_center.y << ", " 
				  << voxel_center.z << ")";
		Logger::instance().log_debug(debug_msg.str());
		
		std::ostringstream bounds_msg;
		bounds_msg << "  Bounds: [" << voxel_min.x << "," 
				  << voxel_max.x << "] [" 
				  << voxel_min.y << "," 
				  << voxel_max.y << "] [" 
				  << voxel_min.z << "," 
				  << voxel_max.z << "]";
		Logger::instance().log_debug(bounds_msg.str());
	}
	
	// Step 1: Compute signed distance from voxel center to mesh surface
	double min_distance = std::numeric_limits<double>::max();
	bool found_surface = false;
	uint32_t closest_layer_id = 0;
	
	for (size_t i = 0; i < layers.size(); ++i) {
		const auto& layer = layers[i];
		double distance = compute_signed_distance_to_mesh(voxel_center, layer.mesh);
		if (std::abs(distance) < std::abs(min_distance)) {
			min_distance = distance;
			closest_layer_id = static_cast<uint32_t>(i);
			found_surface = true;
		}
	}
	
	// Step 2: Classify voxel based on distance to surface
	if (!found_surface) {
		return result; // Empty result (outside all geometry)
	}
	
	// Step 3: Determine which layer(s) contain the voxel and assign correct material
	bool voxel_intersects_geometry = false;
	uint32_t correct_tissue_id = 0;
	uint8_t correct_layer_id = 0;
	bool found_containing_layer = false;
	
	// Check if voxel center is inside any layer - assign material from the containing layer
	// For nested geometries, later layers (inner ones) should override earlier layers (outer ones)
	for (const auto& layer : layers) {
		if (layer.contains_point(voxel_center)) {
			voxel_intersects_geometry = true;
			correct_tissue_id = layer.tissue_id;
			correct_layer_id = layer.id;
			found_containing_layer = true;
			if (is_debug_voxel) {
				std::ostringstream debug_msg;
				debug_msg << "  Center inside layer " << layer.id << " (material " << layer.tissue_id << ")";
				Logger::instance().log_debug(debug_msg.str());
			}
			// Continue checking other layers to allow inner layers to override outer layers
		}
	}
	
	if (is_debug_voxel && !found_containing_layer) {
		Logger::instance().log_debug("  Center not inside any layer");
	}
	
	// Also check if surface passes through voxel (distance < half voxel size)
	if (!voxel_intersects_geometry && std::abs(min_distance) <= SURFACE_THRESHOLD) {
		voxel_intersects_geometry = true;
		// For boundary voxels not containing center, use closest layer material
		if (!found_containing_layer) {
			correct_tissue_id = layers[closest_layer_id].tissue_id;
			correct_layer_id = layers[closest_layer_id].id;
		}
	}
	
	if (!voxel_intersects_geometry) {
		if (is_debug_voxel) {
			std::cout << "  No geometry intersection - returning empty result" << std::endl;
		}
		return result; // Empty result (no intersection)
	}
	
	if (is_debug_voxel) {
		std::ostringstream debug_msg;
		debug_msg << "  Final classification: intersects=" << voxel_intersects_geometry 
				  << ", material=" << correct_tissue_id 
				  << ", boundary=" << (std::abs(min_distance) <= SURFACE_THRESHOLD)
				  << ", min_dist=" << min_distance;
		Logger::instance().log_debug(debug_msg.str());
	}
	
	// Step 4: Compute volume fraction and surface classification for intersecting voxels
	result.is_inside_geometry = true;
	result.dominant_tissue_id = correct_tissue_id; // Use material from containing layer, not closest surface
	result.dominant_layer_id = correct_layer_id;   // Use layer ID for rendering
	
	// IMPROVED Boundary voxel detection: Check if surface passes through voxel volume
	// Method 1: Center-based detection (existing)
	bool center_near_surface = (std::abs(min_distance) <= SURFACE_THRESHOLD);
	
	// Method 2: Volume-based detection - check if surface intersects voxel bounds
	bool surface_intersects_volume = false;
	
	// Sample multiple points within the voxel to detect surface intersection
	const int SAMPLE_DENSITY = 3; // 3x3x3 sampling grid
	double step_size = voxel_size.x / SAMPLE_DENSITY;
	int inside_samples = 0;
	int total_samples = 0;
	
	for (int i = 0; i < SAMPLE_DENSITY; i++) {
		for (int j = 0; j < SAMPLE_DENSITY; j++) {
			for (int k = 0; k < SAMPLE_DENSITY; k++) {
				glm::dvec3 sample_point = voxel_min + glm::dvec3(
					(i + 0.5) * step_size,
					(j + 0.5) * step_size,
					(k + 0.5) * step_size
				);
				
				bool point_inside = false;
				for (const auto& layer : layers) {
					if (layer.contains_point(sample_point)) {
						point_inside = true;
						break;
					}
				}
				
				if (point_inside) inside_samples++;
				total_samples++;
			}
		}
	}
	
	// If not all samples have the same inside/outside state, surface intersects volume
	surface_intersects_volume = (inside_samples > 0 && inside_samples < total_samples);
	
	// Boundary voxel if either center is near surface OR surface intersects volume
	result.is_boundary_voxel = center_near_surface || surface_intersects_volume;
	
	if (is_debug_voxel) {
		std::ostringstream debug_msg;
		debug_msg << "  Boundary detection: center_near=" << center_near_surface 
				  << ", surface_intersects=" << surface_intersects_volume
				  << ", samples=" << inside_samples << "/" << total_samples
				  << ", final_boundary=" << result.is_boundary_voxel;
		Logger::instance().log_debug(debug_msg.str());
	}
	
	if (result.is_boundary_voxel) {
		// For boundary voxels, compute accurate volume fraction
		result.volume_fraction = compute_volume_fraction_sdf(voxel_min, voxel_max, layers, 4);
	} else {
		// For interior voxels, assume full volume
		result.volume_fraction = 1.0;
	}
	
	return result;
}

// Compute accurate volume fraction for surface voxels using SDF sampling
double Volume::compute_volume_fraction_sdf(const glm::dvec3& voxel_min, const glm::dvec3& voxel_max,
										   const std::vector<Layer>& layers, int subdivisions) const {
	
	glm::dvec3 voxel_size = voxel_max - voxel_min;
	
	// Adaptive sampling density - more samples for surface voxels
	int samples_per_axis = std::max(4, std::min(12, subdivisions * 3));
	double sample_spacing = voxel_size.x / samples_per_axis;
	
	int inside_samples = 0;
	int total_samples = 0;
	
	// Regular grid sampling within the voxel
	for (int i = 0; i < samples_per_axis; i++) {
		for (int j = 0; j < samples_per_axis; j++) {
			for (int k = 0; k < samples_per_axis; k++) {
				// Sample point within voxel (offset from edges)
				glm::dvec3 sample_point = voxel_min + glm::dvec3(
					(i + 0.5) * sample_spacing,
					(j + 0.5) * sample_spacing,
					(k + 0.5) * sample_spacing
				);
				
				// Use signed distance to determine if point is inside
				bool point_inside = false;
				for (const auto& layer : layers) {
					double distance = compute_signed_distance_to_mesh(sample_point, layer.mesh);
					if (distance <= 0.0) {  // Negative distance = inside
						point_inside = true;
						break;
					}
				}
				
				if (point_inside) {
					inside_samples++;
				}
				total_samples++;
			}
		}
	}
	
	// Compute raw volume fraction
	double raw_fraction = static_cast<double>(inside_samples) / total_samples;
	
	// Apply minimum threshold to filter out barely-touching voxels
	const double MIN_VOLUME_THRESHOLD = 0.05;  // 5% minimum
	if (raw_fraction < MIN_VOLUME_THRESHOLD) {
		return 0.0;
	}
	
	return raw_fraction;
}

// Compute signed distance from a point to the closest point on a triangle mesh
// Negative distance = inside, positive = outside
double Volume::compute_signed_distance_to_mesh(const glm::dvec3& point, const std::vector<Triangle>& mesh) const {
	
	if (mesh.empty()) {
		return std::numeric_limits<double>::max();
	}
	
	double min_distance_squared = std::numeric_limits<double>::max();
	glm::dvec3 closest_point;
	glm::dvec3 closest_normal;
	
	// Find closest point on the mesh surface
	for (const auto& triangle : mesh) {
		glm::dvec3 closest_on_triangle = closest_point_on_triangle(point, triangle);
		double dist_squared = glm::length2(point - closest_on_triangle);
		
		if (dist_squared < min_distance_squared) {
			min_distance_squared = dist_squared;
			closest_point = closest_on_triangle;
			closest_normal = glm::normalize(glm::cross(
				triangle.v1() - triangle.v0(),
				triangle.v2() - triangle.v0()
			));
		}
	}
	
	double distance = std::sqrt(min_distance_squared);
	
	// Determine sign using normal direction
	glm::dvec3 to_point = point - closest_point;
	double sign = glm::dot(to_point, closest_normal);
	
	return (sign >= 0.0) ? distance : -distance;
}

// Find closest point on a triangle to a given point
glm::dvec3 Volume::closest_point_on_triangle(const glm::dvec3& point, const Triangle& triangle) const {
	
	glm::dvec3 v0 = triangle.v0();
	glm::dvec3 v1 = triangle.v1();
	glm::dvec3 v2 = triangle.v2();
	
	// Vectors from v0 to other vertices
	glm::dvec3 edge0 = v1 - v0;
	glm::dvec3 edge1 = v2 - v0;
	glm::dvec3 v0_to_point = point - v0;
	
	// Compute dot products
	double a = glm::dot(edge0, edge0);
	double b = glm::dot(edge0, edge1);
	double c = glm::dot(edge1, edge1);
	double d = glm::dot(edge0, v0_to_point);
	double e = glm::dot(edge1, v0_to_point);
	
	double det = a * c - b * b;
	double s = b * e - c * d;
	double t = b * d - a * e;
	
	if (s + t < det) {
		if (s < 0.0) {
			if (t < 0.0) {
				// Region 4
				if (d < 0.0) {
					s = std::clamp(-d / a, 0.0, 1.0);
					t = 0.0;
				} else {
					s = 0.0;
					t = std::clamp(-e / c, 0.0, 1.0);
				}
			} else {
				// Region 3
				s = 0.0;
				t = std::clamp(-e / c, 0.0, 1.0);
			}
		} else if (t < 0.0) {
			// Region 5
			s = std::clamp(-d / a, 0.0, 1.0);
			t = 0.0;
		} else {
			// Region 0
			double inv_det = 1.0 / det;
			s *= inv_det;
			t *= inv_det;
		}
	} else {
		if (s < 0.0) {
			// Region 2
			double tmp0 = b + d;
			double tmp1 = c + e;
			if (tmp1 > tmp0) {
				double numer = tmp1 - tmp0;
				double denom = a - 2.0 * b + c;
				s = std::clamp(numer / denom, 0.0, 1.0);
				t = 1.0 - s;
			} else {
				t = std::clamp(-e / c, 0.0, 1.0);
				s = 0.0;
			}
		} else if (t < 0.0) {
			// Region 6
			double tmp0 = b + e;
			double tmp1 = a + d;
			if (tmp1 > tmp0) {
				double numer = tmp1 - tmp0;
				double denom = a - 2.0 * b + c;
				t = std::clamp(numer / denom, 0.0, 1.0);
				s = 1.0 - t;
			} else {
				s = std::clamp(-d / a, 0.0, 1.0);
				t = 0.0;
			}
		} else {
			// Region 1
			double numer = c + e - b - d;
			if (numer <= 0.0) {
				s = 0.0;
			} else {
				double denom = a - 2.0 * b + c;
				s = std::clamp(numer / denom, 0.0, 1.0);
			}
			t = 1.0 - s;
		}
	}
	
	return v0 + s * edge0 + t * edge1;
}
