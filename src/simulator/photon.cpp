#include "photon.hpp"

// Photon constructor implementations
Photon::Photon(uint64_t i, const glm::dvec3& src_origin, const glm::dvec3& src_direction) noexcept 
	: id(i) 
{
	source.origin = src_origin;
	source.direction = src_direction;
}

Photon::Photon(uint64_t i, const glm::dvec3& src_origin, const glm::dvec3& src_direction,
			   const glm::dvec3& spec_direction, const glm::dvec3& src_intersect, 
			   const Triangle& src_triangle) noexcept 
	: id(i) 
{
	source.origin = src_origin;
	source.direction = src_direction;
	source.specular_direction = spec_direction;
	source.intersect = src_intersect;
	source.triangle = src_triangle;
}

// Path management methods (integrated from PhotonPath)
void Photon::add_internal_vertex(std::shared_ptr<Node> vert) noexcept 
{
	if (!path_head) {
		path_head = vert;
		path_last = vert;
	} else {
		path_last->next = vert;
		vert->prev = path_last;
		path_last = vert;
	}
	++num_seg_int;
}

void Photon::add_external_vertex(std::shared_ptr<Node> vert) noexcept 
{
	if (path_last) {
		path_last->emit = vert;
		vert->prev = path_last;
	}
	++num_seg_ext;
}

void Photon::initialize_path(const glm::dvec3& start_pos, double path_weight) 
{
	path_head = std::make_shared<Node>(start_pos, path_weight);
	path_last = path_head;
	num_seg_int = 1;
	num_seg_ext = 1;
}

// Source methods
void Photon::set_source_data(const SourceData& src_data) 
{
	source = src_data;
}

void Photon::set_source_data(uint64_t src_id, const glm::dvec3& origin, const glm::dvec3& src_direction) 
{
	source.id = src_id;
	source.origin = origin;
	source.direction = src_direction;
}