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
void Photon::add_internal_vertex(std::shared_ptr<PhotonNode> vert) noexcept 
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

void Photon::add_external_vertex(std::shared_ptr<PhotonNode> vert) noexcept 
{
	if (path_last) {
		path_last->emit = vert;
		vert->prev = path_last;
	}
	++num_seg_ext;
}

void Photon::initialize_path(const glm::dvec3& start_pos, double path_weight) 
{
	path_head = std::make_shared<PhotonNode>(start_pos, path_weight);
	path_last = path_head;
	num_seg_int = 1;
	num_seg_ext = 1;
}

// Source methods
void Photon::set_source_data(const Source& src_data) 
{
	source = src_data;
}

void Photon::set_source_data(uint64_t src_id, const glm::dvec3& origin, const glm::dvec3& src_direction) 
{
	source.id = src_id;
	source.origin = origin;  
	source.direction = src_direction;
}

// Path analysis methods
glm::dvec3 Photon::get_entrance_position() const noexcept {
	return path_head ? path_head->position : glm::dvec3(0.0);
}

glm::dvec3 Photon::get_entrance_direction() const noexcept {
	return source.direction;
}

glm::dvec3 Photon::get_exit_position() const noexcept {
	// Find the first exit point in the path chain
	if (!path_head) return glm::dvec3(0.0);
	
	auto current = path_head;
	while (current) {
		if (current->emit && current->emit->emitter) {
			return current->emit->emitter->position;
		}
		current = current->next;
	}
	return glm::dvec3(0.0); // No exit found
}

glm::dvec3 Photon::get_exit_direction() const noexcept
{
	// Find the first exit direction in the path chain  
	if (!path_head) return glm::dvec3(0.0);
	
	auto current = path_head;
	while (current) {
		if (current->emit && current->emit->emitter) {
			return current->emit->emitter->direction;
		}
		current = current->next;
	}
	return glm::dvec3(0.0); // No exit found
}

glm::dvec3 Photon::get_termination_position() const noexcept
{
	return position; // Current photon position is termination position
}

glm::dvec3 Photon::get_termination_direction() const noexcept
{
	return direction; // Current photon direction is termination direction
}

uint32_t Photon::get_path_scatter_count() const noexcept
{
	if (!path_head) return 0;
	
	uint32_t count = 0;
	auto current = path_head;
	while (current->next) {
		count++;
		current = current->next;
	}
	return count;
}

double Photon::get_total_path_length() const noexcept
{
	if (!path_head) return 0.0;
	
	double total_length = 0.0;
	auto current = path_head;
	while (current->next) {
		auto next = current->next;
		total_length += glm::length(next->position - current->position);
		current = next;
	}
	return total_length;
}

double Photon::get_total_absorption_deposited() const noexcept
{
	return total_energy_absorbed; // Already tracked in photon
}

bool Photon::has_exit() const noexcept
{
	if (!path_head) return false;
	
	auto current = path_head;
	while (current) {
		if (current->emit && current->emit->emitter) {
			return true;
		}
		current = current->next;
	}
	return false;
}

bool Photon::exited_medium() const noexcept
{
	return exit_type != ExitType::NONE;
}

std::string Photon::get_termination_reason() const noexcept
{
	if (!alive) {
		if (weight <= 0.0) return "absorbed";
		if (exit_type != ExitType::NONE) return "exited";
		return "terminated";
	}
	return "active";
}

std::vector<glm::dvec3> Photon::get_all_exit_positions() const
{
	std::vector<glm::dvec3> exit_positions;
	if (!path_head) return exit_positions;
	
	auto current = path_head;
	while (current) {
		if (current->emit && current->emit->emitter) {
			exit_positions.push_back(current->emit->emitter->position);
		}
		current = current->next;
	}
	return exit_positions;
}