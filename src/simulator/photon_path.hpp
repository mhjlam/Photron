#pragma once

#include <cstdint>
#include <memory>

#include "photon_node.hpp"

struct PhotonPath
{
	uint64_t id;
	std::shared_ptr<PhotonNode> head;
	std::shared_ptr<PhotonNode> last;

	uint64_t num_seg_int {1}; // internal segments
	uint64_t num_seg_ext {1}; // emittant segments

	explicit PhotonPath(uint64_t i, std::shared_ptr<PhotonNode> h) noexcept : id(i), head(h), last(h) {}

	void add_internal_vertex(std::shared_ptr<PhotonNode> vert) noexcept {
		last->next = vert;
		vert->prev = last;
		last = vert;
		++num_seg_int;
	}

	void add_external_vertex(std::shared_ptr<PhotonNode> vert) noexcept {
		last->emit = vert;
		vert->prev = last;
		++num_seg_ext;
	}
};
