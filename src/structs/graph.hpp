#pragma once

#include "vertex.hpp"

#include <cstdint>

struct Graph {
	uint64_t id;
	Vertex* head;
	Vertex* last;

	uint64_t num_seg_int; // internal segments
	uint64_t num_seg_ext; // emittant segments

	Graph(uint64_t i, Vertex* h) {
		id = i;
		head = h;
		last = h;
		num_seg_int = 1;
		num_seg_ext = 1;
	}

	void add_internal_vertex(Vertex* vert) {
		last->next = vert;
		vert->prev = last;
		last = vert;
		++num_seg_int;
	}

	void add_external_vertex(Vertex* vert) {
		last->emit = vert;
		vert->prev = last;
		++num_seg_ext;
	}
};
