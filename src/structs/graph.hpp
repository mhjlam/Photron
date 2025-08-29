#pragma once

#include "vertex.hpp"

struct Graph {
	long id;
	Vertex* head;
	Vertex* last;

	long numsegint; // internal segments
	long numsegext; // emittant segments

	Graph(long i, Vertex* h) {
		id = i;
		head = h;
		last = h;

		numsegint = 1;
		numsegext = 1;
	}

	//	~Graph()
	//	{
	//		Vertex* item = head;
	//
	//		if (head->prev)
	//		{
	//			delete head->prev;
	//			head->prev = NULL;
	//		}
	//
	//		while (item)
	//		{
	//			Vertex* old = item;
	//			item = item->next;
	//
	//			delete old;
	//			old = NULL;
	//		}
	//	}

	void AddInternalVertex(Vertex* vert) {
		last->next = vert;
		vert->prev = last;

		last = vert;

		++numsegint;
	}

	void AddExternalVertex(Vertex* vert) {
		last->emit = vert;
		vert->prev = last;

		++numsegext;
	}
};
