#pragma once

struct Record {
	double at; // total absorption
	double rd; // diffuse reflection
	double rs; // specular reflection
	double td; // diffuse transmission
	double ts; // specular transmission

	Record() : at(0), rd(0), rs(0), td(0), ts(0) {}
};
