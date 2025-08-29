#pragma once

#include "point3.hpp"

struct Camera {
	Point3 eye;
	Point3 center;
	Point3 up;

	float theta;
	float phi;
	float radius;

	Point3 pos;
	Point3 rot;
	float speed;

	Camera() {
		eye = Point3();
		center = Point3();
		up = Point3(0, 1, 0);

		theta = 2.8f;
		phi = 2.0f;
		radius = 3.0f;

		pos = Point3(2, 1.5, 1.5);
		rot = Point3(33, -50, 0);
		speed = 0.01f;
	}
};
