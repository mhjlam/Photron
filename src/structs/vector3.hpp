#pragma once

#include "point3.hpp"

#include <cmath>

struct Vector3 {
	double x, y, z;

	Vector3() { x = y = z = 0; }

	Vector3(double xyz) { x = y = z = xyz; }

	Vector3(double vx, double vy, double vz) {
		x = vx;
		y = vy;
		z = vz;
	}

	Vector3(double vx, double vy, double vz, bool norm) {
		x = vx;
		y = vy;
		z = vz;

		if (norm) {
			normalize();
		}
	}

	Vector3(Point3& point, bool norm = false) {
		x = point.x;
		y = point.y;
		z = point.z;

		if (norm) {
			normalize();
		}
	}

	Vector3(const Vector3& vec, bool norm = false) {
		x = vec.x;
		y = vec.y;
		z = vec.z;

		if (norm) {
			normalize();
		}
	}

	// Assignment operator to fix deprecated copy warnings
	Vector3& operator=(const Vector3& vec) {
		if (this != &vec) {
			x = vec.x;
			y = vec.y;
			z = vec.z;
		}
		return *this;
	}

	Vector3(Point3 origin, Point3 target, bool norm = false) {
		x = target.x - origin.x;
		y = target.y - origin.y;
		z = target.z - origin.z;

		if (norm) {
			normalize();
		}
	}

	double length() { return std::sqrt(x * x + y * y + z * z); }

	void normalize() {
		double l = length();

		x /= l;
		y /= l;
		z /= l;
	}

	void reverse() {
		x = -x;
		y = -y;
		z = -z;
	}

	double dot(const Vector3& w) { return (x * w.x + y * w.y + z * w.z); }
	Vector3 cross(const Vector3& w) { return Vector3(y * w.z - z * w.y, z * w.x - x * w.z, x * w.y - y * w.x); }

	Vector3 operator+(const Vector3& w) { return Vector3(x + w.x, y + w.y, z + w.z); }
	Vector3 operator*(const double& d) { return Vector3(x * d, y * d, z * d); }
	Vector3 operator*(const Vector3& w) { return Vector3(x * w.x, y * w.y, z * w.z); }
	bool operator==(const Vector3& w) const { return (x == w.x && y == w.y && z == w.z); }
};
