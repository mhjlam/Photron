#pragma once

#include "../structs/cuboid.hpp"
#include "../structs/point3.hpp"
#include "../structs/ray.hpp"
#include "../structs/triangle.hpp"
#include "../structs/vector3.hpp"

#include <list>
#include <ostream>
#include <sstream>
#include <string>
#include <vector>

/***********************************************************
 * FUNCTION TEMPLATES
 ***********************************************************/

/***********************************************************
 * Returns the sum of the given vector of numerical values.
 ***********************************************************/
template<typename T>
double sum(std::vector<T> vals) {
	if (!vals[0]) {
		return 0.0;
	}

	if (std::isnan(vals[0])) {
		return 0.0;
	}

	double sum = 0.0;

	for (size_t i = 0; i < vals.size(); ++i) {
		sum += vals[i];
	}

	return sum;
}

/***********************************************************
 * Converts a string to a requested numerical value.
 ***********************************************************/
template<typename T>
T str2num(const std::string& str) {
	T num;
	std::stringstream ss(str);
	ss >> num;

	// isnan
	if (num != num) {
		return 0;
	}
	return num ? num : 0;
}

/***********************************************************
 * Converts a string to a requested numerical value,
 * returning a default return value on error.
 ***********************************************************/
template<typename T>
T str2num(const std::string& str, T def) {
	T num;
	std::stringstream ss(str);
	ss >> num;

	// isnan
	if (num != num) {
		return def;
	}
	return num ? num : def;
}

/***********************************************************
 * ARITHMETICS
 ***********************************************************/

double sq(double n);
double distribution_(Point3& p, Point3& q);

#ifndef _WIN32
double round(double n, int precision);
#endif

double rad(double deg);
double deg(double rad);

/***********************************************************
 * COMPARISONS
 ***********************************************************/

bool isbetween(uint32_t x, uint32_t a, uint32_t b);
bool isbetween(double x, double a, double b);

double min3(double a, double b, double c);
double max3(double a, double b, double c);
double min4(double a, double b, double c, double d);
double max4(double a, double b, double c, double d);

/***********************************************************
 * STREAM MANIPULATION
 ***********************************************************/

std::ostream& tab(std::ostream& os);
bool equals(const std::string& s, const std::string& t);
void trim_comment(std::string& str);
void trim_spaces(std::string& str);

/***********************************************************
 * FILE EXTRACTION
 ***********************************************************/

std::vector<double> split(const std::string& str, char splitchar);
std::vector<std::pair<std::string, std::string> > parameter_values(std::list<std::string>& lines);

/***********************************************************
 * LINEAR ALGEBRA
 ***********************************************************/

double dot(Vector3& v, Vector3& w);
Vector3 cross(Vector3& v, Vector3& w);
Vector3 subtract(Point3& v, Point3& w);

/***********************************************************
 * INTERSECTION TESTS
 ***********************************************************/

bool ray_triangle_intersect(Ray& ray, Triangle& triangle, Point3& intersect);
bool ray_plane_intersect(Ray& ray, Vector3& normal, Point3& point, Point3& intersect);

std::pair<bool, Point3> ray_triangle_intersect(Ray& ray, Triangle& triangle);
void ray_triangles_intersections(Ray& ray, std::vector<Triangle>& triangles, std::vector<Point3>& intersections);

double first_ray_triangle_intersect(Ray& ray, std::vector<Triangle>& triangles, Point3& intersect, Triangle& triangle);
double first_ray_cuboid_intersect_internal(Ray& ray, Cuboid& cuboid, Point3& intersect, Vector3& normal);

/***********************************************************
 * GENERATING RANDOM NUMBERS
 ***********************************************************/

double rand3(int* idum);
double random_num(void);
