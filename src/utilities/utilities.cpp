#include "utilities.hpp"

#include <algorithm>
#include <cfloat>
#include <cstdint>
#include <ctime>
#include <iostream>
#include <limits>
#include <numbers>

/***********************************************************
 * ARITHMETICS
 ***********************************************************/

/***********************************************************
 * Returns the perfect square of a number.
 ***********************************************************/
double sq(double n) {
	return n * n;
}

/***********************************************************
 * Returns the distance between two three dimensional points.
 ***********************************************************/
double distribution_(Point3& p, Point3& q) {
	double xd = q.x - p.x;
	double yd = q.y - p.y;
	double zd = q.z - p.z;

	return sqrt(xd * xd + yd * yd + zd * zd);
}

/***********************************************************
 * Rounds a floating-point number to a given precision.
 ***********************************************************/
#ifndef _WIN32
double round(double n, int precision) {
	if (n < 0.000001 && n > -0.000001)
		return 0;

	double places = pow(10.0, precision);
	return round(n * places) / places;
}
#endif

/***********************************************************
 * Returns the given angular degrees in radians.
 ***********************************************************/
double rad(double deg) {
	return deg * std::numbers::pi / 180.0;
}

/***********************************************************
 * Returns the given angular radian in degrees.
 ***********************************************************/
double deg(double rad) {
	return rad * 180.0 / std::numbers::pi;
}

/***********************************************************
 * COMPARISONS
 ***********************************************************/

/***********************************************************
 * Return true if x is between a and b, inclusive.
 ***********************************************************/
bool isbetween(uint32_t x, uint32_t a, uint32_t b) {
	return (x >= a && x <= b);
}

/***********************************************************
 * Return true if x is between a and b, inclusive.
 ***********************************************************/
bool isbetween(double x, double a, double b) {
	return (x >= a && x <= b);
}

/***********************************************************
 * Return the minimum value of three numbers.
 ***********************************************************/
double min3(double a, double b, double c) {
	return std::min(a, std::min(b, c));
}

/***********************************************************
 * Return the maximum value of three numbers.
 ***********************************************************/
double max3(double a, double b, double c) {
	return std::max(a, std::max(b, c));
}

/***********************************************************
 * Return the minimum value of four numbers.
 ***********************************************************/
double min4(double a, double b, double c, double d) {
	return std::min(a, std::min(b, std::min(c, d)));
}

/***********************************************************
 * Return the maximum value of four numbers.
 ***********************************************************/
double max4(double a, double b, double c, double d) {
	return std::max(a, std::max(b, std::max(c, d)));
}

/***********************************************************
 * STRING & STREAM MANIPULATION
 ***********************************************************/

/***********************************************************
 * Add a tab character to the stream.
 ***********************************************************/
std::ostream& tab(std::ostream& os) {
	return os << '\t';
}

/***********************************************************
 * Return true if the two strings are equal.
 ***********************************************************/
bool equals(const std::string& s, const std::string& t) {
	return (s.compare(t) == 0);
}

/***********************************************************
 * Remove all characters after the first occurrence of '#'
 * in a string.
 ***********************************************************/
void trim_comment(std::string& str) {
	size_t position = str.find_first_of('#');

	if (position != str.npos) {
		str.erase(position);
	}
}

/***********************************************************
 * Remove all whitespace characters from a string.
 ***********************************************************/
void trim_spaces(std::string& str) {
	char const* delims = " \t\r\n"; // space, tab, return, newline

	size_t position = str.find_first_not_of(delims);
	if (position == std::string::npos) {
		str = "";                   // String is all whitespace
		return;
	}

	str.erase(0, position);

	position = str.find_last_not_of(delims);
	if (position != std::string::npos) {
		str.erase(position + 1);
	}
}

/***********************************************************
 * FILE EXTRACTION
 ***********************************************************/

/***********************************************************
 * Split a string on a character and returns the individual
 * numerical elements.
 ***********************************************************/
std::vector<double> split(const std::string& str, char split_char) {
	std::vector<double> dblarray;
	dblarray.clear();

	size_t position = str.find(split_char);
	size_t beg = 0;

	// decompose statement
	while (position != std::string::npos) {
		std::string piece = str.substr(beg, position - beg);
		trim_spaces(piece);
		if (!piece.empty()) {
			double val = str2num<double>(piece);
			dblarray.push_back(val);
		}
		beg = position + 1;
		position = str.find(split_char, beg);
	}

	// add the last piece if there's anything left
	if (beg < str.size()) {
		std::string piece = str.substr(beg);
		trim_spaces(piece);
		if (!piece.empty()) {
			double val = str2num<double>(piece);
			dblarray.push_back(val);
		}
	}

	return dblarray;
}

/***********************************************************
 * Extract param/value pairs from a list of strings.
 ***********************************************************/
std::vector<std::pair<std::string, std::string>> parameter_values(std::list<std::string>& lines) {
	std::vector<std::pair<std::string, std::string>> param_values;

	for (const auto& line : lines) {
		size_t delim = line.find_first_of("=");

		// parameter name
		std::string param = line.substr(0, delim);
		trim_spaces(param);

		// parameter value
		std::string value = line.substr(delim + 1);
		trim_comment(value);
		trim_spaces(value);

		param_values.push_back(std::make_pair(param, value));
	}

	return param_values;
}

/***********************************************************
 * LINEAR ALGEBRA
 ***********************************************************/

/***********************************************************
 * Compute the inner product between two vectors.
 ***********************************************************/
double dot(Vector3& v, Vector3& w) {
	return (v.x * w.x + v.y * w.y + v.z * w.z);
}

/***********************************************************
 * Compute the cross product between two vectors.
 ***********************************************************/
Vector3 cross(Vector3& v, Vector3& w) {
	double x = v.y * w.z - v.z * w.y;
	double y = v.z * w.x - v.x * w.z;
	double z = v.x * w.y - v.y * w.x;

	return Vector3(x, y, z);
}

/***********************************************************
 * Compute the vector between a given start and end point.
 ***********************************************************/
Vector3 subtract(Point3& v, Point3& w) {
	return Vector3(v.x - w.x, v.y - w.y, v.z - w.z);
}

/***********************************************************
 * INTERSECTION TESTS
 ***********************************************************/

/***********************************************************
 * Compute the intersection point of a given ray and a
 * given triangle, described by its three vertices.
 *
 * Returns false if there is no intersection point,
 * or true otherwise.
 *
 * Ray equation: r(t) = p + td
 * Barycentric coordinates: (u, v, 1-(u+v))
 *
 * This ray-triangle intersection algorithm is an
 * implementation of Moller-Trumbore algorithm, available at
 * http://www.graphics.cornell.edu/pubs/1997/MT97.html
 ***********************************************************/
bool ray_triangle_intersect(Ray& ray, Triangle& triangle_out, Point3& intersect_out) {
	// ray origin and direction
	Point3 origin = ray.origin;
	Vector3 direction = ray.direction;
	direction.normalize();

	// triangle vertices
	Point3 vert0 = triangle_out.v0;
	Point3 vert1 = triangle_out.v1;
	Point3 vert2 = triangle_out.v2;

	// vectors for two edges sharing v0
	Vector3 edge1 = subtract(vert1, vert0);
	Vector3 edge2 = subtract(vert2, vert0);

	// compute determinant
	Vector3 pvec = cross(direction, edge2);
	double det = dot(edge1, pvec);

	// if determinant is near zero, ray lies on triangle plane
	if (det > -0.000001 && det < 0.000001) {
		return false;
	}

	// inversed determinant
	double det_inv = 1.0 / det;

	// vector from vert0 to ray origin
	Vector3 tvec = subtract(origin, vert0);

	// compute u and test bounds
	double u = dot(tvec, pvec) * det_inv;

	// invalid barycentric coordinate
	if (u < 0.0 || u > 1.0) {
		return false;
	}

	// compute v and test bounds
	Vector3 qvec = cross(tvec, edge1);
	double v = dot(direction, qvec) * det_inv;

	// check if the ray intersects the triangle
	if (v < 0.0 || (u + v) > 1.0) {
		return false;
	}

	// the parametric location of the intersection
	double t = dot(edge2, qvec) * det_inv;

	if (t < 0) {
		return false;
	}

	// (x,y,z) = origin + (t * direction)
	intersect_out.x = origin.x + (t * direction.x);
	intersect_out.y = origin.y + (t * direction.y);
	intersect_out.z = origin.z + (t * direction.z);

	// round off at the origin
	if (abs(intersect_out.x) <= 1E-10) {
		intersect_out.x = 0;
	}
	if (abs(intersect_out.y) <= 1E-10) {
		intersect_out.y = 0;
	}
	if (abs(intersect_out.z) <= 1E-10) {
		intersect_out.z = 0;
	}

	return true;
}

/***********************************************************
 * Return ray-triangle intersection point as a pair<bool, Point3>,
 * where the boolean indicates whether a collision was detected,
 * and if so, Point3 is the intersection point.
 ***********************************************************/
std::pair<bool, Point3> ray_triangle_intersect(Ray& ray, Triangle& triangle_out) {
	std::pair<bool, Point3> out;

	Point3 intersect_out;
	out.first = ray_triangle_intersect(ray, triangle_out, intersect_out);
	out.second = intersect_out;

	return out;
}

/***********************************************************
 * Compute all intersections of a ray with a given list of
 * triangles and initialize a reference to a vector of
 * intersections.
 ***********************************************************/
void ray_triangles_intersections(Ray& ray, std::vector<Triangle>& triangles, std::vector<Point3>& intersections) {
	for (auto& triangle_out : triangles) {
		std::pair<bool, Point3> intersect_out = ray_triangle_intersect(ray, triangle_out);
		if (intersect_out.first) {
			intersections.push_back(intersect_out.second);
		}
	}
}

/***********************************************************
 * Compute the closest ray-triangle intersection point from
 * a vector of triangles.
 *
 * Returns the distance to the first intersection, or
 * std::numeric_limits<double>::max() if there is no intersection.
 ***********************************************************/
double first_ray_triangle_intersect(Ray& ray, std::vector<Triangle>& triangles, Point3& intersect_out,
									Triangle& triangle_out) {
	uint8_t num_hit = 0;
	double min_dist = std::numeric_limits<double>::max();
	Point3 origin = ray.origin;

	for (auto& triangle : triangles) {
		Point3 intersect;
		if (ray_triangle_intersect(ray, triangle_out, intersect)) {
			num_hit++;
			double distance = distribution_(origin, intersect);

			if (distance < min_dist) {
				min_dist = distance;
				intersect_out = intersect;
				triangle_out = triangle;
			}
		}
	}

	// Grazing hits are not good enough
	if (num_hit < 2) {
		return std::numeric_limits<double>::max();
	}

	return min_dist;
}

/***********************************************************
 * Computes the intersection point of a given ray and a
 * given plane, described by its normal vector and a point
 * on the plane.
 *
 * Returns false if there is no intersection point,
 * or true otherwise.
 ***********************************************************/
bool ray_plane_intersect(Ray& ray, Vector3& normal, Point3& point, Point3& intersect_out) {
	/*
	 * Let r(t) = O + tD be the ray equation, where O is the ray origin and D the direction.
	 * Let N be the plane normal and let Q be a point on the plane.
	 * Then any point x on the plane satisfies the equation:  N . (Q - x) = 0
	 *
	 * Substitute r(t) for x and solve for t:
	 *  N . (Q - r(t)) = 0
	 *  N . (Q - O - tD) = 0
	 *  N . (Q - O) = (N . D) * t
	 *
	 *  t = [(Q - O) . N] / (D . N)
	 */

	Vector3 direction = ray.direction;
	Vector3 position = Vector3(point.x - ray.origin.x, point.y - ray.origin.y, point.z - ray.origin.z);

	if (abs(position.x) <= 1E-10) {
		position.x = 0;
	}
	if (abs(position.y) <= 1E-10) {
		position.y = 0;
	}
	if (abs(position.z) <= 1E-10) {
		position.z = 0;
	}

	double DoN = dot(normal, direction);

	if (DoN == 0) {
		intersect_out = ray.origin;
		return false;
	}

	double OoN = dot(normal, position);
	double t = OoN / DoN;

	if (t < -0) {
		return false;
	}

	// set hitpoint coordinates
	intersect_out.x = ray.origin.x + (t * direction.x);
	intersect_out.y = ray.origin.y + (t * direction.y);
	intersect_out.z = ray.origin.z + (t * direction.z);

	// round off if values are close to zero
	if (abs(intersect_out.x) <= 1E-10) {
		intersect_out.x = 0;
	}
	if (abs(intersect_out.y) <= 1E-10) {
		intersect_out.y = 0;
	}
	if (abs(intersect_out.z) <= 1E-10) {
		intersect_out.z = 0;
	}

	return true;
}

/***********************************************************
 * Compute the coordinates and plane normal of the closest intersection of a given ray and the
 * faces of an axis aligned box, which is described by two vertices.
 *
 * Returns the distance to the first intersection, or double max if no intersection was found.
 * (Which should never happen if the given ray originates from inside the given box).
 *
 * References to the intersection point and to the intersected plane normal are initialized.
 ***********************************************************/
double first_ray_cuboid_intersect_internal(Ray& ray, Cuboid& cuboid, Point3& intersect_out, Vector3& normal) {
	double min_dist = std::numeric_limits<double>::max();

	// possible intersection points
	Point3 inter1, inter2, inter3, inter4, inter5, inter6;

	// voxel plane normals are pointed inwards for internal intersection test
	Vector3 normal1 = Vector3(1, 0, 0);         // left face plane (inwards)
	Vector3 normal2 = Vector3(-1, 0, 0);        // right face plane (inwards)
	Vector3 normal3 = Vector3(0, -1, 0);        // top face plane (inwards)
	Vector3 normal4 = Vector3(0, 1, 0);         // bottom face plane (inwards)
	Vector3 normal5 = Vector3(0, 0, -1);        // front face plane (inwards)
	Vector3 normal6 = Vector3(0, 0, 1);         // rear face plane (inwards)

	Point3 point1 = Point3(cuboid.min_x, 0, 0); // left face plane point
	Point3 point2 = Point3(cuboid.max_x, 0, 0); // right face plane point
	Point3 point3 = Point3(0, cuboid.min_y, 0); // top face plane point
	Point3 point4 = Point3(0, cuboid.max_y, 0); // bottom face plane point
	Point3 point5 = Point3(0, 0, cuboid.min_z); // front face plane point
	Point3 point6 = Point3(0, 0, cuboid.max_z); // rear face plane point

	// temporarily store ray-plane intersection tests
	std::vector<std::pair<Point3, Vector3>> intersections; // intersection point / normal vector

	// compute individual ray-plane intersections
	// left face
	if (ray_plane_intersect(ray, normal1, point1, inter1)) {
		intersections.push_back(std::make_pair(inter1, normal1));
	}
	// right face
	if (ray_plane_intersect(ray, normal2, point2, inter2)) {
		intersections.push_back(std::make_pair(inter2, normal2));
	}
	// top face
	if (ray_plane_intersect(ray, normal3, point3, inter3)) {
		intersections.push_back(std::make_pair(inter3, normal3));
	}
	// bottom face
	if (ray_plane_intersect(ray, normal4, point4, inter4)) {
		intersections.push_back(std::make_pair(inter4, normal4));
	}
	// front face
	if (ray_plane_intersect(ray, normal5, point5, inter5)) {
		intersections.push_back(std::make_pair(inter5, normal5));
	}
	// rear face
	if (ray_plane_intersect(ray, normal6, point6, inter6)) {
		intersections.push_back(std::make_pair(inter6, normal6));
	}

	// find first intersection
	for (auto& intersection : intersections) {
		double distance = distribution_(ray.origin, intersection.first);

		if (distance != 0 && distance < min_dist) {
			min_dist = distance;
			intersect_out = intersection.first;
			normal = intersection.second;
		}
	}

	// first intersection found
	if (min_dist != std::numeric_limits<double>::max()) {
		return min_dist;
	}

	// otherwise, look for a zero-distance intersection
	for (auto& intersection : intersections) {
		if (ray.origin == intersection.first) {
			min_dist = 0;
			intersect_out = intersection.first;
			normal = intersection.second;
			break; // stop on the first find
		}
	}

	// either 0 or std::numeric_limits<double>::max()
	return min_dist;
}

/***********************************************************
 * RANDOMIZATION
 ***********************************************************/

/***********************************************************
 * Return a random double value that is between 0.0 and 1.0.
 ***********************************************************/
double random_num(void) {
	static int idum;
	static bool first_time = true;

	if (first_time) {
		// use 16-bit integer as the seed
		idum = -(int)time(nullptr) % (1 << 15);
		rand3(&idum);

		idum = 1;
		first_time = false;
	}

	return ((double)rand3(&idum));
}

/***********************************************************
 * Randomizer from "Numerical Recipes".
 ***********************************************************/
double rand3(int* idum) {
	static char MZ = 0;
	static long MSEED = 161803398;
	static long MBIG = 1000000000;
	static double FAC = 1.0E-9;

	static int inext, inextp;
	static long ma[56];
	static int iff = 0;
	long mj, mk;
	int i, ii, k;

	if (*idum < 0 || iff == 0) {
		iff = 1;
		mj = MSEED - (*idum < 0 ? -*idum : *idum);
		mj %= MBIG;
		ma[55] = mj;
		mk = 1;

		for (i = 1; i <= 54; i++) {
			ii = (21 * i) % 55;
			ma[ii] = mk;
			mk = mj - mk;

			if (mk < MZ) {
				mk += MBIG;
			}

			mj = ma[ii];
		}

		for (k = 1; k <= 4; k++) {
			for (i = 1; i <= 55; i++) {
				ma[i] -= ma[1 + (i + 30) % 55];

				if (ma[i] < MZ) {
					ma[i] += MBIG;
				}
			}
		}

		inext = 0;
		inextp = 31;
		*idum = 1;
	}

	if (++inext == 56) {
		inext = 1;
	}

	if (++inextp == 56) {
		inextp = 1;
	}

	mj = ma[inext] - ma[inextp];

	if (mj < MZ) {
		mj += MBIG;
	}

	ma[inext] = mj;
	return mj * FAC;
}
