#include "utilities.hpp"

#include <algorithm>
#include <cfloat>
#include <ctime>
#include <iostream>

#define PI (3.1415926535897932384626433832795)

typedef unsigned int uint;
typedef unsigned short ushort;

using namespace std;

using namespace std;

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
double dist(Point3& p, Point3& q) {
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
	return deg * PI / 180.0;
}

/***********************************************************
 * Returns the given angular radian in degrees.
 ***********************************************************/
double deg(double rad) {
	return rad * 180.0 / PI;
}

/***********************************************************
 * COMPARISONS
 ***********************************************************/

/***********************************************************
 * Return true if x is between a and b, inclusive.
 ***********************************************************/
bool isbetween(uint x, uint a, uint b) {
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
	return min(a, min(b, c));
}

/***********************************************************
 * Return the maximum value of three numbers.
 ***********************************************************/
double max3(double a, double b, double c) {
	return max(a, max(b, c));
}

/***********************************************************
 * Return the minimum value of four numbers.
 ***********************************************************/
double min4(double a, double b, double c, double d) {
	return min(a, min(b, min(c, d)));
}

/***********************************************************
 * Return the maximum value of four numbers.
 ***********************************************************/
double max4(double a, double b, double c, double d) {
	return max(a, max(b, max(c, d)));
}

/***********************************************************
 * STRING & STREAM MANIPULATION
 ***********************************************************/

/***********************************************************
 * Add a tab character to the stream.
 ***********************************************************/
ostream& tab(ostream& os) {
	return os << '\t';
}

/***********************************************************
 * Return true if the two strings are equal.
 ***********************************************************/
bool equals(const string& s, const string& t) {
	return (s.compare(t) == 0);
}

/***********************************************************
 * Remove all characters after the first occurrence of '#'
 * in a string.
 ***********************************************************/
void trimcomment(string& str) {
	size_t pos = str.find_first_of('#');

	if (pos != str.npos) {
		str.erase(pos);
	}
}

/***********************************************************
 * Remove all whitespace characters from a string.
 ***********************************************************/
void trimspaces(string& str) {
	char const* delims = " \t\r\n"; // space, tab, return, newline

	size_t pos = str.find_first_not_of(delims);
	if (pos == string::npos) {
		str = ""; // String is all whitespace
		return;
	}
	
	str.erase(0, pos);

	pos = str.find_last_not_of(delims);
	if (pos != string::npos) {
		str.erase(pos + 1);
	}
}

/***********************************************************
 * FILE EXTRACTION
 ***********************************************************/

/***********************************************************
 * Split a string on a character and returns the individual
 * numerical elements.
 ***********************************************************/
vector<double> split(const string& str, char splitchar) {
	vector<double> dblarray;
	dblarray.clear();

	size_t pos = str.find(splitchar);
	size_t beg = 0;

	// decompose statement
	while (pos != string::npos) {
		string piece = str.substr(beg, pos - beg);
		trimspaces(piece);
		if (!piece.empty()) {
			double val = str2num<double>(piece);
			dblarray.push_back(val);
		}
		beg = pos + 1;
		pos = str.find(splitchar, beg);
	}

	// add the last piece if there's anything left
	if (beg < str.size()) {
		string piece = str.substr(beg);
		trimspaces(piece);
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
vector<std::pair<std::string, std::string> > getparvals(list<string>& lines) {
	vector<std::pair<std::string, std::string> > parvals;

	for (list<string>::iterator it = lines.begin(); it != lines.end(); ++it) {
		string line = (*it);
		size_t delim = line.find_first_of("=");

		// parameter name
		string param = line.substr(0, delim);
		trimspaces(param);

		// parameter value
		string value = line.substr(delim + 1);
		trimcomment(value);
		trimspaces(value);

		parvals.push_back(std::make_pair(param, value));
	}

	return parvals;
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
bool ray_triangle_intersect(Ray& ray, Triangle& triangle, Point3& inter) {
	// ray origin and direction
	Point3 orig = ray.origin;
	Vector3 dir = ray.direction;
	dir.normalize();

	// triangle vertices
	Point3 vert0 = triangle.vertex0;
	Point3 vert1 = triangle.vertex1;
	Point3 vert2 = triangle.vertex2;

	// vectors for two edges sharing vertex0
	Vector3 edge1 = subtract(vert1, vert0);
	Vector3 edge2 = subtract(vert2, vert0);

	// compute determinant
	Vector3 pvec = cross(dir, edge2);
	double det = dot(edge1, pvec);

	// if determinant is near zero, ray lies on triangle plane
	if (det > -0.000001 && det < 0.000001) {
		return false;
	}

	// inversed determinant
	double detinv = 1.0 / det;

	// vector from vert0 to ray origin
	Vector3 tvec = subtract(orig, vert0);

	// compute u and test bounds
	double u = dot(tvec, pvec) * detinv;

	// invalid barycentric coordinate
	if (u < 0.0 || u > 1.0) {
		return false;
	}

	// compute v and test bounds
	Vector3 qvec = cross(tvec, edge1);
	double v = dot(dir, qvec) * detinv;

	// check if the ray intersects the triangle
	if (v < 0.0 || (u + v) > 1.0) {
		return false;
	}

	// w coordinate (unused)
	// double w = 1 - (u + v);

	// the parametric location of the intersection
	double t = dot(edge2, qvec) * detinv;

	if (t < 0) {
		return false;
	}

	// (x,y,z) = vert0*u + vert1*v + vert2*w
	// inter.x = vert0.x * u + vert1.x * v + vert2.x * w;
	// inter.y = vert0.y * u + vert1.y * v + vert2.y * w;
	// inter.z = vert0.z * u + vert1.z * v + vert2.z * w;

	// (x,y,z) = orig + (t * dir)
	inter.x = orig.x + (t * dir.x);
	inter.y = orig.y + (t * dir.y);
	inter.z = orig.z + (t * dir.z);

	// round off at the origin
	if (abs(inter.x) <= 1E-10) {
		inter.x = 0;
	}
	if (abs(inter.y) <= 1E-10) {
		inter.y = 0;
	}
	if (abs(inter.z) <= 1E-10) {
		inter.z = 0;
	}

	return true;
}

/***********************************************************
 * Return ray-triangle intersection point as a pair<bool, Point3>,
 * where the boolean indicates whether a collision was detected,
 * and if so, Point3 is the intersection point.
 ***********************************************************/
pair<bool, Point3> ray_triangle_intersect(Ray& ray, Triangle& triangle) {
	pair<bool, Point3> out;

	Point3 inter;
	out.first = ray_triangle_intersect(ray, triangle, inter);
	out.second = inter;

	return out;
}

/***********************************************************
 * Compute all intersections of a ray with a given list of
 * triangles and initialize a reference to a vector of
 * intersections.
 ***********************************************************/
void ray_triangles_intersections(Ray& ray, vector<Triangle>& triangles, vector<Point3>& intersections) {
	for (vector<Triangle>::iterator t = triangles.begin(); t != triangles.end(); ++t) {
		pair<bool, Point3> inter = ray_triangle_intersect(ray, *t);
		if (inter.first) {
			intersections.push_back(inter.second);
		}
	}
}

/***********************************************************
 * Compute the closest ray-triangle intersection point from
 * a vector of triangles.
 *
 * Returns the distance to the first intersection, or
 * DBL_MAX if there is no intersection.
 ***********************************************************/
double first_ray_triangle_intersect(Ray& ray, vector<Triangle>& triangles, Point3& inter, Triangle& triangle) {
	ushort numhit = 0;
	double mindist = DBL_MAX;
	Point3 origin = ray.origin;

	for (vector<Triangle>::iterator t = triangles.begin(); t != triangles.end(); ++t) {
		Point3 isect;
		Triangle tri = *t;

		if (ray_triangle_intersect(ray, tri, isect)) {
			numhit++;
			double distance = dist(origin, isect);

			if (distance < mindist) {
				mindist = distance;
				inter = isect;
				triangle = tri;
			}
		}
	}

	// grazing hits are not good enough
	if (numhit < 2) {
		return DBL_MAX;
	}

	return mindist;
}

/***********************************************************
 * Computes the intersection point of a given ray and a
 * given plane, described by its normal vector and a point
 * on the plane.
 *
 * Returns false if there is no intersection point,
 * or true otherwise.
 ***********************************************************/
bool ray_plane_intersect(Ray& ray, Vector3& normal, Point3& point, Point3& inter) {
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

	Vector3 dir = ray.direction;
	Vector3 pos = Vector3(point.x - ray.origin.x, point.y - ray.origin.y, point.z - ray.origin.z);

	if (abs(pos.x) <= 1E-10) {
		pos.x = 0;
	}
	if (abs(pos.y) <= 1E-10) {
		pos.y = 0;
	}
	if (abs(pos.z) <= 1E-10) {
		pos.z = 0;
	}

	double DoN = dot(normal, dir);

	if (DoN == 0) {
		inter = ray.origin;
		return false;
	}

	double OoN = dot(normal, pos);
	double t = OoN / DoN;

	if (t < -0) {
		return false;
	}

	// set hitpoint coordinates
	inter.x = ray.origin.x + (t * dir.x);
	inter.y = ray.origin.y + (t * dir.y);
	inter.z = ray.origin.z + (t * dir.z);

	// round off if values are close to zero
	if (abs(inter.x) <= 1E-10)
		inter.x = 0;
	if (abs(inter.y) <= 1E-10)
		inter.y = 0;
	if (abs(inter.z) <= 1E-10)
		inter.z = 0;

	return true;
}

/***********************************************************
 * Compute the coordinates and plane normal of the closest
 * intersection of a given ray and the faces of an axis
 * aligned box, which is described by two vertices.
 *
 * Returns the distance to the first intersection, or
 * DBL_MAX if no intersection was found (which should never
 * happen if the given ray originates from inside the given
 * box).
 *
 * References to the intersection point and to the
 * intersected plane normal are initialized.
 ***********************************************************/
double first_ray_cuboid_intersect_internal(Ray& ray, Cuboid& cuboid, Point3& inter, Vector3& normal) {
	double mindist = DBL_MAX;

	// possible intersection points
	Point3 inter1, inter2, inter3, inter4, inter5, inter6;

	// voxel plane normals are pointed inwards for internal intersection test
	Vector3 normal1 = Vector3(1, 0, 0);        // left face plane (inwards)
	Vector3 normal2 = Vector3(-1, 0, 0);       // right face plane (inwards)
	Vector3 normal3 = Vector3(0, -1, 0);       // top face plane (inwards)
	Vector3 normal4 = Vector3(0, 1, 0);        // bottom face plane (inwards)
	Vector3 normal5 = Vector3(0, 0, -1);       // front face plane (inwards)
	Vector3 normal6 = Vector3(0, 0, 1);        // rear face plane (inwards)

	Point3 point1 = Point3(cuboid.minx, 0, 0); // left face plane point
	Point3 point2 = Point3(cuboid.maxx, 0, 0); // right face plane point
	Point3 point3 = Point3(0, cuboid.miny, 0); // top face plane point
	Point3 point4 = Point3(0, cuboid.maxy, 0); // bottom face plane point
	Point3 point5 = Point3(0, 0, cuboid.minz); // front face plane point
	Point3 point6 = Point3(0, 0, cuboid.maxz); // rear face plane point

	// compute individual ray-plane intersections
	bool t1 = ray_plane_intersect(ray, normal1, point1, inter1); // left face
	bool t2 = ray_plane_intersect(ray, normal2, point2, inter2); // right face
	bool t3 = ray_plane_intersect(ray, normal3, point3, inter3); // top face
	bool t4 = ray_plane_intersect(ray, normal4, point4, inter4); // bottom face
	bool t5 = ray_plane_intersect(ray, normal5, point5, inter5); // front face
	bool t6 = ray_plane_intersect(ray, normal6, point6, inter6); // rear face

	// temporarily store ray-plane intersection tests
	vector<pair<Point3, Vector3> > intersections; // intersection point / normal vector
	vector<pair<Point3, Vector3> >::iterator iit; // intersection iterator

	if (t1) {
		intersections.push_back(std::make_pair(inter1, normal1));
	}
	if (t2) {
		intersections.push_back(std::make_pair(inter2, normal2));
	}
	if (t3) {
		intersections.push_back(std::make_pair(inter3, normal3));
	}
	if (t4) {
		intersections.push_back(std::make_pair(inter4, normal4));
	}
	if (t5) {
		intersections.push_back(std::make_pair(inter5, normal5));
	}
	if (t6) {
		intersections.push_back(std::make_pair(inter6, normal6));
	}

	// find first intersection
	for (iit = intersections.begin(); iit != intersections.end(); ++iit) {
		double distance = dist(ray.origin, iit->first);

		if (distance != 0 && distance < mindist) {
			mindist = distance;
			inter = iit->first;
			normal = iit->second;
		}
	}

	// first intersection found
	if (mindist != DBL_MAX) {
		return mindist;
	}

	// otherwise, look for a zero-distance intersection
	for (iit = intersections.begin(); iit != intersections.end(); ++iit) {
		if (ray.origin == iit->first) {
			mindist = 0;
			inter = iit->first;
			normal = iit->second;
			break; // stop on the first find
		}
	}

	// either 0 or DBL_MAX
	return mindist;
}

/***********************************************************
 * RANDOMIZATION
 ***********************************************************/

/***********************************************************
 * Return a random double value that is between 0.0 and 1.0.
 ***********************************************************/
double randomnum(void) {
	static int idum;
	static bool first_time = true;

	if (first_time) {
		// use 16-bit integer as the seed
		idum = -(int)time(NULL) % (1 << 15);
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
