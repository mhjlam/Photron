#pragma once

struct Cuboid
{
	double minx, miny, minz;
	double maxx, maxy, maxz;
	
	
	Cuboid(double x1, double y1, double z1, double x2, double y2, double z2)
	{
		minx = x1;
		miny = y1;
		minz = z1;
		maxx = x2;
		maxy = y2;
		maxz = z2;
	}
};
