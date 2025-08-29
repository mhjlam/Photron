#pragma once

struct Point2
{
	double x, y;
	
	Point2()
	{
		x = y = 0;
	}
	
	bool operator==(Point2& p)
	{
		return (x == p.x && y == p.y);
	}
};
