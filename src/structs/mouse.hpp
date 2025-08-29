#pragma once

struct Mouse
{
	int x;
	int y;
	
	unsigned int lmb;
	unsigned int rmb;
	
	Mouse()
	{
		x = y = 0;
		lmb = rmb = 1;
	}
	
	void Update(int nx, int ny)
	{
		x = nx;
		y = ny;
	}
};
