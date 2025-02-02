#include "Surface.h"
#include "ClosedIan.h"

const double DodecaElements[60][3][3] = {
	{
		{-0.25*(1+SQRT5), 0.5, -0.25*(1-SQRT5)},
		{0.5, -0.25*(1-SQRT5), 0.25*(1+SQRT5)},
		{-0.25*(1-SQRT5), 0.25*(1+SQRT5), -0.5}
	},
	{
		{0, 0, 1},
		{1, 0, 0},
		{0, 1, 0}
	},
	{
		{0.25*(1+SQRT5), -0.5, 0.25*(1-SQRT5)},
		{-0.5, 0.25*(1-SQRT5), -0.25*(1+SQRT5)},
		{-0.25*(1-SQRT5), 0.25*(1+SQRT5), -0.5}
	},
	{
		{1, 0, 0},
		{0, 1, 0},
		{0, 0, 1}
	},
	{
		{-0.25*(1-SQRT5), 0.25*(1+SQRT5), 0.5},
		{-0.25*(1+SQRT5), 0.5, 0.25*(1-SQRT5)},
		{-0.5, 0.25*(1-SQRT5), 0.25*(1+SQRT5)}
	},
	{
		{0, -1, 0},
		{0, 0, 1},
		{-1, 0, 0}
	},
	{
		{0.5, -0.25*(1-SQRT5), -0.25*(1+SQRT5)},
		{-0.25*(1-SQRT5), 0.25*(1+SQRT5), 0.5},
		{0.25*(1+SQRT5), -0.5, -0.25*(1-SQRT5)}
	},
	{
		{-0.25*(1-SQRT5), 0.25*(1+SQRT5), -0.5},
		{-0.25*(1+SQRT5), 0.5, -0.25*(1-SQRT5)},
		{0.5, -0.25*(1-SQRT5), 0.25*(1+SQRT5)}
	},
	{
		{0, 1, 0},
		{0, 0, 1},
		{1, 0, 0}
	},
	{
		{0.5, 0.25*(1-SQRT5), -0.25*(1+SQRT5)},
		{-0.25*(1-SQRT5), -0.25*(1+SQRT5), 0.5},
		{-0.25*(1+SQRT5), -0.5, 0.25*(1-SQRT5)}
	},
	{
		{-0.25*(1+SQRT5), 0.5, -0.25*(1-SQRT5)},
		{-0.5, 0.25*(1-SQRT5), -0.25*(1+SQRT5)},
		{0.25*(1-SQRT5), -0.25*(1+SQRT5), 0.5}
	},
	{
		{0.25*(1-SQRT5), -0.25*(1+SQRT5), -0.5},
		{0.25*(1+SQRT5), -0.5, -0.25*(1-SQRT5)},
		{-0.5, 0.25*(1-SQRT5), 0.25*(1+SQRT5)}
	},
	{
		{0.25*(1+SQRT5), 0.5, -0.25*(1-SQRT5)},
		{-0.5, -0.25*(1-SQRT5), 0.25*(1+SQRT5)},
		{-0.25*(1-SQRT5), -0.25*(1+SQRT5), 0.5}
	},
	{
		{0.25*(1+SQRT5), 0.5, -0.25*(1-SQRT5)},
		{0.5, 0.25*(1-SQRT5), -0.25*(1+SQRT5)},
		{0.25*(1-SQRT5), 0.25*(1+SQRT5), -0.5}
	},
	{
		{0.25*(1+SQRT5), 0.5, 0.25*(1-SQRT5)},
		{-0.5, -0.25*(1-SQRT5), -0.25*(1+SQRT5)},
		{0.25*(1-SQRT5), 0.25*(1+SQRT5), 0.5}
	},
	{
		{0.25*(1+SQRT5), -0.5, -0.25*(1-SQRT5)},
		{-0.5, 0.25*(1-SQRT5), 0.25*(1+SQRT5)},
		{0.25*(1-SQRT5), -0.25*(1+SQRT5), -0.5}
	},
	{
		{-0.5, 0.25*(1-SQRT5), 0.25*(1+SQRT5)},
		{-0.25*(1-SQRT5), 0.25*(1+SQRT5), 0.5},
		{-0.25*(1+SQRT5), 0.5, 0.25*(1-SQRT5)}
	},
	{
		{0.25*(1-SQRT5), 0.25*(1+SQRT5), -0.5},
		{-0.25*(1+SQRT5), -0.5, 0.25*(1-SQRT5)},
		{-0.5, -0.25*(1-SQRT5), 0.25*(1+SQRT5)}
	},
	{
		{-0.25*(1+SQRT5), 0.5, 0.25*(1-SQRT5)},
		{0.5, -0.25*(1-SQRT5), -0.25*(1+SQRT5)},
		{0.25*(1-SQRT5), -0.25*(1+SQRT5), -0.5}
	},
	{
		{-0.25*(1-SQRT5), -0.25*(1+SQRT5), 0.5},
		{0.25*(1+SQRT5), 0.5, -0.25*(1-SQRT5)},
		{-0.5, -0.25*(1-SQRT5), 0.25*(1+SQRT5)}
	},
	{
		{-0.25*(1-SQRT5), 0.25*(1+SQRT5), -0.5},
		{0.25*(1+SQRT5), -0.5, 0.25*(1-SQRT5)},
		{-0.5, 0.25*(1-SQRT5), -0.25*(1+SQRT5)}
	},
	{
		{0.25*(1+SQRT5), -0.5, -0.25*(1-SQRT5)},
		{0.5, -0.25*(1-SQRT5), -0.25*(1+SQRT5)},
		{-0.25*(1-SQRT5), 0.25*(1+SQRT5), 0.5}
	},
	{
		{0.5, -0.25*(1-SQRT5), 0.25*(1+SQRT5)},
		{-0.25*(1-SQRT5), 0.25*(1+SQRT5), -0.5},
		{-0.25*(1+SQRT5), 0.5, -0.25*(1-SQRT5)}
	},
	{
		{-0.5, -0.25*(1-SQRT5), 0.25*(1+SQRT5)},
		{0.25*(1-SQRT5), 0.25*(1+SQRT5), -0.5},
		{-0.25*(1+SQRT5), -0.5, 0.25*(1-SQRT5)}
	},
	{
		{-0.5, -0.25*(1-SQRT5), -0.25*(1+SQRT5)},
		{-0.25*(1-SQRT5), -0.25*(1+SQRT5), -0.5},
		{-0.25*(1+SQRT5), -0.5, -0.25*(1-SQRT5)}
	},
	{
		{0.5, 0.25*(1-SQRT5), -0.25*(1+SQRT5)},
		{0.25*(1-SQRT5), 0.25*(1+SQRT5), -0.5},
		{0.25*(1+SQRT5), 0.5, -0.25*(1-SQRT5)}
	},
	{
		{-0.5, 0.25*(1-SQRT5), -0.25*(1+SQRT5)},
		{-0.25*(1-SQRT5), 0.25*(1+SQRT5), -0.5},
		{0.25*(1+SQRT5), -0.5, 0.25*(1-SQRT5)}
	},
	{
		{-0.25*(1-SQRT5), -0.25*(1+SQRT5), 0.5},
		{-0.25*(1+SQRT5), -0.5, 0.25*(1-SQRT5)},
		{0.5, 0.25*(1-SQRT5), -0.25*(1+SQRT5)}
	},
	{
		{0, 0, 1},
		{-1, 0, 0},
		{0, -1, 0}
	},
	{
		{-0.25*(1-SQRT5), -0.25*(1+SQRT5), -0.5},
		{0.25*(1+SQRT5), 0.5, 0.25*(1-SQRT5)},
		{0.5, 0.25*(1-SQRT5), 0.25*(1+SQRT5)}
	},
	{
		{0.25*(1-SQRT5), 0.25*(1+SQRT5), 0.5},
		{-0.25*(1+SQRT5), -0.5, -0.25*(1-SQRT5)},
		{0.5, 0.25*(1-SQRT5), 0.25*(1+SQRT5)}
	},
	{
		{0.5, 0.25*(1-SQRT5), 0.25*(1+SQRT5)},
		{0.25*(1-SQRT5), 0.25*(1+SQRT5), 0.5},
		{-0.25*(1+SQRT5), -0.5, -0.25*(1-SQRT5)}
	},
	{
		{-0.25*(1+SQRT5), -0.5, 0.25*(1-SQRT5)},
		{0.5, 0.25*(1-SQRT5), -0.25*(1+SQRT5)},
		{-0.25*(1-SQRT5), -0.25*(1+SQRT5), 0.5}
	},
	{
		{-0.5, 0.25*(1-SQRT5), 0.25*(1+SQRT5)},
		{0.25*(1-SQRT5), -0.25*(1+SQRT5), -0.5},
		{0.25*(1+SQRT5), -0.5, -0.25*(1-SQRT5)}
	},
	{
		{0.25*(1-SQRT5), -0.25*(1+SQRT5), 0.5},
		{-0.25*(1+SQRT5), 0.5, -0.25*(1-SQRT5)},
		{-0.5, 0.25*(1-SQRT5), -0.25*(1+SQRT5)}
	},
	{
		{0.25*(1+SQRT5), -0.5, 0.25*(1-SQRT5)},
		{0.5, -0.25*(1-SQRT5), 0.25*(1+SQRT5)},
		{0.25*(1-SQRT5), -0.25*(1+SQRT5), 0.5}
	},
	{
		{-0.25*(1+SQRT5), -0.5, -0.25*(1-SQRT5)},
		{-0.5, -0.25*(1-SQRT5), -0.25*(1+SQRT5)},
		{-0.25*(1-SQRT5), -0.25*(1+SQRT5), -0.5}
	},
	{
		{0.25*(1-SQRT5), 0.25*(1+SQRT5), -0.5},
		{0.25*(1+SQRT5), 0.5, -0.25*(1-SQRT5)},
		{0.5, 0.25*(1-SQRT5), -0.25*(1+SQRT5)}
	},
	{
		{0.25*(1-SQRT5), -0.25*(1+SQRT5), -0.5},
		{-0.25*(1+SQRT5), 0.5, 0.25*(1-SQRT5)},
		{0.5, -0.25*(1-SQRT5), -0.25*(1+SQRT5)}
	},
	{
		{-0.5, -0.25*(1-SQRT5), -0.25*(1+SQRT5)},
		{0.25*(1-SQRT5), 0.25*(1+SQRT5), 0.5},
		{0.25*(1+SQRT5), 0.5, 0.25*(1-SQRT5)}
	},
	{
		{-0.25*(1+SQRT5), 0.5, 0.25*(1-SQRT5)},
		{-0.5, 0.25*(1-SQRT5), 0.25*(1+SQRT5)},
		{-0.25*(1-SQRT5), 0.25*(1+SQRT5), 0.5}
	},
	{
		{-0.5, -0.25*(1-SQRT5), 0.25*(1+SQRT5)},
		{-0.25*(1-SQRT5), -0.25*(1+SQRT5), 0.5},
		{0.25*(1+SQRT5), 0.5, -0.25*(1-SQRT5)}
	},
	{
		{-0.25*(1-SQRT5), 0.25*(1+SQRT5), 0.5},
		{0.25*(1+SQRT5), -0.5, -0.25*(1-SQRT5)},
		{0.5, -0.25*(1-SQRT5), -0.25*(1+SQRT5)}
	},
	{
		{0.5, -0.25*(1-SQRT5), 0.25*(1+SQRT5)},
		{0.25*(1-SQRT5), -0.25*(1+SQRT5), 0.5},
		{0.25*(1+SQRT5), -0.5, 0.25*(1-SQRT5)}
	},
	{
		{0, 0, -1},
		{-1, 0, 0},
		{0, 1, 0}
	},
	{
		{0.25*(1-SQRT5), 0.25*(1+SQRT5), 0.5},
		{0.25*(1+SQRT5), 0.5, 0.25*(1-SQRT5)},
		{-0.5, -0.25*(1-SQRT5), -0.25*(1+SQRT5)}
	},
	{
		{1, 0, 0},
		{0, -1, 0},
		{0, 0, -1}
	},
	{
		{0, -1, 0},
		{0, 0, -1},
		{1, 0, 0}
	},
	{
		{0.5, 0.25*(1-SQRT5), 0.25*(1+SQRT5)},
		{-0.25*(1-SQRT5), -0.25*(1+SQRT5), -0.5},
		{0.25*(1+SQRT5), 0.5, 0.25*(1-SQRT5)}
	},
	{
		{-0.25*(1+SQRT5), -0.5, -0.25*(1-SQRT5)},
		{0.5, 0.25*(1-SQRT5), 0.25*(1+SQRT5)},
		{0.25*(1-SQRT5), 0.25*(1+SQRT5), 0.5}
	},
	{
		{0.25*(1-SQRT5), -0.25*(1+SQRT5), 0.5},
		{0.25*(1+SQRT5), -0.5, 0.25*(1-SQRT5)},
		{0.5, -0.25*(1-SQRT5), 0.25*(1+SQRT5)}
	},
	{
		{0, 1, 0},
		{0, 0, -1},
		{-1, 0, 0}
	},
	{
		{-0.25*(1+SQRT5), -0.5, 0.25*(1-SQRT5)},
		{-0.5, -0.25*(1-SQRT5), 0.25*(1+SQRT5)},
		{0.25*(1-SQRT5), 0.25*(1+SQRT5), -0.5}
	},
	{
		{-1, 0, 0},
		{0, 1, 0},
		{0, 0, -1}
	},
	{
		{-1, 0, 0},
		{0, -1, 0},
		{0, 0, 1}
	},
	{
		{-0.5, 0.25*(1-SQRT5), -0.25*(1+SQRT5)},
		{0.25*(1-SQRT5), -0.25*(1+SQRT5), 0.5},
		{-0.25*(1+SQRT5), 0.5, -0.25*(1-SQRT5)}
	},
	{
		{-0.25*(1-SQRT5), -0.25*(1+SQRT5), -0.5},
		{-0.25*(1+SQRT5), -0.5, -0.25*(1-SQRT5)},
		{-0.5, -0.25*(1-SQRT5), -0.25*(1+SQRT5)}
	},
	{
		{0, 0, -1},
		{1, 0, 0},
		{0, -1, 0}
	},
	{
		{0.25*(1+SQRT5), 0.5, 0.25*(1-SQRT5)},
		{0.5, 0.25*(1-SQRT5), 0.25*(1+SQRT5)},
		{-0.25*(1-SQRT5), -0.25*(1+SQRT5), -0.5}
	},
	{
		{0.5, -0.25*(1-SQRT5), -0.25*(1+SQRT5)},
		{0.25*(1-SQRT5), -0.25*(1+SQRT5), -0.5},
		{-0.25*(1+SQRT5), 0.5, 0.25*(1-SQRT5)}
	}
};
