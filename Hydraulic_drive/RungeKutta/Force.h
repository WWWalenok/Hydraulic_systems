#pragma once
#include <cmath>
#include <iomanip>
#include "RKSolver.h"

struct IForse
{
	virtual double F(double *, void *) = 0;
};

struct Forse_manipulator: IForse
{
	double
		A[2]{ 30, 10 },
		B[2]{ 0,0 },
		H[2]{ A[0] - B[0],A[1] - B[1] },
		h = sqrt(H[0] * H[0] + H[1] * H[1]),
		l = 25,
		d = 10,
		dh = 0,
		dd = 0,
		Ft[2]{ +0,-0 },
		Fn,
		cosa = 0,
		sina = 0;

	double F(double *var, void *_t)
	{
		double
			&t = var[0],
			&v = var[1],
			&x = var[2],
			&p1 = var[3],
			&p2 = var[4];

		//Dummy
		return 0;

		double
			l0 = l + x;
		cosa = (h * h + d * d - l0 * l0) / (2 * h * d);
		sina = sqrt(1 - cosa * cosa);

		double HN[2]
		{
			H[1],
			-H[0]
		};
		if(HN[1] < 0)
		{
			HN[0] *= -1;
			HN[1] *= -1;
		}

		double
			E[2]
		{
			(H[0] * cosa + HN[0] * sina) / h,
			(HN[1] * sina + H[1] * cosa) / h
		},
			L[2]
		{
			E[0] * d,
			E[1] * d
		},
			N[2]
		{
			(L[0] - A[0] + B[0]) / l0,
			(L[1] - A[1] - B[0]) / l0
		},
			D[2]
		{
			E[0] * (d + dd) + HN[0] * (dh),
			E[1] * (d + dd) + HN[1] * (dh)
		};


		double
			ret = (D[0] * Ft[1] - D[1] * Ft[0]) / (L[0] * N[1] - L[1] * N[0]);
		Fn = ret;
		return ret;
	}
};

