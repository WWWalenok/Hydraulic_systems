#pragma once
#include <cmath>
#include <iomanip>
#include "RKSolver.h"
#include "Force.h"

struct PeredMotor
{

	double
		K = 1,
		T = 1;

	double
		*U = 0;
		

	static double DDX(double *var, void *_t)
	{
		double
			&t = var[0],
			&x = var[1],
			&dx = var[2];

		PeredMotor *T = (PeredMotor *)(_t);

		double *U = T->U;

		if(U == 0)
			return 0;

		double u_t = MAX2(-1, MIN2(1, *U));

		return (u_t - x) * T->K / T->T - dx / T->T;
	}

	static double DX(double *var, void *_t)
	{
		double
			&t = var[0],
			&x = var[1],
			&dx = var[2];

		return dx;

	}

	RKSolver<2> solver;
};




