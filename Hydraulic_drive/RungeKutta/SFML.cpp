#include <iostream>
#include <SFML/Graphics.hpp>
#include<thread>
#include<time.h>
#include <ctime>
#include <ratio>
#include <chrono>
#include <fstream>
#include "RKSolver.h"

const double
G = 9.81,
PI = 3.1415926535;

// прямой цилиндр
double
V1 = 0.1,
S1 = 0.1
;
// обратный цилиндр
double
V2 = 0.05,
S2 = 0.1
;
// Свойства рабочей жидкости
double
ro = 860,
E = 1.56E6
;
// Харрактеристики дроссилей
double
mu = 0.62,
S = 0.01
;
// Харрактеристики системы
double
min_x = 0,
max_x = 1,
m = 2,
b_prop = 5,
p_s = 0,
p_i = 1000
;
// Управляющий сигнал
double
U = 1
;
// Геометрия
double
A[2]{ 1,-1 },
B[2]{ 0,0 },
H[2]
{
	A[0] - B[0],
	A[1] - B[1],
},
h = sqrt(H[0] * H[0] + H[1] * H[1]),
l = 0.5,
d = 1,
dh = 2,
dd = -3,
Ft[2]
{
	+0, 
	-0.1 * G
},
Fn
;
double
cosa = 0,
sina = 0;
;

#define sign(x) (x == 0 ? 0 : (x > 0 ? 1 : -1))

double F(double* var)
{
	double
		& t = var[0],
		& v = var[1],
		& x = var[2],
		& p1 = var[3],
		& p2 = var[4];



	double
		l0 = l + x;
	cosa = (h * h + d * d - l0 * l0) / (2 * h * d);
	sina = sqrt(1 - cosa * cosa);

	double
		E[2]
	{
		(H[0] * cosa - H[1] * sina) / h,
		(H[0] * sina + H[1] * cosa) / h
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
		E[0] * (d + dd) + E[1] * (dh),
		E[1] * (d + dd) - E[0] * (dh)
	};


	double
		ret = -(D[0] * Ft[1] - D[1] * Ft[0]) / (L[0] * N[1] - L[1] * N[0]);
	Fn = ret;
	return ret;
}

double DV(double* var)
{
	double
		& t = var[0],
		& v = var[1],
		& x = var[2],
		& p1 = var[3],
		& p2 = var[4];

	double ret = (S1 * p1 - S2 * p2 - b_prop * v + F(var)) / m;

	return ret;
}

double  DX(double* var)
{
	double
		& t = var[0],
		& v = var[1],
		& x = var[2],
		& p1 = var[3],
		& p2 = var[4];

	double ret = v;

	if (x >= max_x && v > 0)
		ret = 0;
	if (x <= min_x && v < 0)
		ret = 0;

	return ret;
}

double  DP1(double* var)
{
	double
		& t = var[0],
		& v = var[1],
		& x = var[2],
		& p1 = var[3],
		& p2 = var[4],
		Q1 = 0;

	double _S = S * fabs(U);

	if (U > 0)
	{
		Q1 = mu * _S * sign(p_i - p1) * sqrt(2 / ro * fabs(p_i - p1));
	}
	else
	{
		Q1 = -mu * _S * sign(p1 - p_s) * sqrt(2 / ro * fabs(p_s - p1));
	}

	double Vt1 = V1 + S1 * (x - min_x);

	double ret = E / Vt1 * (Q1 - S1 * v);

	return ret;
}

double  DP2(double* var)
{

	double
		& t = var[0],
		& v = var[1],
		& x = var[2],
		& p1 = var[3],
		& p2 = var[4],
		Q2 = 0;

	double _S = S * fabs(U);

	if (U > 0)
		Q2 = -mu * _S * sign(p2 - p_s) * sqrt(2 / ro * fabs(p_s - p2));
	else
		Q2 = mu * _S * sign(p_i - p2) * sqrt(2 / ro * fabs(p_i - p2));

	double Vt2 = V2 + S2 * (max_x - x);

	double ret = E / Vt2 * (Q2 + S2 * v);

	return ret;
}

void main()
{
	RKSolver<4> solver;

	solver.funcs[0] = DV;
	solver.funcs[1] = DX;
	solver.funcs[2] = DP1;
	solver.funcs[3] = DP2;
	solver.State[3] = (U > 0 ? p_i : p_s);
	solver.State[4] = (U < 0 ? p_i : p_s);
	solver.h = 0.002;
	F(solver.State);
	double ba = acos(cosa);
	std::ofstream Out("out.txt");
	for (int i = 0; i < 10000; i++)
	{
		F(solver.State);
		if(i % 10 == 0)
		Out
			<< int(solver.State[0] * 1000) << "\t"
			<< solver.State[1] << "\t"
			<< solver.State[2] << "\t"
			<< solver.State[3] * S1 << "\t"
			<< solver.State[4] * S2 << "\t"
			<< Fn << "\t"
			<< (acos(cosa) - ba) * 180 / PI << '\n';
		solver.Calc();
		
		double& x = solver.State[2];
		if (x > max_x)
			x = max_x;
		if (x < min_x)
			x = min_x;
	}

}

