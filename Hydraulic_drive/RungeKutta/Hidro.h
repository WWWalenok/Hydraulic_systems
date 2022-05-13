#pragma once
#include <cmath>
#include <iomanip>
#include "RKSolver.h"
#include "Force.h"

#define sign(x) (x == 0 ? 0 : (x > 0 ? 1 : -1))


const double
G = 9.81,
PI = 3.1415926535;

static inline double U_to_S(double u)
{
	u = 2 * u - 1;
	u = MAX2(-1, MIN2(1, u));
	return
		1 - (acos(u) - u * sqrt(1 - u * u)) / PI;
}

static double U_to_S_dead_zone(double u, double dz)
{

	if(fabs(u) < dz)
	{
		u = 0;
	}
	else
	{
		u = u - (u > 0 ? dz : -dz);
		u /= 1 - dz;
	}

	return U_to_S(u);
}

inline float D2S(float D)
{
	return  std::powf(D * 0.5, 2) * 3.1415926535;
}


static double DV(double *var, void *_t);

static  double  DX(double *var, void *_t);

static  double DP1(double *var, void *_t);

static double DP2(double *var, void *_t);

struct H_System
{

	// прямой цилиндр
	double
		V1 = 0.01,
		S1 = D2S(100 * 1E-3);
	;
	// обратный цилиндр
	double
		V2 = 0.01,
		S2 = D2S(100 * 1E-3) - D2S(50 * 1E-3);
	;
	// Свойства рабочей жидкости
	double
		ro = 860,
		E = 1.56E7
		;
	// Харрактеристики дроссилей
	double
		mu = 0.62,
		S = D2S((14 - 6) * 1E-3)
		;
	// Харрактеристики системы
	double
		min_x = 0,
		max_x = 10,
		m = 10,
		b_prop = 5,
		f_tr_suh = 0,
		p_sliv = 0.1 * 1E6,
		p_input = 10 * 1E6
		;

	// Управляющий сигнал
	double
		U = 1,
		DU_max = 10,
		p1_sliv_stat[2]{-0.01, -1}, 
		p1_input_stat[2]{0.01, 1},
		p2_sliv_stat[2]{0.01, 1},
		p2_input_stat[2]{-0.01, -1}
		;

	// Регулятор

	double
		target_x = 0.5,
		K_p = 1,
		K_d = 0,
		K_i = 0;

	// Дополнительные обекты
	IForse *
		force = (new Forse_manipulator());

	RKSolver<4> solver = RKSolver<4>();

	double
		&t = solver.State[0],
		&v = solver.State[1],
		&x = solver.State[2],
		&p1 = solver.State[3],
		&p2 = solver.State[4],
		&dt = solver.h;

	inline double U1_input()
	{
		float X = (p1_input_stat[0] - U) / (p1_input_stat[0] - p1_input_stat[1]);

		return MAX2(0, MIN2(1, X));
	}

	inline double U1_sliv()
	{
		float X = (p1_sliv_stat[0] - U) / (p1_sliv_stat[0] - p1_sliv_stat[1]);

		return MAX2(0, MIN2(1, X));
	}

	inline double U2_input()
	{
		float X = (p2_input_stat[0] - U) / (p2_input_stat[0] - p2_input_stat[1]);

		return MAX2(0, MIN2(1, X));
	}

	inline double U2_sliv()
	{
		float X = (p2_sliv_stat[0] - U) / (p2_sliv_stat[0] - p2_sliv_stat[1]);

		return MAX2(0, MIN2(1, X));
	}

	void Reset()
	{
		{
			Forse_manipulator *T = dynamic_cast<Forse_manipulator *>(force);
			if(T != 0)
			{
				T->l = std::max(T->l, fabs(T->h - T->d) * 1.01);

				max_x = std::min(max_x, (T->h + T->d) * 0.99 - T->l);

			}
		}


		solver.State[0] = 0;
		solver.State[1] = 0;
		solver.State[2] = 0.5 * (max_x + min_x);
		solver.State[3] = p_input;
		solver.State[4] = p_sliv;

		{

			double pd = (S1 * solver.State[3] - S2 * solver.State[4] + force->F(solver.State, this)) / (S1 + S2);

			solver.State[3] = solver.State[3] - pd;
			solver.State[4] = solver.State[4] + pd;

		}
	}

	H_System()
	{
		solver.funcs[0] = DV;
		solver.funcs[1] = DX;
		solver.funcs[2] = DP1;
		solver.funcs[3] = DP2;
		Reset();

	}

	void Calc()
	{
		double OU = U;
		U = (target_x - solver.State[2]) * K_p;
		double DU = U - OU;

		DU = MAX2(-DU_max * solver.h, MIN2(DU_max * solver.h, DU));
		U = OU + DU;

		U = MAX2(-1, MIN2(1, U));

		solver.Calc(this);

		solver.State[2] = MAX2(min_x, MIN2(max_x, solver.State[2]));

	}

};

static double DV(double *var, void *_t)
{
	double
		&t = var[0],
		&v = var[1],
		&x = var[2],
		&p1 = var[3],
		&p2 = var[4];

	H_System *T = (H_System *)(_t);
	double _F = (T->S1 * p1 - T->S2 * p2 - T->b_prop * v + T->force->F(var, _t));
	double ret = _F / T->m;
	return ret;

	if(-_F / (T->solver.h * T->m * v) > 1)
	{
		if(_F * _F < T->f_tr_suh * T->f_tr_suh)
			ret = -v * T->solver.h;
	}
}

static  double  DX(double *var, void *_t)
{
	double
		&t = var[0],
		&v = var[1],
		&x = var[2],
		&p1 = var[3],
		&p2 = var[4];

	H_System *T = (H_System *)(_t);

	double ret = v;

	if(x >= T->max_x && v > 0)
		ret = 0;
	if(x <= T->min_x && v < 0)
		ret = 0;

	return ret;
}

static  double DP1(double *var, void *_t)
{
	double
		&t = var[0],
		&v = var[1],
		&x = var[2],
		&p1 = var[3],
		&p2 = var[4];

	H_System *T = (H_System *)(_t);

	double Q = 0;

	Q += T->mu * U_to_S(T->U1_input()) * T->S * sign(T->p_input - p1) * sqrt(2 / T->ro * fabs(T->p_input - p1));
	Q += T->mu * U_to_S(T->U1_sliv()) * T->S * sign(T->p_sliv - p1) * sqrt(2 / T->ro * fabs(T->p_sliv - p1));

	double Vt1 = T->V1 + T->S1 * (x - T->min_x);

	double ret = T->E / Vt1 * (Q - T->S1 * v);

	return ret;
}

static double DP2(double *var, void *_t)
{
	double
		&t = var[0],
		&v = var[1],
		&x = var[2],
		&p1 = var[3],
		&p2 = var[4];

	H_System *T = (H_System *)(_t);

	double Q = 0;

	Q += T->mu * U_to_S(T->U2_input()) * T->S * sign(T->p_input - p1) * sqrt(2 / T->ro * fabs(T->p_input - p1));
	Q += T->mu * U_to_S(T->U2_sliv()) * T->S * sign(T->p_sliv - p1) * sqrt(2 / T->ro * fabs(T->p_sliv - p1));

	double Vt2 = T->V2 + T->S2 * (T->max_x - x);

	double ret = T->E / Vt2 * (Q + T->S2 * v);

	return ret;
}

static  double  DU(double *var, void *_t)
{
	double
		&t = var[0],
		&v = var[1],
		&x = var[2],
		&p1 = var[3],
		&p2 = var[4];

	H_System *T = (H_System *)(_t);

	double ret = -v * T->K_p + (T->target_x - x) * T->K_p / T->solver.h;

	ret = MAX2(-T->DU_max, MIN2(T->DU_max, ret));

	return ret;
}
