#pragma once
#include <cmath>
#include <iomanip>
#include "RKSolver.h"
#include "Force.h"

#define sign(x) (x == 0 ? 0 : (x > 0 ? 1 : -1))


const double
G = 9.81,
PI = 3.1415926535;

static inline double U_to_S(double y)
{
	y = 2 * y - 1;
	y = MAX2(-1, MIN2(1, y));
	return
		1 - (acos(y) - y * sqrt(1 - y * y)) / PI;
}

static double U_to_S_dead_zone(double y, double dz)
{

	if(fabs(y) < dz)
	{
		y = 0;
	}
	else
	{
		y = y - (y > 0 ? dz : -dz);
		y /= 1 - dz;
	}

	return U_to_S(y);
}

inline float D2S(float D)
{
	return  std::powf(D * 0.5, 2) * 3.1415926535;
}

static double DV(double *var, void *_t);

static  double  DX(double *var, void *_t);

static  double DP1(double *var, void *_t);

static double DP2(double *var, void *_t);

static double DY(double *var, void *_t);

static double DDY(double *var, void *_t);

inline double Get_U(double stat[2], double U)
{
	float X = (stat[0] - U) / (stat[0] - stat[1]);

	return MAX2(0, MIN2(1, X));
}

struct H_System
{

	// прямой цилиндр
	double
		V1 = 0.01,
		S1 = D2S(80 * 1E-3);
	;  
	// обратный цилиндр
	double
		V2 = 0.01,
		S2 = D2S(80 * 1E-3) - D2S(60 * 1E-3);
	;
	// Свойства рабочей жидкости
	double
		ro = 860,
		E = 2.02E7
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
		m = 1,
		b = 10,
		f_tr_suh = 0,
		p_sliv = 135,
		p_input = 10 * 1E6
		;

	// Золотник
	double
		U = 1,
		U_K = 30,
		U_T = 0.01,
		p1_sliv_stat[2]{	0.0,	-1}, 
		p1_input_stat[2]{	0.02,	1},
		p2_sliv_stat[2]{	-0.0,	1},
		p2_input_stat[2]{	-0.02,	-1}
		;

	// Регулятор

	double
		target_x = 0.5,
		K_p = 10,
		K_d = 0,
		K_i = 0;

	// Дополнительные обекты
	IForse *
		force = (new Forse_manipulator());

	RKSolver<6> solver = RKSolver<6>();

	void Calc()
	{
		//U = 1;
		U = MAX2(-1, MIN2(1, U));

		solver.Calc(this);

		solver.State[2] = MAX2(min_x, MIN2(max_x, solver.State[2]));
		solver.State[5] = MAX2(-1, MIN2(1, solver.State[5]));
	}

	double Get_K()
	{
		double ret = 
			ro * (S1 * S1 * S1 + S2 * S2 * S2) / (2 * mu * mu * S * S);
		return
			ret;
	}

	double Get_F_0()
	{
		double ret =
			p_input * S1 - p_sliv * S2;
		return
			ret;
	}

	double Get_V_r()
	{
		double
			K = Get_K(),
			F_0 = Get_F_0(),
			ret =  b / (2 * K) * (sqrt(1 + 4 * K * abs(F_0) / (b * b)) - 1);

		return
			ret;
	}

	double Get_P1_r()
	{
		double
			v_r = Get_V_r(),
			ret = p_input - v_r * abs(v_r) * ro * S1 * S1 / (S * S * 2 * mu * mu);

		return
			ret;
	}

	double Get_P2_r()
	{
		double
			v_r = Get_V_r(),
			ret = p_sliv + v_r * abs(v_r) * ro * S2 * S2 / (S * S * 2 * mu * mu);

		return
			ret;
	}

	double Get_DP_t(double v, double a)
	{
		double
			F_0 = Get_F_0(),
			ret = (F_0 - v * b - a * m) / (S1 + S2);

		return
			ret;
	}

	double Get_P1_t(double v, double a)
	{
		return
			p_input - Get_DP_t(v, a);
	}

	double Get_P2_t(double v, double a)
	{
		return
			p_sliv + Get_DP_t(v, a);
	}

	double Get_T()
	{
		double
			K = Get_K(),
			v_r = Get_V_r(),
			ret =  m / (2 * K * abs(Get_V_r())+b);

		return
			ret;
	}

	inline double GetQ1(double *var)
	{
		double
			&t = var[0],
			&v = var[1],
			&x = var[2],
			&p1 = var[3],
			&p2 = var[4],
			&y = var[5],
			&dy = var[6];

		double Q = 0;
		Q += mu * Get_U(p1_input_stat, y) * S * sign(p_input - p1) * sqrt(2 / ro * fabs(p_input - p1));
		Q += mu * Get_U(p1_sliv_stat, y) * S * sign(p_sliv - p1) * sqrt(2 / ro * fabs(p_sliv - p1));

		return Q;
	}

	inline double GetQ2(double *var)
	{
		double
			&t = var[0],
			&v = var[1],
			&x = var[2],
			&p1 = var[3],
			&p2 = var[4],
			&y = var[5],
			&dy = var[6];

		double Q = 0;
		Q += mu * Get_U(p2_input_stat, y) * S * sign(p_input - p2) * sqrt(2 / ro * fabs(p_input - p2));
		Q += mu * Get_U(p2_sliv_stat, y) * S * sign(p_sliv - p2) * sqrt(2 / ro * fabs(p_sliv - p2));

		return Q;
	}

	inline double GetQI(double *var)
	{
		double
			&t = var[0],
			&v = var[1],
			&x = var[2],
			&p1 = var[3],
			&p2 = var[4],
			&y = var[5],
			&dy = var[6];

		double Q = 0;
		Q += mu * Get_U(p1_input_stat, y) * S * sign(p_input - p1) * sqrt(2 / ro * fabs(p_input - p1));
		Q += mu * Get_U(p2_input_stat, y) * S * sign(p_input - p2) * sqrt(2 / ro * fabs(p_input - p2));

		return Q;
	}

	inline double GetQS(double *var)
	{
		double
			&t = var[0],
			&v = var[1],
			&x = var[2],
			&p1 = var[3],
			&p2 = var[4],
			&y = var[5],
			&dy = var[6];

		double Q = 0;
		Q += mu * Get_U(p1_sliv_stat, y) * S * sign(p_sliv - p1) * sqrt(2 / ro * fabs(p_sliv - p1));
		Q += mu * Get_U(p2_sliv_stat, y) * S * sign(p_sliv - p2) * sqrt(2 / ro * fabs(p_sliv - p2));

		return Q;
	}

	double GetK_old()
	{
		double
			K = Get_K(),
			F_0 = Get_F_0(),
			v_r = Get_V_r(),
			ret = sign(F_0)*(abs(F_0) - b*b/(4 * K) * std::pow(sqrt(S*S +  (4 * K * abs(F_0)) / (b*b)), 2)) / b;	

		return
			ret;
	}



	void Reset(
		bool base_init = true,
		double v = 0, 
		double x = 0,
		double p1 = 0, 
		double p2 = 0,
		double u = 0)
	{
		{
			Forse_manipulator *T = dynamic_cast<Forse_manipulator *>(force);
			if(T != 0)
			{
				T->l = std::max(T->l, fabs(T->h - T->d) * 1.01);

				max_x = std::min(max_x, (T->h + T->d) * 0.99 - T->l);

			}
		}

		if(base_init)
		{
			solver.State[0] = 0;
			solver.State[1] = 0;
			solver.State[2] = (max_x + min_x) * 0.5;
			solver.State[3] = p_input;
			solver.State[4] = p_sliv;
			solver.State[5] = 0;

			{

				double pd = (S1 * solver.State[3] - S2 * solver.State[4] + force->F(solver.State, this)) / (S1 + S2);

				solver.State[3] = p_input - pd;
				solver.State[4] = p_sliv  + pd;

			}
		}
		else
		{
			solver.State[0] = 0;
			solver.State[1] = v;
			solver.State[2] = x;
			solver.State[3] = p1;
			solver.State[4] = p2;
			solver.State[5] = u;
		}
		//solver.State[6] = 0;

	}

	H_System()
	{
		solver.funcs[0] = DV;
		solver.funcs[1] = DX;
		solver.funcs[2] = DP1;
		solver.funcs[3] = DP2;
		solver.funcs[4] = DY;
		solver.funcs[5] = DDY;
		Reset();
	}
};

static double DV(double *var, void *_t)
{
	double
		&t = var[0],
		&v = var[1],
		&x = var[2],
		&p1 = var[3],
		&p2 = var[4],
		&y = var[5],
		&dy = var[6];

	H_System *T = (H_System *)(_t);
	double _F = (T->S1 * p1 - T->S2 * p2 - T->b * v);
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
		&p2 = var[4],
		&y = var[5],
		&dy = var[6];

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
		&p2 = var[4],
		&y = var[5],
		&dy = var[6];

	H_System *T = (H_System *)(_t);

	double Q = T->GetQ1(var);

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
		&p2 = var[4],
		&y = var[5],
		&dy = var[6];

	H_System *T = (H_System *)(_t);

	double Q = T->GetQ2(var);

	double Vt2 = T->V2 + T->S2 * (T->max_x - x);

	double ret = T->E / Vt2 * (Q + T->S2 * v);

	return ret;
}

static double DY(double *var, void *_t)
{
	double
		&t = var[0],
		&v = var[1],
		&x = var[2],
		&p1 = var[3],
		&p2 = var[4],
		&y = var[5],
		&dy = var[6];

	double ret = dy;

	if(y >= 1 && dy > 0)
		ret = 0;
	if(y <= -1 && dy < 0)
		ret = 0;

	return ret;
}

static double DDY(double *var, void *_t)
{
	double
		&t = var[0],
		&v = var[1],
		&x = var[2],
		&p1 = var[3],
		&p2 = var[4],
		&y = var[5],
		&dy = var[6];

	H_System *T = (H_System *)(_t);

	double u_t = T->U;

	u_t = MAX2(-1, MIN2(1, u_t));

	return u_t * T->U_K / T->U_T - dy / T->U_T -  y * T->U_K / T->U_T;
}