#include <iostream>
#include <iomanip>  
#include <SFML/Graphics.hpp>
#include<thread>
#include<time.h>
#include <ctime>
#include <ratio>
#include <chrono>
#include <fstream>
#include "RKSolver.h"
#include "PoliLine.h"

#include <vector>

const double
G = 9.81,
PI = 3.1415926535;

static double U_to_S(double u)
{
	u = 2 * u - 1;
	u = MAX2(-1, MIN2(1, u));
	return
		1 - (acos(u) - u * sqrt(1 - u * u)) / PI;
}

static double U_to_S_dead_zone(double u, double dz)
{

	if (fabs(u) < dz)
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

struct H_System
{

	// прямой цилиндр
	double
		V1 = 0.1,
		S1 = 0.1
		;
	// обратный цилиндр
	double
		V2 = 0.1,
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
		m = 1,
		b_prop = 0.0,
		f_tr_suh = 1,
		p_s = 0,
		p_i = 1000
		;
	// Управляющий сигнал
	double
		DU_max = 1,
		target_x = 0.5,
		U_deadZone = 0.05,
		K_p = 10

		;
	// Допущения

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
	l =  0.5,
	d = 1,
	dh = 0,
	dd = 0,
	Ft[2]
	{
		+0,
		-0
	},
		Fn
		;
	double
		cosa = 0,
		sina = 0;
	;

#define sign(x) (x == 0 ? 0 : (x > 0 ? 1 : -1))

	static double F(double* var, void* _t)
	{
		double
			& t = var[0],
			& v = var[1],
			& x = var[2],
			& p1 = var[3],
			& p2 = var[4];

		H_System *T = (H_System*)(_t);

		double
			l0 = T->l + x;
		T->cosa = (T->h * T->h + T->d * T->d - l0 * l0) / (2 * T->h * T->d);
		T->sina = sqrt(1 - T->cosa * T->cosa);

		double HN[2]
		{
			T->H[1],
			-T->H[0]
		};
		if (HN[1] < 0)
		{
			HN[0] *= -1;
			HN[1] *= -1;
		}

		double
			E[2]
		{
			(T->H[0] * T->cosa + HN[0] * T->sina) / T->h,
			(HN[1] * T->sina + T->H[1] * T->cosa) / T->h
		},
			L[2]
		{
			E[0] * T->d,
			E[1] * T->d
		},
			N[2]
		{
			(L[0] - T->A[0] + T->B[0]) / l0,
			(L[1] - T->A[1] - T->B[0]) / l0
		},
			D[2]
		{
			E[0] * (T->d + T->dd) + HN[0] * (T->dh),
			E[1] * (T->d + T->dd) + HN[1] * (T->dh)
		};


		double
			ret = (D[0] * T->Ft[1] - D[1] * T->Ft[0]) / (L[0] * N[1] - L[1] * N[0]);
		T->Fn = ret;
		return ret;
	}

	static  double DV(double* var, void* _t)
	{
		double
			& t = var[0],
			& v = var[1],
			& x = var[2],
			& p1 = var[3],
			& p2 = var[4];

		H_System* T = (H_System*)(_t);
		double _F = (T->S1 * p1 - T->S2 * p2 - T->b_prop * v + F(var, _t));
		double ret = _F / T->m;
		return ret;

		if (-_F / (T->solver.h * T->m * v) > 1)
		{			
			if (_F * _F < T->f_tr_suh * T->f_tr_suh)
				ret = -v * T->solver.h;
		}
		


	}

	static  double  DX(double* var, void* _t)
	{
		double
			& t = var[0],
			& v = var[1],
			& x = var[2],
			& p1 = var[3],
			& p2 = var[4];

		H_System* T = (H_System*)(_t);

		double ret = v;

		if (x >= T->max_x && v > 0)
			ret = 0;
		if (x <= T->min_x && v < 0)
			ret = 0;

		return ret;
	}

	static  double DP1(double* var, void* _t)
	{
		double
			& t = var[0],
			& v = var[1],
			& x = var[2],
			& p1 = var[3],
			& p2 = var[4],
			& u = var[5];

		H_System* T = (H_System*)(_t);

		double 
			_S = T->S * U_to_S_dead_zone(fabs(u), T->U_deadZone),
			Q1 = 0;

		if (u > 0)
			Q1 = T->mu * _S * sign(T->p_i - p1) * sqrt(2 / T->ro * fabs(T->p_i - p1));
		else
			Q1 = -T->mu * _S * sign(p1 - T->p_s) * sqrt(2 / T->ro * fabs(T->p_s - p1));

		double Vt1 = T->V1 + T->S1 * (x - T->min_x);

		double ret = T->E / Vt1 * (Q1 - T->S1 * v);

		return ret;
	}

	static double DP2(double* var, void* _t)
	{
		double
			& t = var[0],
			& v = var[1],
			& x = var[2],
			& p1 = var[3],
			& p2 = var[4],
			& u = var[5];

		H_System* T = (H_System*)(_t);

		double 
			_S = T->S * U_to_S_dead_zone(fabs(u), T->U_deadZone),
			Q2 = 0;

		if (u > 0)
			Q2 = -T->mu * _S * sign(p2 - T->p_s) * sqrt(2 / T->ro * fabs(T->p_s - p2));
		else
			Q2 = T->mu * _S * sign(T->p_i - p2) * sqrt(2 / T->ro * fabs(T->p_i - p2));

		double Vt2 = T->V2 + T->S2 * (T->max_x - x);

		double ret = T->E / Vt2 * (Q2 + T->S2 * v);

		return ret;
	}

	static  double  DU(double* var, void* _t)
	{
		double
			& t = var[0],
			& v = var[1],
			& x = var[2],
			& p1 = var[3],
			& p2 = var[4];

		H_System* T = (H_System*)(_t);

		double ret = -v * T->K_p + (T->target_x - x) * 1e20;

		ret = MAX2(-T->DU_max, MIN2(T->DU_max, ret));

		return ret;
	}


	RKSolver<5> solver;

	H_System()
	{
		solver.funcs[0] = DV;
		solver.funcs[1] = DX;
		solver.funcs[2] = DP1;
		solver.funcs[3] = DP2;
		solver.funcs[4] = DU;
		solver.State[2] = 1;
		solver.State[5] = 0.0000000000001 * sign(target_x - solver.State[2]);
		solver.State[3] = (solver.State[5] > 0 ? p_i : p_s);
		solver.State[4] = (solver.State[5] < 0 ? p_i : p_s);

		max_x = std::min(max_x, (h + d) * 0.99 - l);

		F(solver.State, this);

		{

			double pd = (S1 * solver.State[3] - S2 * solver.State[4] + Fn) / (S1 + S2);

			solver.State[3] = solver.State[3] - pd;
			solver.State[4] = solver.State[4] + pd;

		}

	}

	void Calc()
	{
		solver.Calc(this);


		solver.State[5] = MAX2(-1, MIN2(1, solver.State[5]));

		solver.State[2] = MAX2(min_x, MIN2(max_x, solver.State[2]));
	}

};

void main()
{

	H_System HD;

	HD.solver.h = 1e-20;
	HD.F(HD.solver.State, &HD);
	double ba = acos(HD.cosa);
	std::ofstream Out("out.txt");
	Out << std::scientific;
	Out << std::setprecision(15);
	Out
		<< "T" << ","
		<< "V" << ","
		<< "X" << ","
		<< "P1" << ","
		<< "P2" << ","
		<< "F" << ","
		<< "H" << ","
		<< "N" << ","
		<< "A" << ","
		<< "U" << ","
		<< "CU" << ","
		<< "DT" << "\n";
	std::vector<double>
		t,
		x,
		v,
		p1,
		p2,
		n,
		us;

	float
		target = 1,
		dU = 10;
	const float deadZone = 0.0;
	float U = deadZone * 2, OU = deadZone * 2;

	float 
		dt = 0.0001,
		ot = -dt * 2;

	float maxTime = 30;

	int j = 0;

	float I = 0, BI = 0;

	double mh = HD.solver.h;
	//HD.solver.State[2] = target;
	//HD.solver.State[1] = 1;
	std::cout << std::scientific;
	std::cout << std::setprecision(2);
	double &T = HD.solver.State[0];
	HD.solver.State[2] = 1;
	for (int i = 0, j = 0, J = 0, KK = 0; T < maxTime; i++)
	{
		if (KK == 0 && T > 7)
		{
			KK++;
			HD.target_x = 0.1;
		}
		else if (KK == 1 && T > 15)
		{
			KK++;
			HD.target_x = 0.3;
		}
		else if (KK == 2 && T > 20)
		{
			KK++;
			HD.target_x = 0.2;
		}
		else if (KK == 3 && T > 25)
		{
			KK++;
			HD.target_x = 0.25;
		}

		double oh = HD.solver.h;
		if (i % 100 == 0)
		{
			HD.solver.UpateH(&HD, MIN2(1e-6, HD.solver.h * 2), 1e-0);
			//mh = mh * 0.99 + HD.solver.h * 0.01;
			//HD.solver.h = MIN2(mh, HD.solver.h);
		}



		HD.F(HD.solver.State, &HD);
		if (HD.solver.State[0] - ot >= dt)
		{
			
			while (HD.solver.State[0] - ot >= dt)
				ot += dt;
			t.push_back(HD.solver.State[0]);
			x.push_back(HD.solver.State[2]);
			v.push_back(HD.solver.State[1]);
			p1.push_back(HD.solver.State[3]);
			p2.push_back(HD.solver.State[4]);
			n.push_back(HD.Fn);
			us.push_back(HD.solver.State[5]);
			Out
				<< HD.solver.State[0] << ",\t"
				<< HD.solver.State[1] << ",\t"
				<< HD.solver.State[2] << ",\t"
				<< HD.solver.State[3] * HD.S1 << ",\t"
				<< HD.solver.State[4] * HD.S2 << ",\t"
				<< HD.solver.State[3] * HD.S1 - HD.solver.State[4] * HD.S2 + HD.Fn << ",\t"
				<< (HD.solver.State[3] * HD.S1 - HD.solver.State[4] * HD.S2 + HD.Fn - HD.b_prop * HD.solver.State[1]) / HD.m << ",\t"
				<< HD.Fn << ","
				<< (acos(HD.cosa) - ba) * 180 / PI << ","
				<< U_to_S_dead_zone(HD.solver.State[5], HD.U_deadZone) << ","
				<< HD.solver.State[5] << ","
				<< HD.solver.h << "\n";;
			if(HD.solver.State[0] > 1)
				if(mh > 1E-4 && fabsl(HD.solver.State[1]) < 1E-10)
				{
					//break;
				}
			Out.flush();
			if (j % 1000 == 0)
			{
				for (int i = 0; i < 20; i++)
					std::cout << (i / 19.0 >= HD.solver.State[0] / maxTime ? "_" : "X");
				std::cout << HD.solver.h << ' ' << HD.solver.State[0] << '\r';
			}
			j++;
		}
		if (J > 1000000)
		{
			for (int i = 0; i < 20; i++)
				std::cout << (i / 19.0 >= HD.solver.State[0] / maxTime ? "_" : "X");
			std::cout << HD.solver.h << ' ' << HD.solver.State[0] << '\r';
			J = 0;
		}
		J++;
		HD.Calc();
		
		double& _x = HD.solver.State[2];
		if (_x > HD.max_x)
			_x = HD.max_x;
		if (_x < HD.min_x)
			_x = HD.min_x;
	}
	std::cout << '\n';
	Out.flush();
	Out.close();

	system("python main.py");

}

