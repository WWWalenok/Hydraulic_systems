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

inline float D2S(float D)
{
	return  std::powf(D * 0.5, 2) * 3.1415926535;
}

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
		E = 1.56E9
		;
	// Харрактеристики дроссилей
	double
		mu = 0.62,
		S = D2S((14 - 6) * 1E-3)
		;
	// Харрактеристики системы
	double
		min_x = 0,
		max_x = 90,
		m = 10,
		b_prop = 5,
		f_tr_suh = 0,
		p_s = 0.1 * 1E6,
		p_i = 10 * 1E6
		;
	// Управляющий сигнал
	double
		U = 1,
		DU_max = 10,
		target_x = 1000,
		U_deadZone = 0.00,
		K_p = 1

		;
	// Допущения

	// Геометрия
	double
		A[2]{ 200, 50},
		B[2]{ 0,0 },
		H[2]
	{
		A[0] - B[0],
		A[1] - B[1],
	},
	h = sqrt(H[0] * H[0] + H[1] * H[1]),
	l = max_x + 350 * 1E-2,
	d = 50,
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

	float V_r()
	{
		double* var = solver.State;
		double
			& t = var[0],
			& v = var[1],
			& x = var[2],
			& p1 = var[3],
			& p2 = var[4];
		float
			F0 = p_i * S1 - p_s * S2 + F(var, this),
			K = ro * (S1 * S1 * S1 + S2 * S2 * S2) / (2 * mu * mu),
			Sgn = (F0 > 0 ? 1 : -1),
			_S = S * U_to_S_dead_zone(fabs(U), U_deadZone);
		float ret = Sgn * _S * _S / (2 * K) * (sqrt(b_prop * b_prop + 4 * K / (_S * _S) * fabs(F0)) - b_prop);		
		return ret;
	}

#define sign(x) (x == 0 ? 0 : (x > 0 ? 1 : -1))

	static double F(double* var, void* _t)
	{
		double
			& t = var[0],
			& v = var[1],
			& x = var[2],
			& p1 = var[3],
			& p2 = var[4];

		//Dummy
		return 0;

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

	static double DV(double* var, void* _t)
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
			& p2 = var[4];

		H_System* T = (H_System*)(_t);

		double 
			_S = T->S * U_to_S_dead_zone(fabs(T->U), T->U_deadZone),
			Q1 = 0;

		if (T->U > 0)
			Q1 = T->mu * _S * sign(T->p_i - p1) * sqrt(2 / T->ro * fabs(T->p_i - p1));
		else
			Q1 = T->mu * _S * sign( T->p_s - p1) * sqrt(2 / T->ro * fabs(T->p_s - p1));

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
			& p2 = var[4];

		H_System* T = (H_System*)(_t);

		double 
			_S = T->S * U_to_S_dead_zone(fabs(T->U), T->U_deadZone),
			Q2 = 0;

		if (T->U > 0)
			Q2 = T->mu * _S * sign(T->p_s - p2) * sqrt(2 / T->ro * fabs(T->p_s - p2));
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

		double ret = -v * T->K_p + (T->target_x - x) * T->K_p / T->solver.h;

		ret = MAX2(-T->DU_max, MIN2(T->DU_max, ret));

		return ret;
	}


	RKSolver<4> solver;

	void Reset()
	{
		l = std::max(l, fabs(h - d) * 1.01);

		max_x = std::min(max_x, (h + d) * 0.99 - l);


		solver.State[0] = 0;
		solver.State[1] = 0;
		solver.State[2] = 0.5 * (max_x + min_x);
		solver.State[3] = p_i;
		solver.State[4] = p_s;

		F(solver.State, this);


		{

			double pd = (S1 * solver.State[3] - S2 * solver.State[4] + Fn) / (S1 + S2);

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

void Var1()
{
	H_System HD;
	HD.U = 0;
	HD.K_p = 0;
	HD.solver.h = 1e-7;

	float
		dt = 0.001,
		ot = -dt * 2;

	float maxTime = 25;

	double& T = HD.solver.State[0];

	HD.Reset();

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
		<< "CU" << ","
		<< "U" << ","
		<< "DT" << "\n";
	std::vector<double>
		t,
		x,
		v,
		p1,
		p2,
		n,
		us;
	int j = 0;

	float I = 0, BI = 0;

	double mh = HD.solver.h;
	HD.solver.State[2] = HD.max_x * 0.5;
	HD.solver.State[1] = 1;
	std::cout << std::scientific;
	std::cout << std::setprecision(2);
	//HD.solver.State[2] = 0.5;

	int medC = maxTime / HD.solver.h / 100;
	HD.target_x = HD.max_x * 0.9;
	std::cout << HD.V_r() << std::endl; 
	for (int i = 0, j = 0, J = 0, KK = 0; T < maxTime; i++)
	{
		if (KK == 0 && T > 7)
		{
			KK++;
			//HD.target_x = 0.1;
		}
		else if (KK == 1 && T > 15)
		{
			KK++;
			//HD.target_x = 0.3;
		}
		else if (KK == 2 && T > 20)
		{
			KK++;
			//HD.target_x = 0.2;
		}
		else if (KK == 3 && T > 25)
		{
			KK++;
			//HD.target_x = 0.25;
		}

		double oh = HD.solver.h;
		if (i % 100 == -1)
		{
			HD.solver.UpateH(&HD, MIN2(1e-6, HD.solver.h * 2), 1e-0);
			//mh = mh * 0.99 + HD.solver.h * 0.01;
			//HD.solver.h = MIN2(mh, HD.solver.h);
		}



		HD.F(HD.solver.State, &HD);
		if (T - ot >= dt)
		{
			while (T - ot >= dt)
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
				<< U_to_S_dead_zone(fabs(HD.U), HD.U_deadZone) << ","
				<< HD.U << ","
				<< HD.solver.h << "\n";;
			if (HD.solver.State[0] > 1)
				if (mh > 1E-4 && fabsl(HD.solver.State[1]) < 1E-10)
				{
					//break;
				}
			Out.flush();
		}
		if (J > medC)
		{
			for (int i = 0; i < 20; i++)
				std::cout << (i / 19.0 >= HD.solver.State[0] / maxTime ? "_" : "X");
			std::cout << HD.solver.h << ' ' << HD.solver.State[0] << '\r';
			J = 0;
		}
		J++;
		HD.Calc();
	}
	std::cout << '\n';
	Out.flush();
	Out.close();

	system("python main.py");
}

void Var2()
{
	H_System HD;

	HD.solver.h = 1e-6;

	float
		dt = 0.001,
		ot = -dt * 2;

	float maxTime = 0.5;

	double& T = HD.solver.State[0];

	for (int a = -5; a <= 5; a++) for (int b = -5; b <= 5; b++) for (int c = -5; c <= 5; c++)
	{
		HD.Reset();

		std::ofstream fout = std::ofstream(
			"out\\smple_" +
			std::to_string(a) +
			"_" +
			std::to_string(b) +
			"_" +
			std::to_string(c) +
			".csv"
		);
		std::cout
			<< "smple_" +
			std::to_string(a) +
			"_" +
			std::to_string(b) +
			"_" +
			std::to_string(c) +
			'\n';
		fout << std::scientific;
		fout << std::setprecision(10);
		fout << "t,\tv,\tx,\tp1,\tp2,\ta\n";
		HD.solver.State[1] = a * 3.34e-02 * 0.5;

		float minp = MIN2(HD.solver.State[3], HD.solver.State[4]);

		HD.solver.State[3] += b * minp * 0.1;
		HD.solver.State[4] += c * minp * 0.1;
		for (int i = 0, j = 0, J = 0, KK = 0; T <= maxTime; i++)
		{
			J++;
			if (i % 100 == 0)
			{
				fout
					<< HD.solver.State[0] << ",\t"
					<< HD.solver.State[1] << ",\t"
					<< HD.solver.State[2] << ",\t"
					<< HD.solver.State[3] << ",\t"
					<< HD.solver.State[4] << ",\t"
					<< HD.DV(HD.solver.State, &HD) << "\n"
					;
			}
			HD.Calc();
		}
	}

}

void main()
{
	Var1();

}

