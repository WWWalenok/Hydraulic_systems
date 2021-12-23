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
		m = 10,
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
	l =  0.5,
	d = 1,
	dh = 0,
	dd = 0,
	Ft[2]
	{
		+0,
		-5
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

		double ret = (T->S1 * p1 - T->S2 * p2 - T->b_prop * v + F(var, _t)) / T->m;

		return ret;
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

		double _S = T->S * fabs(T->U),
			Q1 = 0;

		if (T->U > 0)
		{
			Q1 = T->mu * _S * sign(T->p_i - p1) * sqrt(2 / T->ro * fabs(T->p_i - p1));
		}
		else
		{
			Q1 = -T->mu * _S * sign(p1 - T->p_s) * sqrt(2 / T->ro * fabs(T->p_s - p1));
		}

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
		double _S = T->S * fabs(T->U),
			Q2 = 0;

		if (T->U > 0)
			Q2 = -T->mu * _S * sign(p2 - T->p_s) * sqrt(2 / T->ro * fabs(T->p_s - p2));
		else
			Q2 = T->mu * _S * sign(T->p_i - p2) * sqrt(2 / T->ro * fabs(T->p_i - p2));

		double Vt2 = T->V2 + T->S2 * (T->max_x - x);

		double ret = T->E / Vt2 * (Q2 + T->S2 * v);

		return ret;
	}

	RKSolver<4> solver;

	H_System()
	{
		solver.funcs[0] = DV;
		solver.funcs[1] = DX;
		solver.funcs[2] = DP1;
		solver.funcs[3] = DP2;
		solver.State[3] = (U > 0 ? p_i : p_s);
		solver.State[4] = (U < 0 ? p_i : p_s);

		max_x = std::min(max_x, (h + d) * 0.99 - l);

		F(solver.State, this);

		{

			double pd = (S1 * solver.State[3] - S2 * solver.State[4] + Fn) / (S1 + S2);

			solver.State[3] = solver.State[3] - pd;
			solver.State[4] = solver.State[4] + pd;

		}

	}

};

void main()
{

	H_System HD;

	HD.solver.h = 0.0000001;
	HD.F(HD.solver.State, &HD);
	double ba = acos(HD.cosa);
	std::ofstream Out("out.txt");
	Out << std::scientific;
	Out << std::setprecision(10);
	Out
		<< "T" << ","
		<< "V" << ","
		<< "X" << ","
		<< "P1" << ","
		<< "P2" << ","
		<< "F" << ","
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
		k_p = 150,
		k_d = 0,
		k_i = 0,
		target = 0.5,
		dU = 10;
	float U = 0, OU = 0;

	float 
		dt = 0.001,
		ot = -dt * 2;

	float maxTime = 100;

	int j = 0;

	float I = 0, BI = 0;

	for (int i = 0; HD.solver.State[0] < 100; i++)
	{

		const float deadZone = 0.01;
		if (i % 1 == 0)
			HD.solver.UpateH(&HD, HD.solver.h * 2, 1e-12);
		U = HD.solver.State[1] * k_d + (target - HD.solver.State[2]) * k_p + I * k_i;
		if (fabs(target - HD.solver.State[2]) < 1e-7)
			U = 0;

		U = fmin(fmax(-1, U), 1);

		if (fabs(OU - U) > dU * HD.solver.h)
		{
			U = OU + ((U - OU) > 0 ? dU * HD.solver.h : -dU * HD.solver.h);
		}
		OU = U;
		float TU = U;

		I = I + (fabs(target - HD.solver.State[2]) < 1e-3 ? (target - HD.solver.State[2]) * HD.solver.h : 0);
		if (I * (target - HD.solver.State[2]) < 0)
			I *= .99;
		if (fabs(TU) < deadZone)
		{
			TU = 0;
		}
		else
		{
			TU = TU - (TU > 0 ? deadZone : -deadZone);
			TU /= 1 - deadZone;
		}

		HD.U = TU;

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
			us.push_back(TU);
			Out
				<< HD.solver.State[0] << ","
				<< HD.solver.State[1] << ","
				<< HD.solver.State[2] << ","
				<< HD.solver.State[3] * HD.S1 << ","
				<< HD.solver.State[4] * HD.S2 << ","
				<< HD.Fn << ","
				<< (acos(HD.cosa) - ba) * 180 / PI << ","
				<< TU << ","
				<< I << ","
				<< HD.solver.h << "\n";;
			if(HD.solver.State[0] > 1)
				if(fabs(v[x.size() - 1]) < 1E-7 && fabs(v[x.size() - 2]) < 1E-7)
				{
					//break;
				}

			if (j % 100 == 0)
			{
				for (int i = 0; i < 20; i++)
					std::cout << (i / 20.0 > HD.solver.State[0] / maxTime ? "_" : "X");
				std::cout << '\r';
			}
			j++;
		}
		HD.solver.Calc(&HD);
		
		double& _x = HD.solver.State[2];
		if (_x > HD.max_x)
			_x = HD.max_x;
		if (_x < HD.min_x)
			_x = HD.min_x;
	}
	std::cout << '\n';
	Out.flush();
	Out.close();
	return;
	uint32_t H = 800, W = 400;

	sf::ContextSettings context_setting(0, 0, 2);
	sf::RenderWindow window(sf::VideoMode(H, W), "SFML window", sf::Style::Default, context_setting);
	sf::CircleShape shape(100.f);
	shape.setFillColor(sf::Color::Green);

	PoliLine 
		X,
		Y,
		S_X,
		S_V,
		S_V0,
		S_P1,
		S_P2,
		S_N,
		S_T
		;

	X.Set(new sf::Vertex[2], 2);
	Y.Set(new sf::Vertex[2], 2);
	S_V0.Set(new sf::Vertex[2], 2);
	S_P1.Set(new sf::Vertex[t.size()], t.size());
	S_P2.Set(new sf::Vertex[t.size()], t.size());
	S_N.Set(new sf::Vertex[t.size()], t.size());
	S_X.Set(new sf::Vertex[t.size()], t.size());
	S_T.Set(new sf::Vertex[t.size()], t.size());
	S_V.Set(new sf::Vertex[t.size()], t.size());

	X.SetColor(sf::Color::White);
	Y.SetColor(sf::Color::White);
	S_P1.SetColor(sf::Color::Red);
	S_P2.SetColor(sf::Color::Blue);
	S_N.SetColor(sf::Color::Green);
	S_X.SetColor(sf::Color::Red);
	S_T.SetColor(sf::Color::Green);
	S_V.SetColor(sf::Color::Green);
	S_V0.SetColor(sf::Color::Red);

	X[0].position = { float(10), float(W - 10) };
	X[1].position = { float(H - 10), float(W - 10) };

	S_V0[0].position = { float(10), float(W - 10) };
	S_V0[1].position = { float(H - 10), float(W - 10) };

	Y[0].position = { float(10), float(W - 10) };
	Y[1].position = { float(10), float(10) };

	double Y_min = 1e20;
	double Y_max = -1e20;

	double
		Vmin = 1e20,
		Vmax = -1e20;

	int Mode = 1;

	for (int i = 0; i < t.size(); i++)
	{
		S_N[i].position = { float(10 + i * (H - 20) / float(t.size() - 1.0)), float(W * 0.5) };
		S_V[i].position.x = S_T[i].position.x = S_X[i].position.x = S_P1[i].position.x = S_P2[i].position.x = S_N[i].position.x;
		Y_min = std::min(p1[i], std::min(p2[i], std::min(n[i], Y_min)));
		Y_max = std::max(p1[i], std::max(p2[i], std::max(n[i], Y_max)));
		Vmin = std::min(v[i], Vmin);
		Vmax = std::max(v[i], Vmax);
	}

	for (int i = 0; i < t.size(); i++)
	{
		S_N[i].position.y = W - 10 - (W - 20) * ((n[i] - Y_min) / (Y_max - Y_min));
		S_P1[i].position.y = W - 10 - (W - 20) * ((p1[i] - Y_min) / (Y_max - Y_min));
		S_P2[i].position.y = W - 10 - (W - 20) * ((p2[i] - Y_min) / (Y_max - Y_min));
	}

	for (int i = 0; i < t.size(); i++)
	{
		S_X[i].position.y = W - 10 - (W - 20) * ((x[i] - HD.min_x) / (HD.max_x - HD.min_x));
		S_T[i].position.y = W - 10 - (W - 20) * ((target - HD.min_x) / (HD.max_x - HD.min_x));
	}

	for (int i = 0; i < t.size(); i++)
	{
		S_V[i].position.y = W - 10 - (W - 20) * ((v[i] - Vmin) / (Vmax - Vmin));
	}

	S_V0[1].position.y = S_V0[0].position.y = W - 10 - (W - 20) * ((0 - Vmin) / (Vmax - Vmin));

	while (window.isOpen())
	{
		sf::Event event;
		while (window.pollEvent(event))
		{
			if (event.type == sf::Event::Closed)
				window.close();
		}

		window.clear();

		window.draw(X);
		window.draw(Y);

		switch (Mode)
		{
		case 0:
			window.draw(S_V);
			window.draw(S_V0);

			break;
		case 1:
			window.draw(S_T);
			window.draw(S_X);

			break;
		default:
			break;
		}

		window.display();
	}


}

