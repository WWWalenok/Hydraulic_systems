﻿#include <iostream>
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
		S = 0.02
		;
	// Харрактеристики системы
	double
		min_x = 0,
		max_x = 2,
		m = 2,
		b_prop = 5,
		p_s = 0,
		p_i = 2000
		;
	// Управляющий сигнал
	double
		U = 1
		;
	// Геометрия
	double
		A[2]{ 1,-2 },
		B[2]{ 0,0 },
		H[2]
	{
		A[0] - B[0],
		A[1] - B[1],
	},
	h = sqrt(H[0] * H[0] + H[1] * H[1]),
	l =  1.5,
	d = 1,
	dh = 0,
	dd = 2,
	Ft[2]
	{
		+0,
		-m * G
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
		solver.funcs[0] = &H_System::DV;
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

	HD.solver.h = 0.002;
	HD.F(HD.solver.State, &HD);
	double ba = acos(HD.cosa);
	std::ofstream Out("out.txt");
	Out << std::scientific;

	std::vector<double>
		t,
		x,
		v,
		p1,
		p2,
		n;

	for (int i = 0; i < 10000; i++)
	{
		HD.F(HD.solver.State, &HD);
		if (i % 10 == 0)
		{
			t.push_back(HD.solver.State[0]);
			x.push_back(HD.solver.State[2]);
			v.push_back(HD.solver.State[1]);
			p1.push_back(HD.solver.State[3]);
			p2.push_back(HD.solver.State[4]);
			n.push_back(HD.Fn);
			Out
				<< HD.solver.State[0] << "\t"
				<< HD.solver.State[1] << "\t"
				<< HD.solver.State[2] << "\t"
				<< HD.solver.State[3] * HD.S1 << "\t"
				<< HD.solver.State[4] * HD.S2 << "\t"
				<< HD.Fn << "\t"
				<< (acos(HD.cosa) - ba) * 180 / PI << '\n';
		}
		HD.solver.Calc(&HD);
		
		double& _x = HD.solver.State[2];
		if (_x > HD.max_x)
			_x = HD.max_x;
		if (_x < HD.min_x)
			_x = HD.min_x;
	}

	Out.flush();
	Out.close();

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
		S_P1,
		S_P2,
		S_N
		;

	X.Set(new sf::Vertex[2], 2);
	Y.Set(new sf::Vertex[2], 2);
	S_P1.Set(new sf::Vertex[t.size()], t.size());
	S_P2.Set(new sf::Vertex[t.size()], t.size());
	S_N.Set(new sf::Vertex[t.size()], t.size());

	X.SetColor(sf::Color::White);
	Y.SetColor(sf::Color::White);
	S_P1.SetColor(sf::Color::Red);
	S_P2.SetColor(sf::Color::Blue);
	S_N.SetColor(sf::Color::Green);

	X[0].position = { float(10), float(W - 10) };
	X[1].position = { float(H - 10), float(W - 10) };

	Y[0].position = { float(10), float(W - 10) };
	Y[1].position = { float(10), float(10) };

	double Y_min = 1e20;
	double Y_max = -1e20;

	for (int i = 0; i < t.size(); i++)
	{
		S_N[i].position = { float(10 + i * (H - 20) / float(t.size() - 1.0)), float(W * 0.5) };
		S_P1[i].position.x = S_P2[i].position.x = S_N[i].position.x;
		Y_min = std::min(p1[i], std::min(p2[i], std::min(n[i], Y_min)));
		Y_max = std::max(p1[i], std::max(p2[i], std::max(n[i], Y_max)));
	}

	for (int i = 0; i < t.size(); i++)
	{
		S_N[i].position.y = W - 10 - (W - 20) * ((n[i] - Y_min) / (Y_max - Y_min));
		S_P1[i].position.y = W - 10 - (W - 20) * ((p1[i] - Y_min) / (Y_max - Y_min));
		S_P2[i].position.y = W - 10 - (W - 20) * ((p2[i] - Y_min) / (Y_max - Y_min));
	}

	while (window.isOpen())
	{
		sf::Event event;
		while (window.pollEvent(event))
		{
			if (event.type == sf::Event::Closed)
				window.close();
		}

		window.clear();
		//window.draw(shape);
		window.draw(X);
		window.draw(Y);
		window.draw(S_N);
		window.draw(S_P1);
		window.draw(S_P2);
		window.display();
	}


}
