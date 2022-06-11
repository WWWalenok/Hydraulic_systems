#include <iostream>
#include <fstream>
#include "Hidro.h"
#include <vector>
#include <string>
#include <strstream>

void PrintHD(std::ofstream *fout, H_System &HD)
{

}

float maxTime = 10 * 1.000;
float dt = 0.001;

double ie = 0;

float
	K_p = 1,
	K_d = 0,
	K_i = 0;

float U = 4;

float GetU(float T, float x, float v, double dt)
{
	double e = U - x;
	ie += e * dt;
	return e * K_p + v * K_d + (ie - e * dt * 0.5) * K_i;

	if(T < 0.1)
		U = 0;
	else if(T < 1.1)
		U = 1;
	else if(T < 2.1)
		U = 0;
	else if(T < 3.1)
		U = -1;
	else if(T < 4)
		U = 0;
	else if(T < 4.5)
		U = 0.5;
	else if(T < 5)
		U = -0.5;
	else if(T < 10)
		U = 0;
	else if(T < 10)
		U = 0;

	return U;
}

void Base(H_System& HD, std::string name)
{
	Forse_manipulator *force = dynamic_cast<Forse_manipulator *>(HD.force);
	HD.solver.h = 1e-7;
	float
		ot = -dt * 2;
	double &T = HD.solver.State[0];
	force->F(HD.solver.State, &HD);
	double ba = acos(force->cosa);
	std::ofstream Out(name);
	Out << std::scientific;
	Out << std::setprecision(15);
	Out
		<< "T" << ","
		<< "V" << ","
		<< "X" << ","
		<< "P1" << ","
		<< "P2" << ","
		<< "A" << ","
		<< "F" << ","
		<< "F1" << ","
		<< "F2" << ","
		<< "Q1" << ","
		<< "Q2" << ","
		<< "QI" << ","
		<< "QS" << ","
		<< "Y" << ","
		<< "U" << ","
		<< "DT" << "\n";
	int j = 0;

	float I = 0, BI = 0;

	double mh = HD.solver.h;

	int medC = maxTime / HD.solver.h / 250;

	for(int i = 0, j = 0, J = medC * 2, KK = 0; T <= maxTime; i++)
	{
		HD.U = GetU(T, HD.solver.State[2], HD.solver.State[1], HD.solver.h);
		HD.U = MAX2(-1, MIN2(1, HD.U));
		force->F(HD.solver.State, &HD);
		if(T - ot >= dt)
		{
			while(T - ot >= dt)
				ot += dt;
			Out
				<< HD.solver.State[0] << ",\t"
				<< HD.solver.State[1] << ",\t"
				<< HD.solver.State[2] << ",\t"
				<< HD.solver.State[3] << ",\t"
				<< HD.solver.State[4] << ",\t"
				<< (HD.solver.State[3] * HD.S1 - HD.solver.State[4] * HD.S2 + force->Fn - HD.solver.State[1] * HD.b) / HD.m << ",\t"
				<< HD.solver.State[3] * HD.S1 - HD.solver.State[4] * HD.S2 + force->Fn << ",\t"
				<< HD.solver.State[3] * HD.S1 << ",\t"
				<< HD.solver.State[4] * HD.S2 << ",\t"
				<< HD.GetQ1(HD.solver.State) << ",\t"
				<< HD.GetQ2(HD.solver.State) << ",\t"
				<< HD.GetQI(HD.solver.State) << ",\t"
				<< HD.GetQS(HD.solver.State) << ",\t"
				<< HD.solver.State[5] << ",\t"
				<< HD.U << ","
				<< HD.solver.h << "\n";
			Out.flush();
		}
		if(J > medC)
		{
			for(int i = 0; i < 20; i++)
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
}

void Linear(H_System& HD, std::string name, float P_K, float P_T)
{
	float
		ot = -dt * 2;

	float dt_solve = 1e-7;

	std::ofstream Out(name);
	Out << std::scientific;
	Out << std::setprecision(15);
	Out
		<< "T" << ","
		<< "X" << ","
		<< "V" << ","
		<< "A" << ","
		<< "P1" << ","
		<< "P2" << ","
		<< "U" << "\n";
	int j = 0;

	float I = 0, BI = 0;

	std::cout << std::scientific;
	std::cout << std::setprecision(2);
	//HD.solver.State[2] = 0.5;

	double 
		T = 0,
		x = 0,
		v = 0,
		a = 0,
		u = 0,
		du = 0,
		ddu = 0,
		p1 = 0,
		p2 = 0;

	int medC = maxTime / dt_solve / 250;

	HD.target_x = 1;

	p1 = HD.Get_P1_t(v, a);
	p2 = HD.Get_P2_t(v, a);

	for(int i = 0, j = 0, J = medC * 2, KK = 0; T <= maxTime; i++)
	{
		HD.U = GetU(T, x, v, dt_solve);
		HD.U = MAX2(-1, MIN2(1, HD.U));
		if(T - ot >= dt)
		{
			while(T - ot >= dt)
				ot += dt;
			Out
				<< T << ",\t"
				<< x << ",\t"
				<< v << ",\t"
				<< a << ",\t"
				<< p1 << ",\t"
				<< p2 << ",\t"
				<< HD.U << "\n";
			Out.flush();
		}
		if(J > medC)
		{
			for(int i = 0; i < 20; i++)
				std::cout << (i / 19.0 > T / maxTime ? "_" : "X");
			std::cout << dt_solve << ' ' << T << '\r';
			J = 0;
		}
		//solve
		{
			double nddu = HD.U * HD.U_K / HD.U_T - du /HD.U_T -  u * HD.U_K / HD.U_T;
			double ndu = du + dt_solve * (nddu);
			double nu = u + dt_solve * (ndu);
			
			double na = (u * P_K - v) / P_T;
			double nv = v + dt_solve * (na  + a) * 0.5;
			double nx = x + dt_solve * (nv + v) * 0.5;

			u = nu;
			du = ndu;

			x = nx;
			v = nv;
			a = na;

			T += dt_solve;

			p1 = HD.Get_P1_t(v, a);
			p2 = HD.Get_P2_t(v, a);
		}
		J++;
	}
	std::cout << '\n';
	Out.flush();
	Out.close();
}

void NonLinear(H_System& HD, std::string name)
{
	float
		ot = -dt * 2;

	float dt_solve = 1e-10;

	std::ofstream Out(name);
	Out << std::scientific;
	Out << std::setprecision(15);
	Out
		<< "T" << ","
		<< "X" << ","
		<< "V" << ","
		<< "A" << ","
		<< "U" << "\n";
	int j = 0;

	float I = 0, BI = 0;

	std::cout << std::scientific;
	std::cout << std::setprecision(2);
	//HD.solver.State[2] = 0.5;

	double 
		T = 0,
		x = 0,
		v = 0,
		a = 0,
		u = 0,
		du = 0,
		ddu = 0;

	int medC = maxTime / dt_solve / 250;

	double KK = HD.Get_K();
	double F_0 = HD.Get_F_0();

	double stat[10];
	for(int i = 0, j = 0, J = medC * 2; T <= maxTime; i++)
	{
		HD.U = GetU(T, x, v, dt_solve);
		HD.U = MAX2(-1, MIN2(1, HD.U));
		if(T - ot >= dt)
		{
			while(T - ot >= dt)
				ot += dt;
			Out
				<< T << ",\t"
				<< x << ",\t"
				<< v << ",\t"
				<< a << ",\t"
				<< u << "\n";
			Out.flush();
		}
		if(J > medC)
		{
			for(int i = 0; i < 20; i++)
				std::cout << (i / 19.0 > T / maxTime ? "_" : "X");
			std::cout << dt_solve << ' ' << T << '\r';
			J = 0;
		}
		//solve


		stat[0] = T;
		stat[1] = v;
		stat[2] = x;
		stat[3] = 0;
		stat[4] = 0;
		stat[5] = u;
		stat[4] = 0;
		{
			double nddu = HD.U * HD.U_K / HD.U_T - du /HD.U_T -  u * HD.U_K / HD.U_T;
			double ndu = du + dt_solve * (nddu);
			double nu = u + dt_solve * (ndu);

			double na = -v * 2 / dt_solve - a;
			double nv = 0;
			double nx = x;

			if(u > 0)
				F_0 = HD.S1 * HD.p_input - HD.S2 * HD.p_sliv + HD.force->F(stat, &HD);
			else
				F_0 = HD.S1 * HD.p_sliv - HD.S2 * HD.p_input + HD.force->F(stat, &HD);

			if(u * u > 1e-10)
			{
				na = (-v * abs(v) * KK / (u * u) - HD.b * v + F_0) / HD.m;
				nv = v + dt_solve * (na + a) * 0.5;
				nx = x + dt_solve * (nv + v) * 0.5;
			}

			u = nu;
			du = ndu;

			x = nx;
			v = nv;
			a = na;

			T += dt_solve;
		}
		J++;
	}
	std::cout << '\n';
	Out.flush();
	Out.close();
}

void LinearStab(H_System& HD, std::string name)
{
	float
		ot = -dt * 2;

	float dt_solve = 1e-7;


	float control = 1;

	std::ofstream Out(name);
	Out << std::scientific;
	Out << std::setprecision(15);
	Out
		<< "T" << ","
		<< "X" << ","
		<< "V" << ","
		<< "A" << ","
		<< "P1" << ","
		<< "P2" << ","
		<< "DP1" << ","
		<< "DP2" << ","
		<< "U" << "\n";
	int j = 0;

	float I = 0, BI = 0;

	std::cout << std::scientific;
	std::cout << std::setprecision(2);
	//HD.solver.State[2] = 0.5;


	int medC = maxTime / dt_solve / 250;

	double
		v_r = HD.Get_V_r(),
		T0 = sqrt(2 / HD.ro) * HD.mu * HD.S / v_r,
		x0 = (HD.max_x - HD.min_x) * 0.5,
		p1_r = HD.Get_P1_r(),
		p2_r = HD.Get_P2_r()
		;

	double 
		T = 0,
		x = x0,
		u = -v_r * 0.25,
		v = v_r * 0.75,
		a = 0,
		y = 1,
		p1 = p1_r,
		p2 = p2_r,
		w1 = 0,
		w2 = 0,
		dw1 = 0,
		dw2 = 0;

	for(int i = 0, j = 0, J = medC * 2, KK = 0; T <= maxTime; i++)
	{
		HD.U = GetU(T, x, v, dt_solve);
		HD.U = MAX2(-1, MIN2(1, HD.U));
		if(T - ot >= dt)
		{
			while(T - ot >= dt)
				ot += dt;
			Out
				<< T << ",\t"
				<< x << ",\t"
				<< v << ",\t"
				<< a << ",\t"
				<< p1 << ",\t"
				<< p2 << ",\t"
				<< dw1 << ",\t"
				<< dw2 << ",\t"
				<< y << "\n";
			Out.flush();
		}
		if(J > medC)
		{
			for(int i = 0; i < 20; i++)
				std::cout << (i / 19.0 > T / maxTime ? "_" : "X");
			std::cout << dt_solve << ' ' << T << '\r';
			J = 0;
		}
		//solve
		{
			double ndw1 = HD.S1 * HD.E / (HD.S1 * x0 + HD.V1) * (-abs(v_r) * T0 * T0 / (HD.S1 * HD.S1) - u);
			double ndw2 = HD.S2 * HD.E / (HD.S2 * x0 + HD.V2) * (-abs(v_r) * T0 * T0 / (HD.S2 * HD.S2) + u);

			w1 = w1 + dt_solve * (ndw1 + dw1) * 0.5;
			w2 = w2 + dt_solve * (ndw2 + dw2) * 0.5;

			double na = (w1 * HD.S1 - w2 * HD.S2 - u * HD.b) / HD.m;
			u = u + dt_solve * (na  + a) * 0.5;

			x = x0;
			v = u + v_r;
			a = na;
			p1 = p1_r + w1;
			p2 = p2_r + w2;

			T += dt_solve;

		}
		J++;
	}
	std::cout << '\n';
	Out.flush();
	Out.close();
}



void main()
{

	std::string
		fnames = "",
		out_name_base = "\"base\"";



	H_System HD;

	HD.target_x = 1;

	float K_T = HD.Get_T();
	float K_K = HD.Get_V_r();

	K_d = (30 * K_T - 1) / K_K;
	K_p = (304 * K_T) / K_K;
	K_i = (1040 * K_T) / K_K;

	if(0)
	{
		std::strstream temp;

		temp << std::scientific << std::setprecision(3)<< "stab" << '\0';

		std::string fname = "out_" + std::string(temp.str()) + ".csv";

		temp.freeze(false);

		LinearStab(HD, fname);
		//HD.E *= powf(10, -1.0);
		fnames += " \"" + std::string(temp.str()) + "\"";
	}
	if(1)
	{
		HD.Reset();
		std::strstream temp;

		temp << std::scientific << std::setprecision(3)<< "real" << '\0';

		std::string fname = "out_" + std::string(temp.str()) + ".csv";

		temp.freeze(false);

		Base(HD, fname);
		//HD.E *= powf(10, -1.0);
		fnames += " \"" + std::string(temp.str()) + "\"";
	}




	if(1)
	{
		std::strstream temp;

		temp << std::scientific << std::setprecision(3)<< "linear 1" << '\0';

		std::string fname = "out_" + std::string(temp.str()) + ".csv";

		temp.freeze(false);

		Linear(HD, fname, HD.Get_V_r(), HD.Get_T());
		//HD.E *= powf(10, -1.0);
		fnames += " \"" + std::string(temp.str()) + "\"";
	}

	if(0)
	{
		std::strstream temp;

		temp << std::scientific << std::setprecision(3)<< "linear 2" << '\0';

		std::string fname = "out_" + std::string(temp.str()) + ".csv";

		temp.freeze(false);

		Linear(HD, fname, HD.GetK_old(), HD.m / HD.b);
		//HD.E *= powf(10, -1.0);
		fnames += " \"" + std::string(temp.str()) + "\"";
	}

	if(0)
	{
		std::strstream temp;

		temp << std::scientific << std::setprecision(3)<< "non linear" << '\0';

		std::string fname = "out_" + std::string(temp.str()) + ".csv";

		temp.freeze(false);

		NonLinear(HD, fname);
		//HD.E *= powf(10, -1.0);
		fnames += " \"" + std::string(temp.str()) + "\"";
	}

	system("python manual_multy.py");
	//system(("python double_1.py " + out_name_base + fnames).c_str());

	return;
}

