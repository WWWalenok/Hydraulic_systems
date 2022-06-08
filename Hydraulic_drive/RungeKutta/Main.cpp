#include <iostream>
#include <fstream>
#include "Hidro.h"
#include <vector>
#include <string>
#include <strstream>

void PrintHD(std::ofstream *fout, H_System &HD)
{

}

float maxTime = 5 * 1.000;
float dt = 0.001;

void Var2()
{
	H_System HD;

	HD.solver.h = 1e-6;

	float
		ot = -dt * 2;

	double &T = HD.solver.State[0];

	for(int a = -5; a <= 5; a++) for(int b = -5; b <= 5; b++) for(int c = -5; c <= 5; c++)
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
		for(int i = 0, j = 0, J = 0, KK = 0; T <= maxTime; i++)
		{
			J++;
			if(i % 100 == 0)
			{
				fout
					<< HD.solver.State[0] << ",\t"
					<< HD.solver.State[1] << ",\t"
					<< HD.solver.State[2] << ",\t"
					<< HD.solver.State[3] << ",\t"
					<< HD.solver.State[4] << ",\t"
					<< DV(HD.solver.State, &HD) << "\n"
					;
			}
			HD.Calc();
		}
	}

}

void Var1()
{
	H_System HD;
	HD.Reset();
	Forse_manipulator *force = dynamic_cast<Forse_manipulator *>(HD.force);
	HD.U = 0;
	HD.K_p = 15;
	HD.solver.h = 1e-7;

	float
		ot = -dt * 2;

	double &T = HD.solver.State[0];

	force->F(HD.solver.State, &HD);
	double ba = acos(force->cosa);
	std::ofstream Out("out.csv");
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
	int j = 0;

	float I = 0, BI = 0;

	double mh = HD.solver.h;
	HD.solver.State[2] = HD.max_x * 0.01;
	HD.solver.State[1] = 0
		;
	HD.target_x = HD.max_x * 0.1;

	std::cout << std::scientific;
	std::cout << std::setprecision(2);
	//HD.solver.State[2] = 0.5;

	std::cout << HD.Get_K() << std::endl;

	int medC = maxTime / HD.solver.h / 250;
	HD.target_x = HD.max_x * 0.1;
	for(int i = 0, j = 0, J = medC * 2, KK = 0; T <= maxTime + dt * 10; i++)
	{
		force->F(HD.solver.State, &HD);
		if(T - ot >= dt)
		{
			while(T - ot >= dt)
				ot += dt;
			Out
				<< HD.solver.State[0] << ",\t"
				<< HD.solver.State[1] << ",\t"
				<< HD.solver.State[2] << ",\t"
				<< HD.solver.State[3] * HD.S1 << ",\t"
				<< HD.solver.State[4] * HD.S2 << ",\t"
				<< HD.solver.State[3] * HD.S1 - HD.solver.State[4] * HD.S2 + force->Fn << ",\t"
				<< DV(HD.solver.State, &HD) << ",\t"
				<< force->Fn << ","
				<< (acos(force->cosa) - ba) * 180 / PI << ","
				<< HD.U << ","
				<< HD.solver.State[5] << ","
				<< HD.solver.h << "\n";;
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

	system("python main.py");
}

float GetU(float T, int& j)
{
	float TT = 0.5 * T;
	int a = int(TT) % 4;
	float b = TT - int(TT);

	float U = 0;
	if(0)
	switch(a)
	{
		case 0:
		return b;
		case 1:
		return 1 - b;
		case 2:
		return - b;
		case 3:
		return -1 + b;
		default:
		break;
	}

	if(T < 0.1)
		U = 0;
	else if(T < 1.1)
		U = 1;
	else if(T < 1.6)
		U = 0;
	else if(T < 2.6)
		U = -1;
	else if(T < 3)
		U = 0;
	else if(T < 3.5)
		U = 0.5;
	else if(T < 4)
		U = -0.5;
	else if(T < 10)
		U = 0;
	else if(T < 10)
		U = 0;

	return U;

	switch(j)
	{
		case 0:
		{
			U = 1;
			if(T > 1)
				j++;
			break;
		}
		case 1:
		{
			U = 0;
			if(T > 1.5)
				j++;
			break;
		}
		case 2:
		{
			U = -1;
			if(T > 2)
				j++;
			break;
		}
		case 3:
		{
			U = 0.5;
			if(T > 3)
				j++;
			break;
		}
		case 4:
		{
			U = -0.5;
			if(T > 4)
				j++;
			break;
		}
		default:
		break;
	}
	return U;
}

void Var1_2(H_System& HD, std::string name)
{
	HD.Reset();
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
	HD.solver.State[2] = 2;
	HD.solver.State[1] = 0;

	int medC = maxTime / HD.solver.h / 250;

	for(int i = 0, j = 0, J = medC * 2, KK = 0; T <= maxTime; i++)
	{
		HD.U = GetU(T, j);
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
				<< (HD.solver.State[3] * HD.S1 - HD.solver.State[4] * HD.S2 + force->Fn - HD.solver.State[1] * HD.b_prop) / HD.m << ",\t"
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

void Var_LLin(H_System& HD, std::string name, float P_K, float P_T)
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
		<< "U" << "\n";
	int j = 0;

	float I = 0, BI = 0;

	std::cout << std::scientific;
	std::cout << std::setprecision(2);
	//HD.solver.State[2] = 0.5;

	double 
		T = 0,
		x = 2,
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
		HD.U = GetU(T, j);
		//HD.U = (HD.target_x - x) * HD.K_p;
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

void Var_Norm(H_System& HD, std::string name)
{
	float
		ot = -dt * 2;

	float dt_solve = 1e-10;


	float control = 1;

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
		t = 0,
		x = 0,
		v = 0,
		a = 0,
		u = 0,
		du = 0,
		ddu = 0;

	int medC = maxTime / dt_solve / 250;

	double KK = HD.Get_K();
	double F_0 = HD.Get_F_0();


	for(int i = 0, j = 0, J = medC * 2; t <= maxTime; i++)
	{
		HD.U = (HD.target_x - x) * HD.K_p;
		HD.U = MAX2(-1, MIN2(1, HD.U));
		if(t - ot >= dt)
		{
			while(t - ot >= dt)
				ot += dt;
			Out
				<< t << ",\t"
				<< x << ",\t"
				<< v << ",\t"
				<< a << ",\t"
				<< u << "\n";
			Out.flush();
		}
		if(J > medC)
		{
			for(int i = 0; i < 20; i++)
				std::cout << (i / 19.0 > t / maxTime ? "_" : "X");
			std::cout << dt_solve << ' ' << t << '\r';
			J = 0;
		}
		//solve
		{
			double nddu = control * HD.U_K / HD.U_T - du /HD.U_T -  u * HD.U_K / HD.U_T;
			double ndu = du + dt_solve * (nddu);
			double nu = u + dt_solve * (ndu);

			double na = -v * 2 / dt_solve - a;
			double nv = 0;
			double nx = x;

			if(u * u > 1e-10)
			{
				na = (-v * abs(v) * KK / (u * u) - HD.b_prop * v + F_0) / HD.m;
				nv = v + dt_solve * (na + a) * 0.5;
				nx = x + dt_solve * (nv + v) * 0.5;
			}

			u = nu;
			du = ndu;

			x = nx;
			v = nv;
			a = na;

			t += dt_solve;
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

	{
		std::strstream temp;

		temp << std::scientific << std::setprecision(3)<< "real" << '\0';

		std::string fname = "out_" + std::string(temp.str()) + ".csv";

		temp.freeze(false);

		Var1_2(HD, fname);
		HD.E *= powf(10, -1.0);
		fnames += " \"" + std::string(temp.str()) + "\"";
	}

	if(false)
	{
		std::strstream temp;

		temp << std::scientific << std::setprecision(3)<< "linear" << '\0';

		std::string fname = "out_" + std::string(temp.str()) + ".csv";

		temp.freeze(false);

		Var_LLin(HD, fname, HD.Get_V_r(), HD.Get_T());
		HD.E *= powf(10, -1.0);
		fnames += " \"" + std::string(temp.str()) + "\"";
	}

	if(false)
	{
		std::strstream temp;

		temp << std::scientific << std::setprecision(3)<< "non linear" << '\0';

		std::string fname = "out_" + std::string(temp.str()) + ".csv";

		temp.freeze(false);

		Var_Norm(HD, fname);
		HD.E *= powf(10, -1.0);
		fnames += " \"" + std::string(temp.str()) + "\"";
	}

	system("python main.py ");
	//system(("python double_1.py " + out_name_base + fnames).c_str());

	return;

	HD.E = 2.062E9;
	for(int i = 0; i < 2; i++)
	{
		std::strstream temp;

		temp << std::scientific << std::setprecision(3)<< "E = " << HD.E << '\0';

		std::string fname = "out_" + std::string(temp.str()) + ".csv";

		temp.freeze(false);

		Var1_2(HD, fname);
		HD.E *= powf(10, -1.0);
		fnames += " \"" + std::string(temp.str()) + "\"";
	}

	system(("python multi.py " + out_name_base + fnames).c_str());
	//system(("python multy_vToy.py " + out_name_second + fnames).c_str());
}

