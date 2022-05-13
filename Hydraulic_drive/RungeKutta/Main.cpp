#include <iostream>
#include <fstream>
#include "Hidro.h"
#include <vector>
#include <string>

void PrintHD(std::ofstream *fout, H_System &HD)
{

}


void Var2()
{
	H_System HD;

	HD.solver.h = 1e-6;

	float
		dt = 0.001,
		ot = -dt * 2;

	float maxTime = 0.5;

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
	HD.K_p = 3;
	HD.solver.h = 1e-7;

	float
		dt = 0.001,
		ot = -dt * 2;

	float maxTime = 1;

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
	HD.solver.State[2] = HD.max_x * 0.1;
	HD.solver.State[1] = 0
		;
	HD.target_x = HD.max_x * 0.9;

	std::cout << std::scientific;
	std::cout << std::setprecision(2);
	//HD.solver.State[2] = 0.5;

	std::cout << HD.target_x << std::endl;

	int medC = maxTime / HD.solver.h / 250;
	HD.target_x = HD.max_x * 0.9;
	for(int i = 0, j = 0, J = medC * 2, KK = 0; T < maxTime; i++)
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
				<< HD.U << ","
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

void Var1_2(H_System& HD, std::string name)
{
	HD.Reset();
	Forse_manipulator *force = dynamic_cast<Forse_manipulator *>(HD.force);
	HD.U = 0;
	HD.K_p = 3;
	HD.solver.h = 1e-7;

	float
		dt = 0.001,
		ot = -dt * 2;

	float maxTime = 1;

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
	HD.solver.State[2] = HD.max_x * 0.1;
	HD.solver.State[1] = 0
		;
	HD.target_x = HD.max_x * 0.9;

	std::cout << std::scientific;
	std::cout << std::setprecision(2);
	//HD.solver.State[2] = 0.5;

	std::cout << HD.E << std::endl;

	int medC = maxTime / HD.solver.h / 250;
	HD.target_x = HD.max_x * 0.9;
	for(int i = 0, j = 0, J = medC * 2, KK = 0; T < maxTime; i++)
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
				<< HD.U << ","
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
}


void main()
{
	H_System HD;

	HD.E = 1.56E9;
	for(int i = 0; i < 4; i++)
	{
		Var1_2(HD, "out_" + std::to_string(i + 1) + ".csv");
		HD.E /= 10;
	}
}

