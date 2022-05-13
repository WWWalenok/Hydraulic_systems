#pragma once
#include <thread>
#include <mutex>



#define MAX2(a,b) (a > b ? a : b)

#define MIN2(a,b) (a < b ? a : b)

template<unsigned char Size>
struct RKSolver
{
	static double BaseFunction(double* A, void* B)
	{
		return 1;
	}

	typedef double (*FuncPtr)(double*, void*);

	FuncPtr funcs[Size];

	double h = 1;

	double State[Size + 1];

	RKSolver()
	{
		for (int i = 0; i < Size; i++)
		{
			funcs[i] = BaseFunction;
		}
		for (int i = 0; i < Size + 1; i++)
		{
			State[i] = 0;
		}
	}

private:

	double K[Size][4];
	double select[Size + 1];

public:

	void Calc(void* _t)
	{
		for (int i = 0; i < Size + 1; i++)
			select[i] = State[i];

		// K1
		for (int i = 0; i < Size + 1; i++)
			select[i] = State[i];

		for (int i = 0; i < Size; i++)
			K[i][0] = h * funcs[i](select, _t);

		// K2
		select[0] = State[0] + h * 0.5;
		for (int i = 1; i < Size + 1; i++)
			select[i] = State[i] + K[i - 1][0] * 0.5;

		for (int i = 0; i < Size; i++)
			K[i][1] = h * funcs[i](select, _t);

		// K3
		select[0] = State[0] + h * 0.5;
		for (int i = 1; i < Size + 1; i++)
			select[i] = State[i] + K[i - 1][1] * 0.5;

		for (int i = 0; i < Size; i++)
			K[i][2] = h * funcs[i](select, _t);

		// K4
		select[0] = State[0] + h;
		for (int i = 1; i < Size + 1; i++)
			select[i] = State[i] + K[i - 1][2];

		for (int i = 0; i < Size; i++)
			K[i][3] = h * funcs[i](select, _t);

		State[0] = State[0] + h;
		for (int i = 1; i < Size + 1; i++)
		{
			State[i] = State[i] + (K[i - 1][0] + 2 * K[i - 1][1] + 2 * K[i - 1][2] + K[i - 1][3]) / 6.0;
		}
	}

	void UpateH(void* _t, float bh = 1, float eps = 1e-7)
	{
		double OState[Size + 1];
		double State1[Size + 1];
		double State2[Size + 1];

		for (int i = 0; i < Size + 1; i++)
			OState[i] = State[i];
		h = bh;
		float err = eps * 2;
		double oh = h;
		while (err > eps)
		{
			oh = h;
			err = 0;

			Calc(_t);

			for (int i = 0; i < Size + 1; i++)
				State1[i] = State[i];

			for (int i = 0; i < Size + 1; i++)
				State[i] = OState[i];

			h = h * 0.5;

			Calc(_t);
			Calc(_t);

			for (int i = 0; i < Size + 1; i++)
				State2[i] = State[i];

			for (int i = 0; i < Size + 1; i++)
				State[i] = OState[i];
			
			for (int i = 1; i < Size + 1; i++)
			{
				err = MAX2(err, fabs(State1[i] - State2[i]));
			}
		}

		h = oh;

	}

};


