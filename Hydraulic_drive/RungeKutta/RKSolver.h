#pragma once



template<unsigned char Size>
struct RKSolver
{
	static double BaseFunction(double* A)
	{
		return 1;
	}

	typedef double (*FuncPtr)(double*);

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

	void Calc()
	{
		double K[Size][4];
		double select[Size + 1];

		// K1
		for (int i = 0; i < Size + 1; i++)
			select[i] = State[i];

		for (int i = 0; i < Size; i++)
			K[i][0] = h * funcs[i](select);

		// K2
		select[0] = State[0] + h * 0.5;
		for (int i = 1; i < Size + 1; i++)
			select[i] = State[i] + K[i - 1][0] * 0.5;

		for (int i = 0; i < Size; i++)
			K[i][1] = h * funcs[i](select);

		// K3
		select[0] = State[0] + h * 0.5;
		for (int i = 1; i < Size + 1; i++)
			select[i] = State[i] + K[i - 1][1] * 0.5;

		for (int i = 0; i < Size; i++)
			K[i][2] = h * funcs[i](select);

		// K4
		select[0] = State[0] + h;
		for (int i = 1; i < Size + 1; i++)
			select[i] = State[i] + K[i - 1][2];

		for (int i = 0; i < Size; i++)
			K[i][3] = h * funcs[i](select);

		State[0] = State[0] + h;
		for (int i = 1; i < Size + 1; i++)
		{
			State[i] = State[i] + (K[i - 1][0] + 2 * K[i - 1][1] + 2 * K[i - 1][2] + K[i - 1][3]) / 6.0;
		}
	}

};


