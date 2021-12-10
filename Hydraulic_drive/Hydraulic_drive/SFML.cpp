#include <iostream>
#include <SFML/Graphics.hpp>
#include <thread>
#include <time.h>
#include <ctime>
#include <ratio>
#include <chrono>

struct PoliLine: public sf::Drawable
{
public:
	sf::Vertex* m_vertices = 0;
	uint32_t count = 0;

	sf::Vertex& operator[](uint32_t i)
	{
		return m_vertices[i % count];
	}

	void Set(sf::Vertex* n_vertices, uint32_t count)
	{
		m_vertices = n_vertices;
		PoliLine::count = count;
	}

	void SetColor(sf::Color color)
	{
		for (int i = 0; i < count; i++)
		{
			m_vertices[i].color = color;
		}
	}

private:
	virtual void draw(sf::RenderTarget& target, sf::RenderStates states) const
	{
		target.draw(m_vertices, count, sf::PrimitiveType::LineStrip, states);
	}
};

#define sign(x) (x == 0 ? 0 : (x > 0 ? 1 : -1))

#define step(x) (x > 0 ? 1 : -1)

#define max(a,b) (a > b ? a : b)
#define min(a,b) (a < b ? a : b)

enum Mode
{
	FRWARD,
	BACK,
	BLOCK,
	SLIV
};

struct DoubleActingHydraulicDrive
{
	// прямой цилиндр
	double
		p1 = 1,
		V1 = 0.05,
		Vt1 = 0.05,
		S1 = 0.1,
		Q1 = 1
		;
	// обратный цилиндр
	double
		p2 = 1,
		V2 = 0.05,
		Vt2 = 0.05,
		S2 = 0.05,
		Q2 = 1
		;
	// Свойства рабочей жидкости
	double
		ro = 860,
		E = 1.56E9
		;
	// Харрактеристики дроссилей
	double
		mu = 0.62,
		S = 0.001
		;
	// Харрактеристики системы
	double
		x = 0,
		min_x = 0,
		max_x = 1,
		v = 0,
		m = 1,
		F = -9.81,
		b_prop = 5,
		b_abs = 0.0,
		c = 0,
		p_s = 0,
		p_i = 10000,
		dv = 0,
		dp1 = 0,
		dp2 = 0
		;


	Mode mode = BLOCK;

	void Step(double dt)
	{
		switch (mode)
		{
		case FRWARD:
			Q1 = mu * S * step(p_i - p1) * sqrt(2 / ro * abs(p_i - p1));
			Q2 = -mu * S * step(p2 - p_s) * sqrt(2 / ro * abs(p_s - p2));
			break;
		case BACK:
			Q2 = mu * S * sign(p_i - p2) * sqrt(2 / ro * abs(p_i - p2));
			Q1 = -mu * S * step(p1 - p_s) * sqrt(2 / ro * abs(p_s - p1));
			break;
		case BLOCK:
			Q1 = 0;
			Q2 = 0;
		case SLIV:
			Q1 = -mu * S * step(p1 - p_s) * sqrt(2 / ro * abs(p_s - p1));
			Q2 = -mu * S * step(p2 - p_s) * sqrt(2 / ro * abs(p_s - p2));
			break;
		default:
			break;
		}
		Vt1 = V1 + S1 * (x - min_x);
		Vt2 = V2 + S2 * (max_x - x);

		dp1 = E / Vt1 * (Q1 - S1 * v);
		dp2 = E / Vt2 * (Q2 + S2 * v);
		dv = (p1 * S1 - p2 * S2 - c * x - b_prop * v + F) / m;

		p1 += dp1 * dt;
		p2 += dp2 * dt;
		//p1 = max(p1, p_s);
		//p2 = max(p2, p_s);
		double tv = v;
		v += dv * dt;
		double _v = b_abs * dt / m;
		v = (_v > abs(v) ? 0 : v - sign(v) * _v);
		dv = (_v > abs(v) ? 0 : dv - sign(v) * _v / dt);
		x += tv * dt + dv * dt * dt * 0.5;
		x = max(min(x, max_x), min_x);
		//v = ((v - sign(v) * b_abs * dt * 0.5 / m) * v > 0 ? v - sign(v) * b_abs / m * 0.5 * dt : 0);
		if (x == min_x && v < 0) v = 0;
		if (x == max_x && v > 0) v = 0;
		//if (dp1 > 2 * abs(p_i - p_s) / dt)
		//    std::cout << "WArNING!\n- dp/dt extreme values : " + std::to_string(dp1) + ".\n- This is " + std::to_string(abs(dp1) / abs(p_i - p_s) * dt) + " times more than (Pi - Ps)/dt.\n";

		//if (dp2 > 2 * abs(p_i - p_s) / dt)
		//    std::cout << "WArNING!\n- dp/dt extreme values : " + std::to_string(dp2) + ".\n- This is " + std::to_string(abs(dp2) / abs(p_i - p_s) * dt) + " times more than (Pi - Ps)/dt.\n";
	}
};


struct DoubleActingHydraulicDrive_NoDp
{
	// прямой цилиндр
	double
		S1 = 0.1
		;
	// обратный цилиндр
	double
		S2 = 0.1
		;
	// Свойства рабочей жидкости
	double
		ro = 860
		;
	// Харрактеристики дроссилей
	double
		mu = 0.62,
		S = 0.01
		;
	// Харрактеристики системы
	double
		x = 0,
		min_x = 0,
		max_x = 1,
		v = 0,
		m = 1,
		F = -0,
		b = 5,
		c = 0,
		p_s = 1,
		p_i = 10000,
		dv = 0
		;
	Mode mode = BLOCK;
	void Step(double dt)
	{
		double p1, p2;
		switch (mode)
		{
		case FRWARD:
			p1 = p_i;
			p2 = p_s;
			break;
		case BACK:
			p2 = p_i;
			p1 = p_s;
			break;
		case BLOCK:
			p1 = 0;
			p2 = 0;
			break;
		case SLIV:
			p2 = p_s;
			p1 = p_s;
			break;
		default:
			break;
		}

		double
			K0 = ro * (S1 * S1 * S1 + S2 * S2 * S2) / (2 * mu * mu * S * S),
			F0 = +p1 * S1 - p2 * S2 + F,
			C0 = sqrt(b*b + 4 * K0 * abs(F0)),
			mv = (sign(F0) * (C0 - b)) / (2 * K0);

		dv = (-v * abs(v) * K0 - b * v + F0) / m;

		double tv = v;
		if (abs(v) > abs(mv) && v*mv > 0)
		{
			v = mv;
			dv = (v - tv) / dt;
		}
		else
			v += dv * dt;
		if (abs(v) > abs(mv) && v*mv > 0)
		{
			v = mv;
			dv = (v - tv) / dt;
		}
		if (mode == BLOCK || S == 0)
		{
			dv = 0;
			v = 0;
		}
		x += tv * dt + dv * dt * dt * 0.5;
		x = max(min(x, max_x), min_x);
		if (x == min_x && v < 0) v = 0;
		if (x == max_x && v > 0) v = 0;
	}
};



void main()
{

	uint32_t H = 800, W = 400;

	sf::ContextSettings context_setting(0, 0, 2);
	sf::RenderWindow window(sf::VideoMode(H, W), "SFML window", sf::Style::Default, context_setting);
	sf::CircleShape shape(100.f);
	shape.setFillColor(sf::Color::Green);

	PoliLine X, Y, F;

	X.Set(new sf::Vertex[2], 2);
	Y.Set(new sf::Vertex[2], 2);
	F.Set(new sf::Vertex[500], 500);

	X.SetColor(sf::Color::White);
	Y.SetColor(sf::Color::White);
	F.SetColor(sf::Color::Red);

	X[0].position = {float(10), float(W - 10)};
	X[1].position = {float(H - 10), float(W - 10)};

	Y[0].position = {float(10), float(W - 10)};
	Y[1].position = {float(10), float(10)};

	for (int i = 0; i < 500; i++)
	{
		F[i].position = {float(10 + i * (H - 20) / 499.0), float(W * 0.5)};
	}

	DoubleActingHydraulicDrive_NoDp A;
	A.mode = FRWARD;
	double dt = 1E-2;
	int I = 0, _I = 1 / dt;
	//while ((A.x - A.min_x) / (A.max_x - A.min_x) < 0.0 - 0.000005)
	//{
	//	A.Step(dt);
	//	if ((I = (I + 1) % _I) == 0)
	//		std::cout << A.x << "      \r";
	//}
	double
		S = A.S
		;
	for (int i = 0, K = 0; i < 500; i++)
	{
		F[i].position.y = (W - 10) - (W - 20) * (((A.x - A.min_x) / (A.max_x - A.min_x) * 1 - 0.5) * 1 + 0.5);
		if (i > 200)
			K++;
		if (K > 0)
		{
			A.S = S * max(0, min(1, (abs(K - 10) / 5.0 - 1)));
		}
		if (K == 10)
		{
			A.mode = Mode::BACK;

			//A.S = 0;
		}
		for (int t = 0; t < 1; t++)
		{
			A.Step(dt);
		}
		std::cout << A.x << "      \r";
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
		window.draw(F);
		window.display();
	}

}
