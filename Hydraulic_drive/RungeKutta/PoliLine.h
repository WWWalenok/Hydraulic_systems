#pragma once
#include <SFML/Graphics.hpp>

struct PoliLine : public sf::Drawable
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
