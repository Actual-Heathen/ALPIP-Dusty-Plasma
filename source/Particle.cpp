#include "../headers/Particle.h"

Particle::Particle()
{
    setX(0);
    setY(0);
    setMass(0);
}

Particle::Particle(float x, float y, float m)
{
    setX(x);
    setY(y);
    setMass(m);
}

void Particle::setX(float x)
{
    x_pos = x;
}

void Particle::setY(float y)
{
    y_pos = y;
}

void Particle::setMass(float m)
{
    mass = m;
}

float Particle::getX()
{
    return x_pos;
}

float Particle::getY()
{
    return y_pos;
}

float Particle::getMass()
{
    return mass;
}