#include "../headers/Particle.h"
#include <math.h>

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
float Particle::getSpeed()
{
    return speed;
}

void Particle::addAcceleration(float x, float y)
{
    accelX += x;
    accelY += y;
    accel = sqrt(pow(x,2) + pow(y,2));
}

void Particle::move(float s)
{
    float dX = velX*s + .5*accelX*s*s;
    float dY = velY*s + .5*accelY*s*s;

    velX = dX;
    velY = dY;

    float xV = dX/s;
    float yV = dY/s;
    speed = sqrt(xV*xV + yV*yV);

    x_pos += dX;
    y_pos += dY;

    resetAcc();
}

void Particle::resetAcc()
{
    accelX = 0;
    accelY = 0;
    accel = 0;
}

void Particle::resetVel()
{
    velX = 0;
    velY = 0;
    speed = 0;
}