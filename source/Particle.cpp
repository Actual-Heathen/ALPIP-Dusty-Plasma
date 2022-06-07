#include "../headers/Particle.h"
#include "../headers/Point.h"
#include <math.h>
#include <vector>
#include <stdlib.h>
#include <iostream>

#define PI 3.14159265
#define G 6.6743*pow(10,-11)

Particle::Particle()
{
    setX(0);
    setY(0);
    setMass(0);
}

Particle::Particle(long double x, long double y, long double m)
{
    setX(x);
    setY(y);
    setMass(m);
}

void Particle::setX(long double x)
{
    x_pos = x;
}

void Particle::setY(long double y)
{
    y_pos = y;
}

void Particle::setMass(long double m)
{
    mass = m;
}

long double Particle::getX()
{
    return x_pos;
}

long double Particle::getY()
{
    return y_pos;
}

long double Particle::getMass()
{
    return mass;
}
long double Particle::getSpeed()
{
    return speed;
}

void Particle::addAcceleration(double spacing, std::vector<std::vector<double>> rho)
{
    int iXm = abs(floor(x_pos/spacing));
    int iXp = iXm + 1;
    if (iXp >= rho.size())
    {
        iXp = iXp-rho.size();
    }


    int iYm = floor(y_pos/spacing);
    int iYp = iYm + 1;
    if (iYp >=rho.size())
    {
        iYp = iYp-rho.size();
    }
    //cout << iXp<< ", "<< wXp << "<-ixp "<< iYp << ", "<< wYp << "<-iyp\n";
    //cout << iXm<< " ,"<< wXm << "<-ixm "<< iYm << ", "<< wYm << "<-iym\n";

    double dX = iXm*spacing-x_pos;
    double dY = iYm*spacing-y_pos;
    double d = sqrt(dX*dX+dY*dY);
    double psi = (4*PI*G*rho[iXm][iYm])*exp(-2*d);
    dX = (dX/d)*psi;
    dY = (dY/d)*psi;

    accelX += dX;
    accelY += dY;

    dX = iXm*spacing-x_pos;
    dY = iYp*spacing-y_pos;
    d = sqrt(dX*dX+dY*dY);
    psi = (4*PI*G*rho[iXm][iYp])*exp(-2*d);
    dX = (dX/d)*psi;
    dY = (dY/d)*psi;

    accelX += dX;
    accelY += dY;

    dX = iXp*spacing-x_pos;
    dY = iYm*spacing-y_pos;
    d = sqrt(dX*dX+dY*dY);
    psi = (4*PI*G*rho[iXp][iYm])*exp(-2/d);
    dX = (dX/d)*psi;
    dY = (dY/d)*psi;

    accelX += dX;
    accelY += dY;

    dX = iXp*spacing-x_pos;
    dY = iYp*spacing-y_pos;
    d = sqrt(dX*dX+dY*dY);
    psi = (4*PI*G*rho[iXp][iYp])*exp(-2/d);
    dX = (dX/d)*psi;
    dY = (dY/d)*psi;

    accelX += dX;
    accelY += dY;
    //std::cout<< accelX<<"<-accel\n";
}

void Particle::move(double s, double size)
{
    long double dX = velX*s + .5*accelX*s*s;
    long double dY = velY*s + .5*accelY*s*s;

    long double xV = dX/s;
    long double yV = dY/s;

    velX = xV;
    velY = yV;
    speed = sqrt(xV*xV + yV*yV);

    x_pos += dX;
    y_pos += dY;
    //std::cout<<accelX<<","<<dX<<"<-dX\n";

    resetAcc();
    accelX = 0;
    accelY = 0;
    accel = 0;

    if(x_pos > size)
    {
        x_pos = x_pos-size;
    }
    if(y_pos > size)
    {
        y_pos = y_pos-size;
    }
    if(x_pos < 0)
    {
        x_pos = x_pos+size;
    }
    if(y_pos < 0)
    {
        y_pos = y_pos+size;
    }
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