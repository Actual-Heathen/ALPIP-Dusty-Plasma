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

Particle::Particle(long double x, long double y, long double m, double vx, double vy)
{
    setX(x);
    setY(y);
    setMass(m);
    velX = vx;
    velY = vy;
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

void Particle::addAcceleration(double spacing, std::vector<std::vector<double>> dpsix, std::vector<std::vector<double>> dpsiy)
{
    int iXm = floor(x_pos/spacing);                 //calculate iXm & iXp
    int iXp = iXm + 1;

    double wXm = 1- abs((x_pos/spacing)-iXm);       //weight calculations
    double wXp = 1- abs((x_pos/spacing)-iXp);

    if (iXp >= dpsix.size())                                 //adjust iXp if in "ghost region"
    {
        iXp = iXp - dpsix.size();
    }

    int iYm = floor(y_pos/spacing);                 //calculate iYm and iYp
    int iYp = iYm + 1;

    double wYm = 1- abs((y_pos/spacing)-iYm);       //wieght calculation
    double wYp = 1- abs((y_pos/spacing)-iYp);

    if (iYp >= dpsiy.size())                                 //ghost region adjust
    {
        iYp = iYp - dpsiy.size();
    }

                //W I P//W I P//W I P// //W I P//W I P//calculate psi and use it to add accleleration will change//W I P//W I P//W I P// //W I P//W I P//
    double dX = wXm*wYm*dpsix[iXm][iYm];
    double dY = wXm*wYm*dpsiy[iXm][iYm];

    accelX += dX;
    accelY += dY;

    dX = wXm*wYp*dpsix[iXm][iYp];
    dY = wXm*wYp*dpsiy[iXm][iYp];

    accelX += dX;
    accelY += dY;


    dX = wXp*wYm*dpsix[iXp][iYm];
    dY = wXp*wYm*dpsiy[iXp][iYm];

    accelX += dX;
    accelY += dY;

    dX = wXp*wYp*dpsix[iXp][iYp];
    dY = wXp*wYp*dpsiy[iXp][iYp];

    accelX += dX;
    accelY += dY;

}

void Particle::move(double s, double size)
{
    // using d= v0t+ 1/2at^2 to calculate distance traveled
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

    if(x_pos > size)        //loop particle around periodically
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
