#include "../headers/Particle.h"
#include "../headers/Point.h"
#include <math.h>

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

void Particle::addAcceleration(double rho[][3], double spacing)
{
        int iXm = abs(floor(x_pos/spacing));
        int iXp = iXm + 1;


        int iYm = floor(y_pos/spacing);
        int iYp = iYm + 1;

        
        //cout << iXp<< ", "<< wXp << "<-ixp "<< iYp << ", "<< wYp << "<-iyp\n";
        //cout << iXm<< " ,"<< wXm << "<-ixm "<< iYm << ", "<< wYm << "<-iym\n";

        double dX = iXm*spacing-x_pos;
        double dY = iYm*spacing-y_pos;
        double d = sqrt(dX*dX+dY*dY);
        double psi = 4*PI*G*rho[iXm][iYm];
        dX = (dX/d)*psi;
        dY = (dY/d)*psi;

        accelX += dX;
        accelY += dY;

        dX = iXm*spacing-x_pos;
        dY = iYp*spacing-y_pos;
        d = sqrt(dX*dX+dY*dY);
        psi = 4*PI*G*rho[iXm][iYp];
        dX = (dX/d)*psi;
        dY = (dY/d)*psi;

        accelX += dX;
        accelY += dY;

        dX = iXp*spacing-x_pos;
        dY = iYm*spacing-y_pos;
        d = sqrt(dX*dX+dY*dY);
        psi = 4*PI*G*rho[iXp][iYm];
        dX = (dX/d)*psi;
        dY = (dY/d)*psi;

        accelX += dX;
        accelY += dY;

        dX = iXp*spacing-x_pos;
        dY = iYp*spacing-y_pos;
        d = sqrt(dX*dX+dY*dY);
        psi = 4*PI*G*rho[iXp][iYp];
        dX = (dX/d)*psi;
        dY = (dY/d)*psi;

        accelX += dX;
        accelY += dY;
}

void Particle::move(long double s)
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

    resetAcc();
    accelX = 0;
    accelY = 0;
    accel = 0;
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