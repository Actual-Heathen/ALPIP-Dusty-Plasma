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

Particle::Particle( double x,  double y,  double m, double vx, double vy)
{
    setX(x);
    setY(y);
    setMass(m);
    velX = vx;
    velY = vy;
    accel = 0;
    accelX = 0;
    accelY = 0;
    charge = 1;
    speed = sqrt(pow(velX,2)+pow(velY,2));
}

void Particle::setX( double x)
{
    x_pos = x;
}

void Particle::setY( double y)
{
    y_pos = y;
}

void Particle::setMass( double m)
{
    mass = m;
}

 double Particle::getX()
{
    return x_pos;
}

 double Particle::getY()
{
    return y_pos;
}

 double Particle::getMass()
{
    return mass;
}
 double Particle::getSpeed()
{
    return speed;
}

void Particle::addAcceleration(double spacing, std::vector<std::vector<double>> dpsix, std::vector<std::vector<double>> dpsiy, std::vector<std::vector<double>> dphix, std::vector<std::vector<double>> dphiy, double E[3], double B[3], double time)
{

    int iXm = floor(x_pos/spacing);                 //calculate iXm & iXp
    int iXp = iXm + 1;

    double wXm = 1- abs((x_pos/spacing)-iXm);       //weight calculations
    double wXp = 1- abs((x_pos/spacing)-iXp);

    if ((double)iXp >= dpsix.size())                                 //adjust iXp if in "ghost region"
    {
        iXp = iXp - dpsix.size();
    }

    int iYm = floor(y_pos/spacing);                 //calculate iYm and iYp
    int iYp = iYm + 1;

    double wYm = 1- abs((y_pos/spacing)-iYm);       //wieght calculation
    double wYp = 1- abs((y_pos/spacing)-iYp);

    //std::cout << dpsix[iXm][iYm] << "wXm\n"; 

    if ((double)iYp >= dpsiy.size())                                 //ghost region adjust
    {
        iYp = iYp - dpsiy.size();
    }

    //std::cout <<iYm << "," <<iXp<<"yxmp\n";

    double dX = wXm*wYm*dpsix[iXm][iYm];
    double dY = wXm*wYm*dpsiy[iXm][iYm];
    double dPX = wXm*wYm*dphix[iXm][iYm];
    double dPY = wXm*wYm*dphiy[iXm][iYm];

    accelX = dX;
    accelY = dY;
    double electricX = -1*dPX;
    double electricY = -1*dPY;
    //std::cout <<"mm\n";

    dX = wXm*wYp*dpsix[iXm][iYp];
    dY = wXm*wYp*dpsiy[iXm][iYp];
    dPX = wXm*wYp*dphix[iXm][iYp];
    dPY = wXm*wYp*dphiy[iXm][iYp];

    accelX += dX;
    accelY += dY;
    electricX -= dPX;
    electricY -= dPY;

    //std::cout <<"mp\n";

    dX = wXp*wYm*dpsix[iXp][iYm];
    dY = wXp*wYm*dpsiy[iXp][iYm];
    dPX = wXp*wYm*dphix[iXp][iYm];
    dPY = wXp*wYm*dphiy[iXp][iYm];

    accelX += dX;
    accelY += dY;
    electricX -= dPX;
    electricY -= dPY;

    //std::cout <<"pm\n";

    dX = wXp*wYp*dpsix[iXp][iYp];
    dY = wXp*wYp*dpsiy[iXp][iYp];
    dPX = wXp*wYp*dphix[iXp][iYp];
    dPY = wXp*wYp*dphiy[iXp][iYp];

    accelX += dX;
    accelY += dY;
    electricX -= dPX;
    electricY -= dPY;

    double vM[3] = {velX + ((charge/mass)*(E[0]+electricX)+accelX)*(time/2), velY + ((charge/mass)*(E[1]+electricY)+accelY)*(time/2),velZ + ((charge/mass)*E[2]+accelZ)*(time/2)};
    double T[3] = {(charge/mass)*B[0]*(time/2),(charge/mass)*B[1]*(time/2),(charge/mass)*B[2]*(time/2)};
    double vS[3] = {vM[0]+vM[1]*T[2]-vM[2]*T[1], vM[1]+vM[2]*T[0]-vM[0]*T[2], vM[2]+vM[0]*T[1]-vM[1]*T[0]};
    double S[3] = {(2*T[0])/(1+pow(T[0],2)+pow(T[1],2)+pow(T[2],2)), (2*T[1])/(1+pow(T[0],2)+pow(T[1],2)+pow(T[2],2)), (2*T[2])/(1+pow(T[0],2)+pow(T[1],2)+pow(T[2],2))};
    double vP[3] = {vM[0]+vS[1]*S[2]-vS[2]*S[1], vM[1]+vS[2]*S[0]-vS[0]*S[2], vM[2]+vS[0]*S[1]-vS[1]*S[0]};

    velX = vP[0] + ((charge/mass)*(electricX+E[0])+accelX)*(time/2);
    velY = vP[1] + ((charge/mass)*(electricY+E[1])+accelY)*(time/2);
    velZ = vP[2] + ((charge/mass)*(E[2])+accelZ)*(time/2);
    

    speed = sqrt(pow(velX,2) + pow(velY, 2) + pow(velZ, 2));
}

void Particle::move(double s, double size)
{


    x_pos += velX*s;
    y_pos += velY*s;
 

    resetAcc();
    accelX = 0.0;
    accelY = 0.0;
    accel = 0.0;

    if(x_pos >= size)        //loop particle around periodically
    {
        x_pos = fmod(x_pos,size);
    }
    if(y_pos >= size)
    {
        y_pos = fmod(y_pos,size);
    }
    if(x_pos < 0)
    {
        x_pos = size - fmod(x_pos*-1.0,size);
    }
    if(y_pos < 0)
    {
        y_pos = size- fmod(y_pos*-1.0,size);
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