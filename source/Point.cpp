#include "../headers/Point.h"
#include "../headers/Particle.h"
#include <math.h>

#define G 6.6743*pow(10,-11)

Point::Point()
{
    x_pos = 0;
    y_pos = 0;
}

Point::Point(long double x, long double y)
{
    x_pos = x;
    y_pos = y;
}

long double Point::getX()
{
    return x_pos;
}

long double Point::getY()
{
    return y_pos;
}

long double Point::getXGV()
{
    return GXVector;
}

long double Point::getYGV()
{
    return GYVector;
}

long double Point::getGVM()
{
    return GVectorMagnitude;
}

void Point::addGravity(Particle p)
{
    long double dX = p.getX() - x_pos;
    long double dY = p.getY() - y_pos;
    long double distance = sqrt(pow(dX,2) + pow(dY,2));
    long double accel = G*p.getMass()/pow(distance,2);

    GXVector += (dX*accel)/distance;
    GYVector += (dY*accel)/distance;
    GVectorMagnitude = sqrt(pow(GXVector,2) + pow(GYVector,2));

}

void Point::resetGVector()
{
    GXVector = 0;
    GYVector = 0;
    GVectorMagnitude = 0;

}