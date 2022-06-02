#include "../headers/Point.h"
#include "../headers/Particle.h"
#include <math.h>

#define G 6.6743*pow(10,-11)

Point::Point()
{
    x_pos = 0;
    y_pos = 0;
}

Point::Point(float x, float y)
{
    x_pos = x;
    y_pos = y;
}

float Point::getX()
{
    return x_pos;
}

float Point::getY()
{
    return y_pos;
}

float Point::getXGV()
{
    return GXVector;
}

float Point::getYGV()
{
    return GYVector;
}

float Point::getGVM()
{
    return GVectorMagnitude;
}

void Point::addGravity(Particle p)
{
    float dX = p.getX() - x_pos;
    float dY = p.getY() - y_pos;
    float distance = sqrt(pow(dX,2) + pow(dY,2));
    float accel = G*p.getMass()/pow(distance,2);

    GXVector += (dX*accel)/distance;
    GYVector += (dY*accel)/distance;
    GVectorMagnitude = sqrt(pow(GXVector,2) + pow(GYVector,2));

}