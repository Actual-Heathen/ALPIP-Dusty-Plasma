#include <stdlib.h>
#include <iostream>
#include <math.h>

#define G 6.6743

using namespace std;

float gravitationalForce(float distance)
{
    return (G*pow(10.0,-11.0))/pow(distance,2.0);
}

float distance(float x_1, float y_1, float x_2, float y_2)
{
    return sqrt(pow(x_2 - x_1, 2.0) + pow(y_2 - y_1, 2.0));
}