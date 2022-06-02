#include <time.h>
#include <stdlib.h>
#include <iostream>
#include <math.h>
#include "basicFunctions.cpp"
#include "../headers/Particle.h"
#include "../headers/Point.h"

using namespace std;


int main()
{
    srand(time(NULL));

    Particle dust[3];
    Point points[4];

    points[0] = Point(0.0f,0.0f);
    points[1] = Point(0.0f,1.0f);
    points[2] = Point(1.0f,0.0f);
    points[3] = Point(1.0f,1.0f);

    for (int i = 0; i < 3; i++)
    {
        Particle temp(float(rand())/float(RAND_MAX), float(rand())/float(RAND_MAX),1.0);

        dust[i] = temp;
    }

    for (int i = 0; i < 3; i++)
    {
        for (int j = 0; j < 4; j++)
        {
            points[j].addGravity(dust[i]);
            cout << "point " << j <<", paticle " << i << ":\n";
            cout <<  "(x,y): (" << dust[i].getX() << ", " << dust[i].getY() << ")\n";
            cout << "gravity: " << points[j].getGVM() << "\n\n";
        }
    }

}