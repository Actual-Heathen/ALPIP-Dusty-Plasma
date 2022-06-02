#include <time.h>
#include <stdlib.h>
#include <iostream>
#include <math.h>
#include "basicFunctions.cpp"
#include "Particle.h"
#include "Point.h"

using namespace std;


struct vector2D
{
    float x = 0;
    float y = 0;
    float magnitude = 0;
};

int main()
{
    srand(time(NULL));

    Particle dust[3];
    vector2D points[4];

    for (int i = 0; i < 3; i++)
    {
        Particle temp(float(rand())/float(RAND_MAX), float(rand())/float(RAND_MAX),1.0);

        dust[i] = temp;
    }

    for (int i = 0; i < 3; i++)
    {
        float tempDist = distance(dust[i].getX(),dust[i].getY(), 0,0);
        float tempG = gravitationalForce(tempDist);
        float tempV = sqrt(pow(dust[i].getX(),2) + pow(dust[i].getY(),2));
        points[0].x += (dust[i].getX() * tempG)/tempV;
        points[0].y += (dust[i].getY() * tempG)/tempV;
        points[0].magnitude += sqrt(pow(points[0].x,2) + pow(points[0].y,2));

        cout << tempDist << "\t";
        cout << tempG << "\t";
        cout << tempV << "\n";

        cout << points[0].x << "\t";
        cout << points[0].y << "\t";
        cout << points[0].magnitude << "\n\n";

        tempDist = distance(dust[i].getX(),dust[i].getY(), 0,1);
        tempG = gravitationalForce(tempDist);
        tempV = sqrt(pow(dust[i].getX(),2) + pow(dust[i].getY(),2));
        points[1].x += (dust[i].getX() * tempG)/tempV;
        points[1].y += (((dust[i].getY()-1) * tempG)/tempV);
        points[1].magnitude += sqrt(pow(points[1].x,2) + pow(points[1].y,2));

        cout << tempDist << "\t";
        cout << tempG << "\t";
        cout << tempV << "\n";

        cout << points[1].x << "\t";
        cout << points[1].y << "\t";
        cout << points[1].magnitude << "\n\n";

        tempDist = distance(dust[i].getX(),dust[i].getY(), 1,0);
        tempG = gravitationalForce(tempDist);
        tempV = sqrt(pow(dust[i].getX(),2) + pow(dust[i].getY(),2));
        points[2].x += (((dust[i].getX()-1) * tempG)/tempV);
        points[2].y += (dust[i].getY() * tempG)/tempV;
        points[2].magnitude += sqrt(pow(points[2].x,2) + pow(points[2].y,2));

        cout << tempDist << "\t";
        cout << tempG << "\t";
        cout << tempV << "\n";

        cout << points[2].x << "\t";
        cout << points[2].y << "\t";
        cout << points[2].magnitude << "\n\n";

        tempDist = distance(dust[i].getX(),dust[i].getY(), 1,1);
        tempG = gravitationalForce(tempDist);
        tempV = sqrt(pow(dust[i].getX(),2) + pow(dust[i].getY(),2));
        points[3].x += (((dust[i].getX()-1) * tempG)/tempV);
        points[3].y += (((dust[i].getY()-1) * tempG)/tempV);
        points[3].magnitude += sqrt(pow(points[3].x,2) + pow(points[3].y,2));

        cout << tempDist << "\t";
        cout << tempG << "\t";
        cout << tempV << "\n";

        cout << points[3].x << "\t";
        cout << points[3].y << "\t";
        cout << points[3].magnitude << "\n\n";
    }
    
    float random = float(rand())/float(RAND_MAX);
    float random2 = float(rand())/float(RAND_MAX);
    cout << random << "\n";
    cout << random2 << "\n";

    cout << gravitationalForce(random) << "\n";
    cout << gravitationalForce(random2) << "\n";

    cout << distance(random,random,random2,random2) << "\n";

}