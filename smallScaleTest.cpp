#include <time.h>
#include <stdlib.h>
#include <iostream>
#include <math.h>
#include "basicFunctions.cpp"

using namespace std;

struct dust
{
    float x_pos;
    float y_pos;
    float mass;
};


int main()
{
    srand(time(NULL));
    
    float random = float(rand())/float(RAND_MAX);
    float random2 = float(rand())/float(RAND_MAX);
    cout << random << "\n";
    cout << random2 << "\n";

    cout << gravitationalForce(random) << "\n";
    cout << gravitationalForce(random2) << "\n";

    cout << distance(random,random,random2,random2) << "\n";

}