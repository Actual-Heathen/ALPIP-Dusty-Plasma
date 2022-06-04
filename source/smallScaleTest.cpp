#include <time.h>
#include <stdlib.h>
#include <iostream>
#include <math.h>
#include "../headers/Particle.h"
#include "../headers/Point.h"
#include <vector>

using namespace std;
#define paticleCount 1
//#define pointCount 2
#define spacing .00001
#define gridDiv 3

int main()
{
    srand(time(NULL));

    vector<Particle> dust;
    //Point points[pointCount];

    double rho[gridDiv][gridDiv];


    for (int i = 0; i < paticleCount; i++) //set particle position
    {
        Particle temp( (double(rand())/double(RAND_MAX)*spacing*(gridDiv-1)), (double(rand())/double(RAND_MAX))*spacing*(gridDiv-1),1.0);

        dust.push_back(temp);
        //cout << dust[i].getX() << ", " << dust[i].getY() << "\n";
    }

    for (int l = 0; l< 600; l++)
    {
        for (int i = 0; i < gridDiv; i++)   //initialize grid
        {
            for (int j = 0; j < gridDiv; j++)
            {
                rho[i][j] = 0;
            }
        }

        for (int i = 0; i < paticleCount; i++) //calculate rho
        {
            int iXm = floor(dust[i].getX()/spacing);
            int iXp = iXm + 1;
            double wXm = 1- abs((dust[i].getX()-iXm*spacing)/spacing);
            double wXp = 1- abs((dust[i].getX()-iXp*spacing)/spacing);

            int iYm = floor(dust[i].getY()/spacing);
            int iYp = iYm + 1;
            double wYm = 1- abs((dust[i].getY()-iYm*spacing)/spacing);
            double wYp = 1- abs((dust[i].getY()-iYp*spacing)/spacing);
            
            //cout << iXp<< ", "<< wXp << "<-ixp "<< iYp << ", "<< wYp << "<-iyp\n";
            //cout << iXm<< " ,"<< wXm << "<-ixm "<< iYm << ", "<< wYm << "<-iym\n";

            rho[iXm][iYm] += wXm*wYm;
            rho[iXm][iYp] += wXm*wYp;
            rho[iXp][iYm] += wXp*wYm;
            rho[iXp][iYp] += wXp*wYp;
        }

        double rhoTemp = 0;

        for (int i = 0; i < gridDiv; i++)
        {
            for (int j = 0; j < gridDiv; j++)
            {
                //cout << rho[i][j] << "\n";
                rhoTemp += rho[i][j];
            }
        }

        //cout << rhoTemp << " rho sum\n";

        for (int i = 0; i < paticleCount; i++) //calculate gravity
        {
            dust[i].addAcceleration(rho,spacing);
            dust[i].move(1);
            if (i == 0)
                cout<< dust[i].getX() << ", " << dust[i].getY() << "\n ";
            //if (i == 1)
            //    cout<< dust[i].getX() << ", " << dust[i].getY() << "\n ";
            // if (i == 2)
            //     cout<< dust[i].getX() << ", " << dust[i].getY() << ", ";
            // if (i == 3)
            //     cout<< dust[i].getX() << ", " << dust[i].getY() << ", ";
            // if (i == 4)
            //     cout<< dust[i].getX() << ", " << dust[i].getY() << "\n";
        }
    }
}