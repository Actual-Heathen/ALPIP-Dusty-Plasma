#include <time.h>
#include <stdlib.h>
#include <iostream>
#include <math.h>
#include "../headers/Particle.h"
#include "../headers/Point.h"
#include <vector>
#include <fstream>

using namespace std;
#define paticleCount 512*512*10
//#define pointCount 2
#define spacing 1.0
#define gridDiv 512

int main()
{
    srand(time(NULL));
    ofstream data;
    data.open("data.d");
    vector<Particle> dust;
    //Point points[pointCount];

    vector<vector<double>> rho(gridDiv, vector<double> (gridDiv));

    cout << "declared\n";
    for (int i = 0; i < paticleCount; i++) //set particle position
    {
        Particle temp( ((double)rand()/(double)RAND_MAX)*spacing*(gridDiv), ((double)rand()/(double)RAND_MAX)*spacing*(gridDiv),1.0);
        //cout << temp.getX()<<temp.getY()<<"\n";
        dust.push_back(temp);
        //cout << dust[i].getX() << ", " << dust[i].getY() << "\n";
    }

    for (int l = 0; l< 1; l++)
    {
        for (int i = 0; i < gridDiv; i++)   //initialize grid
        {
            for (int j = 0; j < gridDiv; j++)
            {
                rho[i][j] = 0;
            }

        }
        cout << "rho reset\n";

        for (int i = 0; i < paticleCount; i++) //calculate rho
        {
            int iXm = floor(dust[i].getX()/spacing);
            int iXp = iXm + 1;
            double wXm = 1- abs((dust[i].getX()/spacing)-iXm);
            double wXp = 1- abs((dust[i].getX()/spacing)-iXp);
            if (iXp >= gridDiv)
            {
                iXp = iXp - gridDiv;
            }
            int iYm = floor(dust[i].getY()/spacing);
            int iYp = iYm + 1;
         
            double wYm = 1- abs((dust[i].getY()/spacing)-iYm);
            double wYp = 1- abs((dust[i].getY()/spacing)-iYp);
             if (iYp >= gridDiv)
            {
                iYp = iYp - gridDiv;
            }            
            //cout << iXp<< ", "<< wXp << "<-ixp "<< iYp << ", "<< wYp << "<-iyp\n";
            //cout << iXm<< " ,"<< wXm << "<-ixm "<< iYm << ", "<< wYm << "<-iym\n";

            rho[iXm][iYm] += wXm*wYm;
            rho[iXm][iYp] += wXm*wYp;
            rho[iXp][iYm] += wXp*wYm;
            rho[iXp][iYp] += wXp*wYp;
        }

        cout << "rho set\n";
        double rhoTemp = 0;

        for (int i = 0; i < gridDiv; i++)
        {
            for (int j = 0; j < gridDiv; j++)
            {
                //cout << rho[i][j] << "\n";
                rhoTemp += rho[i][j];
                data << rho[i][j]<<"\n";
            }
        }
        cout << "rho counted";

        cout << rhoTemp << " rho sum\n";

        int work_done = 0;

        #pragma omp parallel for num_threads(6) schedule(static)
        for (int i = 0; i < paticleCount; i++) //calculate gravity
        {
            dust[i].addAcceleration(spacing, rho);
            dust[i].move(.01);
            if (i == 0)
            {
                //data<< dust[i].getX() << ", " << dust[i].getY() << "\n ";
            }
            //if (i == 1)
            //    cout<< dust[i].getX() << ", " << dust[i].getY() << "\n ";
            // if (i == 2)
            //     cout<< dust[i].getX() << ", " << dust[i].getY() << ", ";
            // if (i == 3)
            //     cout<< dust[i].getX() << ", " << dust[i].getY() << ", ";
            // if (i == 4)
            //     cout<< dust[i].getX() << ", " << dust[i].getY() << "\n";
            #pragma omp atomic
            work_done++;
            
            if ((work_done % 1000) == 0)
                cout <<"number "<< work_done << "counted\n";
        }
        cout << "grav calc\n";
    }
    dust.clear();
    rho.clear();
    data.close();
}