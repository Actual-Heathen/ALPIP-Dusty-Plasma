#include "../headers/Particle.h"
#include "../headers/Point.h"
#include <time.h>
#include <stdlib.h>
#include <iostream>
#include <math.h>
#include <vector>
#include <fstream>
#include <fftw3.h>

using namespace std;
#define paticleCount 5000
#define spacing (1*pow(10,-7))
#define gridDiv 100
#define loopCount 200

int main()
{
    srand(time(NULL));                                          //seed random number generator

    ofstream data;                                              //open data files
    ofstream coor;
    data.open("../data/data.d");
    //coor.open("data.csv");
    
    vector<Particle> dust;                                      //declare dust and point grid
    vector<vector<double>> rho(gridDiv, vector<double> (gridDiv));

    cout << "declared\n";
    
    for (int i = 0; i < paticleCount; i++)                      //set particle position
    {
        Particle temp( ((double)rand()/(double)RAND_MAX)*spacing*(gridDiv), ((double)rand()/(double)RAND_MAX)*spacing*(gridDiv),1.0,0.0,0.0);
        dust.push_back(temp);
        //cout << dust[i].getX() << " " << dust[i].getY() << "\n";
    }
    char fn[50];
    for (int l = 0; l< loopCount; l++)                          //main loop
    {
        fftw_plan p;
        fftw_complex *in, *out;
        in = fftw_alloc_complex(gridDiv*gridDiv);
        out = fftw_alloc_complex(gridDiv*gridDiv);

        snprintf(fn, sizeof fn, "../data/points%05d.d",l);
        ofstream points; points.open(fn);
        snprintf(fn, sizeof fn, "../data/density%05d.d",l);
        ofstream density; density.open(fn);



        for (int i = 0; i < gridDiv; i++)                       //initialize/reset  point grid
        {
            for (int j = 0; j < gridDiv; j++)
            {
                rho[i][j] = 0.0;
            }

        }
        //cout << "rho reset\n";

        for (int i = 0; i < paticleCount; i++)                  //calculate rho
        {
            double sp = spacing;

            int iXm = floor(dust[i].getX()/sp);                 //calculate iXm & iXp
            int iXp = iXm + 1;

            double wXm = 1- abs((dust[i].getX()/sp)-iXm);       //weight calculations
            double wXp = 1- abs((dust[i].getX()/sp)-iXp);
            
            if (iXp >= gridDiv)                                 //adjust iXp if in "ghost region"
            {
                iXp = iXp - gridDiv;
            }

            int iYm = floor(dust[i].getY()/sp);                 //calculate iYm and iYp
            int iYp = iYm + 1;
         
            double wYm = 1- abs((dust[i].getY()/sp)-iYm);       //wieght calculation
            double wYp = 1- abs((dust[i].getY()/sp)-iYp);
            
            if (iYp >= gridDiv)                                 //ghost region adjust
            {
                iYp = iYp - gridDiv;
            }            

            rho[iXm][iYm] += wXm*wYm;                           //add weights to points
            rho[iXm][iYp] += wXm*wYp;
            rho[iXp][iYm] += wXp*wYm;
            rho[iXp][iYp] += wXp*wYp;
        }



        //cout << "rho set\n";
        double rhoTemp = 0;                                     //debugging//

        for (int i = 0; i < gridDiv; i++)
        {
            for (int j = 0; j < gridDiv; j++)
            {
                //cout << rho[i][j] << "\n";
                rhoTemp += rho[i][j];                           //debugging//
                density << rho[j][i]<<"\n";                        //add grid values to file
            }
        }
        //cout << "rho counted";    

        //cout << rhoTemp << " rho sum\n";                      //debugging//

        int work_done = 0;  //in paralell serial counter

        for (int i = 0; i < paticleCount; i++)
        {
           points << dust[i].getX() << " "<<dust[i].getY()<<"\n"; //write particle 0's coordinatess to csv
        }

        #pragma omp parallel for num_threads(6) schedule(static)//define parallel section
        {
            for (int i = 0; i < paticleCount; i++)              //calculate gravity
            {
                dust[i].addAcceleration(spacing, rho);          //add Acceleration based on densities
                dust[i].move(1, gridDiv*spacing);              //move particle

                #pragma omp atomic
                work_done++;
                
                // if ((work_done % 1000) == 0) //debugging couter
                //     cout <<"number "<< work_done << "counted\n";
            }
        }
        //cout << "grav calc\n";
        points.close();
        density.close();
    }
    dust.clear();
    rho.clear();
    data.close();
    coor.close();
}