#include "../headers/Particle.h"
#include "../headers/Point.h"
#include <time.h>
#include <stdlib.h>
#include <iostream>
#include <math.h>
#include <vector>
#include <fstream>
#include <fftw3.h>
#include <random>
#include <unistd.h>

using namespace std;
#define ppc 100
#define gridSize (100)
#define gridDiv 100
#define loopCount 1000
#define threads 16
int main()
{
    int particleCount = (ppc*pow((gridDiv),2));
    //int particleCount = 25;
    double spacing = (double)gridSize/(double)(gridDiv);
    double E[3] = {0,0,0};
    double B[3] = {0,0,1};
    double timeStep = .05;
    srand(time(NULL));                                          //seed random number generator
    double energy = 0;
    double meanRho = 0;
    double lambda = 1;

    ofstream density;                                              //open data files
    //ofstream coor;
    ofstream kEnergyF;
    ofstream gEnergyF;
    ofstream energyF;
    ofstream fou;
    ofstream adj;
    ofstream fp;
    ofstream rhoS;
    density.open("data/density.d");
    kEnergyF.open("data/kEnergy.d");
    gEnergyF.open("data/gEnergy.d");
    energyF.open("data/energy.d");
    fou.open("data/fTransform.d");
    adj.open("data/aFTransform.d");
    fp.open("data/psi.d");
    rhoS.open("data/rhoS.d");

    vector<Particle> dust;                                      //declare dust and point grid
    vector<vector<double>> rho(gridDiv, vector<double> (gridDiv));
    vector<vector<double>> psi(gridDiv, vector<double> (gridDiv));
    vector<vector<double>> phi(gridDiv, vector<double> (gridDiv));
    vector<vector<double>> dpsix(gridDiv, vector<double> (gridDiv));
    vector<vector<double>> dpsiy(gridDiv, vector<double> (gridDiv));
    vector<vector<double>> dphix(gridDiv, vector<double> (gridDiv));
    vector<vector<double>> dphiy(gridDiv, vector<double> (gridDiv));
    cout << particleCount << "\n";
    default_random_engine generator(time(NULL));
    normal_distribution<double> normal(0,1);
    double vth = 1;
    

        for (int i = 0; i < particleCount; i++)                      //set particle position
        {
            double vx=vth*normal(generator);
            double vy=vth*normal(generator);
            Particle temp( ((double)rand()/(double)RAND_MAX)*spacing*(gridDiv), ((double)rand()/(double)RAND_MAX)*spacing*(gridDiv),1,vx,vy);
            //Particle temp( 2.5, 0,1,-0.1,0.0);
            dust.push_back(temp);
            //cout << dust[i].getX() << " " << dust[i].getY() << "\n";
        }
    
    double *Kx;
    double *Ky;
    Kx = (double*) malloc(sizeof(double)*gridDiv);
    Ky = (double*) malloc(sizeof(double)*gridDiv);

    int test = ceil(gridDiv/2.0);
    #pragma omp parallel for num_threads(threads) schedule(static)
        for (int i = 0; i < test ; i++)
        {
            Kx[i] = (2*i*M_PI)/(gridSize);
            Kx[gridDiv-1-i] = -(i+1)*(2*M_PI)/((gridSize));
            Ky[i] = (2*i*M_PI)/(gridSize);
            Ky[gridDiv-1-i] = -(i+1)*(2*M_PI)/((gridSize));
        }
    for (int l = 0; l< loopCount; l++)                          //main loop
    {
        if (l%10 == 0)
            cout << l << "\n";
        
        double gEnergy = 0;
        double kEnergy = 0;
        energy = 0;
        fftw_plan p;
        fftw_plan r;
        fftw_plan phiR;
        fftw_complex *in, *out;
        fftw_complex *phiIn, *phiOut;
        in = fftw_alloc_complex(gridDiv*gridDiv);
        out = fftw_alloc_complex(gridDiv*gridDiv);
        phiIn = fftw_alloc_complex(gridDiv*gridDiv);
        phiOut = fftw_alloc_complex(gridDiv*gridDiv);
        p = fftw_plan_dft_2d(gridDiv,gridDiv,in,out,FFTW_FORWARD,FFTW_ESTIMATE);
        r = fftw_plan_dft_2d(gridDiv,gridDiv,out,in,FFTW_BACKWARD,FFTW_ESTIMATE);
        phiR = fftw_plan_dft_2d(gridDiv, gridDiv, phiIn, phiOut, FFTW_BACKWARD, FFTW_ESTIMATE);

        //cout << "dec\n";

        #pragma omp parallel for num_threads(threads) schedule(static)
            for (int i = 0; i < gridDiv; i++)                       //initialize/reset  point grid
            {
                for (int j = 0; j < gridDiv; j++)
                {
                    rho[i][j] = 0.0;
                }

            }
        //cout << "rho reset\n";

            for (int i = 0; i < particleCount; i++)                  //calculate rho
            {
                double sp = spacing;

                int iXm = floor(dust[i].getX()/sp);                 //calculate iXm & iXp
                if (iXm >= gridDiv)
                {
                    iXm -= gridDiv;
                }
                if (iXm < 0)
                {
                    iXm += gridDiv;
                }
                int iXp = iXm + 1;
                //cout << iXm << "-";
                double wXm = 1- abs((dust[i].getX()/sp)-iXm);       //weight calculations
                double wXp = 1- abs((dust[i].getX()/sp)-iXp);

                if (iXp >= gridDiv)                                 //adjust iXp if in "ghost region"
                {
                    iXp = iXp - gridDiv;
                }
                //cout << iXp <<", ";

                int iYm = floor(dust[i].getY()/sp);                 //calculate iYm and iYp
                if (iYm >= gridDiv)
                {
                    iYm -= gridDiv;
                }
                if (iYm < 0)
                {
                    iYm += gridDiv;
                }
                int iYp = iYm + 1;
                //cout << iYm <<"\n";
                double wYm = 1- abs((dust[i].getY()/sp)-iYm);       //wieght calculation
                double wYp = 1- abs((dust[i].getY()/sp)-iYp);

                if (iYp >= gridDiv)                                 //ghost region adjust
                {
                    iYp = iYp - gridDiv;
                }
                //cout << iYp << "\n";

                //cout << iXm<<"\n";        

                rho[iXm][iYm] += ((wXm*wYm))/ppc;                           //add weights to points
                rho[iXm][iYp] += ((wXm*wYp))/ppc;
                rho[iXp][iYm] += ((wXp*wYm))/ppc;
                rho[iXp][iYp] += ((wXp*wYp))/ppc;
                //cout << i << "\n";

                kEnergy += .5*dust[i].getMass()*pow(dust[i].getSpeed(),2)/ppc;
            }


        //cout << "rho set\n";

        for (int i = 0; i < gridDiv; i++)
        {
            for (int j = 0; j < gridDiv; j++)
            {
                //cout << rho[i][j] << "\n";
                in[j+i*gridDiv][0] = rho[i][j];  //add rho to fftw input
                density<<rho[j][i]<<"\n";
                in[j+i*gridDiv][1] = 0;
                meanRho += pow(rho[i][j], 2);
            }
        }

        meanRho /= pow(gridDiv, 2);
        rhoS << timeStep*l <<" "<<meanRho << "\n";

        fftw_execute(p);        //Execute the FFT
        //calculate Kx, Ky//

        out[0][0] = 0;
        out[0][1] = 0;
        for (int i = 0; i < gridDiv; i++)
        {
            for (int j = 0; j < gridDiv; j++)
            {
                //cout << rho[i][j] << "\n";
                double temp = out[i*gridDiv+j][0];
                //fou << temp << "\n";
                //cout << temp << "\n";
		        double tempC = out[i*gridDiv+j][1];

                double tempPhi = temp;
                double tempPhiC = tempC;

                if (Ky[j] != 0 || Kx[i] != 0)
                {
                    temp = temp/(pow(Ky[j],2)+pow(Kx[i],2));
                    tempC = tempC/(pow(Ky[j],2)+pow(Kx[i],2));

                    tempPhi = tempPhi/(pow(Ky[j],2)+pow(Kx[i],2)+(pow(lambda,-2)));
                    tempPhiC = tempPhiC/(pow(Ky[j],2)+pow(Kx[i],2)+(pow(lambda,-2)));
                }


		        out[i*gridDiv +j][0] = temp;
                phiIn[i*gridDiv +j][0] = tempPhi;
                adj << temp<<"\n";
		        out[i*gridDiv+j][1] = tempC;
                phiIn[i*gridDiv+j][1] = tempPhiC;
		        
            }
        }


        fftw_execute(r);
        fftw_execute(phiR);

        #pragma omp parallel for num_threads(threads) schedule(static)
            for (int i = 0; i < gridDiv;i++)
            {
                for(int j = 0; j < gridDiv; j++)
                {
                    psi[i][j] = in[i*gridDiv+j][0]/pow(gridDiv,2);
                    phi[i][j] = phiOut[i*gridDiv+j][0]/pow(gridDiv,2);
                }
            }

        for (int i = 0; i < gridDiv; i++)
        {
            for (int j = 0; j < gridDiv; j++)
            {
                int xm;
                int xp;
                int ym;
                int yp;
                xm = i - 1;
                xp = i + 1;
                ym = j - 1;
                yp = j + 1;

                if (i == 0)
                    xm = gridDiv-1;
                if(j == 0)
                    ym = gridDiv-1;
                if (i >= gridDiv-1)
                    xp = 0;
                if (j >= gridDiv-1)
                    yp = 0;

                dpsix[i][j] = ((psi[xp][j]-psi[xm][j]))/(2*spacing);
                dpsiy[i][j] = ((psi[i][yp]-psi[i][ym]))/(2*spacing);

                dphix[i][j] = ((phi[xp][j]-phi[xm][j]))/(2*spacing);
                dphiy[i][j] = ((phi[i][yp]-phi[i][ym]))/(2*spacing);
                fp<< psi[j][i]<<"\n";
                fou<<phi[j][i]<<"\n";
                gEnergy += (-0.5*rho[i][j]*psi[i][j]);
            }
	    }

        
        //cout << "rho counted\n";

        //cout << rhoTemp << " rho sum\n";                      //debugging//

        // for (int i = 0; i < particleCount; i++)
        // {
        //    //points << dust[i].getX() << " "<<dust[i].getY()<<"\n"; //write particle 0's coordinatess to csv
        //    //coor << dust[i].getX() << " "<<dust[i].getY()<<"\n";
        //    //pointTime << timeStep*l << " "<<  dust[0].getY() <<"\n";
        // }

        #pragma omp parallel for num_threads(threads) schedule(static)//define parallel section
            for (int i = 0; i < particleCount; i++)              //calculate gravity
            {
                dust[i].addAcceleration(spacing, dpsix, dpsiy, dphix, dphiy, E,B,timeStep);          //add Acceleration based on densities
                //cout<<"accel\n";
                dust[i].move(timeStep, gridSize);              //move particle
                //cout<<"move\n";
		        //cout <<"moved\n";
                // if ((work_done % 1000) == 0) //debugging couter
                //     cout <<"number "<< work_done << "counted\n";
            }

        //points.close();
        //density.close();
        fftw_destroy_plan(p);
        fftw_destroy_plan(r);
        fftw_free(in);
        fftw_free(out);
        energy = kEnergy + gEnergy;
        if (energy > 100000)
            cout << energy<< "e\n";
        //cout << gEnergy<< ", "<< kEnergy << "\n";
        energyF << timeStep*l<<" "<< energy <<"\n";
        gEnergyF << timeStep*l<<" "<< gEnergy<<"\n";
        kEnergyF << timeStep*l<<" "<< kEnergy<<"\n";
        energy = 0;
        meanRho = 0;

    }

    dust.clear();
    rho.clear();
    density.close();
    kEnergyF.close();
    gEnergyF.close();
    energyF.close();
    fou.close();
    adj.close();
    fp.close();
    rhoS.close();
}
