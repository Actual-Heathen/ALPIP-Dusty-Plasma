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
#define ppc 100
#define gridSize (1)
#define gridDiv 50
#define loopCount 2000
#define G (6.6743*pow(10,-11))
#define threads 12
int main()
{
    int particleCount = (ppc*pow((gridDiv),2));
    //int particleCount = 25;
    double spacing = (double)gridSize/(double)(gridDiv);
    double E[3] = {0,0,0};
    double B[3] = {0,0,0};
    double timeStep = .01;
    srand(time(NULL));                                          //seed random number generator
    double energy = 0;
    double meanRho = 0;

    ofstream data;                                              //open data files
    ofstream coor;
    ofstream pointTime;
    ofstream energyPlt;
    ofstream energy2;
    data.open("data/density.d");
    coor.open("data/pointsF.d");
    pointTime.open("data/pointTime.d");
    energyPlt.open("data/energy.d");
    energy2.open("data/energy2.d");
    ofstream fou;                                              //open data files
    ofstream adj;
    fou.open("data/fTransform.d");
    adj.open("data/aFTransform.d");
    ofstream fp;                                              //open data files
    fp.open("data/psi.d");

    vector<Particle> dust;                                      //declare dust and point grid
    vector<vector<double>> rho(gridDiv, vector<double> (gridDiv));
    vector<vector<double>> psi(gridDiv, vector<double> (gridDiv));
    vector<vector<double>> dpsix(gridDiv, vector<double> (gridDiv));
    vector<vector<double>> dpsiy(gridDiv, vector<double> (gridDiv));
    //cout << "declared\n";
    cout << particleCount << "\n";

        for (int i = 0; i < particleCount; i++)                      //set particle position
        {
            Particle temp( ((double)rand()/(double)RAND_MAX)*spacing*(gridDiv), ((double)rand()/(double)RAND_MAX)*spacing*(gridDiv),1,0.0,0.0);
            //Particle temp( 2.5, 0,1,-0.1,0.0);
            dust.push_back(temp);
            //cout << dust[i].getX() << " " << dust[i].getY() << "\n";
        }
    
    for (int l = 0; l< loopCount; l++)                          //main loop
    {
        cout << l << "\n";
        double gEnergy = 0;
        double kEnergy = 0;
        energy = 0;
        fftw_plan p;
        fftw_plan r;
        fftw_complex *in, *out;
        in = fftw_alloc_complex(gridDiv*gridDiv);
        out = fftw_alloc_complex(gridDiv*gridDiv);

        p = fftw_plan_dft_2d(gridDiv,gridDiv,in,out,FFTW_FORWARD,FFTW_ESTIMATE);
        r = fftw_plan_dft_2d(gridDiv,gridDiv,out,in,FFTW_BACKWARD,FFTW_ESTIMATE);

        // snprintf(fn, sizeof fn, "../data/points%05d.d",l);
        // ofstream points; points.open(fn);
        // snprintf(fn, sizeof fn, "../data/density%05d.d",l);
        // ofstream density; density.open(fn);

        //cout << "dec\n";

        #pragma omp parallel for num_threads(threads) schedule(static)
        {
            for (int i = 0; i < gridDiv; i++)                       //initialize/reset  point grid
            {
                for (int j = 0; j < gridDiv; j++)
                {
                    rho[i][j] = 0.0;
                }

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
                data<<rho[j][i]<<"\n";
                in[j+i*gridDiv][1] = 0;
                meanRho += pow(rho[i][j], 2);
                //density << rho[j][i]<<"\n";                        //add grid values to file
            }
        }

        meanRho /= pow(gridDiv, 2);

        fftw_execute(p);        //Execute the FFT
        //calculate Kx, Ky//
        double *Kx;
        double *Ky;
        Kx = (double*) malloc(sizeof(double)*gridDiv);
        Ky = (double*) malloc(sizeof(double)*gridDiv);

        #pragma omp parallel for num_threads(threads) schedule(static)
        {
            for (int i = 0; i < gridDiv/2.0; i++)
            {
                Kx[i] = (2*i*M_PI)/(gridSize);
                Kx[gridDiv-1-i] = -(i+1)*(2*M_PI)/((gridSize));
                Ky[i] = (2*i*M_PI)/(gridSize);
                Ky[gridDiv-1-i] = -(i+1)*(2*M_PI)/((gridSize));
            }
        }

        for (int i = 0; i < gridDiv; i++)
        {
            for (int j = 0; j < gridDiv; j++)
            {
                //cout << rho[i][j] << "\n";
                double temp = out[i*gridDiv+j][0];
                fou << temp << "\n";
                //cout << temp << "\n";
		        double tempC = out[i*gridDiv+j][1];
                if (Ky[j] != 0 || Kx[i] != 0)
                {
                    temp = temp/(pow(Ky[j],2)+pow(Kx[i],2));
                    tempC = tempC/(pow(Ky[j],2)+pow(Kx[i],2));
                }


		        out[i*gridDiv +j][0] = temp;
                adj << temp<<"\n";
		        out[i*gridDiv+j][1] = tempC;
		        
            }
        }
        out[0][0] = 0;
        out[0][1] = 0;
        delete[] Kx;
        delete[] Ky;


        fftw_execute(r);

        #pragma omp parallel for num_threads(threads) schedule(static)
        {
            for (int i = 0; i < gridDiv;i++)
            {
                for(int j = 0; j < gridDiv; j++)
                {
                    psi[i][j] = in[i*gridDiv+j][0]/pow(gridDiv,2);
                    
                }
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
                fp<< psi[j][i]<<"\n";
                gEnergy += (-0.5*rho[i][j]*psi[i][j]);
            }
	    }

        
        //cout << "rho counted\n";

        //cout << rhoTemp << " rho sum\n";                      //debugging//

        int work_done = 0;  //in paralell serial counter

        // for (int i = 0; i < particleCount; i++)
        // {
        //    //points << dust[i].getX() << " "<<dust[i].getY()<<"\n"; //write particle 0's coordinatess to csv
        //    //coor << dust[i].getX() << " "<<dust[i].getY()<<"\n";
        //    //pointTime << timeStep*l << " "<<  dust[0].getY() <<"\n";
        // }

        #pragma omp parallel for num_threads(threads) schedule(static)//define parallel section
        {
            for (int i = 0; i < particleCount; i++)              //calculate gravity
            {
                dust[i].addAcceleration(spacing, dpsix, dpsiy,E,B,timeStep);          //add Acceleration based on densities
                //cout<<"accel\n";
                dust[i].move(timeStep, gridSize);              //move particle
                //cout<<"move\n";
		        //cout <<"moved\n";
                #pragma omp atomic
                work_done++;
                
                // if ((work_done % 1000) == 0) //debugging couter
                //     cout <<"number "<< work_done << "counted\n";
            }
        }
        //cout << "grav calc\n";
        //points.close();
        //density.close();
        fftw_destroy_plan(p);
        fftw_destroy_plan(r);
        fftw_free(in);
        fftw_free(out);
        energy = kEnergy + gEnergy;
        //cout << gEnergy<< ", "<< kEnergy << "\n";
        energyPlt << timeStep*l<<" "<< energy <<"\n";
        pointTime << timeStep*l<<" "<< gEnergy<<"\n";
        energy2 << timeStep*l<<" "<< kEnergy<<"\n";
        energy = 0;
        meanRho = 0;

    }

    dust.clear();
    rho.clear();
    data.close();
    coor.close();
    fou.close();
    fp.close();
    adj.close();
    pointTime.close();
}