#include <fftw3.h>
#include <iostream>
#include <stdlib.h>
#include <math.h>
#include <fstream>

#define PI 3.14159265

int main()
{
    //fftw_complex *out, *in;
    fftw_plan p;
    fftw_plan p2;
    int N = 100;

    std::ofstream data;
    data.open("../data/data.d");

    
    fftw_complex *in, *out, *in2, *out2;

    in = (fftw_complex*) fftw_malloc(sizeof(fftw_complex)*N);
    out = (fftw_complex*) fftw_malloc(sizeof(fftw_complex)*N);
    in2 = (fftw_complex*) fftw_malloc(sizeof(fftw_complex)*N*N);
    out2 = (fftw_complex*) fftw_malloc(sizeof(fftw_complex)*N*N);
    

    //out = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * N);

    p = fftw_plan_dft_1d(N,in,out,FFTW_FORWARD, FFTW_ESTIMATE);
    p2 = fftw_plan_dft_2d(N,N,in2,out2,FFTW_FORWARD,FFTW_ESTIMATE);

    for (int i = 0; i < N; i++)
    {
        for (int j = 0;j <N;j++)
        {
            double x = ((2.0*PI*i)/N);
            double y = ((2.0*PI*j)/N);
            in2[j+(i*N)][0] = sin(25*x)+sin(25*y);
            std::cout << x <<"\n";
            in2[j+(i*N)][1] = 0;
            //data << in2[i][0]<<"\n";     //print sin(X)
        }
    }

    for (int i =0; i< N; i++)
    {
        in[i][0] = sin(100*(2*PI*i)/N);
        in[i][1] = 0;
    }
    fftw_execute(p);
    fftw_execute(p2);

    for (int i = 0; i < N; i++)
    {
        for (int j = 0;j <N;j++)
        {
            float temp = (float)out2[i+(j*N)][0];
            data << temp << "\n"; // print data
        }
    }
    for (int i  = 0; i < N;i++)
    {
        //data<<out[i][0]<<"\n";      //print fft 1d
    }

    //memeory deallocation//
    fftw_destroy_plan(p);
    free(out);
    free(in);
    fftw_destroy_plan(p2);
    fftw_free(out2);
    fftw_free(in2);

}