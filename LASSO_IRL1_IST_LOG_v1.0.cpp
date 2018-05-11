// Written by Sophie Fosson 2018
// Algorithm 4 - Lasso IST IRL1

#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <math.h>
#include <iostream>
#include <fstream>
#include <sstream>
#include <algorithm>
#include <iomanip>
#include <Eigen/Dense>
using namespace std;
using namespace Eigen;

#define pi 3.1415926535897932384626433832795028841971
#define MAX_ITERATIONS 5e5
#define LAMBDA 1e-5
#define ERROR_TOL 1e-10

double sigma, tau, BETA, SNR;
int N, M, h, i, j, k, t, t_out, t_tot, t_save, run;
bool thereisnoise;

//Standard Gaussian r.v. (mean=0; std=sigma)  via Box-Muller Algorithm
double Gaussian_Noise(double std) {
    return std*(sqrt(-2*log( ((double) rand() / (RAND_MAX))  )))*cos(2*pi* ((double) rand() / (RAND_MAX))  );
}

int main (int argc, char *argv[]) {
	if (argc < 5) {
		cerr << "Usage: " << argv[0] << " <signal size N> <sparsity k> <# measurements M> <run>  <optional: SNR>" << endl;
                cerr << "Ex: ./run_irl1_ist 256 30 100 1" << endl;
		return EXIT_FAILURE;
	}
	N = atoi (argv[1]);
	k = atoi (argv[2]);
	M = atoi (argv[3]);
	run = atoi(argv[4]);
  	 if (argc == 6) {
            SNR = atof(argv[5]);
            thereisnoise=1;
        }
        else
            thereisnoise=0;
        
	tau=2.5e-1;

    
	srand(10*run*N);
	VectorXd y, x,  xOriginal, difference, x_old, z, beta, final_error;
	MatrixXd A, A_T, A_T_A, ID, AS; //X= dynamic d=double
    	A.resize(M, N);
 	xOriginal.resize(N);
	x.resize(N);
    	final_error.resize(N);
        difference.resize(N);
   	beta.resize(N);
	y.resize(M);
 	x_old.resize(N);
 	ID.setIdentity(N,N);
	A_T_A.resize(N,N);
 	z.resize (N);

    	for (i = 0; i < M; i++) {
		for (h=0; h< N; h++) {
			A(i,h)=1./sqrt(M)*Gaussian_Noise(1);
		}
	}
	A_T=A.transpose();
    	A_T_A=A_T*A;

	/** 1. GENERATE SIGNAL **/
	VectorXd index_nonzero(N);
	for (j = 0; j < N; j++) {
		index_nonzero(j)=j;
	}

	PermutationMatrix<Dynamic,Dynamic> perm(N);
	perm.setIdentity();
	std::random_shuffle(perm.indices().data(), perm.indices().data()+perm.indices().size());
    	index_nonzero = perm * index_nonzero; // permute rows
	
	for (j = 0; j < k; j++) {
            xOriginal(index_nonzero(j))= Gaussian_Noise(1);
        }
    

	/** 2. GENERATE MEASUREMENTS **/
	y = A*xOriginal;
    	/** noise **/
	if (thereisnoise) {
        	sigma= sqrt(xOriginal.squaredNorm()/(pow(10.0, 0.1*SNR)*M));
        	for (i = 0; i < M; i++) {
            		y(i) =y(i)+ Gaussian_Noise(sigma);
        	}
	}
    	/** 3. RECOVERY **/	
	for (t = 1; t < MAX_ITERATIONS; t++)    {
		for (i = 0; i < N; i++)  {
			if ( abs(x(i))>1e2) {
				cerr << x(i) << "  BUM"  << " " <<i << " " << t <<endl;
				return 0;
			}
		}
		z=x+tau*(A_T*y- A_T_A*x);
		for (i = 0; i < N; i++) {
			beta(i)=1./(abs(x(i))+0.1);
			if ( abs( z(i) ) <   LAMBDA*beta(i))
				x(i) = 0;
 			else if ( z(i) > 0)
				x(i) = z(i)- LAMBDA*beta(i);
			else if ( z(i) < 0)
				x(i) = z(i) + LAMBDA*beta(i);
        		//beta(i)=BETA- abs(x(i));
			
 		}
		difference = x- x_old;

        	//STOP CRITERION
        	if (difference.squaredNorm() < ERROR_TOL) {
            		break;
        	}
        	x_old=x;
	
	}  /** end t **/	

	final_error = x- xOriginal;
	cout << final_error.squaredNorm()/xOriginal.squaredNorm() << "  " << (final_error.lpNorm<Infinity>()<1e-3) << "  " << t << endl;
	return EXIT_SUCCESS;
}


