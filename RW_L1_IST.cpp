/** Written by Sophie M. Fosson **/
/** 7 Oct 2018 **/
/** Code for simulations of paper "A biconvex analysis for Lasso l1 reweighting", IEEE Signal Process. Lett. 2018
/** LASSO IRL1 IST (Algorithm 2)**/

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
#include <ctime>
using namespace std;
using namespace Eigen;

#define pi 3.1415926535897932384626433832795028841971
#define MAX_ITERATIONS 5e5
#define ERROR_TOL 1e-10

double sigma, tau, BETA, SNR, LAMBDA;
int N, M, h, i, j, k, t, run, penalization;
bool thereisnoise;

//Standard Gaussian r.v. (mean=0; std=sigma)  via Box-Muller Algorithm
double Gaussian_Noise(double std) {
    return std*(sqrt(-2*log( ((double) rand() / (RAND_MAX))  )))*cos(2*pi* ((double) rand() / (RAND_MAX))  );
}

int main (int argc, char *argv[]) {
	if (argc < 6) {
		cerr << "Usage: " << argv[0] << " <signal size N> <sparsity k> <# measurements M> <run> < penalization type [1 = log(|x|+eps); 2 = a|x|-0.5|x|^2; 3 = (|x| + eps)^q ] >  <optional: SNR>" << endl;
                cerr << "Ex: ./gorwlist 256 30 100 1" << endl;
		return EXIT_FAILURE;
	}
	N = atoi (argv[1]);
	k = atoi (argv[2]);
	M = atoi (argv[3]);
	run = atoi(argv[4]);
	penalization = atoi(argv[5]);
  	
	if (argc == 7) {
            SNR = atof(argv[6]);
            thereisnoise=1;
        }
        else
            thereisnoise=0;
        
	tau=2.5e-1;

	srand(10*run*N);
	VectorXd y, x,  xOriginal, difference, x_old, z, beta, final_error;
	MatrixXd A, A_T, A_T_A, ID; //X= dynamic d=double
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
	LAMBDA=1e-5;
    	/** noise **/
	if (thereisnoise) {
		LAMBDA=1e-4;
        	sigma= sqrt(xOriginal.squaredNorm()/(pow(10.0, 0.1*SNR)*M));
        	for (i = 0; i < M; i++) {
            		y(i) =y(i)+ Gaussian_Noise(sigma);
        	}
	}
    	/** 3. RECOVERY **/
	time_t tstart, tend; 
 	tstart = time(0);

	for (t = 1; t < MAX_ITERATIONS; t++)    {
		/*for (i = 0; i < N; i++)  {
			if ( abs(x(i))>1e2) {
				cerr << x(i) << "  not stable"  << " " <<i << " " << t <<endl;
				return 0;
			}
		}*/
		z=x+tau*(A_T*y- A_T_A*x);

		if (penalization == 1) {
                        for (i = 0; i < N; i++) {
                                beta(i)=1./(abs(x(i)) +0.1);

				if ( abs( z(i) ) <   LAMBDA*beta(i))
					x(i) = 0;
 				else if ( z(i) > 0)
					x(i) = z(i)- LAMBDA*beta(i);
				else if ( z(i) < 0)
					x(i) = z(i) + LAMBDA*beta(i);

                        }
                } else if (penalization == 2) {
                        for (i = 0; i < N; i++) {
                                beta(i)=( 2-abs(x(i))>0 ? 2-abs(x(i)) : 0 );  
				if ( abs( z(i) ) <   LAMBDA*beta(i))
					x(i) = 0;
 				else if ( z(i) > 0)
					x(i) = z(i)- LAMBDA*beta(i);
				else if ( z(i) < 0)
					x(i) = z(i) + LAMBDA*beta(i);

                        }
                } else if (penalization == 3) {
                        for (i = 0; i < N; i++) {
                                beta(i)=0.5*pow( abs(x(i))+0.1,0.5-1);
				if ( abs( z(i) ) <   LAMBDA*beta(i))
					x(i) = 0;
 				else if ( z(i) > 0)
					x(i) = z(i)- LAMBDA*beta(i);
				else if ( z(i) < 0)
					x(i) = z(i) + LAMBDA*beta(i);

                        }
                }
		difference = x- x_old;

        	/** STOP CRITERION **/
        	if (difference.squaredNorm() < ERROR_TOL) {
            		break;
        	}
        	x_old=x;
	}  /** end t **/	
	tend = time(0);
	final_error = x- xOriginal;
        cout << final_error.squaredNorm()/xOriginal.squaredNorm() << "  " << (final_error.lpNorm<Infinity>()<1e-3) << "  " << t << "  "<<difftime(tend, tstart) << endl;
	return EXIT_SUCCESS;
}


