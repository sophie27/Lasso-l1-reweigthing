// Written by Sophie Fosson 2018
// Algorithm 3 - Lasso IRL1

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
#define MAX_ITERATIONS_OUT 1 // =1 for Classical Lasso; >1 for l1 reweigthing
#define LAMBDA 1e-5
#define ERROR_TOL 1e-10
#define D 1

double sigma, RSE, RSE_save, ROUNDE,  rho2, rho, FUNC_TRUE, FUNC, aux, BETA, SNR;
int N, M, h, i, j, k, t, t_out, t_tot, t_save, run;
bool thereisnoise;

//Standard Gaussian r.v. (mean=0; std=sigma)  via Box-Muller Algorithm
double Gaussian_Noise(double std) {
    return std*(sqrt(-2*log( ((double) rand() / (RAND_MAX))  )))*cos(2*pi* ((double) rand() / (RAND_MAX))  );
}

int main (int argc, char *argv[]) {
	if (argc < 5) {
		cerr << "Usage: " << argv[0] << " <signal size N> <sparsity k> <# measurements M> <run>  <optional: SNR>" << endl;
                cerr << "Ex: ./run_irl1 256 30 100 1" << endl;
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

        rho = 1;
	srand(10*run*N);
	
	VectorXd y, x, xOriginal,  x_old, difference, z, dual, beta, final_error;
	MatrixXd A, A_T, QI, ID, Q;
        A.resize(M, N);
 	xOriginal.resize(N);
	x.resize(N);
        final_error.resize(N);
	difference.resize(N);
        beta.resize(N);
        y.resize(M);
 	x_old.resize(N);
 	Q.resize(N,N);
 	ID.setIdentity(N,N);
	dual.resize (N);
 	z.resize (N);

        for (i = 0; i < M; i++) {
		for (h=0; h< N; h++) {
			A(i,h)=1./sqrt(M)*Gaussian_Noise(1);
		}
	}
	A_T=A.transpose();
    
	/** 1. GENERATE SIGNAL **/
	VectorXd index_nonzero(N);
	for (j = 0; j < N; j++) {
		index_nonzero(j)=j;
	}

	PermutationMatrix<Dynamic,Dynamic> perm(N);
	perm.setIdentity();
	std::random_shuffle(perm.indices().data(), perm.indices().data()+perm.indices().size());
        index_nonzero = perm * index_nonzero; 
	
	for (j = 0; j < k; j++) {
            xOriginal(index_nonzero(j))= Gaussian_Noise(1);
	}
	
	/** 2. GENERATE MEASUREMENTS **/
	y = A*xOriginal;
        
	if (thereisnoise) {
            sigma= sqrt(xOriginal.squaredNorm()/(pow(10.0, 0.1*SNR)*M));
            for (i = 0; i < M; i++) {
                y(i) =y(i)+ Gaussian_Noise(sigma);
            }
	}

        /** 3. RECOVERY VIA ADMM**/	
        Q=(A_T*A+rho*ID); 
        QI=Q.inverse();
        for (i = 0; i < N; i++) {
            x(i)=0;
            z(i)=0;
            dual(i)=0;
        }
    

        t_tot=0;
        for (t_out = 0; t_out < MAX_ITERATIONS_OUT; t_out ++)    {
                for (i = 0; i < N; i++) {
                        beta(i)=1./(abs(x(i)) +0.1);
                        x_old(i)=1e3;
                        dual(i)=0;
                        z(i)=0;
                }
                for (t = 1; t < MAX_ITERATIONS; t++)    {
                        for (i = 0; i < N; i++)  {
                                if ( abs(x(i))>1e2) {
                                        cerr << x(i) << " not stable "  <<i << " " << t <<endl;
                                        return 0;
                                }
                        }
                        x=QI*( A_T*y + rho*z - dual );
                        for (i = 0; i < N; i++) {
                                if ( abs( x(i)+ dual(i)/rho ) <   LAMBDA*beta(i)/rho )
                                        z(i) = 0;
                                else if ( x(i) + dual(i)/rho > 0)
                                        z(i) = x(i) + dual(i)/rho - LAMBDA*beta(i)/rho;
                                else if ( x(i)+ dual(i)/rho < 0)
					z(i) = x(i)+ dual(i)/rho + LAMBDA*beta(i)/rho;
                        }  
                        dual=dual+rho*(x-z);
                        difference = x- x_old;

                        //STOP CRITERION
                        if (difference.squaredNorm() < ERROR_TOL) {
                                t_tot = t_tot+t;
                                break;
                        }
                        x_old=x;
                }/** end t **/
        }
        final_error = x- xOriginal;
        cout << final_error.squaredNorm()/xOriginal.squaredNorm() << "  " << (final_error.lpNorm<Infinity>()<1e-3) << "  " << t_tot << endl;
	return EXIT_SUCCESS;
}
