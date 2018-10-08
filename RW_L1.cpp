/** Written by Sophie M. Fosson **/
/** 7 Oct 2018 **/
/** Code for simulations of paper "A biconvex analysis for Lasso l1 reweighting", IEEE Signal Process. Lett. 2018
/**  LASSO IRL1 (Algorithm 1) **/

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
#define MAX_ITERATIONS_OUT 1 // =1 for Classical Lasso; >1 for l1 reweigthing (=3 in the paper)
#define ERROR_TOL 1e-10

double sigma, RSE,  rho, FUNC_TRUE, FUNC, aux, BETA, SNR, LAMBDA;
int N, M, h, i, j, k, t, t_out, t_tot, run, penalization;
bool thereisnoise;

//Standard Gaussian r.v. (mean=0; std=sigma)  via Box-Muller Algorithm
double Gaussian_Noise(double std) {
    return std*(sqrt(-2*log( ((double) rand() / (RAND_MAX))  )))*cos(2*pi* ((double) rand() / (RAND_MAX))  );
}

int main (int argc, char *argv[]) {
	if (argc < 6) {
		cerr << "Usage: " << argv[0] << " <signal size N> <sparsity k> <# measurements M> <run> < penalization type [1 = log(|x|+eps); 2 = a|x|-0.5|x|^2; 3 = (|x| + eps)^q ] > <optional: SNR>"  << endl;
                cerr << "Ex: ./gorwl 256 30 100 1" << endl;
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
        
	LAMBDA = 1e-5;
	if (thereisnoise) {
	    LAMBDA = 1e-4;
            sigma= sqrt(xOriginal.squaredNorm()/(pow(10.0, 0.1*SNR)*M));
            for (i = 0; i < M; i++) {
                y(i) =y(i)+ Gaussian_Noise(sigma);
            }
	}


        /** 3. RECOVERY VIA ADMM**/	
	time_t tstart, tend; 
 	tstart = time(0);

        Q=(A_T*A+rho*ID); 
        QI=Q.inverse();
        for (i = 0; i < N; i++) {
            x(i)=0;
            z(i)=0;
            dual(i)=0;
        }
    

        t_tot=0;
        for (t_out = 0; t_out < MAX_ITERATIONS_OUT; t_out ++)    {
		if (penalization == 0) {
                	for (i = 0; i < N; i++) {
                        	beta(i)=10; //for lasso
                        	x_old(i)=1e3;
                        	dual(i)=0;
                        	z(i)=0;
                	}
		} else if (penalization == 1) {
                	for (i = 0; i < N; i++) {
                        	beta(i)=1./(abs(x(i)) +0.1); 
                        	x_old(i)=1e3;
                        	dual(i)=0;
                        	z(i)=0;
                	}
		} else if (penalization == 2) {
                	for (i = 0; i < N; i++) {
				beta(i)=( 2-abs(x(i))>0 ? 2-abs(x(i)) : 0 ); 
                        	x_old(i)=1e3;
                        	dual(i)=0;
                        	z(i)=0;
                	} 
		} else if (penalization == 3) {
                	for (i = 0; i < N; i++) {
				beta(i)=0.5*pow( abs(x(i))+0.1,0.5-1); 
                        	x_old(i)=1e3;
                        	dual(i)=0;
                        	z(i)=0;
                	}
		}


                for (t = 1; t < MAX_ITERATIONS; t++)    {
                        /*for (i = 0; i < N; i++)  {
                                if ( abs(x(i))>1e2) {
                                        cerr << x(i) << " not stable "  <<i << " " << t <<endl;
                                        return 0;
                                }
                        }*/
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

                        /*+ STOP CRITERION **/
                        if (difference.squaredNorm() < ERROR_TOL) {
                                t_tot = t_tot+t;
                                break;
                        }
                        x_old=x;
                }/** end t **/
        }
	tend = time(0);
        final_error = x- xOriginal;
        cout << final_error.squaredNorm()/xOriginal.squaredNorm() << "  " << (final_error.lpNorm<Infinity>()<1e-3) << "  " << t_tot << "  "<<difftime(tend, tstart) << endl;
	return EXIT_SUCCESS;
}
