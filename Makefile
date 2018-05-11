CC=g++ 
CFLAGS= -w -Wall -I ./eigen/
all:
	$(CC) $(CFLAGS) LASSO_IRL1_LOG_v1.0.cpp -o run_irl1
	$(CC) $(CFLAGS) LASSO_IRL1_IST_LOG_v1.0.cpp -o run_irl1_ist


