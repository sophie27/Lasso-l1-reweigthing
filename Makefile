CC=g++ 
CFLAGS= -w -Wall -I ./eigen/
all:
	$(CC) $(CFLAGS) 12aug2018_rwl1.cpp -o gorwl
	$(CC) $(CFLAGS) 12aug2018_rwl1ist.cpp -o gorwlist



.PHONY: opt
opt:override CFLAGS += -O2
opt:all

.PHONY: prof
prof:override CFLAGS += -p
prof:all

.PHONY: clean
clean:
	rm go*
