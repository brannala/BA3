P = BA3
INCLUDE = -Iinclude -I/Users/bruce/gsl-2.2.1/gsl-2.2.1/gsl
LIBS = -L/Users/bruce/gsl-2.2.1/gsl-2.2.1/ 
CFLAGS = -Wall -O3
LDLIBS = -lgsl -lcblas
CC=g++
all:
	$(CC) -s src/main.cpp -o BA3 $(INCLUDE) $(LDLIBS)
