P = BA3
INCLUDE = -Iinclude -I/Users/bruce/gsl-2.2.1/gsl-2.2.1/
LIBS = -L/Users/bruce/gsl-2.2.1/gsl-2.2.1/libsX86_64 
CFLAGS = -Wall -O3
LDLIBS = -lgsl -lcblas
CC=g++
all:
	$(CC) $(CFLAGS) src/main.cpp $(INCLUDE) $(LIBS) $(LDLIBS) -o BA3X86_64
