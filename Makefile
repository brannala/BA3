P = BA3
INCLUDE = -Iinclude -I/home/bruce/github/BA3/gsl/include/
LIBS = -L/home/bruce/github/BA3/gsl/lib/ 
CFLAGS = -Wall -O3 --static
LDLIBS = -lgsl -lgslcblas
CC=g++
all:
	$(CC) $(CFLAGS) -s src/main.cpp -o BA3 $(INCLUDE) $(LIBS) $(LDLIBS)
