P = BA3
INCLUDE = -Iinclude -I/Users/bruce/gsl-2.2.1/gsl-2.2.1/
LIBS = -L/Users/bruce/BA3/lib32 
CFLAGS = -Wall -O3 -arch i386 -mmacosx-version-min=10.5
LDLIBS = -lgsl -lcblas
CC=g++
all:
	$(CC) $(CFLAGS) src/main.cpp $(INCLUDE) $(LIBS) $(LDLIBS) -o BA3i386
