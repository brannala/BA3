P = BA3
VPATH = src include
CFLAGS = -O3
INCLFLAGS = -I/usr/local/include/
IFLAGS =  -I include
LIBFLAGS =  -L/usr/local/lib
LDLIBS = -lgsl -lgslcblas -lhts -lm
CC=g++ -std=c++11

# Unified build - supports both SNP and microsatellite data
all: BA3

BA3: main.o
	$(CC) $^ $(CFLAGS) $(LIBFLAGS) $(LDLIBS) -o BA3

main.o: main.cpp BA3.h
	$(CC) $(CFLAGS) $(INCLFLAGS) $(IFLAGS) -c $< -o main.o

.Phony: clean
clean:
	$(RM) main.o
	$(RM) BA3

.Phony: tidy
tidy:
	$(RM) main.o

# Legacy targets for backward compatibility
BA3SNP: BA3
	cp BA3 BA3SNP

BA3MSAT: BA3
	cp BA3 BA3MSAT
