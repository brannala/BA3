P = BA3
VPATH = src include
CFLAGS = -O3
INCLFLAGS = -I/usr/local/include/
IFLAGS =  -I include
LIBFLAGS =  -L/usr/local/lib
# LDLIBS = -lgsl -lcblas # GSL 1.x, pre-2015
# LDLIBS = -lgsl -lgslcblas # GSL 2.0, 2015
LDLIBS := $(shell if ldconfig -p | grep -q libgslcblas; then echo "-lgsl -lgslcblas"; else echo "-lgsl -lcblas"; fi)
CC = g++ -std=c++11

all: BA3SNP BA3MSAT

BA3SNP: mainSNP.o
	$(CC) $^ $(CFLAGS) $(LIBFLAGS) $(LDLIBS) -o BA3SNP

mainSNP.o: main.cpp BA3.h
	$(CC) $(CFLAGS) $(INCLFLAGS) $(IFLAGS) -DSNP -c $< -o mainSNP.o

BA3MSAT: mainMSAT.o
	$(CC) $^ $(CFLAGS) $(LIBFLAGS) $(LDLIBS) -o BA3MSAT

mainMSAT.o: main.cpp BA3.h
	$(CC) $(CFLAGS) $(INCLFLAGS) $(IFLAGS) -DMSAT -c $< -o mainMSAT.o

.PHONY: clean
clean:
	$(RM) mainSNP.o mainMSAT.o
	$(RM) BA3SNP BA3MSAT

.PHONY: tidy
tidy:
	$(RM) mainSNP.o mainMSAT.o

