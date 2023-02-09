P = BA3
VPATH = src include
CFLAGS = --static -O3
INCLFLAGS = -I/usr/local/include/
IFLAGS =  -I include
LIBFLAGS =  -L/usr/local/lib
LDLIBS = -lgsl -lcblas
CC=g++
all: BA3SNP BA3MSAT
BA3SNP: mainSNP.o 
	$(CC) $^ $(CFLAGS) $(LIBFLAGS) $(LDLIBS) -o BA3SNP
mainSNP.o: main.cpp BA3.h
	$(CC) $(CFLAGS) $(INCLFLAGS) $(IFLAGS) -DSNP -c $< -o mainSNP.o
BA3MSAT: mainMSAT.o 
	$(CC) $^ $(CFLAGS) $(LIBFLAGS) $(LDLIBS) -o BA3MSAT
mainMSAT.o: main.cpp BA3.h
	$(CC) $(CFLAGS) $(INCLFLAGS) $(IFLAGS) -DMSAT -c $< -o mainMSAT.o
.Phony: clean
clean:
	$(RM) mainSNP.o mainMSAT.o
	$(RM) BA3SNP BA3MSAT
.Phony: tidy
tidy:
	$(RM) mainSNP.o mainMSAT.o
