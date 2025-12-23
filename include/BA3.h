
/*
    Copyright (C) 2007-2023 Bruce Rannala

    This program is free software: you can redistribute it and/or modify
    it under the terms of the GNU Affero General Public License as
    published by the Free Software Foundation, either version 3 of the
    License, or (at your option) any later version.

    This program is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU Affero General Public License for more details.

    You should have received a copy of the GNU Affero General Public License
    along with this program.  If not, see <http://www.gnu.org/licenses/>.

    Contact: Bruce Rannala <brannala@ucdavis.edu>
    Department of Evolution and Ecology
    University of California Davis
*/

#ifndef BA3_HEADERS
#define BA3_HEADERS
#include <iostream>
#include <iomanip>
#include <fstream>
#include <vector>
#include <algorithm>
#include <regex>
#include <string>
#include <string.h>
#include <sstream>
#include <map>
#include <set>
#include <math.h>
#include <gsl/gsl_sf_gamma.h>
#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>
#include <gsl/gsl_permutation.h>
#include <ctype.h>
#include <stdio.h>
#include <stdlib.h>
#include <getopt.h>
#include <cmath>
#include <chrono>

using namespace std;

const string VERSION="3.4.0";
const string RELEASEDATE="12/23/2025";

// Unified limits - supports both SNP and microsatellite data
// Actual memory allocated dynamically based on dataset size
const int MAXLOCI=100000;
const int MAXALLELE=500;
const int MAXPOPLN=100;
const int MAXINDIV=20000;
const int MAXLINELENGTH=1000000;

// Use int16_t to support up to 500 alleles (values -2 to 499)
typedef int16_t GenotypeType;

// Global variables for actual dataset dimensions (set after reading data)
extern unsigned int gNoLoci;
extern unsigned int gNoIndiv;
extern unsigned int gMaxAlleles;  // Maximum alleles across all loci

struct indiv
{
	GenotypeType (*genotype)[2];  // Dynamically allocated: genotype[noLoci][2]
	vector<int> missingGenotypes;
	unsigned int samplePopln;
	unsigned int migrantPopln;
	unsigned int migrantAge;
	double logL;

	indiv() : genotype(nullptr), samplePopln(0), migrantPopln(0), migrantAge(0), logL(0.0) {}
};

// Memory allocation helpers
void allocateGenotypes(indiv *individuals, unsigned int noIndiv, unsigned int noLoci);
void freeGenotypes(indiv *individuals, unsigned int noIndiv);

struct ancestryProbs
{
	double poplnPostProb[MAXPOPLN][3];
};

// Structure for Savage-Dickey density ratio test statistics
struct SavageDickeyStats
{
	double kernelSum;      // Sum of kernel evaluations at zero
	double sumM;           // Running sum of migration rates (for mean)
	double sumM2;          // Running sum of squared migration rates (for variance)
	long int countNearZero; // Count of samples where m < threshold
	long int nSamples;     // Total number of samples
};

void printBanner(void);
void checkDataSize(unsigned int &outNoIndiv, unsigned int &outNoLoci, unsigned int &outNoPopln, unsigned int &outMaxAlleles);
void readInputFile(indiv *sampleIndiv, unsigned int &noIndiv, unsigned int &noLoci, unsigned int &noPopln, unsigned int *noAlleles);
void getEmpiricalAlleleFreqs(double ***alleleFreqs, indiv *sampleIndiv, unsigned int *noAlleles, unsigned int noPopln, unsigned int noLoci, unsigned int noIndiv);
void fillMigrantCounts(indiv *sampleIndiv, long int ***migrantCounts, unsigned int noIndiv, unsigned int noPopln);
double migCountLogProb(long int ***migrantCounts, double **migrationRates, unsigned int noPopln);
double logLik(const indiv& Indiv, double ***alleleFreqs, double ***logAlleleFreqs, double *FStat, double *log1MinusFStat, unsigned int noLoci);
double oneLocusLogLik(const indiv& Indiv, double ***alleleFreqs, double ***logAlleleFreqs, double *FStat, double *log1MinusFStat, int chosenLocus);
void updateLogAlleleFreqs(double ***logAlleleFreqs, double ***alleleFreqs, unsigned int popln, unsigned int locus, unsigned int noAlleles);
void updateAllLogAlleleFreqs(double ***logAlleleFreqs, double ***alleleFreqs, unsigned int noPopln, unsigned int noLoci, unsigned int *noAlleles);
void updateLogFStat(double *logFStat, double *log1MinusFStat, double *FStat, unsigned int noPopln);
// void proposeMigrantAncDrop(int &migrantPopln, int &migrantAge, int samplePopln, int noPopln, long int ***migrantCounts);
void proposeMigrantAncDrop(unsigned int &migrantPopln, unsigned int &migrantAge, unsigned int samplePopln, int noPopln, long int ***migrantCounts);
void proposeMigrantAncAdd(unsigned int &migrantPopAdd, unsigned int &migrantAgeAdd,unsigned int migrantPopDrop, unsigned int migrantAgeDrop,
						  unsigned int samplePopln, int noPopln);
int countNonEmptyAncestryCategories(long int ***migrantCounts, unsigned int samplePopln, unsigned int noPopln);

// Savage-Dickey density ratio test functions
void initSavageDickeyStats(SavageDickeyStats **sdStats, unsigned int noPopln);
void updateSavageDickeyStats(SavageDickeyStats **sdStats, double **migrationRates,
                              unsigned int noPopln, double bandwidth);
void computeSavageDickeyBayesFactors(SavageDickeyStats **sdStats, unsigned int noPopln,
                                      double priorDensityAtZero, std::ostream &out,
                                      const std::vector<std::string> &poplnNames);
void freeSavageDickeyStats(SavageDickeyStats **sdStats, unsigned int noPopln);

// VCF input functions
void readMetadataFile(const char *metaFileName,
                      std::map<std::string, std::string> &indivToPopln);
void readVCFFile(const char *vcfFileName,
                 const std::map<std::string, std::string> &indivToPopln,
                 indiv *sampleIndiv,
                 unsigned int &noIndiv, unsigned int &noLoci,
                 unsigned int &noPopln, unsigned int *noAlleles);
void checkVCFDataSize(const char *vcfFileName,
                      const std::map<std::string, std::string> &indivToPopln,
                      unsigned int &outNoIndiv, unsigned int &outNoLoci,
                      unsigned int &outNoPopln, unsigned int &outMaxAlleles);

#endif
