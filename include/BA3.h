
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

// #define SNP  compiling for SNP data */
// #define MSAT /* compiling for microsatellite data */ 

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

using namespace std;

#ifdef SNP
const int MAXLOCI=50000;
const int MAXALLELE=4;
const int MAXPOPLN=100;
const int MAXINDIV=5000;
const int MAXLINELENGTH=1000000;
#endif

#ifdef MSAT
const int MAXLOCI=500;
const int MAXALLELE=500;
const int MAXPOPLN=100;
const int MAXINDIV=5000;
const int MAXLINELENGTH=1000000;
#endif

struct indiv
{
	int genotype[MAXLOCI][2];
	vector<int> missingGenotypes;
	unsigned int samplePopln;
	unsigned int migrantPopln;
	unsigned int migrantAge;
	double logL;
};

struct ancestryProbs
{
	double poplnPostProb[MAXPOPLN][3];
};

void printBanner(void);
void checkDataSize(void);
void readInputFile(indiv *sampleIndiv, unsigned int &noIndiv, unsigned int &noLoci, unsigned int &noPopln, unsigned int *noAlleles);
void getEmpiricalAlleleFreqs(double ***alleleFreqs, indiv *sampleIndiv, unsigned int *noAlleles, unsigned int noPopln, unsigned int noLoci, unsigned int noIndiv);
void fillMigrantCounts(indiv *sampleIndiv, long int ***migrantCounts, unsigned int noIndiv, unsigned int noPopln);
double migCountLogProb(long int ***migrantCounts, double **migrationRates, unsigned int noPopln);
double logLik(indiv Indiv, double ***alleleFreqs, double *FStat, unsigned int noLoci);
double oneLocusLogLik(indiv Indiv, double ***alleleFreqs, double *FStat, int chosenLocus);
// void proposeMigrantAncDrop(int &migrantPopln, int &migrantAge, int samplePopln, int noPopln, long int ***migrantCounts);
void proposeMigrantAncDrop(unsigned int &migrantPopln, unsigned int &migrantAge, unsigned int samplePopln, int noPopln, long int ***migrantCounts);
void proposeMigrantAncAdd(unsigned int &migrantPopAdd, unsigned int &migrantAgeAdd,unsigned int migrantPopDrop, unsigned int migrantAgeDrop, 
						  unsigned int samplePopln, int noPopln);


#endif
