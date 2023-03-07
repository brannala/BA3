
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

#include"../include/BA3.h"

// global command line arguments

struct globalArgs {
	int seed;
	int mciter;
	int sampling;
	int burnin;
	char outfileName[100];
	int usingOutfile;
	double deltaM;
	double deltaA;
	double deltaF;
	int verbose;
	int settings;
	int genotypes;
	int trace;
	int debug;
	int nolikelihood;
} gArgs;

static const char *optString = "s:i:n:b:o:m:a:f:vugtdph?";

static const struct option longOpts[] = {
	{ "seed", required_argument, NULL, 's' },
	{ "iterations", required_argument, NULL, 'i' },
	{ "sampling", required_argument, NULL, 'n' },
	{ "burnin", required_argument, NULL, 'b' },
	{ "output", required_argument, NULL, 'o' },
	{ "deltaM", required_argument, NULL, 'm' },
	{ "deltaA", required_argument, NULL, 'a' },
	{ "deltaF", required_argument, NULL, 'f' },
	{ "verbose", no_argument, NULL, 'v' },
	{ "settings", no_argument, NULL, 'u' },
	{ "genotypes", no_argument, NULL, 'g' },
	{ "debug", no_argument, NULL, 'd' },
	{ "help", no_argument, NULL, 'h' },
	{ "trace", no_argument, NULL, 't' },
	{ "nolikelihood", no_argument, NULL, 'p'},
	{ NULL, no_argument, NULL, 0 }
};

// global variables

gsl_rng * r;
char infileName[100];
std::ofstream mcmcout;
std::ofstream tracefile;
std::ofstream indivout;
std::ifstream mcmcin;

typedef std::map<string, unsigned int> IndivMap;

// map between indivID and name in input file
IndivMap indivIDMap;

// map between poplnID and popln name in input file
IndivMap poplnIDMap;

// map between locusID and locus name in input file
IndivMap locusIDMap;

// map between alleleID and allele name at each locus
IndivMap alleleIDMap[MAXLOCI];

// debugging variables
bool NOMIGRATEMCMC=false;
bool NOANCMCMC=false;
bool NOALLELEMCMC=false;
bool NOFSTATMCMC=false;
bool NOMISSINGDATA=false;
bool NOLIKELIHOOD=false;

std::vector<std::string> split(const std::string str,
                               const std::string regex_str) {
    std::regex regexz(regex_str);
    return {std::sregex_token_iterator(str.begin(), str.end(), regexz, -1),
      std::sregex_token_iterator()};
}

int main( int argc, char *argv[] )
{

	unsigned int noIndiv=0;
	unsigned int noLoci=0;
	unsigned int noPopln=0;
	unsigned int noMissingGenotypes=0;
	unsigned int noAlleles[MAXLOCI];

	/* use gnu getopt code to parse command line options */

	int opt = 0;
	int longIndex = 0;

	/* Initialize gArgs */

	gArgs.seed = 10;
	gArgs.mciter = 500000;
	gArgs.sampling = 100;
	gArgs.burnin = 50000;
	gArgs.deltaM = 0.10;
	gArgs.deltaA = 0.10;
	gArgs.deltaF = 0.10;
	gArgs.verbose = 0;
	gArgs.settings = 0;
	gArgs.genotypes = 0;
	gArgs.trace = 0;
	gArgs.debug = 0;
	gArgs.nolikelihood=0;
	gArgs.usingOutfile = 1;

	strcpy(gArgs.outfileName,"BA3out.txt");

	/* parse command line options */

	opt = getopt_long( argc, argv, optString, longOpts, &longIndex );
	while( opt != -1 ) {
		switch( opt ) {
			case 's':
				gArgs.seed = atoi(optarg);	/* true */
				break;

			case 'i':
				gArgs.mciter = atoi(optarg);
				break;

			case 'n':
				gArgs.sampling = atoi(optarg);
				break;

			case 'b':
				gArgs.burnin = atoi(optarg);
				break;

			case 'o':
				{ strcpy(gArgs.outfileName,optarg); gArgs.usingOutfile = 1; }
				break;

			case 'm':
				gArgs.deltaM = atof(optarg);
				break;

			case 'a':
				gArgs.deltaA = atof(optarg);
				break;

			case 'f':
				gArgs.deltaF = atof(optarg);
				break;

			case 'v':
				gArgs.verbose++;
				break;

			case 'u':
				gArgs.settings++;
				break;

			case 'g':
				gArgs.genotypes++;
				break;

			case 't':
				gArgs.trace++;
				break;

			case 'p':
				gArgs.nolikelihood++;
				break;

			case 'd':
				gArgs.debug++;
				break;

			case 'h':
				std::cout << "usage: BA3 [-sinbomaf] [-vhugpt] file ..." << "\n" << "    mcmc analysis of recent migration rates\n";
				exit(0);

			case '?':
				std::cout << "usage: BA3 [-sinbomaf] [-vhugpt] file ..." << "\n" << "    mcmc analysis of recent migration rates\n";
				exit(0);

			default:
				/* You won't actually get here. */
				break;
		}

		opt = getopt_long( argc, argv, optString, longOpts, &longIndex );
	}

	/* get input file name */

	if (optind < argc)
	{
		strcpy(infileName,argv[optind]);
	}
	else
	{
		std::cout << "error: no input file specified...\n";
		exit(1);
	}

	printBanner();

	unsigned long int mciter=gArgs.mciter;
	unsigned long int aSeed=gArgs.seed;

	if (gArgs.usingOutfile)
	{
		mcmcout.open(gArgs.outfileName, std::ios::out);
	}

	if (gArgs.trace)
	{
		tracefile.open("BA3trace.txt", std::ios::out);
	}

	if (gArgs.genotypes)
	{
		indivout.open("BA3indiv.txt", std::ios::out);
	}

	if (gArgs.nolikelihood)
	{
		NOLIKELIHOOD=true;
	}

	mcmcin.open(infileName, std::ios::in);
	if (!mcmcin)
	{
		std::cerr << "\nerror: cannot open file: "  << infileName << " quitting...\n";
		exit(1);
	}


	// initialize gsl random number generator and specify a seed
	r = gsl_rng_alloc (gsl_rng_taus);
	gsl_rng_set(r,aSeed);

	// initialize dynamic structures for data and parameters

	indiv *sampleIndiv = new indiv[MAXINDIV];
	if(sampleIndiv==NULL) cerr << "ran out of memory";

	// check that proposal step parameters are sane (e.g., <0 and <=1)
	  if((gArgs.deltaM<=0)||(gArgs.deltaM>1))
	    { std::cerr << "\nerror: deltaM=" << gArgs.deltaM << " not in interval (0,1]. "  << infileName << " quitting...\n";
	      exit(1); }
	  if((gArgs.deltaA<=0)||(gArgs.deltaA>1))
	    { std::cerr << "\nerror: deltaA=" << gArgs.deltaA << " not in interval (0,1]. "  << infileName << " quitting...\n";
	      exit(1); }
	  if((gArgs.deltaF<=0)||(gArgs.deltaF>1))
	    { std::cerr << "\nerror: deltaF=" << gArgs.deltaF << " not in interval (0,1]. "  << infileName << " quitting...\n";
	      exit(1); }

	vector<int> missingData;


	checkDataSize();
	readInputFile(sampleIndiv, noIndiv, noLoci, noPopln, noAlleles);

	double ***ancP;
	if((ancP = new double**[noIndiv])==0) cerr << "ran out of memory";
	for(unsigned int i = 0; i < noIndiv; i++)
	{
		ancP[i] = new double*[noPopln];
		for(unsigned int j = 0; j < noPopln; j++)
			ancP[i][j] = new double[3];
	}


	for (unsigned int j=0; j<noIndiv; j++)
	{
		for(unsigned int l=0;l<noPopln;l++)
		{
			ancP[j][l][0]=0.0;
			ancP[j][l][1]=0.0;
			ancP[j][l][2]=0.0;
		}
	}



	// identify missing data and initialize imputed genotypes
	for (unsigned int i=0; i<noIndiv; i++)
	{
		bool hasMissing=false;
		for (unsigned int j=0; j<noLoci; j++)
		{
		  if((sampleIndiv[i].genotype[j][0] == -1)||(sampleIndiv[i].genotype[j][1] == -1))
		     {

		       if (noAlleles[j]==0)
			 {
			   std::cerr << "\nerror: missing all data at locus! quitting...\n";
			   exit(1);
			 }

		       
		       noMissingGenotypes+=1;
		       sampleIndiv[i].missingGenotypes.push_back(j);
		       sampleIndiv[i].genotype[j][0] = gsl_rng_uniform_int(r, noAlleles[j]);
		       sampleIndiv[i].genotype[j][1] = gsl_rng_uniform_int(r, noAlleles[j]);
		       if (!hasMissing)
			 {
			   missingData.push_back(i);
			   hasMissing=true;
			 }
		     }
		}
	}

	size_t N=noIndiv;
	gsl_permutation * p = gsl_permutation_alloc (N);

	if (gArgs.verbose)
	{
		cout << "\nInput file: " << infileName;
		cout << "\nOutput file: " << gArgs.outfileName << "\n";
		cout << "Individuals: " << noIndiv << " Populations: " << noPopln << " Loci: " << noLoci;
		cout << " Missing genotypes: " << noMissingGenotypes << "\n\n";
		cout << "Locus:(Number of Alleles)\n";
		for(unsigned int l=0; l<noLoci; l++)
		{
			string locusName;
			IndivMap::iterator iterLocus = locusIDMap.begin();
			while (iterLocus != locusIDMap.end())
			{
				if(iterLocus->second == l)
					locusName=iterLocus->first;
				iterLocus++;
			}
			if((l % 10)==0) cout << "\n";
			cout << locusName << ":" << noAlleles[l] << " ";
		}
		cout << "\n\n";
	}


	double ***alleleFreqs;
	alleleFreqs = new double**[noPopln];
	for(unsigned int i = 0; i < noPopln; i++)
	{
		alleleFreqs[i] = new double*[noLoci];
		for(unsigned int j = 0; j < noLoci; j++)
			alleleFreqs[i][j] = new double[MAXALLELE];
	}

	double ***avgAlleleFreqs;
	avgAlleleFreqs = new double**[noPopln];
	for(unsigned int i = 0; i < noPopln; i++)
	{
		avgAlleleFreqs[i] = new double*[noLoci];
		for(unsigned int j = 0; j < noLoci; j++)
			avgAlleleFreqs[i][j] = new double[MAXALLELE];
	}

	for(unsigned int l = 0; l < noPopln; l++)
		for(unsigned int i = 0; i < noLoci; i++)
			for(unsigned int j = 0; j < noAlleles[i]; j++)
				avgAlleleFreqs[l][i][j]=0.0;

	double ***varAlleleFreqs;
	varAlleleFreqs = new double**[noPopln];
	for(unsigned int i = 0; i < noPopln; i++)
	{
		varAlleleFreqs[i] = new double*[noLoci];
		for(unsigned int j = 0; j < noLoci; j++)
			varAlleleFreqs[i][j] = new double[MAXALLELE];
	}

	for(unsigned int l = 0; l < noPopln; l++)
		for(unsigned int i = 0; i < noLoci; i++)
			for(unsigned int j = 0; j < noAlleles[i]; j++)
				varAlleleFreqs[l][i][j]=0.0;

	double **migrationRates;
	migrationRates = new double*[noPopln];
	for(unsigned int i = 0; i < noPopln; i++)
		migrationRates[i] = new double[noPopln+1];

	double **avgMigrationRates;
	avgMigrationRates = new double*[noPopln];
	for(unsigned int i = 0; i < noPopln; i++)
		avgMigrationRates[i] = new double[noPopln+1];

	for(unsigned int i = 0; i < noPopln; i++)
		for (unsigned int j = 0; j < noPopln; j++)
			avgMigrationRates[i][j]=0.0;

	double **varMigrationRates;
	varMigrationRates = new double*[noPopln];
	for(unsigned int i = 0; i < noPopln; i++)
		varMigrationRates[i] = new double[noPopln+1];

	for(unsigned int i = 0; i < noPopln; i++)
		for (unsigned int j = 0; j < noPopln; j++)
			varMigrationRates[i][j]=0.0;

	long int ***migrantCounts;
	migrantCounts = new long int**[noPopln];
	for(unsigned int i = 0; i < noPopln; i++)
	{
		migrantCounts[i] = new long int*[noPopln];
		for(unsigned int j = 0; j < noPopln; j++)
			migrantCounts[i][j] = new long int[3];
	}

	for(unsigned int i = 0; i < noPopln; i++)
		for(unsigned int j = 0; j < noPopln; j++)
			for(int k = 0; k < 3; k++)
				migrantCounts[i][j][k]=0;


	double *FStat = new double[noPopln];
	double *avgFStat = new double[noPopln];
	double *varFStat = new double[noPopln];
	for(unsigned int i=0; i< noPopln; i++)
	{
		FStat[i]=0.0;
		avgFStat[i]=0.0;
		varFStat[i]=0.0;
	}

	for(unsigned int l=0; l<noPopln; l++)
	{
		for(unsigned int k=0; k<=noPopln; k++)
			if((l!=k)&&(k!=noPopln))
				migrationRates[l][k] = (1.0/3.0)*(1.0/noPopln);
		migrationRates[l][l] = 1.0-(1.0/3.0)*((noPopln-1.0)/noPopln);
		migrationRates[l][noPopln] = (1.0/3.0)*((noPopln-1.0)/noPopln);
	}

	getEmpiricalAlleleFreqs(alleleFreqs,sampleIndiv,noAlleles,noPopln,noLoci,noIndiv);

	/* uncomment for uniform distn of identical freqs in all populations */
    /*	for(int i=0; i < noPopln; i++)
		for(int j=0; j < noLoci; j++)
			for(int l=0; l < noAlleles[j]; l++)
				alleleFreqs[i][j][l] = 1.0/noAlleles[j]; */

	if(gArgs.debug)
	{
		for(unsigned int i=0; i < noPopln; i++)
			for(unsigned int j=0; j < noLoci; j++)
			{
				cout << "aF[" << i << "][" << j << "]";
				for(unsigned int l=0; l < noAlleles[j]; l++)
				{
					std::cout.setf(std::ios::fixed, std::ios::floatfield);
					cout << setprecision(3) << "[" << l << "]:" << alleleFreqs[i][j][l] << " ";
				}
				cout << "\n";
			}
	}

	for(unsigned int l = 0; l < noIndiv; l++)
	{	sampleIndiv[l].migrantAge=0; sampleIndiv[l].migrantPopln=sampleIndiv[l].samplePopln; }
	fillMigrantCounts(sampleIndiv,migrantCounts,noIndiv,noPopln);

	for(unsigned int i = 0; i < noIndiv; i++)
		sampleIndiv[i].logL = logLik(sampleIndiv[i],alleleFreqs,FStat,noLoci);

	if(gArgs.debug)
	{
		std::cout << "mCLp: " << migCountLogProb(migrantCounts,migrationRates,noPopln) << "\n";
		for(unsigned int i = 0; i < noPopln; i++)
			for(unsigned int j = 0; j < noPopln; j++)
				for(int k = 0; k < 3; k++)
				{
					std::cout << "migrantCounts[" << i << "][" << j << "][" << k << "]:";
					std::cout << (int) migrantCounts[i][j][k] << "\n";
				}
	}



	/* Metropolis-Hastings algorithm */
	/* Loop for mciter */

	/* debugging variables */

	unsigned long int iter=0;
	double migrationAcceptRate=0.0;
	double ancestryAcceptRate=0.0;
	double allelefreqAcceptRate=0.0;
	double FStatAcceptRate=0.0;
	double GenotypeAcceptRate=0.0;

	gsl_permutation_init (p);

	for(unsigned int i = 1; i <= mciter; i++)

	{

		indiv tempIndiv;
		long int chosenIndiv;
		unsigned int migrantPopln, migrantAge, samplePopln;
		unsigned int migrantPopAdd, migrantAgeAdd;
		double alpha;
		double logPrMHR;
		double logLprop;
		double dtLogL;



if(!NOANCMCMC)
{

		/* propose modified migrant ancestry for a random individual */


		samplePopln = gsl_rng_uniform_int(r, noPopln);

		proposeMigrantAncDrop(migrantPopln, migrantAge, samplePopln, noPopln, migrantCounts);
		gsl_ran_shuffle (r, p->data, N, sizeof(size_t));
		bool foundIndiv=false;
		int k=0;
		while (!foundIndiv)
		{
			if((sampleIndiv[gsl_permutation_get(p,k)].migrantPopln == migrantPopln)&&
			   (sampleIndiv[gsl_permutation_get(p,k)].migrantAge == migrantAge)&&
			   (sampleIndiv[gsl_permutation_get(p,k)].samplePopln == samplePopln))
			{ chosenIndiv=gsl_permutation_get(p,k); foundIndiv = true; }
			else { k++; }
		}
		proposeMigrantAncAdd(migrantPopAdd, migrantAgeAdd,migrantPopln, migrantAge, samplePopln, noPopln);

		tempIndiv.samplePopln = sampleIndiv[chosenIndiv].samplePopln;
		tempIndiv.migrantPopln = migrantPopAdd;
		tempIndiv.migrantAge = migrantAgeAdd;
		for(unsigned int j = 0; j < noLoci; j++)
		{
			tempIndiv.genotype[j][0] = sampleIndiv[chosenIndiv].genotype[j][0];
			tempIndiv.genotype[j][1] = sampleIndiv[chosenIndiv].genotype[j][1];
		}

		// calculate change of logL for genetic data with new migrant ancestry
		if (!NOLIKELIHOOD)
		{
			logLprop = logLik(tempIndiv,alleleFreqs,FStat,noLoci);
			dtLogL = logLprop - sampleIndiv[chosenIndiv].logL;
		}

		// calculate change of logPr for migrant counts with new migrant ancestry
		double dtLogPrCount=0.0;

		if ((tempIndiv.migrantPopln != sampleIndiv[chosenIndiv].migrantPopln)||
			(tempIndiv.migrantAge != sampleIndiv[chosenIndiv].migrantAge))
		{
		if(tempIndiv.migrantAge == 0)
		{
			dtLogPrCount += (log(1.0-3.0*migrationRates[tempIndiv.samplePopln][noPopln])-
			log(migrantCounts[tempIndiv.samplePopln][tempIndiv.migrantPopln][0]+1.0));

		}
		else
		if (tempIndiv.migrantAge == 1)
		{
			dtLogPrCount += (log(migrationRates[tempIndiv.samplePopln][tempIndiv.migrantPopln])-
			log(migrantCounts[tempIndiv.samplePopln][tempIndiv.migrantPopln][1]+1.0));

		}
		else
		if(tempIndiv.migrantAge == 2)
		{
			dtLogPrCount += (log(2.0*migrationRates[tempIndiv.samplePopln][tempIndiv.migrantPopln])-
			log(migrantCounts[tempIndiv.samplePopln][tempIndiv.migrantPopln][2]+1.0));
		}

		if(sampleIndiv[chosenIndiv].migrantAge == 0)
		{
			dtLogPrCount -= (log(1.0-3.0*migrationRates[sampleIndiv[chosenIndiv].samplePopln][noPopln]) -
			log(migrantCounts[sampleIndiv[chosenIndiv].samplePopln][sampleIndiv[chosenIndiv].migrantPopln][0]));
		}
		else
			if (sampleIndiv[chosenIndiv].migrantAge == 1)
			{
				dtLogPrCount -= (log(migrationRates[sampleIndiv[chosenIndiv].samplePopln][sampleIndiv[chosenIndiv].migrantPopln]) -
				log(migrantCounts[sampleIndiv[chosenIndiv].samplePopln][sampleIndiv[chosenIndiv].migrantPopln][1]));
			}
			else
				if(sampleIndiv[chosenIndiv].migrantAge == 2)
				{
					dtLogPrCount -= (log(2.0*migrationRates[sampleIndiv[chosenIndiv].samplePopln][sampleIndiv[chosenIndiv].migrantPopln]) -
					log(migrantCounts[sampleIndiv[chosenIndiv].samplePopln][sampleIndiv[chosenIndiv].migrantPopln][2]));
				}
		}
		// Acceptance-rejection step
		alpha = gsl_rng_uniform(r);
		if(!NOLIKELIHOOD)
			logPrMHR = dtLogPrCount + dtLogL;
		else
			logPrMHR = dtLogPrCount;

		if(alpha <= exp(logPrMHR))
		{
			sampleIndiv[chosenIndiv].migrantAge = tempIndiv.migrantAge;
			sampleIndiv[chosenIndiv].migrantPopln = tempIndiv.migrantPopln;
			sampleIndiv[chosenIndiv].logL = logLprop;
			fillMigrantCounts(sampleIndiv,migrantCounts,noIndiv,noPopln);
			ancestryAcceptRate = (1.0/i)+((i-1.0)/i)*ancestryAcceptRate;
		}
		else ancestryAcceptRate = ((i-1.0)/i)*ancestryAcceptRate;
}

if(!NOMIGRATEMCMC)
{


	/* propose a change to migration matrix */

	double propMigrationRates[MAXPOPLN];
	unsigned int sourcePopln;
	sourcePopln = gsl_rng_uniform_int(r,noPopln);
	migrantPopln = gsl_rng_uniform_int(r,noPopln-1);
	if(migrantPopln >= sourcePopln) migrantPopln += 1;
	for (unsigned int j=0; j < noPopln; j++)
	{
		if(j != migrantPopln)
			propMigrationRates[j] = migrationRates[sourcePopln][j];
		else
		{
			propMigrationRates[j] = gArgs.deltaM*(gsl_rng_uniform(r)-0.5) + migrationRates[sourcePopln][j];
			while ((propMigrationRates[j]<0)||(propMigrationRates[j]>(1.0/3.0-migrationRates[sourcePopln][noPopln]+migrationRates[sourcePopln][j])))
			{
				if (propMigrationRates[j]<0) {
					propMigrationRates[j]=abs(propMigrationRates[j]);
				}
				if (propMigrationRates[j]>(1.0/3.0-migrationRates[sourcePopln][noPopln]+migrationRates[sourcePopln][j]))
					propMigrationRates[j]=2.0*(1.0/3.0-migrationRates[sourcePopln][noPopln]+migrationRates[sourcePopln][j])-propMigrationRates[j];
			}
		}
	}
	propMigrationRates[noPopln] = migrationRates[sourcePopln][noPopln] - migrationRates[sourcePopln][migrantPopln] + propMigrationRates[migrantPopln];
		propMigrationRates[sourcePopln] = 1.0 - propMigrationRates[noPopln];
	double logPrCurr=0.0, logPrProp=0.0;
	for(unsigned int l = 0; l < noPopln; l++)
		if(sourcePopln != l)
			{
				logPrCurr += migrantCounts[sourcePopln][l][1]*log(migrationRates[sourcePopln][l])-gsl_sf_lnfact(migrantCounts[sourcePopln][l][1]);
				logPrCurr += migrantCounts[sourcePopln][l][2]*log(2.0*migrationRates[sourcePopln][l])-gsl_sf_lnfact(migrantCounts[sourcePopln][l][2]);
			}
	logPrCurr += migrantCounts[sourcePopln][sourcePopln][0]*log(1.0-3*migrationRates[sourcePopln][noPopln]) - gsl_sf_lnfact(migrantCounts[sourcePopln][sourcePopln][0]);
	for(unsigned int l = 0; l < noPopln; l++)
		if(sourcePopln != l)
		{
			logPrProp += migrantCounts[sourcePopln][l][1]*log(propMigrationRates[l])-gsl_sf_lnfact(migrantCounts[sourcePopln][l][1]);
			logPrProp += migrantCounts[sourcePopln][l][2]*log(2.0*propMigrationRates[l])-gsl_sf_lnfact(migrantCounts[sourcePopln][l][2]);
		}
	logPrProp += migrantCounts[sourcePopln][sourcePopln][0]*log(1.0-3*propMigrationRates[noPopln]) - gsl_sf_lnfact(migrantCounts[sourcePopln][sourcePopln][0]);

	// Acceptance-rejection step
	alpha = gsl_rng_uniform(r);

	// debugging
	// alpha = 0;


	logPrMHR = logPrProp - logPrCurr;
	if(alpha <= exp(logPrMHR))
	{
		for (unsigned int k=0; k<=noPopln; k++)
			migrationRates[sourcePopln][k] = propMigrationRates[k];
		migrationAcceptRate = (1.0/i)+((i-1.0)/i)*migrationAcceptRate;
	}
	else migrationAcceptRate = ((i-1.0)/i)*migrationAcceptRate;

}

if(!NOALLELEMCMC)
{

	/* propose a change to a population allele frequency */

	double propAlleleFreq[MAXALLELE];
	double origAlleleFreq[MAXALLELE];
	unsigned int chosenPopln = gsl_rng_uniform_int(r,noPopln);
	unsigned int chosenLocus = gsl_rng_uniform_int(r,noLoci);
	unsigned int chosenAllele = gsl_rng_uniform_int(r,noAlleles[chosenLocus]);
	double logLOrig[MAXINDIV], logLProposed[MAXINDIV];
	if (!NOLIKELIHOOD)
	{
		for (unsigned int l=0; l<noIndiv; l++)
			logLOrig[l] = oneLocusLogLik(sampleIndiv[l],alleleFreqs,FStat,chosenLocus);
	}

	propAlleleFreq[chosenAllele] = 	abs(alleleFreqs[chosenPopln][chosenLocus][chosenAllele]+(gsl_rng_uniform(r)-0.5)*gArgs.deltaA);

	// added July 20, 2011 (further modified Nov 7, 2011)
	if (propAlleleFreq[chosenAllele] > 1.0) {
		propAlleleFreq[chosenAllele] = abs(1.0 - propAlleleFreq[chosenAllele]);
	}

	for (unsigned int l=0; l<noAlleles[chosenLocus]; l++)
		if(l!=chosenAllele)
		{
				propAlleleFreq[l] = exp(log(alleleFreqs[chosenPopln][chosenLocus][l])+
				  log(1.0-propAlleleFreq[chosenAllele])-log(1.0-alleleFreqs[chosenPopln][chosenLocus][chosenAllele]));
		}

	// check that allele freqs sum to one and correct if needed.
	double sum=0.0;
	for (unsigned int l=0; l<noAlleles[chosenLocus]; l++)
		sum+=propAlleleFreq[l];
	if (sum>1.0001||sum<0.9999) {
	  for(unsigned int l=0; l<noAlleles[chosenLocus];l++)
	    propAlleleFreq[l]=propAlleleFreq[l]/sum;


	  /*		double sum2=0.0;
		for (int l=0; l<noAlleles[chosenLocus]-1; l++)
	       	sum2+=propAlleleFreq[l];
			propAlleleFreq[noAlleles[chosenLocus]] = 1.0 - sum2; */
	}

	for (unsigned int l=0; l<noAlleles[chosenLocus]; l++)
	{
		origAlleleFreq[l] = alleleFreqs[chosenPopln][chosenLocus][l];
		alleleFreqs[chosenPopln][chosenLocus][l] = propAlleleFreq[l];
	}
	if (!NOLIKELIHOOD)
	{
		for (unsigned int l=0; l<noIndiv; l++)
			logLProposed[l] = oneLocusLogLik(sampleIndiv[l],alleleFreqs,FStat,chosenLocus);
		logLprop=0; dtLogL=0;
		for (unsigned int l=0; l<noIndiv; l++)
		{
			dtLogL +=  (logLProposed[l] - logLOrig[l]);
		}
	}
	// Acceptance-rejection step
	alpha = gsl_rng_uniform(r);
	if (!NOLIKELIHOOD)
	{
		logPrMHR = dtLogL;
	}
	else
	{
		logPrMHR = 1;
	}

	if(alpha <= exp(logPrMHR))
	{
		if (!NOLIKELIHOOD)
		{
			for (unsigned int k=0; k<noIndiv; k++)
				sampleIndiv[k].logL = sampleIndiv[k].logL - logLOrig[k] + logLProposed[k];
		}
		allelefreqAcceptRate = (1.0/i)+((i-1.0)/i)*allelefreqAcceptRate;
	}
	else
	{
		for (unsigned int k=0; k<noAlleles[chosenLocus]; k++)
			alleleFreqs[chosenPopln][chosenLocus][k] = origAlleleFreq[k];
		allelefreqAcceptRate = ((i-1.0)/i)*allelefreqAcceptRate;
	}

}

if(!NOFSTATMCMC)
{
	/* propose a change to a population inbreeding coefficient */

	double origLik=0.0,propLik=0.0;
	double logLVec[MAXINDIV];
	double propFStat[MAXPOPLN];
	for (unsigned int l=0; l<noPopln; l++)
		propFStat[l]=FStat[l];
	unsigned int chosenPopln = gsl_rng_uniform_int(r,noPopln);
	double prop=abs(FStat[chosenPopln]+(gsl_rng_uniform(r)-0.5)*gArgs.deltaF);
	if (prop<=1)
		propFStat[chosenPopln] = prop;
	else
		propFStat[chosenPopln] = 2.0-prop;
	if (!NOLIKELIHOOD)
	{
		for(unsigned int l=0; l<noIndiv; l++)
		{
			origLik+=sampleIndiv[l].logL;
			if((sampleIndiv[l].samplePopln==chosenPopln)||(sampleIndiv[l].migrantPopln==chosenPopln))
			{
				logLVec[l]=logLik(sampleIndiv[l],alleleFreqs,propFStat,noLoci);
			}
			else
			{
				logLVec[l]=sampleIndiv[l].logL;
			}
			propLik+=logLVec[l];
		}
	}

		// Acceptance-rejection step
	alpha = gsl_rng_uniform(r);
	if (!NOLIKELIHOOD)
	{
		logPrMHR = propLik - origLik;
	}
	else
	{
		logPrMHR = 1;
	}
	if(alpha <= exp(logPrMHR))
	{
		FStat[chosenPopln] = propFStat[chosenPopln];
		if (!NOLIKELIHOOD)
		{
			for(unsigned int l=0; l<noIndiv; l++)
			{
				sampleIndiv[l].logL = logLVec[l];
			}
		}
		FStatAcceptRate = (1.0/i)+((i-1.0)/i)*FStatAcceptRate;
	}
	else
	FStatAcceptRate = ((i-1.0)/i)*FStatAcceptRate;
}

if(!NOMISSINGDATA)
{
	/* propose a change to a missing genotype */
	if(noMissingGenotypes>0)
	{
		double propLogL=0.0, origLogL=0.0;
		int chooseIndiv = missingData[gsl_rng_uniform_int(r, missingData.size())];
		int chosenLocus = sampleIndiv[chooseIndiv].missingGenotypes[gsl_rng_uniform_int(r, sampleIndiv[chooseIndiv].missingGenotypes.size())];
		origLogL=oneLocusLogLik(sampleIndiv[chooseIndiv],alleleFreqs,FStat,chosenLocus);
		int origAllele1=sampleIndiv[chooseIndiv].genotype[chosenLocus][0];
		int origAllele2=sampleIndiv[chooseIndiv].genotype[chosenLocus][1];
		sampleIndiv[chooseIndiv].genotype[chosenLocus][0] = gsl_rng_uniform_int(r, noAlleles[chosenLocus]);
		sampleIndiv[chooseIndiv].genotype[chosenLocus][1] = gsl_rng_uniform_int(r, noAlleles[chosenLocus]);
		propLogL=oneLocusLogLik(sampleIndiv[chooseIndiv],alleleFreqs,FStat,chosenLocus);

		// Acceptance-rejection step
		alpha = gsl_rng_uniform(r);
		if (!NOLIKELIHOOD)
		{
			logPrMHR = propLogL - origLogL;
		}
		else
		{
			logPrMHR = 1;
		}

		if(alpha <= exp(logPrMHR))
		{
			GenotypeAcceptRate = (1.0/i)+((i-1.0)/i)*GenotypeAcceptRate;
			sampleIndiv[chooseIndiv].logL = sampleIndiv[chooseIndiv].logL - origLogL + propLogL;
		}
		else
		{
			sampleIndiv[chooseIndiv].genotype[chosenLocus][0] = origAllele1;
			sampleIndiv[chooseIndiv].genotype[chosenLocus][1] = origAllele2;
			GenotypeAcceptRate = ((i-1.0)/i)*GenotypeAcceptRate;
		}
	}
}

	// Print logL to trace file
		if (i==1) {
			tracefile << "State\t" << "LogProb\t";
		}
		if(gArgs.trace && ((i % gArgs.sampling)==0))
		{
			double logLG=0.0, logLM=0.0;
			for (unsigned int m=0; m < noIndiv; m++) { logLG += sampleIndiv[m].logL; }
			logLM = migCountLogProb(migrantCounts,migrationRates,noPopln);
			tracefile << i << "\t" << logLM + logLG << "\t";
		}



	// Summarize mcmc samples and print likelihoods to screen

		if(gArgs.verbose)
		{
			if((i % 10000)==0)
			{
				double logLG=0.0, logLM=0.0;
				if (!NOLIKELIHOOD)
				{
					for(unsigned int m=0; m < noIndiv; m++) { logLG += sampleIndiv[m].logL; }
				}
				else
				{
					logLG=0;
				}

				logLM = migCountLogProb(migrantCounts,migrationRates,noPopln);
				std::cout.setf(std::ios::fixed, std::ios::floatfield);
				std::cout << std::setprecision(2) << "logP(M): " << logLM << " logL(G): ";
				std::cout << logLG << " logL: " << logLM + logLG << " \% done: " << std::flush;
				if(i < gArgs.burnin)
					std::cout << "[" << i/(mciter*1.0) << "]" << std::flush;
				else
					std::cout << "(" << i/(mciter*1.0) << ")" << std::flush;
				std::cout << " \% accepted: (" << migrationAcceptRate << ", " << ancestryAcceptRate << ", " << allelefreqAcceptRate << ", ";
				std::cout << FStatAcceptRate << ", " << GenotypeAcceptRate << ")" << "\r" << std::flush;
			}
			if((i % (mciter/10))==0)
			{
				double logLG=0.0, logLM=0.0;
				if (!NOLIKELIHOOD)
				{
					for(unsigned int m=0; m < noIndiv; m++) { logLG += sampleIndiv[m].logL; }
				}
				else
				{
					logLG=0;
				}

				logLM = migCountLogProb(migrantCounts,migrationRates,noPopln);
				std::cout.setf(std::ios::fixed, std::ios::floatfield);
				std::cout << std::setprecision(2) << "logP(M): " << logLM << " logL(G): " << logLG << " logL: ";
				std::cout << logLM + logLG << " \% done: "<< std::flush;
				if(i < gArgs.burnin)
					std::cout << "[" << i/(mciter*1.0) << "]" << std::flush;
				else
					std::cout << "(" << i/(mciter*1.0) << ")" << std::flush;
				std::cout << " \% accepted: (" << migrationAcceptRate << ", " << ancestryAcceptRate << ", " << allelefreqAcceptRate << ", ";
				std::cout << FStatAcceptRate <<  ", " << GenotypeAcceptRate << ")" << "\n" << std::flush;
			}

		}
		else
		{
			std::cout.setf(std::ios::fixed, std::ios::floatfield);
			std::cout << std::setprecision(2);
			if((i % (mciter/100))==0)
			{
				int l=0;
				cout << "\r \% done: " << i/(mciter*1.0) << " | " << flush;
				while(l <= (50*(i/(mciter*1.0))))
				{
					cout << "=" << flush;
					l++;
				}
			}
		}

		// Print migration rates to trace file
		if ((i==1) && gArgs.trace)
		{
			for (unsigned int l = 0; l < noPopln; l++)
				for (unsigned int k = 0; k < noPopln; k++)
					tracefile << "m[" << l << "][" << k << "]\t";
			tracefile << "\n";
		}
		if(gArgs.trace && ((i % gArgs.sampling)==0))
		{
			for (unsigned int l = 0; l < noPopln; l++)
				for (unsigned int k = 0; k < noPopln; k++)
					tracefile << migrationRates[l][k] << "\t" ;
			tracefile << "\n";
		}





		if(((i % gArgs.sampling)==0)&&(i > gArgs.burnin))
		{
			double sqrDiffMean=0.0;

			for (unsigned int l = 0; l < noPopln; l++)
				for (unsigned int k = 0; k < noPopln; k++)
				{
					if(iter > 1)
					{
						sqrDiffMean=(migrationRates[l][k]-avgMigrationRates[l][k])*(migrationRates[l][k]-avgMigrationRates[l][k])/(iter+1.0);
						varMigrationRates[l][k] = ((iter-1.0)/iter)*varMigrationRates[l][k]+sqrDiffMean;
					}
						avgMigrationRates[l][k] = avgMigrationRates[l][k]+(migrationRates[l][k]-avgMigrationRates[l][k])/(1.0+iter);
				}

			for (unsigned int l = 0; l < noPopln; l++)
				for (unsigned int k = 0; k < noLoci; k++)
					for(unsigned int m = 0; m < noAlleles[k]; m++)
					{
						if(iter > 1)
						{
							sqrDiffMean=(alleleFreqs[l][k][m]-avgAlleleFreqs[l][k][m])*(alleleFreqs[l][k][m]-avgAlleleFreqs[l][k][m])/(iter+1.0);
							varAlleleFreqs[l][k][m] = ((iter-1.0)/iter)*varAlleleFreqs[l][k][m]+sqrDiffMean;
						}
						avgAlleleFreqs[l][k][m] = avgAlleleFreqs[l][k][m]+(alleleFreqs[l][k][m]-avgAlleleFreqs[l][k][m])/(1.0+iter);
					}

			for (unsigned int l=0; l < noPopln; l++)
			{
				if(iter > 1)
				{
					sqrDiffMean=(FStat[l]-avgFStat[l])*(FStat[l]-avgFStat[l])/(iter+1.0);
					varFStat[l] = ((iter-1.0)/iter)*varFStat[l]+sqrDiffMean;
				}
				avgFStat[l] = avgFStat[l]+(FStat[l]-avgFStat[l])/(1.0+iter);
			}

			for (unsigned int l=0; l < noIndiv; l++)
			{
				if (sampleIndiv[l].migrantAge == 0)
				{
					ancP[l][sampleIndiv[l].samplePopln][0]=
					ancP[l][sampleIndiv[l].samplePopln][0]*(iter/(iter+1.0))+(1.0/(iter+1.0));
					for(unsigned int k=0; k<noPopln; k++)
						for(int b=0; b<3; b++)
							if (!((k==sampleIndiv[l].samplePopln)&&(b==0)))
							{
								ancP[l][k][b]=ancP[l][k][b]*(iter/(iter+1.0));
							}
				}
				else
					if (sampleIndiv[l].migrantAge == 1)
					{
						ancP[l][sampleIndiv[l].migrantPopln][1]=
						ancP[l][sampleIndiv[l].migrantPopln][1]*(iter/(iter+1.0))+(1.0/(iter+1.0));
						for(unsigned int k=0; k<noPopln; k++)
							for(unsigned int b=0; b<3; b++)
								if (!((k==sampleIndiv[l].migrantPopln)&&(b==1)))
								{
									ancP[l][k][b]=ancP[l][k][b]*(iter/(iter+1.0));
								}
					}
					else
						if(sampleIndiv[l].migrantAge == 2)
						{
							ancP[l][sampleIndiv[l].migrantPopln][2]=
							ancP[l][sampleIndiv[l].migrantPopln][2]*(iter/(iter+1.0))+(1.0/(iter+1.0));
							for(unsigned int k=0; k<noPopln; k++)
								for(unsigned int b=0; b<3; b++)
									if (!((k==sampleIndiv[l].migrantPopln)&&(b==2)))
									{
										ancP[l][k][b]=ancP[l][k][b]*(iter/(iter+1.0));
									}
						}

				}


		iter+=1;
		}
}

	if (!gArgs.verbose)
		cout << " | \n";

if(gArgs.debug)
	{
	for(unsigned int l = 0; l < noPopln; l++)
		for(unsigned int j = 0; j < noPopln; j++)
			for(int k = 0; k < 3; k++)
			{
				std::cout << "\nmigrantCounts[" << l << "][" << j << "][" << k << "]:";
				std::cout << migrantCounts[l][j][k];
			}
	for (unsigned int l = 0; l < noIndiv; l++)
	{
		std::cout << "\nIndivID: " << l << " migrantPop: " << sampleIndiv[l].migrantPopln << " migrantAge: " << sampleIndiv[l].migrantAge;
		std::cout << " samplePop: " << sampleIndiv[l].samplePopln;
		for (unsigned int k=0; k < noLoci; k++)
		{
			string locusName;
			IndivMap::iterator iterLocus = locusIDMap.begin();
			while (iterLocus != locusIDMap.end())
			{
				if(iterLocus->second == k)
					locusName=iterLocus->first;
				iterLocus++;
			}
			cout << " locus: " << locusName << " " << sampleIndiv[l].genotype[k][0] << " " << sampleIndiv[l].genotype[k][1];
		}
	}
	}


	/* print results to output file */

	mcmcout << "\nInput file: " << infileName << "\n";
	if (gArgs.settings)
	{
		mcmcout << "Random seed=" << gArgs.seed << " MCMC iterations=" << gArgs.mciter << " Burn-in=" << gArgs.burnin << " Sampling interval=" << gArgs.sampling << "\n";
		mcmcout << "Mixing parameters: (dM=" << gArgs.deltaM << ",dA=" << gArgs.deltaA << ",dF=" << gArgs.deltaF << ")" << " Output file=" << gArgs.outfileName << "\n";
	}
	mcmcout << "Individuals: " << noIndiv << " Populations: " << noPopln << " Loci: " << noLoci << "\n\n";
	mcmcout << "Locus:(Number of Alleles)\n";
for(unsigned int l=0; l<noLoci; l++)
  {
    string locusName;
    IndivMap::iterator iterLocus = locusIDMap.begin();
    while (iterLocus != locusIDMap.end())
      {
	if(iterLocus->second == l)
	  locusName=iterLocus->first;
	iterLocus++;
      }
    if((l % 10)==0) mcmcout << "\n";
    mcmcout	 << locusName << ":" << noAlleles[l] << " ";
  }
mcmcout << "\n\n";
mcmcout << "\n Population Index -> Population Label:\n\n";

	IndivMap::iterator iterPopln = poplnIDMap.begin();
	while (iterPopln != poplnIDMap.end())
	{
		mcmcout << " " << iterPopln->second << "->" << iterPopln->first;
		iterPopln++;
	}
	mcmcout << "\n\n Migration Rates:\n";
	mcmcout.setf(std::ios::fixed, std::ios::floatfield);
	for (unsigned int l = 0; l < noPopln; l++)
	{
		mcmcout << "\n";
		for(unsigned int k=0; k < noPopln; k++)
		{
			mcmcout << std::setprecision(4) <<  " m[" << l << "][" << k << "]: " << avgMigrationRates[l][k] << "(";
		    mcmcout << std::setprecision(4) << sqrt(varMigrationRates[l][k]) << ")" << flush;
		}
	}

	mcmcout << "\n\n Inbreeding Coefficients:\n";
	mcmcout.setf(std::ios::fixed, std::ios::floatfield);
	for (unsigned int l = 0; l < noPopln; l++)
	{
		string poplnName;
		IndivMap::iterator iterPopln = poplnIDMap.begin();
		while (iterPopln != poplnIDMap.end())
		{
			if(iterPopln->second == l)
				poplnName = iterPopln->first;
			iterPopln++;
		}
		mcmcout << "\n " << poplnName << " Fstat: " << avgFStat[l] << "(" << sqrt(varFStat[l]) << ")";
	}

	mcmcout << "\n\n Allele Frequencies:";
	for (unsigned int l = 0; l < noPopln; l++)
	{
		string poplnName;
		IndivMap::iterator iterPopln = poplnIDMap.begin();
		while (iterPopln != poplnIDMap.end())
		{
			if(iterPopln->second == l)
				poplnName = iterPopln->first;
			iterPopln++;
		}
		mcmcout << "\n\n " << poplnName;
		for (unsigned int k=0; k < noLoci; k++)
		{
			string locusName;
			IndivMap::iterator iterLocus = locusIDMap.begin();
			while (iterLocus != locusIDMap.end())
			{
				if(iterLocus->second == k)
					locusName=iterLocus->first;
				iterLocus++;
			}
			mcmcout << "\n " << locusName << ">>\n" << " ";
			for (unsigned int j=0; j < noAlleles[k]; j++)
			{
				string alleleName;
				IndivMap::iterator iterAllele = alleleIDMap[k].begin();
				while (iterAllele != alleleIDMap[k].end())
				{
					if(iterAllele->second == j)
						alleleName=iterAllele->first;
					iterAllele++;
				}

				mcmcout << std::setprecision(3) <<  alleleName << ":" << avgAlleleFreqs[l][k][j] << "(";
				mcmcout << std::setprecision(3) << sqrt(varAlleleFreqs[l][k][j]) << ") " << flush;
			}
		}
	}

	/* print out individual genotypes and migrant ancestries to BA3indiv.txt */

	if (gArgs.genotypes)
	{
	indivout.setf(std::ios::fixed, std::ios::floatfield);
	for (unsigned int l=0; l<noIndiv; l++)
	{
		string indivName;
		IndivMap::iterator iterIndiv = indivIDMap.begin();
		while (iterIndiv != indivIDMap.end())
		{
			if(iterIndiv->second == l)
				indivName = iterIndiv->first;
			iterIndiv++;
		}
		indivout << "\n\n Individual: " << indivName << " Source Popln: " << sampleIndiv[l].samplePopln << "\n" << " Genotypes>>";
		for (unsigned int j=0; j<noLoci; j++)
		{
			if((j % 10)==0) indivout << "\n";
			string locusName;
			IndivMap::iterator iterLocus = locusIDMap.begin();
			while (iterLocus != locusIDMap.end())
			{
				if(iterLocus->second == j)
					locusName=iterLocus->first;
				iterLocus++;
			}
			indivout << " " << locusName << ":";
			for (unsigned int m=0; m <= 1; m++)
			{
				vector<int>::iterator it1,it2;
				char targ[]={static_cast<char>(l)}, targ2[]={static_cast<char>(j)};
				it1 = search (missingData.begin(), missingData.end(),targ,targ+1);
				it2 = search (sampleIndiv[l].missingGenotypes.begin(), sampleIndiv[l].missingGenotypes.end(),targ2,targ2+1);

				if((it1 != missingData.end())&&(it2 != sampleIndiv[l].missingGenotypes.end()))
				{
					if(m==0)
						indivout << "?/?";
				}
				else
				{
					string alleleName;
					IndivMap::iterator iterAllele = alleleIDMap[j].begin();
					while (iterAllele != alleleIDMap[j].end())
					{
						if(iterAllele->second == sampleIndiv[l].genotype[j][m])
                        {
							if(m==0) { indivout << iterAllele->first << "/"; }
							else { indivout << iterAllele->first; }
                        }
						iterAllele++;
					}
				}
			}
		}
		indivout << "\n Migrant ancestry>>";
		for(int j=0; j<3; j++)
		{
			indivout << "\n";
			for(unsigned int k=0; k<noPopln; k++)
				indivout << std::setprecision(3) << " [" << k << "," << j << "]:" << ancP[l][k][j];
		}


	}

	/* more to go here */

	}


	std::cout << "\n\n" << " MCMC run completed. Output written to " << gArgs.outfileName << "\n\n";
	delete []sampleIndiv;
	sampleIndiv=NULL;
	mcmcout.close();
	return 0;
}

void printBanner(void)
{
	std::cout << "\n\n";
	std::cout << "                    BayesAss Edition " << VERSION << " (BA3)                    \n";
	std::cout << "                        Released: " << RELEASEDATE << "                        \n";
	std::cout << "                            Bruce Rannala                           \n";
	std::cout << "           Department of Evolution and Ecology at UC Davis          \n";
	std::cout << "\n\n";
}

void checkDataSize(void)
{
  std::string aline;
  vector<string> indNames;
  vector<string> popNames;
  vector<string> locusNames;
  map<string,set<string>> alleleLabels;
  
  while(std::getline(mcmcin,aline))
    {
      std::string regex_str = "\\s+"; 
      std::vector<std::string> list1 = split(aline, regex_str);
      unsigned int s1 = list1.size();
      if(s1 != 5)
	{
	  std::cerr << "\nerror the line: \"" << aline << "\" has an incorrect number of entries. quitting...\n";
	  exit(1);
	}
      indNames.push_back(list1[0]);
      popNames.push_back(list1[1]);
      locusNames.push_back(list1[2]);
    }
  mcmcin.clear();
  mcmcin.seekg( 0, std::ios::beg );
  std::set<string> namesUniq;
  std::set<string> popsUniq;
  std::set<string> lociUniq;
  for( unsigned i = 0; i < indNames.size(); ++i ) namesUniq.insert( indNames[i] );
  for( unsigned i = 0; i < popNames.size(); ++i ) popsUniq.insert( popNames[i] );
  for( unsigned i = 0; i < locusNames.size(); ++i ) lociUniq.insert( locusNames[i] );
  if(namesUniq.size() > MAXINDIV)
    {
      std::cerr << "\nerror: number of individuals:" << namesUniq.size() << " exceeds maximum:" << MAXINDIV << " quitting...\n";
      exit(1);
    }
  if(popsUniq.size() > MAXPOPLN)
    {
      std::cerr << "\nerror: number of populations:" << popsUniq.size() << " exceeds maximum:" << MAXPOPLN << " quitting...\n";
      exit(1);
    }
  if(lociUniq.size() > MAXLOCI)
    {
      std::cerr << "\nerror: number of loci:" << lociUniq.size() << " exceeds maximum:" << MAXLOCI << " quitting...\n";
      exit(1);
    }
  while(std::getline(mcmcin,aline))
	{
	  std::string regex_str = "\\s+"; 
	  std::vector<std::string> list1 = split(aline, regex_str);
	  alleleLabels[list1[2]].insert(list1[3]);
	  alleleLabels[list1[2]].insert(list1[4]);
	}
  for (auto all = alleleLabels.begin(); all != alleleLabels.end(); all++) {
    if(all->second.size() > MAXALLELE)
      {
	std::cerr << "\nerror: number of alleles at locus " << all->first << " is:" << all->second.size() << " exceeds maximum:" << MAXALLELE << " quitting...\n";
	exit(1);
      }
  }
  mcmcin.clear();
  mcmcin.seekg( 0, std::ios::beg );
}

void readInputFile(indiv *sampleIndiv, unsigned int &noIndiv, unsigned int &noLoci, unsigned int &noPopln, unsigned int *noAlleles)
{
	struct oneLine {
		std::string indiv;
		std::string samplePop;
		std::string locus;
		std::string allele1;
		std::string allele2;
	} aline;
	aline.indiv=" ";
	aline.samplePop=" ";
	aline.locus=" ";
	aline.allele1=" ";
	aline.allele2=" ";
	unsigned int indivIter=0, popIter=0, locIter=0, alleleIter[MAXLOCI];
	unsigned int currIndivID, currPoplnID, currLocusID, currAllele1, currAllele2;
    std::string inputLine=" ";

	for(int k = 0; k < MAXINDIV-1; k++)
	  for(int l = 1; l < MAXLOCI-1; l++)
	    {
	      sampleIndiv[k].genotype[l][0] = -2;
	      sampleIndiv[k].genotype[l][1] = -2;
	      }

	for(int l = 0; l < MAXLOCI; l++)
		alleleIter[l] = 0;

	while(std::getline(mcmcin,inputLine))
	{

	int numspaces=0;
	char lastChar='a', nextChar;

	// checks each character in the string
	bool firstChar=false;
	for (int ii=0; ii<int(inputLine.length()); ii++)
	{
	  if(!firstChar)
	    {
	      if(!isspace(inputLine.at(ii)))
		firstChar=true;
	    }
	  else
	    {
	      nextChar = inputLine.at(ii); // gets a character
	      if ((isspace(inputLine[ii]))&&!(isspace(inputLine[ii-1])))
		numspaces++;
	    }
	}
	if(!((numspaces==4)||(numspaces==5)))
	  {
	    cerr << "Error: Incorrect number of entries in line>> " << inputLine << " >>of input file " << infileName << "\n\n"; exit(1);
	  }

		if (inputLine.size()>1)
		{
        std::istringstream ss(inputLine);
		ss >> aline.indiv;
		ss >> aline.samplePop;
		ss >> aline.locus;
		ss >> aline.allele1;
		ss >> aline.allele2;
		if (!ss)
		{
			cerr << "input file error: " << inputLine;
			exit(1);
		}
		if (indivIter == 0)
		{
			indivIDMap.insert(std::pair<string, unsigned int>(aline.indiv, indivIter));
			poplnIDMap.insert(std::pair<string, unsigned int>(aline.samplePop, popIter));
			locusIDMap.insert(std::pair<string, unsigned int>(aline.locus, locIter));
			currLocusID = locIter;
			currIndivID = indivIter;
			currPoplnID = popIter;

			if(aline.allele1 != "0")
			{
				alleleIDMap[locIter].insert(std::pair<string, unsigned int>(aline.allele1, alleleIter[locIter]));
				currAllele1 = alleleIter[currLocusID];
				alleleIter[currLocusID]++;
			}
			else currAllele1 = -1;

			if(aline.allele2 != "0")
			{
				IndivMap::iterator iterAllele = alleleIDMap[currLocusID].begin();
				iterAllele = alleleIDMap[currLocusID].find(aline.allele2);
				if (iterAllele != alleleIDMap[currLocusID].end())
					currAllele2 = iterAllele->second;
				else
				{
					alleleIDMap[currLocusID].insert(std::pair<string, unsigned int>(aline.allele2, alleleIter[currLocusID]));
					currAllele2 = alleleIter[currLocusID];
					alleleIter[currLocusID]++;
				}
			}
			else currAllele2 = -1;

			sampleIndiv[currIndivID].samplePopln = currPoplnID;
			indivIter++;
			popIter++;
			locIter++;
		}
		else
		{
			IndivMap::iterator iterIndiv = indivIDMap.begin();
			iterIndiv = indivIDMap.find(aline.indiv);
			if (iterIndiv != indivIDMap.end())
			{
				currIndivID = iterIndiv->second;
				IndivMap::iterator iterLocus = locusIDMap.begin();
				iterLocus = locusIDMap.find(aline.locus);
				if (iterLocus != locusIDMap.end() )
					currLocusID = iterLocus->second;
				else
				{
					locusIDMap.insert(std::pair<string, unsigned int>(aline.locus,locIter));
					currLocusID = locIter;
					locIter++;
				}
				if(aline.allele1 != "0")
				{
					IndivMap::iterator iterAllele = alleleIDMap[currLocusID].begin();
					iterAllele = alleleIDMap[currLocusID].find(aline.allele1);
					if (iterAllele != alleleIDMap[currLocusID].end())
						currAllele1 = iterAllele->second;
					else
					{
						alleleIDMap[currLocusID].insert(std::pair<string, unsigned int>(aline.allele1, alleleIter[currLocusID]));
						currAllele1 = alleleIter[currLocusID];
						alleleIter[currLocusID]++;
					}
				}
				else currAllele1 = -1;
				if(aline.allele2 != "0")
				{
					IndivMap::iterator iterAllele = alleleIDMap[currLocusID].begin();
					iterAllele = alleleIDMap[currLocusID].find(aline.allele2);
					if (iterAllele != alleleIDMap[currLocusID].end())
						currAllele2 = iterAllele->second;
					else
					{
						alleleIDMap[currLocusID].insert(std::pair<string, unsigned int>(aline.allele2, alleleIter[currLocusID]));
						currAllele2 = alleleIter[currLocusID];
						alleleIter[currLocusID]++;
					}
				}
				else currAllele2 = -1;
			}
			else
			{
				indivIDMap.insert(std::pair<string, unsigned int>(aline.indiv,indivIter));
				currIndivID = indivIter;
				sampleIndiv[currIndivID].samplePopln = currPoplnID;
				indivIter++;
				IndivMap::iterator iterPopln = poplnIDMap.begin();
				iterPopln = poplnIDMap.find(aline.samplePop);
				if (iterPopln != poplnIDMap.end() )
					currPoplnID = iterPopln->second;
				else
				{
					poplnIDMap.insert(std::pair<string, unsigned int>(aline.samplePop,popIter));
					currPoplnID = popIter;
					sampleIndiv[currIndivID].samplePopln = currPoplnID;
					popIter++;
				}
				IndivMap::iterator iterLocus = locusIDMap.begin();
				iterLocus = locusIDMap.find(aline.locus);
				if (iterLocus != locusIDMap.end() )
					currLocusID = iterLocus->second;
				else
				{
					locusIDMap.insert(std::pair<string, unsigned int>(aline.locus,locIter));
					currLocusID = locIter;
					locIter++;
				}
				if(aline.allele1 != "0")
				{
					IndivMap::iterator iterAllele = alleleIDMap[currLocusID].begin();
					iterAllele = alleleIDMap[currLocusID].find(aline.allele1);
					if (iterAllele != alleleIDMap[currLocusID].end())
						currAllele1 = iterAllele->second;
					else
					{
						alleleIDMap[currLocusID].insert(std::pair<string, unsigned int>(aline.allele1, alleleIter[currLocusID]));
						currAllele1 = alleleIter[currLocusID];
						alleleIter[currLocusID]++;
					}
				}
				else currAllele1 = -1;
				if(aline.allele2 != "0")
				{
					IndivMap::iterator iterAllele = alleleIDMap[currLocusID].begin();
					iterAllele = alleleIDMap[currLocusID].find(aline.allele2);
					if (iterAllele != alleleIDMap[currLocusID].end())
						currAllele2 = iterAllele->second;
					else
					{
						alleleIDMap[currLocusID].insert(std::pair<string, unsigned int>(aline.allele2, alleleIter[currLocusID]));
						currAllele2 = alleleIter[currLocusID];
						alleleIter[currLocusID]++;
					}
				}
				else currAllele2 = -1;
			}

		}
	}
	sampleIndiv[currIndivID].genotype[currLocusID][0] = currAllele1;
	sampleIndiv[currIndivID].genotype[currLocusID][1] = currAllele2;
	}
	noIndiv = indivIter--;
	noPopln = popIter--;
	noLoci = locIter--;
	for (unsigned int l = 0; l < noLoci; l++)
		noAlleles[l] = alleleIter[l];
	// check that each individual has an entry for every locus
	for(unsigned int k = 0; k < noIndiv; k++)
	  for(unsigned int l = 0; l < noLoci; l++)
	    if((sampleIndiv[k].genotype[l][0]==-2)||(sampleIndiv[k].genotype[l][1]==-2))
	      {
		string locusName;
		IndivMap::iterator iterLocus = locusIDMap.begin();
		while (iterLocus != locusIDMap.end())
		{
			if(iterLocus->second == l)
				locusName=iterLocus->first;
			iterLocus++;
		}
		string indivName;
		IndivMap::iterator iterIndiv = indivIDMap.begin();
		while (iterIndiv != indivIDMap.end())
		{
			if(iterIndiv->second == k)
				indivName=iterIndiv->first;
			iterIndiv++;
		}

		cerr << "Missing genotype entry for individual: " << indivName << " at locus: " <<  locusName  << "!\n\n"; exit(1);

	      }


}

void getEmpiricalAlleleFreqs(double ***alleleFreqs, indiv *sampleIndiv, unsigned int *noAlleles, unsigned int noPopln, unsigned int noLoci, unsigned int noIndiv)
{
	double epsilon=0.0001; // minimum allele frequency in any population
	int indivPerPopln[noPopln];
	long int ***poplnAlleleCounts;
	poplnAlleleCounts = new long int**[noPopln];
	for(unsigned int i = 0; i < noPopln; i++)
	{
		poplnAlleleCounts[i] = new long int*[noLoci];
		for(unsigned int j = 0; j < noLoci; j++)
			poplnAlleleCounts[i][j] = new long int[noAlleles[j]];
	}

	for (unsigned int l=0; l<noPopln; l++) { indivPerPopln[l] = 0; }

	for (unsigned int l=0; l<noPopln; l++)
		for (unsigned int j=0; j<noLoci; j++)
			for (unsigned int k=0; k<noAlleles[j]; k++)
			{
				int alleleCount=0;
				for (unsigned int m=0; m<noIndiv; m++)
				{
					if(sampleIndiv[m].samplePopln == l)
					{
						if (sampleIndiv[m].genotype[j][0] == k)
							alleleCount++;
						if (sampleIndiv[m].genotype[j][1] == k)
						{ alleleCount++; if(j==0) indivPerPopln[l]++; }
					}
				}
				poplnAlleleCounts[l][j][k] = alleleCount;
			}
	for (unsigned int l=0; l<noPopln; l++)
		for (unsigned int j=0; j<noLoci; j++)
			for (unsigned int k=0; k<noAlleles[j]; k++)
				alleleFreqs[l][j][k] = (poplnAlleleCounts[l][j][k]/(2.0*indivPerPopln[l])+epsilon)/(1.0+noAlleles[j]*epsilon);

	if(gArgs.debug)
	for (unsigned int l=0; l<noPopln; l++)
		for (unsigned int j=0; j<noLoci; j++)
		{
			double sum=0;
			for (unsigned int k=0; k<noAlleles[j]; k++)
				sum+=alleleFreqs[l][j][k];
			cout << "\nSum of initial allele freqs at pop:" << l << " locus:" << j << " = " << sum;
		}

	delete []poplnAlleleCounts;
	poplnAlleleCounts=NULL;
}









void fillMigrantCounts(indiv *sampleIndiv, long int ***migrantCounts, unsigned int noIndiv, unsigned int noPopln)
{
	for(unsigned int i = 0; i < noPopln; i++)
		for(unsigned int j = 0; j < noPopln; j++)
			for(int k = 0; k < 3; k++)
				migrantCounts[i][j][k]=0;
	for(unsigned int i = 0; i < noIndiv; i++)
		migrantCounts[sampleIndiv[i].samplePopln][sampleIndiv[i].migrantPopln][sampleIndiv[i].migrantAge] += 1;
}

double migCountLogProb(long int ***migrantCounts, double **migrationRates, unsigned int noPopln)
{
	double logPr=0.0;
	for(unsigned int l = 0; l < noPopln; l++)
	{
		for(unsigned int k = 0; k < noPopln; k++)
			if(l != k)
			{
				logPr+=migrantCounts[l][k][1]*log(migrationRates[l][k])-gsl_sf_lnfact(migrantCounts[l][k][1]);
				logPr+=migrantCounts[l][k][2]*log(2.0*migrationRates[l][k])-gsl_sf_lnfact(migrantCounts[l][k][2]);
			}
		logPr+=migrantCounts[l][l][0]*log(1.0-3*migrationRates[l][noPopln]) - gsl_sf_lnfact(migrantCounts[l][l][0]);
	}
	return(logPr);
}


double logLik(indiv Indiv, double ***alleleFreqs, double *FStat, unsigned int noLoci)
{
	double logPr=0;
	if(Indiv.migrantAge == 0)
	{
		for(unsigned int i = 0; i < noLoci; i++)
		{
			if(Indiv.genotype[i][0] == Indiv.genotype[i][1])
				logPr+=log((1-FStat[Indiv.samplePopln])*alleleFreqs[Indiv.samplePopln][i][Indiv.genotype[i][0]]*
				alleleFreqs[Indiv.samplePopln][i][Indiv.genotype[i][0]]+FStat[Indiv.samplePopln]*
				alleleFreqs[Indiv.samplePopln][i][Indiv.genotype[i][0]]);
			else
				logPr+=log(2.0*(1-FStat[Indiv.samplePopln])*alleleFreqs[Indiv.samplePopln][i][Indiv.genotype[i][0]]*
				alleleFreqs[Indiv.samplePopln][i][Indiv.genotype[i][1]]);
		}
	}
	else if(Indiv.migrantAge == 1)
	{
		for(unsigned int i = 0; i < noLoci; i++)
		{
			if(Indiv.genotype[i][0] == Indiv.genotype[i][1])
				logPr+=log((1-FStat[Indiv.migrantPopln])*alleleFreqs[Indiv.migrantPopln][i][Indiv.genotype[i][0]]*
						   alleleFreqs[Indiv.migrantPopln][i][Indiv.genotype[i][0]]+FStat[Indiv.migrantPopln]*
						   alleleFreqs[Indiv.migrantPopln][i][Indiv.genotype[i][0]]);
			else
				logPr+=log(2.0*(1-FStat[Indiv.migrantPopln])*alleleFreqs[Indiv.migrantPopln][i][Indiv.genotype[i][0]]*
						   alleleFreqs[Indiv.migrantPopln][i][Indiv.genotype[i][1]]);
		}
	}
	else
	{
		for(unsigned int i = 0; i < noLoci; i++)
		{
			if(Indiv.genotype[i][0] == Indiv.genotype[i][1])
				logPr += log(alleleFreqs[Indiv.migrantPopln][i][Indiv.genotype[i][0]]*alleleFreqs[Indiv.samplePopln][i][Indiv.genotype[i][1]]);
			else
				logPr += log(alleleFreqs[Indiv.migrantPopln][i][Indiv.genotype[i][0]]*alleleFreqs[Indiv.samplePopln][i][Indiv.genotype[i][1]]+
				alleleFreqs[Indiv.migrantPopln][i][Indiv.genotype[i][1]]*alleleFreqs[Indiv.samplePopln][i][Indiv.genotype[i][0]]);
		}
	}
	return logPr;
}

double oneLocusLogLik(indiv Indiv, double ***alleleFreqs, double *FStat, int chosenLocus)
{
	double logPr=0;
	if(Indiv.migrantAge == 0)
	{
			if(Indiv.genotype[chosenLocus][0] == Indiv.genotype[chosenLocus][1])
				logPr+=log((1-FStat[Indiv.samplePopln])*alleleFreqs[Indiv.samplePopln][chosenLocus][Indiv.genotype[chosenLocus][0]]*
						   alleleFreqs[Indiv.samplePopln][chosenLocus][Indiv.genotype[chosenLocus][0]]+FStat[Indiv.samplePopln]*
						   alleleFreqs[Indiv.samplePopln][chosenLocus][Indiv.genotype[chosenLocus][0]]);
			else
				logPr+=log(2.0*(1-FStat[Indiv.samplePopln])*alleleFreqs[Indiv.samplePopln][chosenLocus][Indiv.genotype[chosenLocus][0]]*
						   alleleFreqs[Indiv.samplePopln][chosenLocus][Indiv.genotype[chosenLocus][1]]);
	}
	else if(Indiv.migrantAge == 1)
	{
		if(Indiv.genotype[chosenLocus][0] == Indiv.genotype[chosenLocus][1])
			logPr+=log((1-FStat[Indiv.migrantPopln])*alleleFreqs[Indiv.migrantPopln][chosenLocus][Indiv.genotype[chosenLocus][0]]*
						   alleleFreqs[Indiv.migrantPopln][chosenLocus][Indiv.genotype[chosenLocus][0]]+FStat[Indiv.migrantPopln]*
						   alleleFreqs[Indiv.migrantPopln][chosenLocus][Indiv.genotype[chosenLocus][0]]);
		else
			logPr+=log(2.0*(1-FStat[Indiv.migrantPopln])*alleleFreqs[Indiv.migrantPopln][chosenLocus][Indiv.genotype[chosenLocus][0]]*
						   alleleFreqs[Indiv.migrantPopln][chosenLocus][Indiv.genotype[chosenLocus][1]]);
	}
	else
	{
		if(Indiv.genotype[chosenLocus][0] == Indiv.genotype[chosenLocus][1])
			logPr += log(alleleFreqs[Indiv.migrantPopln][chosenLocus][Indiv.genotype[chosenLocus][0]]*alleleFreqs[Indiv.samplePopln][chosenLocus][Indiv.genotype[chosenLocus][1]]);
		else
			logPr += log(alleleFreqs[Indiv.migrantPopln][chosenLocus][Indiv.genotype[chosenLocus][0]]*alleleFreqs[Indiv.samplePopln][chosenLocus][Indiv.genotype[chosenLocus][1]]+
							 alleleFreqs[Indiv.migrantPopln][chosenLocus][Indiv.genotype[chosenLocus][1]]*alleleFreqs[Indiv.samplePopln][chosenLocus][Indiv.genotype[chosenLocus][0]]);
	}
	return logPr;
}

void proposeMigrantAncDrop(unsigned int &migrantPopln, unsigned int &migrantAge, unsigned int samplePopln, int noPopln, long int ***migrantCounts)
{

	bool emptyMigrantCounts=true;
	unsigned int proposedPopln=0, proposedAge=0;
	while(emptyMigrantCounts)
	{
		int proposedIndex = gsl_rng_uniform_int(r,2*noPopln-1);
		if(proposedIndex==0)
		{ proposedPopln = samplePopln; proposedAge=0; }
		else
		{
			int j=1; proposedPopln=0; proposedAge=1;
			if (proposedPopln == samplePopln)
				proposedPopln += 1;
			while (j < proposedIndex)
			{
				if (proposedAge == 2)
				{
					proposedPopln += 1;
					proposedAge = 1;
				}
				else
				{
					proposedAge = 2;
				}
				if (proposedPopln == samplePopln)	proposedPopln += 1;
				j += 1;
			}
		}
			if (migrantCounts[samplePopln][proposedPopln][proposedAge] > 0)
				emptyMigrantCounts = false;
	}
	migrantPopln = proposedPopln;
	migrantAge = proposedAge;
}


void proposeMigrantAncAdd(unsigned int &migrantPopAdd, unsigned int &migrantAgeAdd, unsigned int migrantPopDrop, unsigned int migrantAgeDrop,
						  unsigned int samplePopln, int noPopln)
{

	bool IsDroppedMigrant=true;
	unsigned int proposedPopln=0, proposedAge=0;
	while(IsDroppedMigrant)
	{
		int proposedIndex = gsl_rng_uniform_int(r,2*noPopln-1);
		if(proposedIndex==0)
		{ proposedPopln = samplePopln; proposedAge=0; }
		else
		{
			int j=1; proposedPopln=0; proposedAge=1;
			if (proposedPopln == samplePopln)
				proposedPopln += 1;
			while (j < proposedIndex)
			{
			/*	if (proposedPopln == samplePopln)
					proposedPopln += 1; */
				if (proposedAge == 2)
				{
					proposedPopln += 1;
					proposedAge = 1;
				}
				else
				{
					proposedAge = 2;
				}
				if(proposedPopln==samplePopln) proposedPopln += 1;
				j += 1;
			}
		}
		if ((proposedAge!=migrantAgeDrop)||(proposedPopln!=migrantPopDrop))
			IsDroppedMigrant = false;
	}
	migrantPopAdd = proposedPopln;
	migrantAgeAdd = proposedAge;
	if((migrantPopAdd==samplePopln)&&(migrantAgeAdd!=0)) cerr << "\n proposing own popln as migrant popln!\n";
}
