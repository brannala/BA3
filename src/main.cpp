
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
#include <htslib/vcf.h>
#include <htslib/synced_bcf_reader.h>
#include <htslib/hts_log.h>

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
	int autotune;
	int useVCF;              // Use VCF input format
	char vcfFileName[256];   // VCF file path
	char metaFileName[256];  // Metadata file path (INDIV POPLN)
	int usingFreqFile;       // Output allele frequencies to separate file
	char freqFileName[256];  // Allele frequencies file path
} gArgs;

static const char *optString = "s:i:n:b:o:m:a:f:V:M:F:vugtdphTN?";

static const struct option longOpts[] = {
	{ "seed", required_argument, NULL, 's' },
	{ "iterations", required_argument, NULL, 'i' },
	{ "sampling", required_argument, NULL, 'n' },
	{ "burnin", required_argument, NULL, 'b' },
	{ "output", required_argument, NULL, 'o' },
	{ "deltaM", required_argument, NULL, 'm' },
	{ "deltaA", required_argument, NULL, 'a' },
	{ "deltaF", required_argument, NULL, 'f' },
	{ "vcf", required_argument, NULL, 'V' },
	{ "meta", required_argument, NULL, 'M' },
	{ "freqfile", required_argument, NULL, 'F' },
	{ "verbose", no_argument, NULL, 'v' },
	{ "settings", no_argument, NULL, 'u' },
	{ "genotypes", no_argument, NULL, 'g' },
	{ "debug", no_argument, NULL, 'd' },
	{ "help", no_argument, NULL, 'h' },
	{ "trace", no_argument, NULL, 't' },
	{ "nolikelihood", no_argument, NULL, 'p'},
	{ "autotune", no_argument, NULL, 'T'},
	{ "noautotune", no_argument, NULL, 'N'},
	{ NULL, no_argument, NULL, 0 }
};

// global variables

gsl_rng * r;
char infileName[100];
std::ofstream mcmcout;
std::ofstream tracefile;
std::ofstream indivout;
std::ofstream freqout;
std::ifstream mcmcin;

typedef std::map<string, unsigned int> IndivMap;

// map between indivID and name in input file
IndivMap indivIDMap;

// map between poplnID and popln name in input file
IndivMap poplnIDMap;

// map between locusID and locus name in input file
IndivMap locusIDMap;

// map between alleleID and allele name at each locus (dynamically allocated)
IndivMap *alleleIDMap = nullptr;

// Global dataset dimensions (set after reading input file)
unsigned int gNoLoci = 0;

// Forward declarations for progress bar and screen output functions
std::string formatTime(double seconds);
std::string formatIterCount(unsigned long int n);
void printProgress(unsigned long int current, unsigned long int total,
                   double elapsedSecs, bool inBurnin);
unsigned int gNoIndiv = 0;
unsigned int gMaxAlleles = 0;

// Memory allocation helpers for genotypes
void allocateGenotypes(indiv *individuals, unsigned int noIndiv, unsigned int noLoci) {
    for (unsigned int i = 0; i < noIndiv; i++) {
        individuals[i].genotype = new GenotypeType[noLoci][2];
        // Initialize to -2 (unset/missing marker)
        for (unsigned int j = 0; j < noLoci; j++) {
            individuals[i].genotype[j][0] = -2;
            individuals[i].genotype[j][1] = -2;
        }
    }
}

void freeGenotypes(indiv *individuals, unsigned int noIndiv) {
    for (unsigned int i = 0; i < noIndiv; i++) {
        if (individuals[i].genotype != nullptr) {
            delete[] individuals[i].genotype;
            individuals[i].genotype = nullptr;
        }
    }
}

// debugging variables
bool NOMIGRATEMCMC=false;
bool NOANCMCMC=false;
bool NOALLELEMCMC=false;
bool NOFSTATMCMC=false;
bool NOMISSINGDATA=false;
bool NOLIKELIHOOD=false;

// autotune constants
const double AUTOTUNE_TARGET_RATE = 0.30;      // target acceptance rate
const double AUTOTUNE_LOWER_BOUND = 0.20;      // adjust up if below this
const double AUTOTUNE_UPPER_BOUND = 0.40;      // adjust down if above this
const int AUTOTUNE_INTERVAL = 100;             // check every N iterations during burn-in
const double AUTOTUNE_ADJUST_FACTOR = 1.1;     // multiply/divide delta by this factor
const double AUTOTUNE_DELTA_MIN = 0.001;       // minimum delta value
const double AUTOTUNE_DELTA_MAX = 0.99;        // maximum delta value

// Fast whitespace-based string splitting (replaces slow regex version)
std::vector<std::string> split(const std::string& str,
                               const std::string& regex_str) {
    std::vector<std::string> tokens;
    std::istringstream iss(str);
    std::string token;
    while (iss >> token) {
        tokens.push_back(token);
    }
    return tokens;
}

//=============================================================================
// VCF Input Functions
//=============================================================================

// Read metadata file mapping individuals to populations
// Format: INDIV_ID POPLN_ID (whitespace separated, one per line)
void readMetadataFile(const char *metaFileName,
                      std::map<std::string, std::string> &indivToPopln) {
    std::ifstream metaFile(metaFileName);
    if (!metaFile) {
        std::cerr << "\nerror: cannot open metadata file: " << metaFileName << "\n";
        exit(1);
    }

    std::string line;
    int lineNum = 0;
    while (std::getline(metaFile, line)) {
        lineNum++;
        // Skip empty lines and comments
        if (line.empty() || line[0] == '#') continue;

        std::istringstream iss(line);
        std::string indivID, poplnID;
        if (!(iss >> indivID >> poplnID)) {
            std::cerr << "\nerror: invalid format in metadata file at line " << lineNum
                      << ": expected 'INDIV POPLN'\n";
            exit(1);
        }
        indivToPopln[indivID] = poplnID;
    }
    metaFile.close();

    if (indivToPopln.empty()) {
        std::cerr << "\nerror: no individuals found in metadata file: " << metaFileName << "\n";
        exit(1);
    }
}

// Check VCF file dimensions and validate against metadata
void checkVCFDataSize(const char *vcfFileName,
                      const std::map<std::string, std::string> &indivToPopln,
                      unsigned int &outNoIndiv, unsigned int &outNoLoci,
                      unsigned int &outNoPopln, unsigned int &outMaxAlleles) {

    // Suppress htslib warnings (they go to stderr by default)
    hts_set_log_level(HTS_LOG_ERROR);

    // Open VCF file
    htsFile *vcfFile = bcf_open(vcfFileName, "r");
    if (!vcfFile) {
        std::cerr << "\nerror: cannot open VCF file: " << vcfFileName << "\n";
        exit(1);
    }

    bcf_hdr_t *hdr = bcf_hdr_read(vcfFile);
    if (!hdr) {
        std::cerr << "\nerror: cannot read VCF header from: " << vcfFileName << "\n";
        bcf_close(vcfFile);
        exit(1);
    }

    // Get sample names from VCF and match with metadata
    int nSamples = bcf_hdr_nsamples(hdr);
    std::set<std::string> poplnSet;
    std::vector<std::string> validSamples;

    for (int i = 0; i < nSamples; i++) {
        std::string sampleName = hdr->samples[i];
        auto it = indivToPopln.find(sampleName);
        if (it != indivToPopln.end()) {
            validSamples.push_back(sampleName);
            poplnSet.insert(it->second);
        }
    }

    if (validSamples.empty()) {
        std::cerr << "\nerror: no samples in VCF match metadata file\n";
        std::cerr << "VCF samples: ";
        for (int i = 0; i < std::min(5, nSamples); i++) {
            std::cerr << hdr->samples[i] << " ";
        }
        if (nSamples > 5) std::cerr << "...";
        std::cerr << "\n";
        bcf_hdr_destroy(hdr);
        bcf_close(vcfFile);
        exit(1);
    }

    // Count variants and find max alleles
    bcf1_t *rec = bcf_init();
    unsigned int nLoci = 0;
    unsigned int maxAlleles = 0;

    while (bcf_read(vcfFile, hdr, rec) == 0) {
        bcf_unpack(rec, BCF_UN_ALL);
        nLoci++;
        if ((unsigned int)rec->n_allele > maxAlleles) {
            maxAlleles = rec->n_allele;
        }
    }

    bcf_destroy(rec);
    bcf_hdr_destroy(hdr);
    bcf_close(vcfFile);

    // Check limits
    if (validSamples.size() > MAXINDIV) {
        std::cerr << "\nerror: number of individuals: " << validSamples.size()
                  << " exceeds maximum: " << MAXINDIV << "\n";
        exit(1);
    }
    if (poplnSet.size() > MAXPOPLN) {
        std::cerr << "\nerror: number of populations: " << poplnSet.size()
                  << " exceeds maximum: " << MAXPOPLN << "\n";
        exit(1);
    }
    if (nLoci > MAXLOCI) {
        std::cerr << "\nerror: number of loci: " << nLoci
                  << " exceeds maximum: " << MAXLOCI << "\n";
        exit(1);
    }
    if (maxAlleles > MAXALLELE) {
        std::cerr << "\nerror: maximum alleles: " << maxAlleles
                  << " exceeds maximum: " << MAXALLELE << "\n";
        exit(1);
    }

    outNoIndiv = validSamples.size();
    outNoLoci = nLoci;
    outNoPopln = poplnSet.size();
    outMaxAlleles = maxAlleles;

    // Set global dimensions
    gNoIndiv = outNoIndiv;
    gNoLoci = outNoLoci;
    gMaxAlleles = outMaxAlleles;
}

// Read VCF file and populate genotype data
void readVCFFile(const char *vcfFileName,
                 const std::map<std::string, std::string> &indivToPopln,
                 indiv *sampleIndiv,
                 unsigned int &noIndiv, unsigned int &noLoci,
                 unsigned int &noPopln, unsigned int *noAlleles) {

    // Suppress htslib warnings (they go to stderr by default)
    hts_set_log_level(HTS_LOG_ERROR);

    // Open VCF file
    htsFile *vcfFile = bcf_open(vcfFileName, "r");
    if (!vcfFile) {
        std::cerr << "\nerror: cannot open VCF file: " << vcfFileName << "\n";
        exit(1);
    }

    bcf_hdr_t *hdr = bcf_hdr_read(vcfFile);
    if (!hdr) {
        std::cerr << "\nerror: cannot read VCF header\n";
        bcf_close(vcfFile);
        exit(1);
    }

    // Build sample index mapping (VCF sample index -> our individual index)
    int nSamples = bcf_hdr_nsamples(hdr);
    std::vector<int> sampleToIndiv(nSamples, -1);  // -1 means skip this sample
    std::map<std::string, unsigned int> poplnNameToID;
    unsigned int indivIdx = 0;
    unsigned int poplnIdx = 0;

    for (int i = 0; i < nSamples; i++) {
        std::string sampleName = hdr->samples[i];
        auto it = indivToPopln.find(sampleName);
        if (it != indivToPopln.end()) {
            sampleToIndiv[i] = indivIdx;

            // Add to indivIDMap
            indivIDMap[sampleName] = indivIdx;

            // Get or create population ID
            std::string poplnName = it->second;
            auto pIt = poplnNameToID.find(poplnName);
            unsigned int pID;
            if (pIt == poplnNameToID.end()) {
                pID = poplnIdx++;
                poplnNameToID[poplnName] = pID;
                poplnIDMap[poplnName] = pID;
            } else {
                pID = pIt->second;
            }

            sampleIndiv[indivIdx].samplePopln = pID;
            sampleIndiv[indivIdx].migrantPopln = pID;
            sampleIndiv[indivIdx].migrantAge = 0;
            indivIdx++;
        }
    }

    noIndiv = indivIdx;
    noPopln = poplnIdx;

    // Read variants
    bcf1_t *rec = bcf_init();
    unsigned int locusIdx = 0;
    int32_t *gt = NULL;
    int ngt = 0;

    while (bcf_read(vcfFile, hdr, rec) == 0) {
        bcf_unpack(rec, BCF_UN_ALL);

        // Create locus name from CHROM:POS or ID
        std::string locusName;
        if (rec->d.id && strcmp(rec->d.id, ".") != 0) {
            locusName = rec->d.id;
        } else {
            locusName = std::string(bcf_hdr_id2name(hdr, rec->rid)) + ":" + std::to_string(rec->pos + 1);
        }
        locusIDMap[locusName] = locusIdx;

        // Number of alleles at this locus (REF + ALTs)
        noAlleles[locusIdx] = rec->n_allele;

        // Initialize alleleIDMap for this locus
        for (int a = 0; a < rec->n_allele; a++) {
            std::string alleleName = rec->d.allele[a];
            alleleIDMap[locusIdx][alleleName] = a;
        }

        // Get genotypes
        int ngt_ret = bcf_get_genotypes(hdr, rec, &gt, &ngt);
        if (ngt_ret <= 0) {
            // No genotype data for this variant - mark all as missing
            for (unsigned int i = 0; i < noIndiv; i++) {
                sampleIndiv[i].genotype[locusIdx][0] = -1;
                sampleIndiv[i].genotype[locusIdx][1] = -1;
            }
        } else {
            int ploidy = ngt_ret / nSamples;

            for (int s = 0; s < nSamples; s++) {
                int indIdx = sampleToIndiv[s];
                if (indIdx < 0) continue;  // Skip samples not in metadata

                int32_t *ptr = gt + s * ploidy;

                // Handle diploid genotypes
                if (ploidy >= 2) {
                    if (bcf_gt_is_missing(ptr[0]) || bcf_gt_is_missing(ptr[1])) {
                        // Missing genotype
                        sampleIndiv[indIdx].genotype[locusIdx][0] = -1;
                        sampleIndiv[indIdx].genotype[locusIdx][1] = -1;
                    } else {
                        // Allele indices (0=REF, 1=first ALT, etc.)
                        int a0 = bcf_gt_allele(ptr[0]);
                        int a1 = bcf_gt_allele(ptr[1]);
                        sampleIndiv[indIdx].genotype[locusIdx][0] = a0;
                        sampleIndiv[indIdx].genotype[locusIdx][1] = a1;
                    }
                } else {
                    // Haploid - treat as homozygous
                    if (bcf_gt_is_missing(ptr[0])) {
                        sampleIndiv[indIdx].genotype[locusIdx][0] = -1;
                        sampleIndiv[indIdx].genotype[locusIdx][1] = -1;
                    } else {
                        int a0 = bcf_gt_allele(ptr[0]);
                        sampleIndiv[indIdx].genotype[locusIdx][0] = a0;
                        sampleIndiv[indIdx].genotype[locusIdx][1] = a0;
                    }
                }
            }
        }

        locusIdx++;
    }

    noLoci = locusIdx;

    if (gt) free(gt);
    bcf_destroy(rec);
    bcf_hdr_destroy(hdr);
    bcf_close(vcfFile);
}

//=============================================================================

int main( int argc, char *argv[] )
{

	unsigned int noIndiv=0;
	unsigned int noLoci=0;
	unsigned int noPopln=0;
	unsigned int noMissingGenotypes=0;
	unsigned int maxAlleles=0;
	unsigned int *noAlleles = nullptr;  // Dynamically allocated after checkDataSize

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
	gArgs.autotune = 1;  // autotune enabled by default
	gArgs.usingOutfile = 1;
	gArgs.useVCF = 0;
	gArgs.vcfFileName[0] = '\0';
	gArgs.metaFileName[0] = '\0';
	gArgs.usingFreqFile = 0;
	gArgs.freqFileName[0] = '\0';

	strncpy(gArgs.outfileName, "BA3out.txt", sizeof(gArgs.outfileName) - 1);
	gArgs.outfileName[sizeof(gArgs.outfileName) - 1] = '\0';

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
				{ strncpy(gArgs.outfileName, optarg, sizeof(gArgs.outfileName) - 1);
				  gArgs.outfileName[sizeof(gArgs.outfileName) - 1] = '\0';
				  gArgs.usingOutfile = 1; }
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

			case 'T':
				gArgs.autotune = 1;
				break;

			case 'N':
				gArgs.autotune = 0;
				break;

			case 'V':
				strncpy(gArgs.vcfFileName, optarg, sizeof(gArgs.vcfFileName) - 1);
				gArgs.vcfFileName[sizeof(gArgs.vcfFileName) - 1] = '\0';
				gArgs.useVCF = 1;
				break;

			case 'M':
				strncpy(gArgs.metaFileName, optarg, sizeof(gArgs.metaFileName) - 1);
				gArgs.metaFileName[sizeof(gArgs.metaFileName) - 1] = '\0';
				break;

			case 'F':
				strncpy(gArgs.freqFileName, optarg, sizeof(gArgs.freqFileName) - 1);
				gArgs.freqFileName[sizeof(gArgs.freqFileName) - 1] = '\0';
				gArgs.usingFreqFile = 1;
				break;

			case 'd':
				gArgs.debug++;
				break;

			case 'h':
				std::cout << "usage: BA3 [-sinbomaf] [-vhugptTN] file ..." << "\n";
				std::cout << "       BA3 -V vcf_file -M meta_file [-sinbomaf] [-vhugptTN]" << "\n";
				std::cout << "    mcmc analysis of recent migration rates\n\n";
				std::cout << "  Input options:\n";
				std::cout << "    file              BA3 format input file (default)\n";
				std::cout << "    -V, --vcf FILE    VCF/BCF input file\n";
				std::cout << "    -M, --meta FILE   Metadata file with INDIV POPLN columns (required with -V)\n\n";
				std::cout << "  MCMC options:\n";
				std::cout << "    -T, --autotune    auto-tune mixing parameters during burn-in (default)\n";
				std::cout << "    -N, --noautotune  disable auto-tuning of mixing parameters\n";
				exit(0);

			case '?':
				std::cout << "usage: BA3 [-sinbomaf] [-vhugptTN] file ..." << "\n" << "    mcmc analysis of recent migration rates\n";
				exit(0);

			default:
				/* You won't actually get here. */
				break;
		}

		opt = getopt_long( argc, argv, optString, longOpts, &longIndex );
	}

	/* get input file name or validate VCF options */
	std::map<std::string, std::string> indivToPopln;  // For VCF mode

	if (gArgs.useVCF) {
		// VCF mode - validate required options
		if (gArgs.metaFileName[0] == '\0') {
			std::cerr << "error: metadata file (-M) required when using VCF input (-V)\n";
			exit(1);
		}
		strncpy(infileName, gArgs.vcfFileName, sizeof(infileName) - 1);
		infileName[sizeof(infileName) - 1] = '\0';
	} else {
		// Standard BA3 format
		if (optind < argc) {
			strncpy(infileName, argv[optind], sizeof(infileName) - 1);
			infileName[sizeof(infileName) - 1] = '\0';
		} else {
			std::cout << "error: no input file specified...\n";
			std::cout << "usage: BA3 [-sinbomaf] [-vhugptTN] file ...\n";
			std::cout << "   or: BA3 -V vcf_file -M meta_file [options]\n";
			exit(1);
		}
	}

	printBanner();

	unsigned long int mciter=gArgs.mciter;
	unsigned long int aSeed=gArgs.seed;

	if (gArgs.usingOutfile)
	{
		mcmcout.open(gArgs.outfileName, std::ios::out);
	}

	if (gArgs.usingFreqFile)
	{
		freqout.open(gArgs.freqFileName, std::ios::out);
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

	// Initialize GSL random number generator
	r = gsl_rng_alloc (gsl_rng_taus);
	gsl_rng_set(r,aSeed);

	// check that proposal step parameters are sane (e.g., <0 and <=1)
	if((gArgs.deltaM<=0)||(gArgs.deltaM>1))
	  { std::cerr << "\nerror: deltaM=" << gArgs.deltaM << " not in interval (0,1]. quitting...\n";
	    exit(1); }
	if((gArgs.deltaA<=0)||(gArgs.deltaA>1))
	  { std::cerr << "\nerror: deltaA=" << gArgs.deltaA << " not in interval (0,1]. quitting...\n";
	    exit(1); }
	if((gArgs.deltaF<=0)||(gArgs.deltaF>1))
	  { std::cerr << "\nerror: deltaF=" << gArgs.deltaF << " not in interval (0,1]. quitting...\n";
	    exit(1); }

	vector<int> missingData;
	indiv *sampleIndiv = nullptr;

	// Read input data - either VCF or standard BA3 format
	if (gArgs.useVCF) {
		// VCF input mode
		readMetadataFile(gArgs.metaFileName, indivToPopln);

		// First pass: determine dataset dimensions
		checkVCFDataSize(gArgs.vcfFileName, indivToPopln, noIndiv, noLoci, noPopln, maxAlleles);

		// Allocate memory based on actual dataset size
		noAlleles = new unsigned int[noLoci];
		for (unsigned int l = 0; l < noLoci; l++) noAlleles[l] = 0;

		alleleIDMap = new IndivMap[noLoci];

		sampleIndiv = new indiv[noIndiv];
		if(sampleIndiv==NULL) { cerr << "ran out of memory"; exit(1); }

		// Allocate genotypes for each individual
		allocateGenotypes(sampleIndiv, noIndiv, noLoci);

		// Second pass: read the VCF data
		readVCFFile(gArgs.vcfFileName, indivToPopln, sampleIndiv, noIndiv, noLoci, noPopln, noAlleles);

		// Update maxAlleles from actual data
		for (unsigned int l = 0; l < noLoci; l++) {
			if (noAlleles[l] > maxAlleles) maxAlleles = noAlleles[l];
		}
		gMaxAlleles = maxAlleles;

		// Continue with common processing (goto is ugly but avoids major restructuring)
		goto common_processing;
	}

	// Standard BA3 format input
	mcmcin.open(infileName, std::ios::in);
	if (!mcmcin)
	{
		std::cerr << "\nerror: cannot open file: "  << infileName << " quitting...\n";
		exit(1);
	}

	// First pass: determine dataset dimensions
	checkDataSize(noIndiv, noLoci, noPopln, maxAlleles);

	// Allocate memory based on actual dataset size
	noAlleles = new unsigned int[noLoci];
	for (unsigned int l = 0; l < noLoci; l++) noAlleles[l] = 0;

	alleleIDMap = new IndivMap[noLoci];

	sampleIndiv = new indiv[noIndiv];
	if(sampleIndiv==NULL) { cerr << "ran out of memory"; exit(1); }

	// Allocate genotypes for each individual
	allocateGenotypes(sampleIndiv, noIndiv, noLoci);

	// Second pass: read the data
	readInputFile(sampleIndiv, noIndiv, noLoci, noPopln, noAlleles);

	// Update maxAlleles from actual data (might differ slightly due to missing data)
	for (unsigned int l = 0; l < noLoci; l++) {
		if (noAlleles[l] > maxAlleles) maxAlleles = noAlleles[l];
	}
	gMaxAlleles = maxAlleles;

common_processing:

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

	// Check for minimum number of populations
	if (noPopln < 2) {
		std::cerr << "\nerror: At least 2 populations are required to estimate migration rates.\n";
		std::cerr << "       Found only " << noPopln << " population(s) in input data.\n\n";
		exit(1);
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
			alleleFreqs[i][j] = new double[noAlleles[j] > 0 ? noAlleles[j] : 1];
	}

	// Pre-computed log allele frequencies for optimization
	double ***logAlleleFreqs;
	logAlleleFreqs = new double**[noPopln];
	for(unsigned int i = 0; i < noPopln; i++)
	{
		logAlleleFreqs[i] = new double*[noLoci];
		for(unsigned int j = 0; j < noLoci; j++)
			logAlleleFreqs[i][j] = new double[noAlleles[j] > 0 ? noAlleles[j] : 1];
	}

	double ***avgAlleleFreqs;
	avgAlleleFreqs = new double**[noPopln];
	for(unsigned int i = 0; i < noPopln; i++)
	{
		avgAlleleFreqs[i] = new double*[noLoci];
		for(unsigned int j = 0; j < noLoci; j++)
			avgAlleleFreqs[i][j] = new double[noAlleles[j] > 0 ? noAlleles[j] : 1];
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
			varAlleleFreqs[i][j] = new double[noAlleles[j] > 0 ? noAlleles[j] : 1];
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
	// Pre-computed log values for FStat optimization
	double *logFStat = new double[noPopln];
	double *log1MinusFStat = new double[noPopln];
	for(unsigned int i=0; i< noPopln; i++)
	{
		FStat[i]=0.0;
		avgFStat[i]=0.0;
		varFStat[i]=0.0;
		logFStat[i] = log(1e-10);  // log(0) approximation for initial FStat=0
		log1MinusFStat[i] = log(1.0);  // log(1-0) = 0
	}

	// Allocate and initialize Savage-Dickey statistics for migration rate testing
	SavageDickeyStats **sdStats = new SavageDickeyStats*[noPopln];
	for (unsigned int i = 0; i < noPopln; i++) {
		sdStats[i] = new SavageDickeyStats[noPopln];
	}
	initSavageDickeyStats(sdStats, noPopln);
	const double SD_BANDWIDTH = 0.02;  // Reference bandwidth for kernel density estimation
	const double PRIOR_DENSITY_AT_ZERO = 3.0;  // Prior is Uniform(0, 1/3), so density = 3

	for(unsigned int l=0; l<noPopln; l++)
	{
		for(unsigned int k=0; k<=noPopln; k++)
			if((l!=k)&&(k!=noPopln))
				migrationRates[l][k] = (1.0/3.0)*(1.0/noPopln);
		migrationRates[l][l] = 1.0-(1.0/3.0)*((noPopln-1.0)/noPopln);
		migrationRates[l][noPopln] = (1.0/3.0)*((noPopln-1.0)/noPopln);
	}

	getEmpiricalAlleleFreqs(alleleFreqs,sampleIndiv,noAlleles,noPopln,noLoci,noIndiv);

	// Initialize log allele frequencies
	updateAllLogAlleleFreqs(logAlleleFreqs, alleleFreqs, noPopln, noLoci, noAlleles);

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
		sampleIndiv[i].logL = logLik(sampleIndiv[i],alleleFreqs,logAlleleFreqs,FStat,log1MinusFStat,noLoci);

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

	// Autotune tracking variables (acceptance counts during tuning window)
	int tuneWindowMigAccept = 0, tuneWindowMigTotal = 0;
	int tuneWindowAlleleAccept = 0, tuneWindowAlleleTotal = 0;
	int tuneWindowFStatAccept = 0, tuneWindowFStatTotal = 0;

	// Heap-allocated arrays for MCMC (avoid stack overflow with large MAXINDIV)
	double *logLOrig = new double[noIndiv];
	double *logLProposed = new double[noIndiv];
	double *logLVec = new double[noIndiv];

	// Allocate tempIndiv for MCMC proposals (reused each iteration)
	indiv tempIndiv;
	tempIndiv.genotype = new GenotypeType[noLoci][2];

	gsl_permutation_init (p);

	// Print pre-run summary
	std::cout << "  Run Configuration:\n";
	std::cout << "    Input:       " << infileName << "\n";
	std::cout << "    Output:      " << gArgs.outfileName << "\n";
	if (gArgs.usingFreqFile)
		std::cout << "    Freq file:   " << gArgs.freqFileName << "\n";
	std::cout << "    Individuals: " << noIndiv << "\n";
	std::cout << "    Populations: " << noPopln << "\n";
	std::cout << "    Loci:        " << noLoci << "\n";
	std::cout << "\n";
	std::cout << "  MCMC Settings:\n";
	std::cout << "    Iterations:  " << formatIterCount(mciter) << "\n";
	std::cout << "    Burn-in:     " << formatIterCount(gArgs.burnin) << "\n";
	std::cout << "    Sampling:    " << gArgs.sampling << "\n";
	std::cout << "    Seed:        " << gArgs.seed << "\n";
	std::cout << "    Mixing:      dM=" << gArgs.deltaM << " dA=" << gArgs.deltaA << " dF=" << gArgs.deltaF;
	if (gArgs.autotune)
		std::cout << " (initial, autotune on)";
	std::cout << "\n\n";

	// Start timing for progress bar
	auto startTime = std::chrono::steady_clock::now();

	for(unsigned int i = 1; i <= mciter; i++)

	{
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
			logLprop = logLik(tempIndiv,alleleFreqs,logAlleleFreqs,FStat,log1MinusFStat,noLoci);
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

		// Compute Hastings correction for asymmetric proposal
		// The proposal samples from non-empty categories only, so we need to correct
		// for the different number of non-empty categories in states X and Y.
		// Hastings ratio = q(Y->X) / q(X->Y) = nonEmptyX / nonEmptyY
		double logHastings = 0.0;
		if ((tempIndiv.migrantPopln != sampleIndiv[chosenIndiv].migrantPopln)||
			(tempIndiv.migrantAge != sampleIndiv[chosenIndiv].migrantAge))
		{
			int nonEmptyX = countNonEmptyAncestryCategories(migrantCounts, samplePopln, noPopln);
			int nonEmptyY = nonEmptyX;

			// If dropping from a category with count=1, it becomes empty
			if (migrantCounts[samplePopln][sampleIndiv[chosenIndiv].migrantPopln][sampleIndiv[chosenIndiv].migrantAge] == 1)
				nonEmptyY--;

			// If adding to an empty category, it becomes non-empty
			if (migrantCounts[samplePopln][tempIndiv.migrantPopln][tempIndiv.migrantAge] == 0)
				nonEmptyY++;

			logHastings = log((double)nonEmptyX) - log((double)nonEmptyY);
		}

		// Acceptance-rejection step
		alpha = gsl_rng_uniform(r);
		if(!NOLIKELIHOOD)
			logPrMHR = dtLogPrCount + dtLogL + logHastings;
		else
			logPrMHR = dtLogPrCount + logHastings;

		if(alpha <= exp(logPrMHR))
		{
			sampleIndiv[chosenIndiv].migrantAge = tempIndiv.migrantAge;
			sampleIndiv[chosenIndiv].migrantPopln = tempIndiv.migrantPopln;
			if (!NOLIKELIHOOD)
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
					propMigrationRates[j]=std::fabs(propMigrationRates[j]);
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
		if (gArgs.autotune && i <= (unsigned int)gArgs.burnin) tuneWindowMigAccept++;
	}
	else migrationAcceptRate = ((i-1.0)/i)*migrationAcceptRate;
	if (gArgs.autotune && i <= (unsigned int)gArgs.burnin) tuneWindowMigTotal++;

}

if(!NOALLELEMCMC)
{

	/* propose a change to a population allele frequency */

	double propAlleleFreq[MAXALLELE];
	double origAlleleFreq[MAXALLELE];
	unsigned int chosenPopln = gsl_rng_uniform_int(r,noPopln);
	unsigned int chosenLocus = gsl_rng_uniform_int(r,noLoci);
	unsigned int chosenAllele = gsl_rng_uniform_int(r,noAlleles[chosenLocus]);
	// logLOrig and logLProposed are now heap-allocated before the loop
	double origLogAlleleFreq[MAXALLELE];  // Store original log allele frequencies
	if (!NOLIKELIHOOD)
	{
		for (unsigned int l=0; l<noIndiv; l++)
			logLOrig[l] = oneLocusLogLik(sampleIndiv[l],alleleFreqs,logAlleleFreqs,FStat,log1MinusFStat,chosenLocus);
	}

	propAlleleFreq[chosenAllele] = 	std::fabs(alleleFreqs[chosenPopln][chosenLocus][chosenAllele]+(gsl_rng_uniform(r)-0.5)*gArgs.deltaA);

	// added July 20, 2011 (further modified Nov 7, 2011)
	if (propAlleleFreq[chosenAllele] > 1.0) {
		propAlleleFreq[chosenAllele] = std::fabs(1.0 - propAlleleFreq[chosenAllele]);
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
		origLogAlleleFreq[l] = logAlleleFreqs[chosenPopln][chosenLocus][l];
		alleleFreqs[chosenPopln][chosenLocus][l] = propAlleleFreq[l];
		logAlleleFreqs[chosenPopln][chosenLocus][l] = log(propAlleleFreq[l]);
	}
	if (!NOLIKELIHOOD)
	{
		for (unsigned int l=0; l<noIndiv; l++)
			logLProposed[l] = oneLocusLogLik(sampleIndiv[l],alleleFreqs,logAlleleFreqs,FStat,log1MinusFStat,chosenLocus);
		logLprop=0; dtLogL=0;
		for (unsigned int l=0; l<noIndiv; l++)
		{
			dtLogL +=  (logLProposed[l] - logLOrig[l]);
		}
	}

	// Acceptance-rejection step
	// Note: No Jacobian correction needed here because the rescaling proposal is symmetric
	// (the forward and reverse transformations are exact inverses)
	alpha = gsl_rng_uniform(r);
	if (!NOLIKELIHOOD)
	{
		logPrMHR = dtLogL;
	}
	else
	{
		logPrMHR = 0.0;
	}

	if(alpha <= exp(logPrMHR))
	{
		if (!NOLIKELIHOOD)
		{
			for (unsigned int k=0; k<noIndiv; k++)
				sampleIndiv[k].logL = sampleIndiv[k].logL - logLOrig[k] + logLProposed[k];
		}
		allelefreqAcceptRate = (1.0/i)+((i-1.0)/i)*allelefreqAcceptRate;
		if (gArgs.autotune && i <= (unsigned int)gArgs.burnin) tuneWindowAlleleAccept++;
	}
	else
	{
		for (unsigned int k=0; k<noAlleles[chosenLocus]; k++)
		{
			alleleFreqs[chosenPopln][chosenLocus][k] = origAlleleFreq[k];
			logAlleleFreqs[chosenPopln][chosenLocus][k] = origLogAlleleFreq[k];
		}
		allelefreqAcceptRate = ((i-1.0)/i)*allelefreqAcceptRate;
	}
	if (gArgs.autotune && i <= (unsigned int)gArgs.burnin) tuneWindowAlleleTotal++;

}

if(!NOFSTATMCMC)
{
	/* propose a change to a population inbreeding coefficient */

	double origLik=0.0,propLik=0.0;
	// logLVec is now heap-allocated before the loop
	double propFStat[MAXPOPLN];
	double propLogFStat[MAXPOPLN];
	double propLog1MinusFStat[MAXPOPLN];
	for (unsigned int l=0; l<noPopln; l++)
	{
		propFStat[l]=FStat[l];
		propLogFStat[l]=logFStat[l];
		propLog1MinusFStat[l]=log1MinusFStat[l];
	}
	unsigned int chosenPopln = gsl_rng_uniform_int(r,noPopln);
	double prop=std::fabs(FStat[chosenPopln]+(gsl_rng_uniform(r)-0.5)*gArgs.deltaF);
	if (prop<=1)
		propFStat[chosenPopln] = prop;
	else
		propFStat[chosenPopln] = 2.0-prop;
	// Update log values for proposed FStat
	propLogFStat[chosenPopln] = (propFStat[chosenPopln] > 1e-15) ? log(propFStat[chosenPopln]) : log(1e-15);
	propLog1MinusFStat[chosenPopln] = (propFStat[chosenPopln] < 1.0 - 1e-15) ? log(1.0 - propFStat[chosenPopln]) : log(1e-15);
	if (!NOLIKELIHOOD)
	{
		for(unsigned int l=0; l<noIndiv; l++)
		{
			origLik+=sampleIndiv[l].logL;
			if((sampleIndiv[l].samplePopln==chosenPopln)||(sampleIndiv[l].migrantPopln==chosenPopln))
			{
				logLVec[l]=logLik(sampleIndiv[l],alleleFreqs,logAlleleFreqs,propFStat,propLog1MinusFStat,noLoci);
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
		logFStat[chosenPopln] = propLogFStat[chosenPopln];
		log1MinusFStat[chosenPopln] = propLog1MinusFStat[chosenPopln];
		if (!NOLIKELIHOOD)
		{
			for(unsigned int l=0; l<noIndiv; l++)
			{
				sampleIndiv[l].logL = logLVec[l];
			}
		}
		FStatAcceptRate = (1.0/i)+((i-1.0)/i)*FStatAcceptRate;
		if (gArgs.autotune && i <= (unsigned int)gArgs.burnin) tuneWindowFStatAccept++;
	}
	else
		FStatAcceptRate = ((i-1.0)/i)*FStatAcceptRate;
	if (gArgs.autotune && i <= (unsigned int)gArgs.burnin) tuneWindowFStatTotal++;
}

if(!NOMISSINGDATA)
{
	/* propose a change to a missing genotype */
	if(noMissingGenotypes>0)
	{
		double propLogL=0.0, origLogL=0.0;
		int chooseIndiv = missingData[gsl_rng_uniform_int(r, missingData.size())];
		int chosenLocus = sampleIndiv[chooseIndiv].missingGenotypes[gsl_rng_uniform_int(r, sampleIndiv[chooseIndiv].missingGenotypes.size())];
		origLogL=oneLocusLogLik(sampleIndiv[chooseIndiv],alleleFreqs,logAlleleFreqs,FStat,log1MinusFStat,chosenLocus);
		int origAllele1=sampleIndiv[chooseIndiv].genotype[chosenLocus][0];
		int origAllele2=sampleIndiv[chooseIndiv].genotype[chosenLocus][1];
		sampleIndiv[chooseIndiv].genotype[chosenLocus][0] = gsl_rng_uniform_int(r, noAlleles[chosenLocus]);
		sampleIndiv[chooseIndiv].genotype[chosenLocus][1] = gsl_rng_uniform_int(r, noAlleles[chosenLocus]);
		propLogL=oneLocusLogLik(sampleIndiv[chooseIndiv],alleleFreqs,logAlleleFreqs,FStat,log1MinusFStat,chosenLocus);

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

// Autotune: adjust delta values during burn-in to achieve target acceptance rate
if (gArgs.autotune && i <= (unsigned int)gArgs.burnin && (i % AUTOTUNE_INTERVAL) == 0 && i > 0)
{
	// Adjust deltaM for migration rate proposals
	if (tuneWindowMigTotal > 0)
	{
		double rate = (double)tuneWindowMigAccept / tuneWindowMigTotal;
		if (rate < AUTOTUNE_LOWER_BOUND)
		{
			gArgs.deltaM /= AUTOTUNE_ADJUST_FACTOR;
			if (gArgs.deltaM < AUTOTUNE_DELTA_MIN) gArgs.deltaM = AUTOTUNE_DELTA_MIN;
		}
		else if (rate > AUTOTUNE_UPPER_BOUND)
		{
			gArgs.deltaM *= AUTOTUNE_ADJUST_FACTOR;
			if (gArgs.deltaM > AUTOTUNE_DELTA_MAX) gArgs.deltaM = AUTOTUNE_DELTA_MAX;
		}
	}

	// Adjust deltaA for allele frequency proposals
	if (tuneWindowAlleleTotal > 0)
	{
		double rate = (double)tuneWindowAlleleAccept / tuneWindowAlleleTotal;
		if (rate < AUTOTUNE_LOWER_BOUND)
		{
			gArgs.deltaA /= AUTOTUNE_ADJUST_FACTOR;
			if (gArgs.deltaA < AUTOTUNE_DELTA_MIN) gArgs.deltaA = AUTOTUNE_DELTA_MIN;
		}
		else if (rate > AUTOTUNE_UPPER_BOUND)
		{
			gArgs.deltaA *= AUTOTUNE_ADJUST_FACTOR;
			if (gArgs.deltaA > AUTOTUNE_DELTA_MAX) gArgs.deltaA = AUTOTUNE_DELTA_MAX;
		}
	}

	// Adjust deltaF for F-statistic proposals
	if (tuneWindowFStatTotal > 0)
	{
		double rate = (double)tuneWindowFStatAccept / tuneWindowFStatTotal;
		if (rate < AUTOTUNE_LOWER_BOUND)
		{
			gArgs.deltaF /= AUTOTUNE_ADJUST_FACTOR;
			if (gArgs.deltaF < AUTOTUNE_DELTA_MIN) gArgs.deltaF = AUTOTUNE_DELTA_MIN;
		}
		else if (rate > AUTOTUNE_UPPER_BOUND)
		{
			gArgs.deltaF *= AUTOTUNE_ADJUST_FACTOR;
			if (gArgs.deltaF > AUTOTUNE_DELTA_MAX) gArgs.deltaF = AUTOTUNE_DELTA_MAX;
		}
	}

	// Reset tuning window counters
	tuneWindowMigAccept = tuneWindowMigTotal = 0;
	tuneWindowAlleleAccept = tuneWindowAlleleTotal = 0;
	tuneWindowFStatAccept = tuneWindowFStatTotal = 0;

	// Print tuning progress in verbose mode
	if (gArgs.verbose && (i % 10000) == 0)
	{
		std::cout << "\n[Autotune @ " << i << "] deltaM=" << gArgs.deltaM
		          << " deltaA=" << gArgs.deltaA << " deltaF=" << gArgs.deltaF << std::flush;
	}
}

	// Print trace file header (once at start, before any data)
		if ((i==1) && gArgs.trace)
		{
			tracefile << "State\t" << "LogProb\t";
			for (unsigned int l = 0; l < noPopln; l++)
				for (unsigned int k = 0; k < noPopln; k++)
					tracefile << "m[" << l << "][" << k << "]\t";
			tracefile << "\n";
		}

	// Print logL to trace file
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
			// Update progress bar every 1% (or more frequently for short runs)
			unsigned long int updateInterval = mciter / 100;
			if (updateInterval < 1) updateInterval = 1;

			if ((i % updateInterval) == 0 || i == mciter)
			{
				auto now = std::chrono::steady_clock::now();
				double elapsed = std::chrono::duration<double>(now - startTime).count();
				bool inBurnin = (i < gArgs.burnin);
				printProgress(i, mciter, elapsed, inBurnin);
			}
		}

		// Print migration rates to trace file
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

			// Update Savage-Dickey statistics for migration rate hypothesis testing
			updateSavageDickeyStats(sdStats, migrationRates, noPopln, SD_BANDWIDTH);

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

	// Free heap-allocated MCMC arrays
	delete[] logLOrig;
	delete[] logLProposed;
	delete[] logLVec;
	delete[] tempIndiv.genotype;

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
		mcmcout << "Mixing parameters" << (gArgs.autotune ? " (autotuned)" : "") << ": (dM=" << gArgs.deltaM << ",dA=" << gArgs.deltaA << ",dF=" << gArgs.deltaF << ")" << " Output file=" << gArgs.outfileName << "\n";
	}
	mcmcout << "Individuals: " << noIndiv << " Populations: " << noPopln << " Loci: " << noLoci << "\n";
mcmcout << "\n Population Labels:\n";

	// Build a vector of population names indexed by their numeric index
	std::vector<std::string> poplnNames(noPopln);
	IndivMap::iterator iterPopln = poplnIDMap.begin();
	while (iterPopln != poplnIDMap.end())
	{
		poplnNames[iterPopln->second] = iterPopln->first;
		iterPopln++;
	}

	for (unsigned int l = 0; l < noPopln; l++)
	{
		mcmcout << "  [" << l << "] " << poplnNames[l] << "\n";
	}

	mcmcout << "\n Migration Rate Matrix m[source][dest]:\n";
	mcmcout << " Mean(SD)\n";
	mcmcout.setf(std::ios::fixed, std::ios::floatfield);

	// Print column headers
	mcmcout << "       ";
	for (unsigned int k = 0; k < noPopln; k++)
	{
		std::ostringstream hdr;
		hdr << "[" << k << "]";
		mcmcout << std::setw(15) << hdr.str();
	}
	mcmcout << "\n";

	// Print each row
	for (unsigned int l = 0; l < noPopln; l++)
	{
		std::ostringstream rowHdr;
		rowHdr << "[" << l << "]";
		mcmcout << " " << std::setw(4) << rowHdr.str() << " ";
		for(unsigned int k=0; k < noPopln; k++)
		{
			std::ostringstream entry;
			entry << std::fixed << std::setprecision(4) << avgMigrationRates[l][k]
			      << "(" << std::setprecision(4) << sqrt(varMigrationRates[l][k]) << ")";
			mcmcout << std::setw(15) << entry.str();
		}
		mcmcout << "\n";
	}

	// Output Savage-Dickey test results for zero migration hypotheses
	computeSavageDickeyBayesFactors(sdStats, noPopln, PRIOR_DENSITY_AT_ZERO, mcmcout, poplnNames);

	mcmcout << "\n Inbreeding Coefficients:\n";
	mcmcout << " Index  Population                     F(SD)\n";
	mcmcout << " -----  ----------                     -----\n";
	mcmcout.setf(std::ios::fixed, std::ios::floatfield);
	for (unsigned int l = 0; l < noPopln; l++)
	{
		std::ostringstream fstat;
		fstat << std::fixed << std::setprecision(4) << avgFStat[l]
		      << "(" << std::setprecision(4) << sqrt(varFStat[l]) << ")";
		mcmcout << " [" << std::setw(2) << l << "]  "
		        << std::left << std::setw(25) << poplnNames[l]
		        << std::right << std::setw(14) << fstat.str() << "\n";
	}

	// Write allele frequencies to separate file if -F specified
	if (gArgs.usingFreqFile)
	{
		// Build locus names vector
		std::vector<std::string> locusNames(noLoci);
		for (unsigned int k = 0; k < noLoci; k++)
		{
			IndivMap::iterator iterLocus = locusIDMap.begin();
			while (iterLocus != locusIDMap.end())
			{
				if (iterLocus->second == k)
					locusNames[k] = iterLocus->first;
				iterLocus++;
			}
		}

		// Write TSV header
		freqout << "Population\tLocus\tAllele\tFrequency\tSD\n";

		for (unsigned int l = 0; l < noPopln; l++)
		{
			for (unsigned int k = 0; k < noLoci; k++)
			{
				for (unsigned int j = 0; j < noAlleles[k]; j++)
				{
					string alleleName;
					IndivMap::iterator iterAllele = alleleIDMap[k].begin();
					while (iterAllele != alleleIDMap[k].end())
					{
						if (iterAllele->second == j)
							alleleName = iterAllele->first;
						iterAllele++;
					}

					freqout << poplnNames[l] << "\t"
					        << locusNames[k] << "\t"
					        << alleleName << "\t"
					        << std::fixed << std::setprecision(6) << avgAlleleFreqs[l][k][j] << "\t"
					        << std::setprecision(6) << sqrt(varAlleleFreqs[l][k][j]) << "\n";
				}
			}
		}
		freqout.close();
		mcmcout << "\n Allele frequencies written to: " << gArgs.freqFileName << "\n";
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
				it1 = find(missingData.begin(), missingData.end(), static_cast<int>(l));
				it2 = find(sampleIndiv[l].missingGenotypes.begin(), sampleIndiv[l].missingGenotypes.end(), static_cast<int>(j));

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


	// Calculate total elapsed time
	auto endTime = std::chrono::steady_clock::now();
	double totalTime = std::chrono::duration<double>(endTime - startTime).count();

	// Print final newline after progress bar, then summary
	std::cout << "\n\n";
	std::cout << "  ┌──────────────────────────────────────┐\n";
	std::cout << "  │            Run Complete              │\n";
	std::cout << "  └──────────────────────────────────────┘\n";
	std::cout << "  Elapsed time: " << formatTime(totalTime) << "\n";
	if (gArgs.autotune) {
		std::cout << "  Final mixing: dM=" << std::fixed << std::setprecision(3) << gArgs.deltaM
		          << " dA=" << gArgs.deltaA << " dF=" << gArgs.deltaF << "\n";
	}
	std::cout << "\n";

	std::cout << "  Output written to: " << gArgs.outfileName << "\n";
	if (gArgs.usingFreqFile)
		std::cout << "  Allele frequencies: " << gArgs.freqFileName << "\n";
	std::cout << "\n";

	// Free genotypes before freeing sampleIndiv
	freeGenotypes(sampleIndiv, noIndiv);
	delete[] sampleIndiv;
	sampleIndiv = NULL;

	// Free noAlleles and alleleIDMap
	delete[] noAlleles;
	delete[] alleleIDMap;
	alleleIDMap = nullptr;

	// Free ancP (3D array: noIndiv x noPopln x 3)
	for(unsigned int i = 0; i < noIndiv; i++)
	{
		for(unsigned int j = 0; j < noPopln; j++)
			delete[] ancP[i][j];
		delete[] ancP[i];
	}
	delete[] ancP;

	// Free alleleFreqs, logAlleleFreqs, avgAlleleFreqs, varAlleleFreqs (3D arrays: noPopln x noLoci x noAlleles[j])
	for(unsigned int i = 0; i < noPopln; i++)
	{
		for(unsigned int j = 0; j < noLoci; j++)
		{
			delete[] alleleFreqs[i][j];
			delete[] logAlleleFreqs[i][j];
			delete[] avgAlleleFreqs[i][j];
			delete[] varAlleleFreqs[i][j];
		}
		delete[] alleleFreqs[i];
		delete[] logAlleleFreqs[i];
		delete[] avgAlleleFreqs[i];
		delete[] varAlleleFreqs[i];
	}
	delete[] alleleFreqs;
	delete[] logAlleleFreqs;
	delete[] avgAlleleFreqs;
	delete[] varAlleleFreqs;

	// Free migrationRates, avgMigrationRates, varMigrationRates (2D arrays: noPopln x noPopln+1)
	for(unsigned int i = 0; i < noPopln; i++)
	{
		delete[] migrationRates[i];
		delete[] avgMigrationRates[i];
		delete[] varMigrationRates[i];
	}
	delete[] migrationRates;
	delete[] avgMigrationRates;
	delete[] varMigrationRates;

	// Free migrantCounts (3D array: noPopln x noPopln x 3)
	for(unsigned int i = 0; i < noPopln; i++)
	{
		for(unsigned int j = 0; j < noPopln; j++)
			delete[] migrantCounts[i][j];
		delete[] migrantCounts[i];
	}
	delete[] migrantCounts;

	// Free FStat arrays
	delete[] FStat;
	delete[] avgFStat;
	delete[] varFStat;

	// Free Savage-Dickey statistics
	freeSavageDickeyStats(sdStats, noPopln);

	// Free GSL objects
	gsl_permutation_free(p);
	gsl_rng_free(r);

	mcmcout.close();
	return 0;
}

void printBanner(void)
{
	const int width = 54;  // Inner width of the box

	// Build centered strings
	std::ostringstream line1, line2;
	line1 << "BayesAss Edition " << VERSION << " (BA3)";
	line2 << "Released: " << RELEASEDATE;

	std::string title = line1.str();
	std::string released = line2.str();
	std::string author = "Bruce Rannala";
	std::string dept = "Department of Evolution and Ecology at UC Davis";

	// Calculate padding for centering
	int pad1 = (width - title.length()) / 2;
	int pad2 = (width - released.length()) / 2;
	int pad3 = (width - author.length()) / 2;
	int pad4 = (width - dept.length()) / 2;

	// Build horizontal line (─ is a multi-byte UTF-8 character)
	std::string hline;
	for (int i = 0; i < width; i++) hline += "─";

	std::cout << "\n";
	std::cout << "  ┌" << hline << "┐\n";
	std::cout << "  │" << std::string(pad1, ' ') << title << std::string(width - pad1 - title.length(), ' ') << "│\n";
	std::cout << "  │" << std::string(pad2, ' ') << released << std::string(width - pad2 - released.length(), ' ') << "│\n";
	std::cout << "  │" << std::string(pad3, ' ') << author << std::string(width - pad3 - author.length(), ' ') << "│\n";
	std::cout << "  │" << std::string(pad4, ' ') << dept << std::string(width - pad4 - dept.length(), ' ') << "│\n";
	std::cout << "  └" << hline << "┘\n";
	std::cout << "\n";
}

// Format seconds into human-readable time string
std::string formatTime(double seconds)
{
	if (seconds < 0) return "--:--";
	int hrs = (int)(seconds / 3600);
	int mins = (int)((seconds - hrs * 3600) / 60);
	int secs = (int)(seconds - hrs * 3600 - mins * 60);

	std::ostringstream oss;
	if (hrs > 0)
		oss << hrs << "h " << std::setw(2) << std::setfill('0') << mins << "m";
	else if (mins > 0)
		oss << mins << "m " << std::setw(2) << std::setfill('0') << secs << "s";
	else
		oss << secs << "s";
	return oss.str();
}

// Format large numbers with K/M suffix
std::string formatIterCount(unsigned long int n)
{
	std::ostringstream oss;
	if (n >= 1000000)
		oss << std::fixed << std::setprecision(1) << (n / 1000000.0) << "M";
	else if (n >= 1000)
		oss << std::fixed << std::setprecision(0) << (n / 1000.0) << "K";
	else
		oss << n;
	return oss.str();
}

// Print progress bar with ETA
void printProgress(unsigned long int current, unsigned long int total,
                   double elapsedSecs, bool inBurnin)
{
	const int barWidth = 30;
	double fraction = (double)current / total;
	int filled = (int)(fraction * barWidth);

	// Calculate ETA
	double eta = -1;
	if (current > 0 && elapsedSecs > 0) {
		eta = (elapsedSecs / current) * (total - current);
	}

	// Build progress bar
	std::cout << "\r  ";
	if (inBurnin)
		std::cout << "Burn-in:  [";
	else
		std::cout << "Sampling: [";

	for (int i = 0; i < barWidth; i++) {
		if (i < filled)
			std::cout << "█";
		else
			std::cout << "░";
	}

	std::cout << "] " << std::setw(3) << (int)(fraction * 100) << "%  "
	          << formatIterCount(current) << "/" << formatIterCount(total);

	if (eta >= 0)
		std::cout << "  ETA: " << formatTime(eta);
	else
		std::cout << "  ETA: --:--";

	std::cout << "     " << std::flush;  // Extra spaces to clear previous longer output
}

void checkDataSize(unsigned int &outNoIndiv, unsigned int &outNoLoci, unsigned int &outNoPopln, unsigned int &outMaxAlleles)
{
  std::string aline;
  std::set<string> namesUniq;
  std::set<string> popsUniq;
  std::set<string> lociUniq;
  map<string,set<string> > alleleLabels;

  // Single pass through file - collect all info at once
  while(std::getline(mcmcin,aline))
    {
      std::istringstream iss(aline);
      std::string indName, popName, locusName, allele1, allele2;
      iss >> indName >> popName >> locusName >> allele1 >> allele2;

      if (!iss || indName.empty())
	{
	  std::cerr << "\nerror the line: \"" << aline << "\" has an incorrect number of entries. quitting...\n";
	  exit(1);
	}

      namesUniq.insert(indName);
      popsUniq.insert(popName);
      lociUniq.insert(locusName);

      /* bug fix Apr 4 2024 don't count missing data symbol "0" as an allele */
      if(allele1 != "0")
	    alleleLabels[locusName].insert(allele1);
      if(allele2 != "0")
	    alleleLabels[locusName].insert(allele2);
    }

  // Check limits
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

  // Find maximum alleles across all loci
  unsigned int maxAlleles = 0;
  for (auto all = alleleLabels.begin(); all != alleleLabels.end(); all++) {
    if(all->second.size() > MAXALLELE)
      {
	std::cerr << "\nerror: number of alleles at locus " << all->first << " is:" << all->second.size() << " exceeds maximum:" << MAXALLELE << " quitting...\n";
	exit(1);
      }
    if(all->second.size() > maxAlleles)
      maxAlleles = all->second.size();
  }

  // Return the counts
  outNoIndiv = namesUniq.size();
  outNoLoci = lociUniq.size();
  outNoPopln = popsUniq.size();
  outMaxAlleles = maxAlleles;

  // Set global dimensions
  gNoIndiv = outNoIndiv;
  gNoLoci = outNoLoci;
  gMaxAlleles = outMaxAlleles;

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

	// Use dynamically allocated array for alleleIter
	unsigned int *alleleIter = new unsigned int[gNoLoci];
	unsigned int indivIter=0, popIter=0, locIter=0;
	unsigned int currIndivID, currPoplnID, currLocusID, currAllele1, currAllele2;
    std::string inputLine=" ";

	// Genotypes are already initialized to -2 by allocateGenotypes()

	for(unsigned int l = 0; l < gNoLoci; l++)
		alleleIter[l] = 0;

	while(std::getline(mcmcin,inputLine))
	{

	int numspaces=0;

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
	// if a locus is missing for any individual, treat it as missing data and warn
	std::set<unsigned int> incompleteLoci;
	for(unsigned int k = 0; k < noIndiv; k++)
	  for(unsigned int l = 0; l < noLoci; l++)
	    if((sampleIndiv[k].genotype[l][0]==-2)||(sampleIndiv[k].genotype[l][1]==-2))
	      {
		// Mark as missing data instead of failing
		sampleIndiv[k].genotype[l][0] = -1;
		sampleIndiv[k].genotype[l][1] = -1;
		incompleteLoci.insert(l);
	      }

	// Print warnings for incomplete loci
	if (!incompleteLoci.empty()) {
		std::cerr << "\n  Warning: " << incompleteLoci.size() << " locus/loci missing for some individuals (treated as missing data):\n";
		for (unsigned int l : incompleteLoci) {
			// Find locus name
			std::string locusName = "unknown";
			for (auto it = locusIDMap.begin(); it != locusIDMap.end(); ++it) {
				if (it->second == l) {
					locusName = it->first;
					break;
				}
			}
			// Count how many individuals are missing this locus
			int missingCount = 0;
			for (unsigned int k = 0; k < noIndiv; k++) {
				if (sampleIndiv[k].genotype[l][0] == -1 && sampleIndiv[k].genotype[l][1] == -1) {
					missingCount++;
				}
			}
			std::cerr << "    " << locusName << " (missing in " << missingCount << " individuals)\n";
		}
		std::cerr << "\n";
	}

	delete[] alleleIter;
}

void getEmpiricalAlleleFreqs(double ***alleleFreqs, indiv *sampleIndiv, unsigned int *noAlleles, unsigned int noPopln, unsigned int noLoci, unsigned int noIndiv)
{
	double epsilon=0.0001; // minimum allele frequency in any population
	std::vector<int> indivPerPopln(noPopln);
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

	// Properly free all levels of the 3D array
	for(unsigned int i = 0; i < noPopln; i++)
	{
		for(unsigned int j = 0; j < noLoci; j++)
			delete[] poplnAlleleCounts[i][j];
		delete[] poplnAlleleCounts[i];
	}
	delete[] poplnAlleleCounts;
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


// Helper functions for log allele frequency optimization
void updateLogAlleleFreqs(double ***logAlleleFreqs, double ***alleleFreqs, unsigned int popln, unsigned int locus, unsigned int nAlleles)
{
	for (unsigned int a = 0; a < nAlleles; a++)
	{
		logAlleleFreqs[popln][locus][a] = log(alleleFreqs[popln][locus][a]);
	}
}

void updateAllLogAlleleFreqs(double ***logAlleleFreqs, double ***alleleFreqs, unsigned int noPopln, unsigned int noLoci, unsigned int *noAlleles)
{
	for (unsigned int p = 0; p < noPopln; p++)
	{
		for (unsigned int l = 0; l < noLoci; l++)
		{
			for (unsigned int a = 0; a < noAlleles[l]; a++)
			{
				logAlleleFreqs[p][l][a] = log(alleleFreqs[p][l][a]);
			}
		}
	}
}

void updateLogFStat(double *logFStat, double *log1MinusFStat, double *FStat, unsigned int noPopln)
{
	for (unsigned int i = 0; i < noPopln; i++)
	{
		// Handle edge case where FStat is 0 or 1
		logFStat[i] = (FStat[i] > 1e-15) ? log(FStat[i]) : log(1e-15);
		log1MinusFStat[i] = (FStat[i] < 1.0 - 1e-15) ? log(1.0 - FStat[i]) : log(1e-15);
	}
}

// Constant log(2) for heterozygote calculations
static const double LOG2 = log(2.0);

double logLik(const indiv& Indiv, double ***alleleFreqs, double ***logAlleleFreqs, double *FStat, double *log1MinusFStat, unsigned int noLoci)
{
	double logPr = 0.0;

	if (Indiv.migrantAge == 0)
	{
		const unsigned int pop = Indiv.samplePopln;
		const double F = FStat[pop];
		const double oneMinusF = 1.0 - F;
		const double log1MF = log1MinusFStat[pop];
		double **freqPop = alleleFreqs[pop];
		double **logFreqPop = logAlleleFreqs[pop];
		const GenotypeType (*geno)[2] = Indiv.genotype;

		for (unsigned int i = 0; i < noLoci; i++)
		{
			const int a0 = geno[i][0];
			const int a1 = geno[i][1];

			if (a0 == a1)
			{
				// Homozygote: use original formula with direct allele freqs
				const double p = freqPop[i][a0];
				logPr += log(oneMinusF * p * p + F * p);
			}
			else
			{
				// Heterozygote: log(2*(1-F)*p0*p1) = log(2) + log(1-F) + log(p0) + log(p1)
				logPr += LOG2 + log1MF + logFreqPop[i][a0] + logFreqPop[i][a1];
			}
		}
	}
	else if (Indiv.migrantAge == 1)
	{
		const unsigned int pop = Indiv.migrantPopln;
		const double F = FStat[pop];
		const double oneMinusF = 1.0 - F;
		const double log1MF = log1MinusFStat[pop];
		double **freqPop = alleleFreqs[pop];
		double **logFreqPop = logAlleleFreqs[pop];
		const GenotypeType (*geno)[2] = Indiv.genotype;

		for (unsigned int i = 0; i < noLoci; i++)
		{
			const int a0 = geno[i][0];
			const int a1 = geno[i][1];

			if (a0 == a1)
			{
				const double p = freqPop[i][a0];
				logPr += log(oneMinusF * p * p + F * p);
			}
			else
			{
				logPr += LOG2 + log1MF + logFreqPop[i][a0] + logFreqPop[i][a1];
			}
		}
	}
	else  // migrantAge == 2 (hybrid)
	{
		const unsigned int migPop = Indiv.migrantPopln;
		const unsigned int samPop = Indiv.samplePopln;
		double **logFreqMig = logAlleleFreqs[migPop];
		double **logFreqSam = logAlleleFreqs[samPop];
		const GenotypeType (*geno)[2] = Indiv.genotype;

		for (unsigned int i = 0; i < noLoci; i++)
		{
			const int a0 = geno[i][0];
			const int a1 = geno[i][1];

			if (a0 == a1)
			{
				// Homozygote: log(p_mig * p_sample) = log(p_mig) + log(p_sample)
				logPr += logFreqMig[i][a0] + logFreqSam[i][a1];
			}
			else
			{
				// Heterozygote: log(p0_mig*p1_sam + p1_mig*p0_sam)
				// Use log-sum-exp for numerical stability
				const double term1 = logFreqMig[i][a0] + logFreqSam[i][a1];
				const double term2 = logFreqMig[i][a1] + logFreqSam[i][a0];
				const double maxTerm = (term1 > term2) ? term1 : term2;
				logPr += maxTerm + log(exp(term1 - maxTerm) + exp(term2 - maxTerm));
			}
		}
	}
	return logPr;
}

double oneLocusLogLik(const indiv& Indiv, double ***alleleFreqs, double ***logAlleleFreqs, double *FStat, double *log1MinusFStat, int chosenLocus)
{
	double logPr = 0.0;
	const int a0 = Indiv.genotype[chosenLocus][0];
	const int a1 = Indiv.genotype[chosenLocus][1];

	if (Indiv.migrantAge == 0)
	{
		const unsigned int pop = Indiv.samplePopln;
		const double F = FStat[pop];
		const double oneMinusF = 1.0 - F;
		const double log1MF = log1MinusFStat[pop];

		if (a0 == a1)
		{
			const double p = alleleFreqs[pop][chosenLocus][a0];
			logPr = log(oneMinusF * p * p + F * p);
		}
		else
		{
			logPr = LOG2 + log1MF + logAlleleFreqs[pop][chosenLocus][a0] + logAlleleFreqs[pop][chosenLocus][a1];
		}
	}
	else if (Indiv.migrantAge == 1)
	{
		const unsigned int pop = Indiv.migrantPopln;
		const double F = FStat[pop];
		const double oneMinusF = 1.0 - F;
		const double log1MF = log1MinusFStat[pop];

		if (a0 == a1)
		{
			const double p = alleleFreqs[pop][chosenLocus][a0];
			logPr = log(oneMinusF * p * p + F * p);
		}
		else
		{
			logPr = LOG2 + log1MF + logAlleleFreqs[pop][chosenLocus][a0] + logAlleleFreqs[pop][chosenLocus][a1];
		}
	}
	else  // migrantAge == 2
	{
		const unsigned int migPop = Indiv.migrantPopln;
		const unsigned int samPop = Indiv.samplePopln;

		if (a0 == a1)
		{
			logPr = logAlleleFreqs[migPop][chosenLocus][a0] + logAlleleFreqs[samPop][chosenLocus][a1];
		}
		else
		{
			const double term1 = logAlleleFreqs[migPop][chosenLocus][a0] + logAlleleFreqs[samPop][chosenLocus][a1];
			const double term2 = logAlleleFreqs[migPop][chosenLocus][a1] + logAlleleFreqs[samPop][chosenLocus][a0];
			const double maxTerm = (term1 > term2) ? term1 : term2;
			logPr = maxTerm + log(exp(term1 - maxTerm) + exp(term2 - maxTerm));
		}
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

/*
 * Count the number of non-empty ancestry categories for a given sample population.
 * This is needed for computing the Hastings correction ratio in the ancestry MCMC.
 *
 * Categories are:
 * - (samplePopln, age=0): non-migrant
 * - (otherPopln, age=1): first-generation migrant from otherPopln
 * - (otherPopln, age=2): second-generation migrant from otherPopln
 *
 * Total possible categories = 1 + 2*(noPopln-1) = 2*noPopln - 1
 */
int countNonEmptyAncestryCategories(long int ***migrantCounts, unsigned int samplePopln, unsigned int noPopln)
{
	int count = 0;

	// Check non-migrant category (samplePopln, age=0)
	if (migrantCounts[samplePopln][samplePopln][0] > 0) count++;

	// Check migrant categories (other populations, ages 1 and 2)
	for (unsigned int p = 0; p < noPopln; p++) {
		if (p != samplePopln) {
			if (migrantCounts[samplePopln][p][1] > 0) count++;
			if (migrantCounts[samplePopln][p][2] > 0) count++;
		}
	}

	return count;
}

/*
 * Savage-Dickey Density Ratio Test Functions
 *
 * These functions implement the Savage-Dickey density ratio test for
 * testing the null hypothesis H0: m_ij = 0 (no migration from pop j to pop i)
 * against the alternative H1: m_ij > 0.
 *
 * The Bayes factor is: BF_01 = p(m_ij=0|data) / p(m_ij=0|prior)
 *
 * We use a boundary-corrected kernel density estimate (reflection method)
 * to estimate the posterior density at zero.
 */

// Threshold for counting samples "near zero"
const double SD_NEAR_ZERO_THRESHOLD = 0.005;

// Initialize Savage-Dickey statistics structure
void initSavageDickeyStats(SavageDickeyStats **sdStats, unsigned int noPopln)
{
	for (unsigned int i = 0; i < noPopln; i++) {
		for (unsigned int j = 0; j < noPopln; j++) {
			sdStats[i][j].kernelSum = 0.0;
			sdStats[i][j].sumM = 0.0;
			sdStats[i][j].sumM2 = 0.0;
			sdStats[i][j].countNearZero = 0;
			sdStats[i][j].nSamples = 0;
		}
	}
}

// Update Savage-Dickey statistics with current migration rate samples
// Uses Gaussian kernel: K(x) = exp(-x^2 / (2*h^2))
void updateSavageDickeyStats(SavageDickeyStats **sdStats, double **migrationRates,
                              unsigned int noPopln, double bandwidth)
{
	for (unsigned int i = 0; i < noPopln; i++) {
		for (unsigned int j = 0; j < noPopln; j++) {
			if (i == j) continue;  // Skip diagonal (non-migration)

			double m = migrationRates[i][j];

			// Update running sums for mean and variance
			sdStats[i][j].sumM += m;
			sdStats[i][j].sumM2 += m * m;
			sdStats[i][j].nSamples++;

			// Update count near zero
			if (m < SD_NEAR_ZERO_THRESHOLD) {
				sdStats[i][j].countNearZero++;
			}

			// Update kernel sum for density estimation at zero
			// Using Gaussian kernel: exp(-m^2 / (2*h^2))
			double h2 = bandwidth * bandwidth;
			sdStats[i][j].kernelSum += exp(-m * m / (2.0 * h2));
		}
	}
}

// Compute and output Savage-Dickey Bayes factors
void computeSavageDickeyBayesFactors(SavageDickeyStats **sdStats, unsigned int noPopln,
                                      double priorDensityAtZero, std::ostream &out,
                                      const std::vector<std::string> &poplnNames)
{
	const double PI = 3.14159265358979323846;

	out << "\n\n Savage-Dickey Density Ratio Test for Zero Migration:\n";
	out << " (Tests H0: m_ij = 0 vs H1: m_ij > 0)\n\n";
	out << " Source  Dest     Mean(SD)       BF_01  log10(BF)  KL(bits)  Interpretation\n";
	out << " ---------------------------------------------------------------------------\n";

	// Prior standard deviation for uniform on [0, 1/3]
	const double priorSD = (1.0/3.0) / sqrt(12.0);  // ≈ 0.0962
	const double LOG2E = 1.4426950408889634;  // log2(e)

	for (unsigned int i = 0; i < noPopln; i++) {
		for (unsigned int j = 0; j < noPopln; j++) {
			if (i == j) continue;  // Skip diagonal

			long int N = sdStats[i][j].nSamples;
			if (N < 100) {
				out << "  [" << std::setw(2) << j << "]  [" << std::setw(2) << i << "]   Insufficient samples\n";
				continue;
			}

			// Compute mean and standard deviation
			double mean = sdStats[i][j].sumM / N;
			double var = (sdStats[i][j].sumM2 / N) - (mean * mean);
			double sd = sqrt(var > 0 ? var : 0);

			// Compute optimal bandwidth using Silverman's rule
			double bandwidth = 1.06 * sd * pow((double)N, -0.2);
			if (bandwidth < 1e-6) bandwidth = 0.01;  // Minimum bandwidth

			// Compute posterior density at zero using reflection method
			// p(0|D) = 2 * (1/(N*h*sqrt(2*pi))) * sum(K(m_i/h))
			double postDensityAtZero = 2.0 * sdStats[i][j].kernelSum / (N * bandwidth * sqrt(2.0 * PI));

			// Alternative: histogram-based estimate
			double histDensity = (double)sdStats[i][j].countNearZero / (N * SD_NEAR_ZERO_THRESHOLD);

			// Use the KDE estimate (more accurate for smooth distributions)
			// but check against histogram for sanity
			if (postDensityAtZero < 0.01 * histDensity || postDensityAtZero > 100 * histDensity) {
				// Large discrepancy - use histogram estimate as fallback
				postDensityAtZero = histDensity;
			}

			// Compute Bayes factor
			double BF01 = postDensityAtZero / priorDensityAtZero;
			double log10BF = log10(BF01 > 0 ? BF01 : 1e-10);

			// Compute KL divergence (information gain) in bits
			// KL(posterior || prior) ≈ log(σ_prior / σ_posterior) for Gaussian approx
			// In bits: KL = log2(σ_prior / σ_posterior)
			double klBits = 0.0;
			if (sd > 0) {
				klBits = LOG2E * log(priorSD / sd);
				if (klBits < 0) klBits = 0;  // Can't have negative info gain
			}

			// Interpretation
			std::string interp;
			if (log10BF > 2) interp = "Decisive for H0";
			else if (log10BF > 1) interp = "Strong for H0";
			else if (log10BF > 0.5) interp = "Substantial for H0";
			else if (log10BF > 0) interp = "Weak for H0";
			else if (log10BF > -0.5) interp = "Weak for H1";
			else if (log10BF > -1) interp = "Substantial for H1";
			else if (log10BF > -2) interp = "Strong for H1";
			else interp = "Decisive for H1";

			// Format Mean(SD) string with 2 decimal places
			std::ostringstream meanStr;
			meanStr << std::fixed << std::setprecision(2) << mean << "(" << std::setprecision(2) << sd << ")";

			out << "  [" << std::setw(2) << j << "]  [" << std::setw(2) << i << "]  ";
			out << std::right << std::setw(12) << meanStr.str();
			out << std::fixed << std::setw(10) << std::setprecision(2) << BF01;
			out << std::setw(9) << std::setprecision(2) << std::showpos << log10BF << std::noshowpos;
			out << std::setw(9) << std::setprecision(2) << klBits;
			out << "  " << interp << "\n";
		}
	}

	out << " ---------------------------------------------------------------------------\n";
	out << " Note: BF_01 > 1 supports H0 (no migration); BF_01 < 1 supports H1 (migration)\n";
	out << " KL(bits) = information gain from prior to posterior\n";
}

// Free Savage-Dickey statistics memory
void freeSavageDickeyStats(SavageDickeyStats **sdStats, unsigned int noPopln)
{
	for (unsigned int i = 0; i < noPopln; i++) {
		delete[] sdStats[i];
	}
	delete[] sdStats;
}
