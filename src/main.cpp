
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
	char outfileName[256];
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
	int collapse;            // integrate out allele frequencies (collapsed sampler); default off
	int fixGamma;            // fix gamma = 1 (original BayesAss: tau=2, migration cap 1/3); default off
	int useVCF;              // Use VCF input format
	char vcfFileName[256];   // VCF file path
	char metaFileName[256];  // Metadata file path (INDIV POPLN)
	int usingFreqFile;       // Output allele frequencies to separate file
	char freqFileName[256];  // Allele frequencies file path
} gArgs;

static const char *optString = "s:i:n:b:o:m:a:f:V:M:F:cGvugtdphTN?";

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
	{ "collapse", no_argument, NULL, 'c'},
	{ "fixgamma", no_argument, NULL, 'G'},
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

// per-locus inheritance type (0=autosomal, 1=X-linked); sized to noLoci at read time
std::vector<int> gLocusType;

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

// Sex-biased dispersal: Beta(a,b) prior on the global female fraction phi.
// a = b = 1 is uniform on (0,1); phi = 1/2 corresponds to no sex bias.
const double PHI_PRIOR_A = 1.0;
const double PHI_PRIOR_B = 1.0;

// Migrant breeding-success multiplier gamma (tau = 2*gamma = age-2:age-1 rate
// ratio). Weakly-informative log-normal prior centered at gamma = 1 (no change
// in migrant breeding success); GAMMA_PROP_SD scales the log-scale random walk.
const double GAMMA_PRIOR_LOGSD = 1.0;
const double GAMMA_PROP_SD = 0.15;

// Collapsed sampler: symmetric Dirichlet concentration on population allele
// frequencies when they are integrated out analytically (Dirichlet-multinomial).
// alpha = 1 reproduces BA3's current implicit uniform prior. See
// doc/collapsed_allele_freqs_plan.md.
const double ALLELE_PRIOR_ALPHA = 1.0;

// Integral over [0,1] of the product of two Beta densities (midpoint rule).
// Used for the Rao-Blackwellized Savage-Dickey test of H0: phi = rho, since the
// posterior density of (phi - rho) at 0 equals int_0^1 f_phi(t) f_rho(t) dt.
static double betaProductIntegral(double a1, double b1, double a2, double b2) {
	const int N = 1000;
	const double h = 1.0 / N;
	double s = 0.0;
	for (int i = 0; i < N; i++) {
		double t = (i + 0.5) * h;
		s += gsl_ran_beta_pdf(t, a1, b1) * gsl_ran_beta_pdf(t, a2, b2);
	}
	return s * h;
}

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

// Read metadata file mapping individuals to populations (and optionally sex).
// Format: INDIV_ID POPLN_ID [SEX] (whitespace separated, one per line).
// SEX is optional; accepts F/M or female/male (case-insensitive). Individuals
// with no sex column are left SEX_UNKNOWN, so autosomal-only, sexless datasets
// continue to work unchanged.
void readMetadataFile(const char *metaFileName,
                      std::map<std::string, std::string> &indivToPopln,
                      std::map<std::string, int> &indivToSex) {
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
        std::string indivID, poplnID, sexStr;
        if (!(iss >> indivID >> poplnID)) {
            std::cerr << "\nerror: invalid format in metadata file at line " << lineNum
                      << ": expected 'INDIV POPLN [SEX]'\n";
            exit(1);
        }
        indivToPopln[indivID] = poplnID;

        // Optional third column: sex
        if (iss >> sexStr) {
            std::string s = sexStr;
            for (char &c : s) c = (char)tolower((unsigned char)c);
            if (s == "f" || s == "female") {
                indivToSex[indivID] = SEX_FEMALE;
            } else if (s == "m" || s == "male") {
                indivToSex[indivID] = SEX_MALE;
            } else {
                std::cerr << "\nerror: unrecognized sex '" << sexStr << "' in metadata file at line "
                          << lineNum << ": expected F/M (or female/male)\n";
                exit(1);
            }
        }
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

// Identify X-linked loci from the VCF CHROM field. Matches "X" or "chrX"
// (case-insensitive); everything else is treated as autosomal.
static bool isXLinkedChrom(const char *chrom) {
    if (!chrom) return false;
    std::string c(chrom);
    for (char &ch : c) ch = (char)tolower((unsigned char)ch);
    if (c.rfind("chr", 0) == 0) c = c.substr(3);  // strip a leading "chr"
    return (c == "x");
}

// Read VCF file and populate genotype data
void readVCFFile(const char *vcfFileName,
                 const std::map<std::string, std::string> &indivToPopln,
                 const std::map<std::string, int> &indivToSex,
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

            // Sex from metadata (SEX_UNKNOWN if not provided)
            auto sIt = indivToSex.find(sampleName);
            sampleIndiv[indivIdx].sex = (sIt != indivToSex.end())
                                        ? (unsigned int)sIt->second : SEX_UNKNOWN;

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

    // Per-locus inheritance type, filled as loci are read (noLoci from the
    // sizing pass in checkVCFDataSize).
    gLocusType.assign(noLoci, LOCUS_AUTOSOMAL);

    while (bcf_read(vcfFile, hdr, rec) == 0) {
        bcf_unpack(rec, BCF_UN_ALL);

        const char *chromName = bcf_hdr_id2name(hdr, rec->rid);

        // Create locus name from CHROM:POS or ID
        std::string locusName;
        if (rec->d.id && strcmp(rec->d.id, ".") != 0) {
            locusName = rec->d.id;
        } else {
            locusName = std::string(chromName) + ":" + std::to_string(rec->pos + 1);
        }
        locusIDMap[locusName] = locusIdx;

        // Classify locus as X-linked or autosomal from the CHROM field
        bool locusIsX = isXLinkedChrom(chromName);
        if (locusIdx < gLocusType.size())
            gLocusType[locusIdx] = locusIsX ? LOCUS_XLINKED : LOCUS_AUTOSOMAL;

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
            // Max ploidy across samples for this record; individual samples may
            // be haploid, in which case their second slot is a vector-end marker
            // (e.g. hemizygous males at an X-linked locus).
            int maxPloidy = ngt_ret / nSamples;

            for (int s = 0; s < nSamples; s++) {
                int indIdx = sampleToIndiv[s];
                if (indIdx < 0) continue;  // Skip samples not in metadata

                int32_t *ptr = gt + s * maxPloidy;

                // Is a real second gene copy present for this sample?
                bool diploidCall = (maxPloidy >= 2) &&
                                   (ptr[1] != bcf_int32_vector_end);

                if (diploidCall) {
                    if (bcf_gt_is_missing(ptr[0]) || bcf_gt_is_missing(ptr[1])) {
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
                    // Single gene copy for this sample at this locus.
                    if (bcf_gt_is_missing(ptr[0])) {
                        sampleIndiv[indIdx].genotype[locusIdx][0] = -1;
                        sampleIndiv[indIdx].genotype[locusIdx][1] = -1;
                    } else {
                        int a0 = bcf_gt_allele(ptr[0]);
                        sampleIndiv[indIdx].genotype[locusIdx][0] = a0;
                        // Hemizygous on X (true single copy); otherwise a haploid
                        // call on a non-X locus is treated as homozygous, as before.
                        sampleIndiv[indIdx].genotype[locusIdx][1] =
                            locusIsX ? HEMIZYGOUS : (GenotypeType)a0;
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
// Collapsed sampler (integrating out allele frequencies), Phase 1: count tables
// and the Dirichlet-multinomial marginal. See doc/collapsed_allele_freqs_plan.md.
//=============================================================================

// Tally gene copies from the current ALL-NATIVE state (every individual's copies
// assigned to its sampled population) into cnt[p][l][a] and cntN[p][l]. A diploid
// genotype contributes two copies; a hemizygous male X (second slot HEMIZYGOUS)
// or a missing copy (value < 0) contributes one/none.
static void initCountsNative(long ***cnt, long **cntN, indiv *ind, unsigned int noIndiv,
                             unsigned int noLoci, unsigned int *noAlleles, unsigned int noPopln)
{
	for (unsigned int p = 0; p < noPopln; p++)
		for (unsigned int l = 0; l < noLoci; l++)
		{
			cntN[p][l] = 0;
			unsigned int A = noAlleles[l] ? noAlleles[l] : 1;
			for (unsigned int a = 0; a < A; a++) cnt[p][l][a] = 0;
		}
	for (unsigned int i = 0; i < noIndiv; i++)
	{
		unsigned int p = ind[i].samplePopln;   // initial state: all native
		const GenotypeType (*g)[2] = ind[i].genotype;
		for (unsigned int l = 0; l < noLoci; l++)
		{
			int a0 = g[l][0], a1 = g[l][1];
			if (a0 >= 0) { cnt[p][l][a0]++; cntN[p][l]++; }
			if (a1 >= 0) { cnt[p][l][a1]++; cntN[p][l]++; }
		}
	}
}

// Dirichlet-multinomial (Polya) log-marginal of the genotype counts, integrating
// out the population allele frequencies under a symmetric Dirichlet(alpha) prior.
static double collapsedLogLik(long ***cnt, long **cntN, unsigned int noPopln,
                              unsigned int noLoci, unsigned int *noAlleles, double alpha)
{
	double lp = 0.0;
	for (unsigned int p = 0; p < noPopln; p++)
		for (unsigned int l = 0; l < noLoci; l++)
		{
			unsigned int A = noAlleles[l];
			if (A == 0) continue;
			long N = cntN[p][l];
			lp += lgamma(A * alpha) - lgamma(A * alpha + (double) N);
			for (unsigned int a = 0; a < A; a++)
				lp += lgamma(alpha + (double) cnt[p][l][a]) - lgamma(alpha);
		}
	return lp;
}

// --- Incremental engine (Phase 2). O(1) per gene copy, no lgamma. ---
// Probability factor of adding one copy of allele a to cell (p,l): the amount by
// which the marginal changes is log of this (evaluated on the pre-add counts).
static inline double addRatio(long ***cnt, long **cntN, unsigned int p,
                              unsigned int l, int a, unsigned int A, double alpha)
{
	return (cnt[p][l][a] + alpha) / (cntN[p][l] + A * alpha);
}
static inline void addCopy(long ***cnt, long **cntN, unsigned int p, unsigned int l, int a)
{ cnt[p][l][a]++; cntN[p][l]++; }
static inline void removeCopy(long ***cnt, long **cntN, unsigned int p, unsigned int l, int a)
{ cnt[p][l][a]--; cntN[p][l]--; }

// Population an individual's gene copies are drawn from for age 0 (native) and
// age 1 (migrant); both copies come from this one population.
static inline unsigned int copyPop(const indiv& ind)
{ return (ind.migrantAge == 1) ? ind.migrantPopln : ind.samplePopln; }

// Per-locus phase/assignment codes stored for age-2 individuals so a later
// remove can undo exactly what an add committed:
//   0 = none (age 0/1, or age-2 homozygote: deterministic)
//   1 = age-2 het, a0 -> source (migrantPopln), a1 -> native (samplePopln)
//   2 = age-2 het, a1 -> source,                a0 -> native
//   3 = age-2 hemizygous male X, single copy -> source
//   4 = age-2 hemizygous male X, single copy -> native

// z-marginalized predictive of a same-population diploid genotype (a0,a1) drawn
// from pop p at locus l, given inbreeding F. Read-only (analytic; equals the
// Phase-3 two-copy count-marginal when F = 0). Homozygote: F/(1-F) mixes a single
// IBD copy vs two independent copies; het: always outbred.
static double genoPred(long ***cnt, long **cntN, unsigned int p, unsigned int l,
                       int a0, int a1, unsigned int A, double alpha, double F)
{
	double N   = cntN[p][l] + A * alpha;
	double Np1 = cntN[p][l] + 1 + A * alpha;
	if (a0 == a1)
	{
		double r  = (cnt[p][l][a0] + alpha) / N;
		double r2 = (cnt[p][l][a0] + 1 + alpha) / Np1;
		return F * r + (1.0 - F) * r * r2;
	}
	return (1.0 - F) * (cnt[p][l][a0] + alpha) / N * (cnt[p][l][a1] + alpha) / Np1;
}

// Precomputed table of natural logs of small non-negative integers (gLogTab[k] =
// log(k); gLogTab[0] unused). With the Dirichlet concentration fixed at
// ALLELE_PRIOR_ALPHA = 1, every count ratio in the collapsed hot path is a ratio
// of integers, so its log is a difference of table entries -- eliminating the bulk
// of the libm log() calls that dominate the collapsed sampler's per-locus cost
// (profiling: ~36% of samples). Sized to cover 2*noIndiv + maxAlleles + slack.
static double *gLogTab = nullptr;
static long gLogTabSize = 0;
static void initLogTab(long maxArg)
{
	gLogTabSize = maxArg + 1;
	gLogTab = new double[gLogTabSize];
	gLogTab[0] = 0.0;   // never indexed (all arguments are >= 1)
	for (long k = 1; k < gLogTabSize; k++) gLogTab[k] = log((double) k);
}

// Memoized homozygote IBD-mixture log, log(F + (1-F)*(n+2)/D), keyed by
// (population, homozygous count n, denominator D = N+1+A). This is the one real
// log() left in the collapsed hot path (F > 0); F is constant between inbreeding
// updates, so the value depends only on (n, D) within an epoch and is shared by
// every locus whose homozygous count state matches -- pooling the log across loci.
// Filled lazily; an O(1) epoch bump invalidates the whole table when any F changes.
static double **gHomLog = nullptr;
static long   **gHomEpoch = nullptr;
static long gHomDDim = 0, gHomEpochCur = 1;
static void initHomLog(unsigned int noPopln, long nDim, long dDim)
{
	gHomDDim = dDim;
	gHomLog   = new double*[noPopln];
	gHomEpoch = new long*[noPopln];
	for (unsigned int p = 0; p < noPopln; p++)
	{
		gHomLog[p]   = new double[nDim * dDim];
		gHomEpoch[p] = new long[nDim * dDim];
		for (long i = 0; i < nDim * dDim; i++) gHomEpoch[p][i] = 0;   // != gHomEpochCur
	}
}

// log(addRatio) via the integer table (assumes alpha == 1): log((n+1)/(N+A)).
static inline double logAddRatioT(long ***cnt, long **cntN, unsigned int p,
                                  unsigned int l, int a, unsigned int A)
{ return gLogTab[cnt[p][l][a] + 1] - gLogTab[cntN[p][l] + A]; }

// log(genoPred) via the table (assumes alpha == 1). The heterozygote is fully
// tabular; the homozygote keeps one real log for the IBD mixture F + (1-F)*r2,
// which vanishes when F == 0 (then it too is tabular).
static inline double logGenoPredT(long ***cnt, long **cntN, unsigned int p, unsigned int l,
                                  int a0, int a1, unsigned int A, double F, double log1mF)
{
	long N = cntN[p][l];
	if (a0 == a1)
	{
		long n = cnt[p][l][a0];
		double logr = gLogTab[n + 1] - gLogTab[N + A];
		if (F == 0.0)
			return logr + gLogTab[n + 2] - gLogTab[N + 1 + A];
		// IBD-mixture log, memoized by (n, D): shared across loci with the same
		// homozygous count state; recomputed only when the epoch (F) changed.
		long D = N + 1 + (long) A;
		long idx = n * gHomDDim + D;
		if (gHomEpoch[p][idx] != gHomEpochCur)
		{
			gHomLog[p][idx] = log(F + (1.0 - F) * (double)(n + 2) / (double) D);
			gHomEpoch[p][idx] = gHomEpochCur;
		}
		return logr + gHomLog[p][idx];
	}
	return log1mF + gLogTab[cnt[p][l][a0] + 1] - gLogTab[N + A]
	              + gLogTab[cnt[p][l][a1] + 1] - gLogTab[N + 1 + A];
}

// Change in the collapsed log-marginal from ADDING individual i's gene copies
// under its ancestry, given the current (without-i) counts. Read-only. Age 0/1
// diploid genotypes marginalize the IBD indicator (genoPred, inbreeding F);
// age-2 marginalizes the latent phase (a priori 1/2 each).
static double computeAddLogProb(long ***cnt, long **cntN, const indiv& ind,
                                unsigned int noLoci, unsigned int *noAlleles, double alpha,
                                double *FStat)
{
	const GenotypeType (*g)[2] = ind.genotype;
	double dl = 0.0;
	if (ind.migrantAge != 2)
	{
		unsigned int p = copyPop(ind);
		double F = FStat[p];
		double log1mF = (F > 0.0) ? log(1.0 - F) : 0.0;   // once per call (het factor)
		for (unsigned int l = 0; l < noLoci; l++)
		{
			unsigned int A = noAlleles[l]; if (A == 0) continue;
			int a0 = g[l][0], a1 = g[l][1];
			if (a1 == HEMIZYGOUS || (a0 >= 0 && a1 < 0))       // hemizygous male X: single copy
			{
				if (a0 >= 0) dl += logAddRatioT(cnt, cntN, p, l, a0, A);
			}
			else if (a0 >= 0 && a1 >= 0)
				dl += logGenoPredT(cnt, cntN, p, l, a0, a1, A, F, log1mF);
		}
	}
	else   // age 2: one copy from source (mig), one from native (nat); distinct cells
	{
		unsigned int mig = ind.migrantPopln, nat = ind.samplePopln;
		for (unsigned int l = 0; l < noLoci; l++)
		{
			unsigned int A = noAlleles[l]; if (A == 0) continue;
			int a0 = g[l][0], a1 = g[l][1];
			if (a1 == HEMIZYGOUS)            // hemizygous male X: single copy
			{
				if (a0 < 0) continue;
				if (ind.migrantSex == SEX_FEMALE)      dl += logAddRatioT(cnt, cntN, mig, l, a0, A);
				else if (ind.migrantSex == SEX_MALE)   dl += logAddRatioT(cnt, cntN, nat, l, a0, A);
				else
				{
					double rm = addRatio(cnt, cntN, mig, l, a0, A, alpha);
					double rn = addRatio(cnt, cntN, nat, l, a0, A, alpha);
					dl += log(0.5 * (rm + rn));      // marginal over parent sex: log of a sum
				}
			}
			else if (a0 >= 0 && a1 >= 0)
			{
				if (a0 == a1)
					dl += logAddRatioT(cnt, cntN, mig, l, a0, A)
					    + logAddRatioT(cnt, cntN, nat, l, a0, A);
				else
				{
					double w0 = addRatio(cnt, cntN, mig, l, a0, A, alpha) * addRatio(cnt, cntN, nat, l, a1, A, alpha);
					double w1 = addRatio(cnt, cntN, mig, l, a1, A, alpha) * addRatio(cnt, cntN, nat, l, a0, A, alpha);
					dl += log(0.5 * (w0 + w1));       // phase marginal: log of a sum
				}
			}
			else if (a0 >= 0)
				dl += log(0.5 * (addRatio(cnt, cntN, mig, l, a0, A, alpha)
				               + addRatio(cnt, cntN, nat, l, a0, A, alpha)));
		}
	}
	return dl;
}

// Permanently ADD individual i's gene copies, sampling the latent per-locus
// codes (stored in assign[]) so removeIndividual can undo it: for age 0/1 a
// homozygous diploid locus draws its IBD indicator (assign=1 if IBD, one copy;
// assign=0 if outbred, two copies); for age 2, the phase (codes 1-4).
static double addIndividual(long ***cnt, long **cntN, const indiv& ind, unsigned char *assign,
                            unsigned int noLoci, unsigned int *noAlleles, double alpha,
                            double *FStat)
{
	const GenotypeType (*g)[2] = ind.genotype;
	// Commit-only: this samples the latent per-locus codes and updates the counts;
	// the add-log-prob is computed separately by computeAddLogProb, so no log() work
	// is done here (the return value is unused at both call sites). Ratios are still
	// formed where a sampling decision needs them.
	if (ind.migrantAge != 2)
	{
		unsigned int p = copyPop(ind);
		double F = FStat[p];
		for (unsigned int l = 0; l < noLoci; l++)
		{
			assign[l] = 0;
			unsigned int A = noAlleles[l]; if (A == 0) continue;
			int a0 = g[l][0], a1 = g[l][1];
			if (a1 == HEMIZYGOUS || (a0 >= 0 && a1 < 0))       // hemizygous male X: single copy
			{
				if (a0 >= 0) addCopy(cnt, cntN, p, l, a0);
			}
			else if (a0 >= 0 && a1 >= 0)
			{
				if (a0 == a1)   // homozygote: sample the IBD indicator
				{
					double rr  = addRatio(cnt, cntN, p, l, a0, A, alpha);
					double rr2 = (cnt[p][l][a0] + 1 + alpha) / (cntN[p][l] + 1 + A * alpha);
					double pz1 = F * rr, pz0 = (1.0 - F) * rr * rr2;
					if (gsl_rng_uniform(r) < pz1 / (pz1 + pz0))
					{ addCopy(cnt, cntN, p, l, a0); assign[l] = 1; }          // IBD: one copy
					else
					{ addCopy(cnt, cntN, p, l, a0); addCopy(cnt, cntN, p, l, a0); }  // outbred: two copies
				}
				else            // het: always outbred, two copies
				{
					addCopy(cnt, cntN, p, l, a0);
					addCopy(cnt, cntN, p, l, a1);
				}
			}
		}
	}
	else
	{
		unsigned int mig = ind.migrantPopln, nat = ind.samplePopln;
		for (unsigned int l = 0; l < noLoci; l++)
		{
			assign[l] = 0;
			unsigned int A = noAlleles[l]; if (A == 0) continue;
			int a0 = g[l][0], a1 = g[l][1];
			if (a1 == HEMIZYGOUS || (a0 >= 0 && a1 < 0))     // single copy (hemizygous male X)
			{
				if (a0 < 0) continue;
				bool toSrc;
				if (ind.migrantSex == SEX_FEMALE)      toSrc = true;
				else if (ind.migrantSex == SEX_MALE)   toSrc = false;
				else
				{
					double rm = addRatio(cnt, cntN, mig, l, a0, A, alpha);
					double rn = addRatio(cnt, cntN, nat, l, a0, A, alpha);
					toSrc = (gsl_rng_uniform(r) < rm / (rm + rn));
				}
				if (toSrc) { addCopy(cnt, cntN, mig, l, a0); assign[l] = 3; }
				else       { addCopy(cnt, cntN, nat, l, a0); assign[l] = 4; }
			}
			else if (a0 >= 0 && a1 >= 0)
			{
				if (a0 == a1)
				{
					addCopy(cnt, cntN, mig, l, a0);
					addCopy(cnt, cntN, nat, l, a0);
				}
				else
				{
					double w0 = addRatio(cnt, cntN, mig, l, a0, A, alpha) * addRatio(cnt, cntN, nat, l, a1, A, alpha);
					double w1 = addRatio(cnt, cntN, mig, l, a1, A, alpha) * addRatio(cnt, cntN, nat, l, a0, A, alpha);
					if (gsl_rng_uniform(r) < w0 / (w0 + w1))
					{ addCopy(cnt, cntN, mig, l, a0); addCopy(cnt, cntN, nat, l, a1); assign[l] = 1; }
					else
					{ addCopy(cnt, cntN, mig, l, a1); addCopy(cnt, cntN, nat, l, a0); assign[l] = 2; }
				}
			}
		}
	}
	return 0.0;
}

// Permanently REMOVE individual i's gene copies, undoing exactly what
// addIndividual committed (using the stored age-2 phase codes).
static void removeIndividual(long ***cnt, long **cntN, const indiv& ind,
                             unsigned char *assign, unsigned int noLoci)
{
	const GenotypeType (*g)[2] = ind.genotype;
	if (ind.migrantAge != 2)
	{
		unsigned int p = copyPop(ind);
		for (unsigned int l = 0; l < noLoci; l++)
		{
			int a0 = g[l][0], a1 = g[l][1];
			if (a1 == HEMIZYGOUS || (a0 >= 0 && a1 < 0))       // hemizygous male X: one copy
			{
				if (a0 >= 0) removeCopy(cnt, cntN, p, l, a0);
			}
			else if (a0 >= 0 && a1 >= 0)
			{
				if (a0 == a1)                    // homozygote: 1 copy if IBD, else 2
				{
					removeCopy(cnt, cntN, p, l, a0);
					if (assign[l] != 1) removeCopy(cnt, cntN, p, l, a0);
				}
				else { removeCopy(cnt, cntN, p, l, a0); removeCopy(cnt, cntN, p, l, a1); }
			}
		}
	}
	else
	{
		unsigned int mig = ind.migrantPopln, nat = ind.samplePopln;
		for (unsigned int l = 0; l < noLoci; l++)
		{
			int a0 = g[l][0], a1 = g[l][1];
			switch (assign[l])
			{
				case 1: removeCopy(cnt, cntN, mig, l, a0); removeCopy(cnt, cntN, nat, l, a1); break;
				case 2: removeCopy(cnt, cntN, mig, l, a1); removeCopy(cnt, cntN, nat, l, a0); break;
				case 3: removeCopy(cnt, cntN, mig, l, a0); break;
				case 4: removeCopy(cnt, cntN, nat, l, a0); break;
				default:  // age-2 homozygote: one copy to each population
					if (a1 != HEMIZYGOUS && a0 >= 0 && a1 >= 0)
					{ removeCopy(cnt, cntN, mig, l, a0); removeCopy(cnt, cntN, nat, l, a0); }
					break;
			}
		}
	}
}

// Inbreeding update (collapsed sampler). Gibbs-resample the IBD indicator of
// every same-population diploid homozygous genotype (age 0/1), keeping the count
// tables consistent, then draw each population's F from its Beta full-conditional
// F_p ~ Beta(1 + #IBD, 1 + #diploid - #IBD). Heterozygotes are outbred (z fixed 0)
// but still count as trials. O(noIndiv * noLoci).
static void inbreedingUpdate(long ***cnt, long **cntN, indiv *ind, unsigned char **assign,
                             unsigned int noIndiv, unsigned int noLoci, unsigned int *noAlleles,
                             unsigned int noPopln, double alpha, double *FStat)
{
	std::vector<long> zSum(noPopln, 0), zTrials(noPopln, 0);
	for (unsigned int i = 0; i < noIndiv; i++)
	{
		if (ind[i].migrantAge == 2) continue;      // hybrids are outbred (no F)
		unsigned int p = copyPop(ind[i]);
		double F = FStat[p];
		const GenotypeType (*g)[2] = ind[i].genotype;
		for (unsigned int l = 0; l < noLoci; l++)
		{
			unsigned int A = noAlleles[l]; if (A == 0) continue;
			int a0 = g[l][0], a1 = g[l][1];
			if (a1 == HEMIZYGOUS || a0 < 0 || a1 < 0) continue;   // not a diploid genotype
			zTrials[p]++;
			if (a0 != a1) { assign[i][l] = 0; continue; }          // het: outbred
			// homozygote: remove its current contribution, resample IBD, re-add
			removeCopy(cnt, cntN, p, l, a0);
			if (assign[i][l] != 1) removeCopy(cnt, cntN, p, l, a0);
			double rr  = (cnt[p][l][a0] + alpha) / (cntN[p][l] + A * alpha);
			double rr2 = (cnt[p][l][a0] + 1 + alpha) / (cntN[p][l] + 1 + A * alpha);
			double pz1 = F * rr, pz0 = (1.0 - F) * rr * rr2;
			if (gsl_rng_uniform(r) < pz1 / (pz1 + pz0))
			{ addCopy(cnt, cntN, p, l, a0); assign[i][l] = 1; zSum[p]++; }
			else
			{ addCopy(cnt, cntN, p, l, a0); addCopy(cnt, cntN, p, l, a0); assign[i][l] = 0; }
		}
	}
	for (unsigned int p = 0; p < noPopln; p++)
		FStat[p] = gsl_ran_beta(r, 1.0 + (double) zSum[p], 1.0 + (double)(zTrials[p] - zSum[p]));
}

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
	gArgs.collapse = 0;  // collapsed (integrated-frequency) sampler off by default
	gArgs.fixGamma = 0;  // estimate gamma by default; --fixgamma pins gamma=1
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

			case 'c':
				gArgs.collapse = 1;
				break;

			case 'G':
				gArgs.fixGamma = 1;
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
				std::cout << "    -c, --collapse    integrate out allele frequencies (collapsed sampler; faster)\n";
				std::cout << "    -G, --fixgamma    fix gamma=1 (original BayesAss: tau=2, migration cap 1/3)\n";
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

	// The collapsed (integrated-allele-frequency) sampler analytically marginalizes
	// the population allele frequencies (Dirichlet-multinomial); see
	// doc/collapsed_allele_freqs_plan.md. It supports the base migration model, the
	// inbreeding (F) and sex-biased dispersal extensions, and missing data.
	if (gArgs.collapse)
		std::cerr << "\nnote: --collapse active; allele frequencies are integrated out "
		             "(collapsed sampler).\n";

	/* get input file name or validate VCF options */
	std::map<std::string, std::string> indivToPopln;  // For VCF mode
	std::map<std::string, int> indivToSex;            // For VCF mode (optional sex column)

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
		readMetadataFile(gArgs.metaFileName, indivToPopln, indivToSex);

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
		readVCFFile(gArgs.vcfFileName, indivToPopln, indivToSex, sampleIndiv, noIndiv, noLoci, noPopln, noAlleles);

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

	// Ensure per-locus type is sized (native text input has no X annotation, so
	// all loci are autosomal; the VCF path has already filled this).
	if (gLocusType.size() != noLoci)
		gLocusType.assign(noLoci, LOCUS_AUTOSOMAL);

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
		       // The standard sampler imputes missing genotypes (data-augmentation
		       // MH below). The collapsed sampler instead leaves them as -1 so the
		       // count tables simply omit the missing gene copies -- the exact
		       // Dirichlet-multinomial marginalization -- since its missing-data MH
		       // is disabled and initCountsNative/computeAddLogProb skip a < 0.
		       if (!gArgs.collapse)
			 {
			   sampleIndiv[i].genotype[j][0] = gsl_rng_uniform_int(r, noAlleles[j]);
			   sampleIndiv[i].genotype[j][1] = gsl_rng_uniform_int(r, noAlleles[j]);
			 }
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

	// Sex-biased dispersal model setup. Active only when sex metadata was
	// supplied (indivToSex non-empty). Requires complete sex information.
	// The sex-biased model composes with the collapsed sampler (Phase 6): the
	// count-table machinery is sex-aware (age-2 male-X source/native copy
	// assignment via migrantSex), and the sigma Gibbs has a collapsed branch.
	bool sexBiasModel = !indivToSex.empty();
	// Sex ratios (female fractions): movement dispersal (age-1 migrants), effective
	// gene flow (age-2 migrant parents = migrants that also bred), and residents.
	double phiMove = 0.5;             // female fraction of age-1 migrants (movement)
	double phiBreed = 0.5;            // female fraction of age-2 migrant parents (gene flow)
	double rho = 0.5;                 // female fraction of residents (non-dispersers)
	double gamma = 1.0;               // migrant breeding-success multiplier (tau = 2*gamma)
	std::vector<unsigned int> xLoci;  // indices of X-linked loci (for phiBreed/sigma updates)
	if (sexBiasModel)
	{
		for (unsigned int i = 0; i < noIndiv; i++)
			if (sampleIndiv[i].sex != SEX_FEMALE && sampleIndiv[i].sex != SEX_MALE)
			{
				std::cerr << "\nerror: sex-biased model enabled but individual " << i
				          << " has no sex in the metadata file.\n"
				          << "       Provide a sex column (F/M) for every individual, or none.\n";
				exit(1);
			}
		for (unsigned int l = 0; l < noLoci; l++)
			if (l < gLocusType.size() && gLocusType[l] == LOCUS_XLINKED)
				xLoci.push_back(l);
		if (gArgs.verbose)
			std::cout << "Sex-biased dispersal model: ON (phi = female fraction of migrants), "
			          << xLoci.size() << " X-linked loci\n\n";
	}


	if (gArgs.verbose)
	{
		cout << "\nInput file: " << infileName;
		cout << "\nOutput file: " << gArgs.outfileName << "\n";
		cout << "Individuals: " << noIndiv << " Populations: " << noPopln << " Loci: " << noLoci;
		cout << " Missing genotypes: " << noMissingGenotypes << "\n";

		// Report parsed sex and X-linked / hemizygous data (sex-biased dispersal)
		{
			unsigned int nFemale = 0, nMale = 0, nUnknownSex = 0;
			for (unsigned int i = 0; i < noIndiv; i++)
			{
				if (sampleIndiv[i].sex == SEX_FEMALE) nFemale++;
				else if (sampleIndiv[i].sex == SEX_MALE) nMale++;
				else nUnknownSex++;
			}
			unsigned int nXLoci = 0;
			for (unsigned int l = 0; l < noLoci; l++)
				if (l < gLocusType.size() && gLocusType[l] == LOCUS_XLINKED) nXLoci++;
			long int nHemizygous = 0;
			for (unsigned int i = 0; i < noIndiv; i++)
				for (unsigned int l = 0; l < noLoci; l++)
					if (sampleIndiv[i].genotype[l][1] == HEMIZYGOUS) nHemizygous++;

			cout << "Sex: " << nFemale << " female, " << nMale << " male";
			if (nUnknownSex) cout << ", " << nUnknownSex << " unspecified";
			cout << "\n";
			cout << "Loci: " << (noLoci - nXLoci) << " autosomal, " << nXLoci
			     << " X-linked; hemizygous male-X genotypes: " << nHemizygous << "\n";
		}
		cout << "\n";
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

	// Posterior accumulators for the movement (phiMove), gene-flow (phiBreed), and
	// resident (rho) female fractions, plus the Savage-Dickey tests.
	double avgMove = 0.0, varMove = 0.0;
	double avgBreed = 0.0, varBreed = 0.0;
	double avgRho = 0.0, varRho = 0.0;
	long int moveGtRhoCount = 0, breedGtRhoCount = 0;   // P(phi > rho) counters
	// Latest sex counts (set by the Gibbs update each iteration).
	long int moveCountF = 0, moveCountM = 0;    // age-1 migrants (movement)
	long int breedCountF = 0, breedCountM = 0;  // age-2 migrant parents (gene flow)
	long int rhoCountF = 0, rhoCountM = 0;      // residents (non-dispersers)
	// RB Savage-Dickey overlap accumulators for the three equality nulls.
	double sdMoveEqRhoSum = 0.0;    // H0: phiMove  = rho
	double sdBreedEqRhoSum = 0.0;   // H0: phiBreed = rho
	double sdMoveEqBreedSum = 0.0;  // H0: phiMove  = phiBreed
	long int sdPhiNSamples = 0;
	// Migrant breeding-success multiplier gamma (posterior mean/var; P(gamma<1)).
	double avgGamma = 0.0, varGamma = 0.0;
	long int gammaLtOneCount = 0, gammaNSamples = 0;

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
				migrationRates[l][k] = (1.0/(1.0+2.0*gamma))*(1.0/noPopln);
		migrationRates[l][l] = 1.0-(1.0/(1.0+2.0*gamma))*((noPopln-1.0)/noPopln);
		migrationRates[l][noPopln] = (1.0/(1.0+2.0*gamma))*((noPopln-1.0)/noPopln);
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
	{	sampleIndiv[l].migrantAge=0; sampleIndiv[l].migrantPopln=sampleIndiv[l].samplePopln;
		sampleIndiv[l].migrantSex=SEX_UNKNOWN; }
	fillMigrantCounts(sampleIndiv,migrantCounts,noIndiv,noPopln);

	// Per-category index lists for O(1) uniform selection of an individual in a
	// (samplePopln, migrantPopln, migrantAge) category, replacing the old per-
	// iteration O(N) gsl_ran_shuffle + scan. catOf/catPos give O(1) removal
	// (swap-with-last); migrantCounts is then maintained incrementally alongside,
	// retiring the O(N) fillMigrantCounts after each accepted ancestry move.
	auto catIndex = [noPopln](unsigned int sp, unsigned int mp, unsigned int age) -> int
	{ return (int)((sp * (unsigned int)noPopln + mp) * 3 + age); };
	std::vector<std::vector<int> > catList((size_t)noPopln * noPopln * 3);
	std::vector<int> catOf(noIndiv), catPos(noIndiv);
	for (unsigned int i = 0; i < noIndiv; i++)
	{
		int c = catIndex(sampleIndiv[i].samplePopln, sampleIndiv[i].migrantPopln, sampleIndiv[i].migrantAge);
		catOf[i] = c; catPos[i] = (int)catList[c].size(); catList[c].push_back((int)i);
	}
	auto catRemove = [&](int i)
	{
		int c = catOf[i], pos = catPos[i], last = catList[c].back();
		catList[c][pos] = last; catPos[last] = pos; catList[c].pop_back();
	};
	auto catAdd = [&](int i, int c)
	{ catOf[i] = c; catPos[i] = (int)catList[c].size(); catList[c].push_back(i); };

	for(unsigned int i = 0; i < noIndiv; i++)
		sampleIndiv[i].logL = logLik(sampleIndiv[i],alleleFreqs,logAlleleFreqs,FStat,log1MinusFStat,noLoci);

	// Collapsed sampler: persistent count tables (gCount/gCountN) and per-age-2
	// phase store (gAssign), initialised from the all-native state. These are
	// maintained through the MCMC (Phase 3). Also runs the Phase 1/2 validation
	// gates. All null / unused unless --collapse.
	long ***gCount = nullptr;
	long **gCountN = nullptr;
	unsigned char **gAssign = nullptr;
	if(gArgs.collapse)
	{
		gCount = new long**[noPopln];
		gCountN = new long*[noPopln];
		for(unsigned int p = 0; p < noPopln; p++)
		{
			gCount[p] = new long*[noLoci];
			gCountN[p] = new long[noLoci];
			for(unsigned int l = 0; l < noLoci; l++)
				gCount[p][l] = new long[noAlleles[l] > 0 ? noAlleles[l] : 1];
		}
		gAssign = new unsigned char*[noIndiv];
		for(unsigned int i = 0; i < noIndiv; i++)
		{
			gAssign[i] = new unsigned char[noLoci];
			for(unsigned int l = 0; l < noLoci; l++) gAssign[i][l] = 0;
		}
		// log table for the collapsed hot path: max index is cntN(<=2*noIndiv) + A
		// + 1, plus a homozygote's cnt+2; +8 slack covers all count-ratio arguments.
		initLogTab(2L * (long)noIndiv + (long)maxAlleles + 8);
		// homozygote IBD-mixture memo: n in [0, 2*noIndiv], D = N+1+A in
		// [0, 2*noIndiv + maxAlleles + 1]; +3 slack on each dimension.
		initHomLog(noPopln, 2L * (long)noIndiv + 3, 2L * (long)noIndiv + (long)maxAlleles + 3);
		initCountsNative(gCount, gCountN, sampleIndiv, noIndiv, noLoci, noAlleles, noPopln);
		double cll = collapsedLogLik(gCount, gCountN, noPopln, noLoci, noAlleles, ALLELE_PRIOR_ALPHA);
		std::cout << "collapsed init logL (Dirichlet-multinomial, alpha="
		          << ALLELE_PRIOR_ALPHA << "): " << std::setprecision(8) << cll << "\n";

		// Phase-2 gate (init state is all-native, so gAssign is unused here).
		{
			double base = cll, maxDiff = 0.0;
			unsigned int cap = noIndiv < 200 ? noIndiv : 200;   // bound the O(P*L) recomputes
			for(unsigned int i = 0; i < cap; i++)
			{
				removeIndividual(gCount, gCountN, sampleIndiv[i], gAssign[i], noLoci);
				double without = collapsedLogLik(gCount, gCountN, noPopln, noLoci, noAlleles, ALLELE_PRIOR_ALPHA);
				double dadd = computeAddLogProb(gCount, gCountN, sampleIndiv[i], noLoci, noAlleles, ALLELE_PRIOR_ALPHA, FStat);
				double d = std::fabs(base - (without + dadd));
				if(d > maxDiff) maxDiff = d;
				addIndividual(gCount, gCountN, sampleIndiv[i], gAssign[i], noLoci, noAlleles, ALLELE_PRIOR_ALPHA, FStat);
			}
			double restored = collapsedLogLik(gCount, gCountN, noPopln, noLoci, noAlleles, ALLELE_PRIOR_ALPHA);
			std::cout << "collapsed Phase-2 gate: max|addLogProb - recompute| = "
			          << std::scientific << maxDiff << ", restored logL diff = "
			          << std::fabs(restored - base) << std::fixed
			          << "  -> " << ((maxDiff < 1e-6 && std::fabs(restored - base) < 1e-6) ? "PASS" : "FAIL")
			          << "\n";
		}
	}

	if(gArgs.debug)
	{
		std::cout << "mCLp: " << migCountLogProb(migrantCounts,migrationRates,noPopln,gamma) << "\n";
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
		// O(1) uniform pick among individuals in the proposed drop category
		// (proposeMigrantAncDrop guarantees a non-empty category, so size > 0).
		{
			const std::vector<int>& cl = catList[catIndex(samplePopln, migrantPopln, migrantAge)];
			chosenIndiv = cl[gsl_rng_uniform_int(r, cl.size())];
		}
		proposeMigrantAncAdd(migrantPopAdd, migrantAgeAdd,migrantPopln, migrantAge, samplePopln, noPopln);

		tempIndiv.samplePopln = sampleIndiv[chosenIndiv].samplePopln;
		tempIndiv.migrantPopln = migrantPopAdd;
		tempIndiv.migrantAge = migrantAgeAdd;
		tempIndiv.sex = sampleIndiv[chosenIndiv].sex;
		for(unsigned int j = 0; j < noLoci; j++)
		{
			tempIndiv.genotype[j][0] = sampleIndiv[chosenIndiv].genotype[j][0];
			tempIndiv.genotype[j][1] = sampleIndiv[chosenIndiv].genotype[j][1];
		}

		// Sex of the first-generation migrant for the proposed ancestry:
		//   age 1 -> the individual is the migrant, so its observed sex;
		//   age 2 -> latent migrant-parent sex, proposed from its prior
		//            Bernoulli(phi). Drawing sigma from the prior makes the
		//            sigma prior and proposal densities cancel in the MH ratio.
		double dtLogSex = 0.0;
		if (sexBiasModel)
		{
			if (tempIndiv.migrantAge == 1)
				tempIndiv.migrantSex = tempIndiv.sex;
			else if (tempIndiv.migrantAge == 2)
				tempIndiv.migrantSex = (gsl_rng_uniform(r) < phiBreed) ? SEX_FEMALE : SEX_MALE;
			else
				tempIndiv.migrantSex = SEX_UNKNOWN;

			// After the sigma prior/proposal cancellation, the residual per-
			// individual sex term is: age-1 own sex ~ phiMove (disperser/mover);
			// age-0/age-2 own sex ~ rho (resident). The age-2 latent sigma factor
			// (drawn from phiBreed) cancels.
			bool obsFemale = (sampleIndiv[chosenIndiv].sex == SEX_FEMALE);
			auto sexTerm = [&](unsigned int age)->double {
				if (age == 1)
					return obsFemale ? log(phiMove) : log(1.0 - phiMove);
				return obsFemale ? log(rho) : log(1.0 - rho);
			};
			dtLogSex = sexTerm(tempIndiv.migrantAge) - sexTerm(sampleIndiv[chosenIndiv].migrantAge);
		}
		else
			tempIndiv.migrantSex = SEX_UNKNOWN;

		// calculate change of logL for genetic data with new migrant ancestry
		if (!NOLIKELIHOOD)
		{
			if (gArgs.collapse)
			{
				// Collapsed: remove the chosen individual from the count tables,
				// then the genotype-likelihood ratio is the difference of the
				// (phase-marginalized) add-log-probs for the proposed vs current
				// ancestry on the without-i counts. Re-added below (accept/reject).
				removeIndividual(gCount, gCountN, sampleIndiv[chosenIndiv], gAssign[chosenIndiv], noLoci);
				double dcur  = computeAddLogProb(gCount, gCountN, sampleIndiv[chosenIndiv], noLoci, noAlleles, ALLELE_PRIOR_ALPHA, FStat);
				double dprop = computeAddLogProb(gCount, gCountN, tempIndiv, noLoci, noAlleles, ALLELE_PRIOR_ALPHA, FStat);
				dtLogL = dprop - dcur;
			}
			else
			{
				logLprop = logLik(tempIndiv,alleleFreqs,logAlleleFreqs,FStat,log1MinusFStat,noLoci);
				dtLogL = logLprop - sampleIndiv[chosenIndiv].logL;
			}
		}

		// calculate change of logPr for migrant counts with new migrant ancestry
		double dtLogPrCount=0.0;

		if ((tempIndiv.migrantPopln != sampleIndiv[chosenIndiv].migrantPopln)||
			(tempIndiv.migrantAge != sampleIndiv[chosenIndiv].migrantAge))
		{
		if(tempIndiv.migrantAge == 0)
		{
			dtLogPrCount += (log(1.0-(1.0+2.0*gamma)*migrationRates[tempIndiv.samplePopln][noPopln])-
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
			dtLogPrCount += (log(2.0*gamma*migrationRates[tempIndiv.samplePopln][tempIndiv.migrantPopln])-
			log(migrantCounts[tempIndiv.samplePopln][tempIndiv.migrantPopln][2]+1.0));
		}

		if(sampleIndiv[chosenIndiv].migrantAge == 0)
		{
			dtLogPrCount -= (log(1.0-(1.0+2.0*gamma)*migrationRates[sampleIndiv[chosenIndiv].samplePopln][noPopln]) -
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
					dtLogPrCount -= (log(2.0*gamma*migrationRates[sampleIndiv[chosenIndiv].samplePopln][sampleIndiv[chosenIndiv].migrantPopln]) -
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
			logPrMHR = dtLogPrCount + dtLogL + logHastings + dtLogSex;
		else
			logPrMHR = dtLogPrCount + logHastings + dtLogSex;

		if(alpha <= exp(logPrMHR))
		{
			unsigned int oldMigPopln = sampleIndiv[chosenIndiv].migrantPopln;
			unsigned int oldMigAge   = sampleIndiv[chosenIndiv].migrantAge;
			sampleIndiv[chosenIndiv].migrantAge = tempIndiv.migrantAge;
			sampleIndiv[chosenIndiv].migrantPopln = tempIndiv.migrantPopln;
			sampleIndiv[chosenIndiv].migrantSex = tempIndiv.migrantSex;
			if (!NOLIKELIHOOD && !gArgs.collapse)
				sampleIndiv[chosenIndiv].logL = logLprop;
			// O(1) incremental update of migrant counts and category lists
			// (replaces the O(N) fillMigrantCounts re-tally). samplePopln equals
			// sampleIndiv[chosenIndiv].samplePopln (selection was from that category).
			migrantCounts[samplePopln][oldMigPopln][oldMigAge]--;
			migrantCounts[samplePopln][tempIndiv.migrantPopln][tempIndiv.migrantAge]++;
			catRemove((int)chosenIndiv);
			catAdd((int)chosenIndiv, catIndex(samplePopln, tempIndiv.migrantPopln, tempIndiv.migrantAge));
			ancestryAcceptRate = (1.0/i)+((i-1.0)/i)*ancestryAcceptRate;
		}
		else ancestryAcceptRate = ((i-1.0)/i)*ancestryAcceptRate;

		// Collapsed: re-add the chosen individual under its final ancestry
		// (proposed if accepted, current if rejected), sampling the age-2 phase.
		if (gArgs.collapse && !NOLIKELIHOOD)
			addIndividual(gCount, gCountN, sampleIndiv[chosenIndiv], gAssign[chosenIndiv],
			              noLoci, noAlleles, ALLELE_PRIOR_ALPHA, FStat);
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
			while ((propMigrationRates[j]<0)||(propMigrationRates[j]>(1.0/(1.0+2.0*gamma)-migrationRates[sourcePopln][noPopln]+migrationRates[sourcePopln][j])))
			{
				if (propMigrationRates[j]<0) {
					propMigrationRates[j]=std::fabs(propMigrationRates[j]);
				}
				if (propMigrationRates[j]>(1.0/(1.0+2.0*gamma)-migrationRates[sourcePopln][noPopln]+migrationRates[sourcePopln][j]))
					propMigrationRates[j]=2.0*(1.0/(1.0+2.0*gamma)-migrationRates[sourcePopln][noPopln]+migrationRates[sourcePopln][j])-propMigrationRates[j];
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
				logPrCurr += migrantCounts[sourcePopln][l][2]*log(2.0*gamma*migrationRates[sourcePopln][l])-gsl_sf_lnfact(migrantCounts[sourcePopln][l][2]);
			}
	logPrCurr += migrantCounts[sourcePopln][sourcePopln][0]*log(1.0-(1.0+2.0*gamma)*migrationRates[sourcePopln][noPopln]) - gsl_sf_lnfact(migrantCounts[sourcePopln][sourcePopln][0]);
	for(unsigned int l = 0; l < noPopln; l++)
		if(sourcePopln != l)
		{
			logPrProp += migrantCounts[sourcePopln][l][1]*log(propMigrationRates[l])-gsl_sf_lnfact(migrantCounts[sourcePopln][l][1]);
			logPrProp += migrantCounts[sourcePopln][l][2]*log(2.0*gamma*propMigrationRates[l])-gsl_sf_lnfact(migrantCounts[sourcePopln][l][2]);
		}
	logPrProp += migrantCounts[sourcePopln][sourcePopln][0]*log(1.0-(1.0+2.0*gamma)*propMigrationRates[noPopln]) - gsl_sf_lnfact(migrantCounts[sourcePopln][sourcePopln][0]);

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

if(!NOALLELEMCMC && !gArgs.collapse)
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

if(!NOFSTATMCMC && !gArgs.collapse)
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

// Collapsed sampler: inbreeding via the IBD-indicator Gibbs (replaces the FStat
// MH), keeping the count tables consistent and drawing F from its conjugate Beta.
// The full sweep is O(noIndiv*noLoci), so it is run once every noIndiv iterations
// (amortized O(noLoci), matching the ancestry move); F is a strong all-data Gibbs
// so this periodic cadence mixes well. Individuals whose ancestry changes have
// their IBD indicators resampled by addIndividual in the meantime.
if(gArgs.collapse && !NOLIKELIHOOD && !NOFSTATMCMC
   && (i % (noIndiv > 1 ? noIndiv : 1) == 0))
{
	inbreedingUpdate(gCount, gCountN, sampleIndiv, gAssign, noIndiv, noLoci,
	                 noAlleles, noPopln, ALLELE_PRIOR_ALPHA, FStat);
	gHomEpochCur++;   // F changed: invalidate the homozygote-mixture memo (O(1))
}

if(!NOMISSINGDATA && !gArgs.collapse)
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

// Metropolis-Hastings update of the migrant breeding-success multiplier gamma
// (tau = 2*gamma = age-2:age-1 rate ratio). Random walk on log(gamma); the age
// counts inform it. Requires (1+tau)*sum_m < 1 in every population.
// gamma is bounded to (0, 1]: gamma=1 is migrants breeding as well as residents
// (the neutral/stationary maximum), gamma<1 is reduced migrant breeding success;
// gamma>1 would mean migrants out-breed residents, impossible under a stationary
// neutral model. The proposal is reflected at log(gamma)=0 (symmetric, so no
// Hastings term). With --fixgamma, gamma stays fixed at 1 and this is skipped.
if (!NOMIGRATEMCMC && !gArgs.fixGamma)
{
	double propLogGamma = log(gamma) + gsl_ran_gaussian(r, GAMMA_PROP_SD);
	if (propLogGamma > 0.0) propLogGamma = -propLogGamma;   // reflect at gamma=1: bound to (0,1]
	double propGamma = exp(propLogGamma);
	double tProp = 2.0 * propGamma, tCur = 2.0 * gamma;
	bool feasible = true;
	for (unsigned int s = 0; s < noPopln && feasible; s++)
		if ((1.0 + tProp) * migrationRates[s][noPopln] >= 1.0) feasible = false;
	if (feasible)
	{
		double curLL = 0.0, propLL = 0.0;
		for (unsigned int s = 0; s < noPopln; s++)
		{
			long int n0 = migrantCounts[s][s][0];
			curLL  += n0 * log(1.0 - (1.0 + tCur)  * migrationRates[s][noPopln]);
			propLL += n0 * log(1.0 - (1.0 + tProp) * migrationRates[s][noPopln]);
			long int n2s = 0;
			for (unsigned int b = 0; b < noPopln; b++)
				if (b != s) n2s += migrantCounts[s][b][2];
			curLL  += n2s * log(tCur);
			propLL += n2s * log(tProp);
		}
		double lsd2 = GAMMA_PRIOR_LOGSD * GAMMA_PRIOR_LOGSD;
		double logPriorCur  = -0.5 * log(gamma) * log(gamma) / lsd2;   // log(gamma) ~ N(0, sd)
		double logPriorProp = -0.5 * propLogGamma * propLogGamma / lsd2;
		if (log(gsl_rng_uniform(r)) < (propLL + logPriorProp) - (curLL + logPriorCur))
			gamma = propGamma;
	}
}

// Sex-biased dispersal: Gibbs-update the migrant-parent sex (sigma) of every
// second-generation migrant, then the migrant (phi) and resident (rho) female
// fractions. phi counts dispersers (age-1 own sex + age-2 migrant-parent sex);
// rho counts non-dispersers (own sex of age-0 and age-2 individuals).
if (sexBiasModel)
{
	long int moveF = 0, moveM = 0;    // age-1 migrants (movement)
	long int breedF = 0, breedM = 0;  // age-2 migrant parents (effective gene flow)
	long int rF = 0, rM = 0;          // residents (non-dispersers)
	for (unsigned int m = 0; m < noIndiv; m++)
	{
		unsigned int age = sampleIndiv[m].migrantAge;
		if (age == 0)
		{
			// Non-migrant resident: own sex informs the resident sex ratio rho.
			if (sampleIndiv[m].sex == SEX_FEMALE) rF++; else rM++;
		}
		else if (age == 1)
		{
			// First-generation migrant (a mover): its observed sex informs phiMove.
			sampleIndiv[m].migrantSex = sampleIndiv[m].sex;
			if (sampleIndiv[m].sex == SEX_FEMALE) moveF++; else moveM++;
		}
		else if (age == 2)
		{
			// A 2nd-gen individual is itself locally born (own sex -> rho); its
			// migrant parent bred successfully, so the parent's sex (sigma) informs
			// phiBreed (effective gene flow).
			if (sampleIndiv[m].sex == SEX_FEMALE) rF++; else rM++;
			unsigned int newSigma;
			if (sampleIndiv[m].sex == SEX_MALE && !xLoci.empty())
			{
				// A male's single X is from the source population if the migrant
				// parent was the mother (sigma = F) and from the native population
				// if the father (sigma = M). Gibbs from the posterior log-odds.
				const unsigned int src = sampleIndiv[m].migrantPopln;
				const unsigned int nat = sampleIndiv[m].samplePopln;
				if (gArgs.collapse)
				{
					// Collapsed sigma Gibbs: the single X copy currently sits in src
					// (sigma = F) or nat (sigma = M). Remove it, form the posterior
					// log-odds from the count-based frequencies on the without-copy
					// tables, resample, and re-add to the chosen population.
					unsigned int oldPop = (sampleIndiv[m].migrantSex == SEX_FEMALE) ? src : nat;
					for (size_t xi = 0; xi < xLoci.size(); xi++)
					{
						unsigned int l = xLoci[xi];
						if (sampleIndiv[m].genotype[l][1] != HEMIZYGOUS) continue;
						int a0 = sampleIndiv[m].genotype[l][0];
						if (a0 < 0) continue;
						removeCopy(gCount, gCountN, oldPop, l, a0);
					}
					double logOdds = log(phiBreed) - log(1.0 - phiBreed);  // prior log-odds F:M
					for (size_t xi = 0; xi < xLoci.size(); xi++)
					{
						unsigned int l = xLoci[xi];
						if (sampleIndiv[m].genotype[l][1] != HEMIZYGOUS) continue;
						int a0 = sampleIndiv[m].genotype[l][0];
						if (a0 < 0) continue;
						logOdds += log(addRatio(gCount, gCountN, src, l, a0, noAlleles[l], ALLELE_PRIOR_ALPHA))
						         - log(addRatio(gCount, gCountN, nat, l, a0, noAlleles[l], ALLELE_PRIOR_ALPHA));
					}
					double pF = 1.0 / (1.0 + exp(-logOdds));
					newSigma = (gsl_rng_uniform(r) < pF) ? SEX_FEMALE : SEX_MALE;
					unsigned int newPop = (newSigma == SEX_FEMALE) ? src : nat;
					unsigned char code = (newSigma == SEX_FEMALE) ? 3 : 4;   // gAssign X-copy code
					for (size_t xi = 0; xi < xLoci.size(); xi++)
					{
						unsigned int l = xLoci[xi];
						if (sampleIndiv[m].genotype[l][1] != HEMIZYGOUS) continue;
						int a0 = sampleIndiv[m].genotype[l][0];
						if (a0 < 0) continue;
						addCopy(gCount, gCountN, newPop, l, a0);
						gAssign[m][l] = code;
					}
				}
				else
				{
					double logOdds = log(phiBreed) - log(1.0 - phiBreed);  // prior log-odds F:M
					for (size_t xi = 0; xi < xLoci.size(); xi++)
					{
						unsigned int l = xLoci[xi];
						if (sampleIndiv[m].genotype[l][1] != HEMIZYGOUS) continue;
						int a0 = sampleIndiv[m].genotype[l][0];
						if (a0 < 0) continue;
						logOdds += logAlleleFreqs[src][l][a0] - logAlleleFreqs[nat][l][a0];
					}
					double pF = 1.0 / (1.0 + exp(-logOdds));
					newSigma = (gsl_rng_uniform(r) < pF) ? SEX_FEMALE : SEX_MALE;
				}
			}
			else
			{
				// Female (or no X loci): sigma has no likelihood effect -> prior.
				newSigma = (gsl_rng_uniform(r) < phiBreed) ? SEX_FEMALE : SEX_MALE;
			}
			if (newSigma != sampleIndiv[m].migrantSex)
			{
				sampleIndiv[m].migrantSex = newSigma;
				// Only a male's genotype likelihood depends on sigma. In collapse
				// mode the count tables were already updated above (no cached logL).
				if (sampleIndiv[m].sex == SEX_MALE && !NOLIKELIHOOD && !gArgs.collapse)
					sampleIndiv[m].logL = logLik(sampleIndiv[m], alleleFreqs, logAlleleFreqs,
					                             FStat, log1MinusFStat, noLoci);
			}
			if (sampleIndiv[m].migrantSex == SEX_FEMALE) breedF++; else breedM++;
		}
	}
	// Beta-Bernoulli conjugate updates of phiMove, phiBreed, and rho; clamp.
	phiMove  = gsl_ran_beta(r, PHI_PRIOR_A + (double)moveF,  PHI_PRIOR_B + (double)moveM);
	phiBreed = gsl_ran_beta(r, PHI_PRIOR_A + (double)breedF, PHI_PRIOR_B + (double)breedM);
	rho      = gsl_ran_beta(r, PHI_PRIOR_A + (double)rF,     PHI_PRIOR_B + (double)rM);
	if (phiMove  < 1e-6) phiMove  = 1e-6;  if (phiMove  > 1.0-1e-6) phiMove  = 1.0-1e-6;
	if (phiBreed < 1e-6) phiBreed = 1e-6;  if (phiBreed > 1.0-1e-6) phiBreed = 1.0-1e-6;
	if (rho      < 1e-6) rho      = 1e-6;  if (rho      > 1.0-1e-6) rho      = 1.0-1e-6;
	// Retain counts for the Rao-Blackwellized Savage-Dickey tests.
	moveCountF = moveF;   moveCountM = moveM;
	breedCountF = breedF; breedCountM = breedM;
	rhoCountF = rF;       rhoCountM = rM;
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
			logLM = migCountLogProb(migrantCounts,migrationRates,noPopln,gamma);
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

				logLM = migCountLogProb(migrantCounts,migrationRates,noPopln,gamma);
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

				logLM = migCountLogProb(migrantCounts,migrationRates,noPopln,gamma);
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

			// Accumulate the migrant breeding-success multiplier gamma.
			if (gammaNSamples > 1)
				varGamma = ((gammaNSamples-1.0)/gammaNSamples)*varGamma
				           + (gamma - avgGamma)*(gamma - avgGamma)/(gammaNSamples+1.0);
			avgGamma = avgGamma + (gamma - avgGamma)/(1.0 + gammaNSamples);
			if (gamma < 1.0) gammaLtOneCount++;
			gammaNSamples++;

			// Accumulate posteriors for phiMove (movement), phiBreed (gene flow),
			// rho (residents), and the three Savage-Dickey equality tests.
			if (sexBiasModel)
			{
				if (iter > 1)
				{
					varMove  = ((iter-1.0)/iter)*varMove  + (phiMove  - avgMove )*(phiMove  - avgMove )/(iter+1.0);
					varBreed = ((iter-1.0)/iter)*varBreed + (phiBreed - avgBreed)*(phiBreed - avgBreed)/(iter+1.0);
					varRho   = ((iter-1.0)/iter)*varRho   + (rho      - avgRho  )*(rho      - avgRho  )/(iter+1.0);
				}
				avgMove  = avgMove  + (phiMove  - avgMove )/(1.0+iter);
				avgBreed = avgBreed + (phiBreed - avgBreed)/(1.0+iter);
				avgRho   = avgRho   + (rho      - avgRho  )/(1.0+iter);
				if (phiMove  > rho) moveGtRhoCount++;
				if (phiBreed > rho) breedGtRhoCount++;

				// Rao-Blackwellized Savage-Dickey: density of a difference at 0 equals
				// the overlap integral of the two conjugate Beta full-conditionals.
				double aMv = PHI_PRIOR_A + (double)moveCountF,  bMv = PHI_PRIOR_B + (double)moveCountM;
				double aBr = PHI_PRIOR_A + (double)breedCountF, bBr = PHI_PRIOR_B + (double)breedCountM;
				double aRh = PHI_PRIOR_A + (double)rhoCountF,   bRh = PHI_PRIOR_B + (double)rhoCountM;
				sdMoveEqRhoSum   += betaProductIntegral(aMv, bMv, aRh, bRh);
				sdBreedEqRhoSum  += betaProductIntegral(aBr, bBr, aRh, bRh);
				sdMoveEqBreedSum += betaProductIntegral(aMv, bMv, aBr, bBr);
				sdPhiNSamples++;
			}

			double dirAlpha[MAXALLELE];
			double dirDraw[MAXALLELE];
			for (unsigned int l = 0; l < noPopln; l++)
				for (unsigned int k = 0; k < noLoci; k++)
				{
					// In collapse mode the frequencies are integrated out, so alleleFreqs
					// is never updated. For output, draw f from its exact Dirichlet
					// full-conditional given the current count table (posterior parameters
					// n[.]+alpha); the resulting mean/SD accumulators then carry the same
					// meaning as the standard sampler's sampled-frequency output.
					if(gArgs.collapse)
					{
						for(unsigned int m = 0; m < noAlleles[k]; m++)
							dirAlpha[m] = (double)gCount[l][k][m] + ALLELE_PRIOR_ALPHA;
						gsl_ran_dirichlet(r, noAlleles[k], dirAlpha, dirDraw);
					}
					for(unsigned int m = 0; m < noAlleles[k]; m++)
					{
						double freqVal = gArgs.collapse ? dirDraw[m] : alleleFreqs[l][k][m];
						if(iter > 1)
						{
							sqrDiffMean=(freqVal-avgAlleleFreqs[l][k][m])*(freqVal-avgAlleleFreqs[l][k][m])/(iter+1.0);
							varAlleleFreqs[l][k][m] = ((iter-1.0)/iter)*varAlleleFreqs[l][k][m]+sqrDiffMean;
						}
						avgAlleleFreqs[l][k][m] = avgAlleleFreqs[l][k][m]+(freqVal-avgAlleleFreqs[l][k][m])/(1.0+iter);
					}
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

	mcmcout << "\n Migration Rate Matrix m[i][j] (fraction of pop i from pop j):\n";
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

	// Migrant breeding success: gamma = relative breeding success of migrants
	// (tau = 2*gamma is the age-2:age-1 rate ratio). gamma = 1 means migrants breed
	// at the resident rate (BA3's original assumption); gamma < 1 means reduced
	// migrant breeding success. The effective (gene-flow) migration rate is ~gamma
	// times the apparent (movement) rate m.
	{
		double pReduced = (gammaNSamples > 0) ? (double)gammaLtOneCount/gammaNSamples : 0.0;
		mcmcout.setf(std::ios::fixed, std::ios::floatfield);
		mcmcout << "\n Migrant breeding success:\n";
		mcmcout << "   gamma (migrant breeding success rel. to residents) = "
		        << std::setprecision(4) << avgGamma << "(" << sqrt(varGamma) << ")\n";
		mcmcout << "   age-2:age-1 rate ratio tau = 2*gamma = " << std::setprecision(4) << 2.0*avgGamma
		        << ";  P(gamma < 1, reduced migrant breeding) = " << std::setprecision(3) << pReduced << "\n";
		mcmcout << "   (gamma = 1: migrants breed like residents; effective gene-flow rate ~ gamma * movement rate m)\n";
	}

	// Sex-biased dispersal: report the female fractions for movement (phiMove,
	// age-1 migrants), effective gene flow (phiBreed, age-2 migrant parents), and
	// residents (rho), plus Savage-Dickey tests of each bias against rho and of
	// movement vs gene flow. Bias is measured against rho (the population sex
	// ratio), not 1/2. The X-linked data specifically informs phiBreed.
	if (sexBiasModel)
	{
		auto interpret = [](double log10BF, const char *h0lab, const char *h1lab) {
			std::ostringstream o;
			const char *lab; const char *strength;
			if (log10BF > 0) { lab = h0lab;
				strength = (log10BF>2)?"Decisive":(log10BF>1)?"Strong":(log10BF>0.5)?"Substantial":"Weak"; }
			else { lab = h1lab;
				strength = (log10BF<-2)?"Decisive":(log10BF<-1)?"Strong":(log10BF<-0.5)?"Substantial":"Weak"; }
			o << strength << " for " << lab; return o.str();
		};
		double priorEq = betaProductIntegral(PHI_PRIOR_A, PHI_PRIOR_B, PHI_PRIOR_A, PHI_PRIOR_B);
		double n = (sdPhiNSamples > 0) ? (double)sdPhiNSamples : 1.0;
		double bfMoveRho  = (sdMoveEqRhoSum  / n) / priorEq;   // H0: phiMove  = rho
		double bfBreedRho = (sdBreedEqRhoSum / n) / priorEq;   // H0: phiBreed = rho
		double bfMoveBreed= (sdMoveEqBreedSum/ n) / priorEq;   // H0: phiMove  = phiBreed
		double l10MoveRho  = log10(bfMoveRho  > 0 ? bfMoveRho  : 1e-10);
		double l10BreedRho = log10(bfBreedRho > 0 ? bfBreedRho : 1e-10);
		double l10MoveBreed= log10(bfMoveBreed> 0 ? bfMoveBreed: 1e-10);

		mcmcout.setf(std::ios::fixed, std::ios::floatfield);
		mcmcout << "\n Sex-biased dispersal (female fractions; compared to residents rho):\n";
		mcmcout << "   rho      residents (non-dispersers)     = " << std::setprecision(4) << avgRho
		        << "(" << sqrt(varRho) << ")\n";
		mcmcout << "   phiMove  age-1 migrants (movement)      = " << std::setprecision(4) << avgMove
		        << "(" << sqrt(varMove) << ")   P(>rho) = " << std::setprecision(3) << (double)moveGtRhoCount/n << "\n";
		mcmcout << "   phiBreed age-2 migrant parents (gene flow) = " << std::setprecision(4) << avgBreed
		        << "(" << sqrt(varBreed) << ")   P(>rho) = " << std::setprecision(3) << (double)breedGtRhoCount/n << "\n";
		mcmcout << "   (phi > rho: female-biased; phi < rho: male-biased. The X data informs phiBreed.)\n";
		mcmcout << "   Savage-Dickey tests (BF_01; <1 supports the alternative):\n";
		mcmcout << "     H0: phiMove  = rho (no sex-biased movement)  : BF_01 = " << std::setprecision(3) << bfMoveRho
		        << "  log10 = " << std::showpos << l10MoveRho << std::noshowpos
		        << "   " << interpret(l10MoveRho, "H0", "sex-biased movement") << "\n";
		mcmcout << "     H0: phiBreed = rho (no sex-biased gene flow) : BF_01 = " << std::setprecision(3) << bfBreedRho
		        << "  log10 = " << std::showpos << l10BreedRho << std::noshowpos
		        << "   " << interpret(l10BreedRho, "H0", "sex-biased gene flow") << "\n";
		mcmcout << "     H0: phiMove  = phiBreed (movement = gene flow): BF_01 = " << std::setprecision(3) << bfMoveBreed
		        << "  log10 = " << std::showpos << l10MoveBreed << std::noshowpos
		        << "   " << interpret(l10MoveBreed, "H0", "movement differs from gene flow") << "\n";
	}

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

	// Free collapsed-sampler count tables and phase store
	if(gArgs.collapse && gCount != nullptr)
	{
		for(unsigned int p = 0; p < noPopln; p++)
		{
			for(unsigned int l = 0; l < noLoci; l++) delete[] gCount[p][l];
			delete[] gCount[p]; delete[] gCountN[p];
		}
		delete[] gCount; delete[] gCountN;
		for(unsigned int i = 0; i < noIndiv; i++) delete[] gAssign[i];
		delete[] gAssign;
		delete[] gLogTab; gLogTab = nullptr;
		for(unsigned int p = 0; p < noPopln; p++) { delete[] gHomLog[p]; delete[] gHomEpoch[p]; }
		delete[] gHomLog;   gHomLog = nullptr;
		delete[] gHomEpoch; gHomEpoch = nullptr;
	}

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

// gamma = migrant relative breeding success; tau = 2*gamma is the age-2:age-1
// rate ratio (P(age2)=tau*m, P(age0)=1-(1+tau)*sum m). gamma=1 (tau=2) is BA3's
// original assumption that migrants breed at the resident rate.
double migCountLogProb(long int ***migrantCounts, double **migrationRates, unsigned int noPopln, double gamma)
{
	const double tau = 2.0 * gamma;
	double logPr=0.0;
	for(unsigned int l = 0; l < noPopln; l++)
	{
		for(unsigned int k = 0; k < noPopln; k++)
			if(l != k)
			{
				logPr+=migrantCounts[l][k][1]*log(migrationRates[l][k])-gsl_sf_lnfact(migrantCounts[l][k][1]);
				logPr+=migrantCounts[l][k][2]*log(tau*migrationRates[l][k])-gsl_sf_lnfact(migrantCounts[l][k][2]);
			}
		logPr+=migrantCounts[l][l][0]*log(1.0-(1.0+tau)*migrationRates[l][noPopln]) - gsl_sf_lnfact(migrantCounts[l][l][0]);
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

			if (a1 == HEMIZYGOUS)
			{
				// Hemizygous male X (single gene copy) from the sampled population
				logPr += logFreqPop[i][a0];
			}
			else if (a0 == a1)
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

			if (a1 == HEMIZYGOUS)
			{
				// Hemizygous male X (single gene copy) from the source population
				logPr += logFreqPop[i][a0];
			}
			else if (a0 == a1)
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

			if (a1 == HEMIZYGOUS)
			{
				// Hemizygous male X (single gene copy) in a 2nd-generation migrant.
				// A male's single X comes from his mother, so the copy is from the
				// source population if the migrant parent was the mother (migrantSex
				// = female) and from the native population if it was the father.
				if (Indiv.migrantSex == SEX_FEMALE)
					logPr += logFreqMig[i][a0];        // X from source (mother migrated)
				else if (Indiv.migrantSex == SEX_MALE)
					logPr += logFreqSam[i][a0];        // X from native (father migrated)
				else
				{
					// Unknown migrant-parent sex: marginalize 1/2 (source + native)
					const double t1 = logFreqMig[i][a0];
					const double t2 = logFreqSam[i][a0];
					const double mx = (t1 > t2) ? t1 : t2;
					logPr += mx + log(exp(t1 - mx) + exp(t2 - mx)) - LOG2;
				}
			}
			else if (a0 == a1)
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

	if (a1 == HEMIZYGOUS)
	{
		// Hemizygous male X (single gene copy)
		if (Indiv.migrantAge == 0)
			logPr = logAlleleFreqs[Indiv.samplePopln][chosenLocus][a0];
		else if (Indiv.migrantAge == 1)
			logPr = logAlleleFreqs[Indiv.migrantPopln][chosenLocus][a0];
		else
		{
			// age 2 male X: single copy from source if migrant parent is the
			// mother (migrantSex = female), else from native; unknown -> 1/2 each.
			if (Indiv.migrantSex == SEX_FEMALE)
				logPr = logAlleleFreqs[Indiv.migrantPopln][chosenLocus][a0];
			else if (Indiv.migrantSex == SEX_MALE)
				logPr = logAlleleFreqs[Indiv.samplePopln][chosenLocus][a0];
			else
			{
				const double t1 = logAlleleFreqs[Indiv.migrantPopln][chosenLocus][a0];
				const double t2 = logAlleleFreqs[Indiv.samplePopln][chosenLocus][a0];
				const double mx = (t1 > t2) ? t1 : t2;
				logPr = mx + log(exp(t1 - mx) + exp(t2 - mx)) - LOG2;
			}
		}
		return logPr;
	}

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

			// Store sample for KL divergence computation
			sdStats[i][j].samples.push_back(m);

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

// Compute KDE with reflection at boundary for a single point
// Uses Gaussian kernel with reflection method for boundary correction at m=0
double computeKDEatPoint(double x, const std::vector<double> &samples, double bandwidth)
{
	const double PI = 3.14159265358979323846;
	const double SQRT_2PI = sqrt(2.0 * PI);
	long int N = samples.size();
	double h = bandwidth;
	double density = 0.0;

	for (long int t = 0; t < N; t++) {
		double m = samples[t];
		// Standard kernel contribution
		double u1 = (x - m) / h;
		density += exp(-0.5 * u1 * u1);
		// Reflection contribution (for boundary at 0)
		double u2 = (x + m) / h;
		density += exp(-0.5 * u2 * u2);
	}

	return density / (N * h * SQRT_2PI);
}

// Compute KL divergence using numerical integration with KDE posterior and uniform prior
// KL(posterior || prior) = integral of p(m|D) * log(p(m|D) / 3) dm over [0, 1/3]
double computeKLDivergence(const std::vector<double> &samples, double bandwidth)
{
	const int N_GRID = 1000;  // Number of grid points
	const double UPPER_BOUND = 1.0 / 3.0;  // Upper limit of support
	const double PRIOR_DENSITY = 3.0;  // Uniform(0, 1/3) has density = 3
	const double DENSITY_FLOOR = 1e-10;  // Minimum density to avoid log(0)
	const double LOG2E = 1.4426950408889634;  // log2(e) for conversion to bits

	double deltaM = UPPER_BOUND / N_GRID;
	double klNats = 0.0;

	for (int i = 0; i < N_GRID; i++) {
		double m = (i + 0.5) * deltaM;  // Midpoint of interval
		double pPost = computeKDEatPoint(m, samples, bandwidth);

		if (pPost > DENSITY_FLOOR) {
			// KL integrand: p(m|D) * log(p(m|D) / p(m))
			// where p(m) = 3 for uniform prior on [0, 1/3]
			klNats += pPost * log(pPost / PRIOR_DENSITY) * deltaM;
		}
	}

	// Convert from nats to bits and ensure non-negative
	double klBits = klNats * LOG2E;
	return (klBits > 0) ? klBits : 0.0;
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

			// Compute KL divergence using numerical integration
			// KL(posterior || prior) with exact uniform prior and KDE posterior
			double klBits = computeKLDivergence(sdStats[i][j].samples, bandwidth);

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
