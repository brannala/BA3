# BayesAss Edition 3.4 User's Manual

Bruce Rannala

&copy;2007-2025 University of California Davis

Last updated 23 Dec 2025

## 1 Installation

The latest version of BA3 (version 3.4.0) can be downloaded
[here](https://github.com/brannala/BA3/releases). Unzip the archive by
double clicking the downloaded file. A folder will be created in your current directory
containing the source code and/or precompiled binary files for
various computer operating systems. There is also an examples directory containing
example input files with samples from either 2 or 3 populations. If you have a compiler
and are adventurous you can try compiling the source code (see below), otherwise refer
to the instructions below to use a precompiled binary file for your specific computer
operating system.

**Note:** As of version 3.4.0, there is a single unified executable `BA3` that handles both
SNP and microsatellite data. The program supports up to 500 alleles per locus, 100,000 loci,
100 populations, and 20,000 individuals. Memory is allocated dynamically based on your dataset size.

### 1.1 Mac OS X
#### Using Homebrew package manager

The recommended way to install BA3 on Mac (or Linux) is using the [homebrew](https://brew.sh/)
package manager. Once you have homebrew installed on your machine you can install BA3 using
the following commands (executed in the Terminal):
```
brew tap brannala/ba3
brew install ba3
```
The program will be installed and the command `BA3` will be available in
the Terminal. You will want to download the documentation and example files (BA3-docs-examples.zip) from
[here](https://github.com/brannala/BA3/releases).

#### Using precompiled executables

Download the zip archived file with the latest version of the software
[here](https://github.com/brannala/BA3/releases). Unzip the
archive by double clicking the downloaded file. A folder will be created in your current
directory containing the program executable file `BA3`, example data files (in subfolder
examples), and this manual (in subfolder docs). Note that there are two binaries available:
"M" for users of newer M1 or M2 CPU Macs and "I" for users of older Intel CPU Macs. The Intel executable
will run on M1/M2 Macs but this is through emulation and is not recommended.

### 1.2 Windows

Download the zip archived file with the latest version of the software [here](https://github.com/brannala/BA3/releases).
Choose either the file `BA3Windows32.zip` or `BA3Windows64.zip` depending on whether you are running a 32 bit (older)
or 64 bit (newer) version of Windows. Unzip the archive by double
clicking the downloaded file. A folder will be created in your current directory containing
the program executable file `BA3.exe`, example data files (in subfolder examples), and
this manual (in subfolder docs).

### 1.3 Linux

Download the zip or tar.gz archived file with the latest version of the software
[here](https://github.com/brannala/BA3/releases). Unzip the
archive by double clicking the downloaded file. A folder will be created in your current
directory containing the program executable file `BA3`, example data files (in subfolder
examples), and this manual (in subfolder docs).

### 1.4 Compiling the program

You do not need to compile the program if you have been successful following either
steps 1.1 or 1.2 above. The source code for the program is found in the source tarball
distribution file named BA3-*.*.*.tar.gz where * indicate the version numbers. The
program uses routines from the gnu scientific library (gsl) and this library (and header
files) must be installed prior to compiling. The gsl library can be found [here](https://www.gnu.org/software/gsl/).
For VCF file support, the htslib library is also required ([htslib](https://github.com/samtools/htslib)).
It is recommended that you install these libraries using a package manager such
as apt in Ubuntu linux or homebrew on Macs. If using a
command line C++ compiler (e.g., g++, c++, etc), with gsl and htslib installed in the standard location, simply execute the
following terminal commands in the directory that contains the source tarball:
```
tar -xvzf BA3-*.tar.gz
cd BA3-*
make
```
This will create the executable file `BA3` in the current directory and typing ./BA3 at the command
prompt will then execute the program.

## 2 Running the program

The BA3 program is a command line program. If you are familiar with the unix terminal
you will find it straightforward to use as it adheres to standard unix conventions for
command line options, etc. If you have never used a terminal (command line) program
you can find a beginners guide [here](https://people.ischool.berkeley.edu/~kevin/unix-tutorial/toc.html).
Detailed instructions for running the program on
the Mac OS X (or other unix-based) operating system are provided below.

### 2.1 Getting BA3 up and running on Mac OS X or Unix

#### 2.1.1 Running BA3 installed with Homebrew
Download the examples and documentation `BA3-docs-examples.zip` and uncompress. If you uncompress
the files on your desktop then in the Terminal type:
```
cd /Users/<login>/Desktop/BA3-docs-examples
```
to run a 3 population example dataset type:
```
BA3 examples/3pop.txt
```

#### 2.1.2 Running BA3 installed from binaries on Mac or Linux
To run the program, you will first need to start the terminal application which can be
found in the Applications/Utilities folder on a Mac. The following description assumes that you have unzipped
the BA3 distribution file on the Desktop. If you have placed it elsewhere you will need to
change the commands to indicate the correct file path. Once you open terminal you will
see a command line prompt. On a Mac, at the prompt type:
```
cd /Users/<login>/Desktop/BA3.*/binary/macosx
```
where <login> is replaced with your account login name. On another unix computer you
will specify the path to the binary that you unpacked (or compiled from source code).
The unix command cd is short for "change directory" and the above command changes
the current working directory from the user's home directory (the default) to the directory
where the BA3 program binary resides. To run the program using an example data
file with 3 populations (contained in the subdirectory examples) type the following command:
```
./BA3 examples/3pop.txt
```
The prefix ./ means "current directory" and tells the operating system to look for the
program file named `BA3` in the current working directory.

#### 2.1.3 Program screen output
You should see output similar to the following:
```
  ┌──────────────────────────────────────────────────────┐
  │             BayesAss Edition 3.4.0 (BA3)             │
  │                 Released: 12/23/2025                 │
  │                    Bruce Rannala                     │
  │   Department of Evolution and Ecology at UC Davis    │
  └──────────────────────────────────────────────────────┘

  Run Configuration:
    Input:       examples/3pop.txt
    Output:      BA3out.txt
    Individuals: 400
    Populations: 3
    Loci:        7

  MCMC Settings:
    Iterations:  5M
    Burn-in:     2.5M
    Sampling:    100
    Seed:        10
    Mixing:      dM=0.1 dA=0.1 dF=0.1 (initial, autotune on)

  Burn-in:  [████████████████░░░░░░░░░░░░░░]  52%  2.6M/5M  ETA: 45s
```

The program displays a pre-run summary showing the input/output files, dataset dimensions,
and MCMC settings. During the run, an in-place progress bar shows completion percentage
and estimated time remaining (ETA). When complete, a post-run summary is displayed:

```
  ┌──────────────────────────────────────┐
  │            Run Complete              │
  └──────────────────────────────────────┘
  Elapsed time: 1m 32s
  Final mixing: dM=0.095 dA=0.312 dF=0.487

  Output written to: BA3out.txt
```

The "Final mixing" values show the auto-tuned mixing parameters used during sampling.

### 2.2 Getting BA3 up and running on Windows

You will first need to run the Windows "Command Prompt" program which will open a
console that you can use to run the `BA3.exe` program. At the command prompt use the
"cd" command to move to the directory where the file `BA3.exe` is found. For example,
```
C:\Users\bruce>cd Desktop\BA3Windows64
C:\Users\bruce\Desktop\BA3Windows64
```
To run the example file 3pop.txt use the command:
```
C:\Users\bruce\Desktop\BA3Windows64>BA3.exe examples\3pop.txt
```

## 3 Data file format

### 3.1 BA3 native format

The BA3 program uses an input file format that is identical to that of earlier BayesAss
releases. The input file should be in a plain text format. DO NOT use a word processor
such as Word to create the input file without explicitly converting it to a text file format
before use. One possible approach is to input the data into a spreadsheet program such as
Excel and then save the file as a "space-delimited text file." Another approach is to install
one of the many available free text editors such as emacs or vi on your computer. Each
line of the input file should have the following format:
```
indivID popID locID allele1 allele2
```
where `indivID` is a unique identifier for the individual, `popID` is a unique identifier of the
individual's source population, `locID` is a unique identifier for the locus, and `allele1`
and `allele2` are the allele labels for each allele of the individual's genotype. The order of
the alleles on the line is arbitrary. Missing alleles are represented using a 0. If there are n
individuals and L loci there will be n × L lines in the input file. See the example data files
distributed with the program.

**Note:** If a locus entry is completely missing for some individuals (not coded as 0/0 but simply
absent from the file), BA3 will print a warning and treat the missing entries as missing data:
```
  Warning: 1 locus/loci missing for some individuals (treated as missing data):
    LG1-1886280 (missing in 2 individuals)
```

### 3.2 VCF file format

As of version 3.4.0, BA3 can read genotype data directly from VCF (Variant Call Format) files.
This is useful for genomic SNP data. To use VCF input, you need two files:

1. A VCF file containing genotype data
2. A metadata file mapping individual IDs to population IDs

The metadata file should have two tab-separated columns:
```
individual_id	population_id
sample1	pop1
sample2	pop1
sample3	pop2
```

To run BA3 with VCF input:
```
./BA3 -V mydata.vcf -M metadata.txt -o output.txt
```

The VCF file should have a proper header including contig definitions. Missing genotypes
(./.) in the VCF are treated as missing data.

## 4 Command line options

The BA3 program has command line options that allow you to control
the way the program runs and the level of detail in the output that it produces. The
command line options are given after the program name and before the input file name.
For example,
```
./BA3 -i10000000 -o myout.txt myin.txt
```
executes the program for 10 million iterations, writing the output to
the file myout.txt and using the input file myin.txt. Some options such as the option
specifying the number of iterations, `-i`, take parameter values while others such as `-v`
do not. Parameter values should follow the option specifier and may, or may not, be
separated from the option specifier by a space.

Table 1 lists all the command line options with a brief description of their parameters and
effects.

| Option            | Values                  | Effect                                      |
|:------------------|:------------------------|:--------------------------------------------|
| -a --deltaA       | $0 < \Delta_A \leq 1.0$ | Mixing parameter for allele frequencies     |
| -b --burnin       | Positive integer        | Number of iterations to discard as burnin   |
| -f --deltaF       | $0 < \Delta_F \leq 1.0$ | Mixing parameter for inbreeding coefficients|
| -g --genotypes    | None                    | Output genotypes and migrant ancestries     |
| -i --iterations   | Positive integer        | Number of iterations for MCMC               |
| -m --deltaM       | $0 < \Delta_M \leq 1.0$ | Mixing parameter for migration rates        |
| -n --sampling     | Positive integer        | Interval between samples for MCMC           |
| -o --output       | String                  | Output file name                            |
| -s --seed         | Positive integer        | Seed for random number generator            |
| -p --nolikelihood | None                    | Fix likelihood to 1 and generate priors     |
| -t --trace        | None                    | Create a trace file to monitor convergence  |
| -u --settings     | None                    | Output options and parameter settings       |
| -v --verbose      | None                    | Use verbose screen output                   |
| -T --autotune     | None                    | Auto-tune mixing parameters (default: on)   |
| -N --noautotune   | None                    | Disable auto-tuning of mixing parameters    |
| -V --vcf          | String                  | VCF input file (requires -M)                |
| -M --meta         | String                  | Metadata file for VCF input                 |
| -F --freqfile     | String                  | Output allele frequencies to separate file  |

**Table 1:** Options available for BA3 program

### 4.1 Random number generator seed

The option -s (--seed) is used to specify a positive integer used to "seed" the random
number generator algorithm. Separate runs of the program started using
same seed will produce exactly the same outcome. To test whether the program is converging
it is important to carry out several independent runs initiated with different seeds.
To start the program using 10456 as the random number seed use the following command:
```
./BA3 -s10456 myin.txt
```
If no seed is specified the default seed is 10.

### 4.2 MCMC iterations, burn-in and sampling interval

The command line option -i (--iterations) specifies the number of iterations for the
Markov chain Monte Carlo (MCMC) analysis. By default the program uses 5,000,000
iterations. The value of the number of iterations should be a positive
integer. For example,
```
./BA3 -i10000000 test.txt
```
will execute the program using the data file test.txt and carry out 10 million iterations. The
option -b (--burnin) is used to specify a positive integer that is the number of iterations
of the MCMC that are discarded before sampling begins. By default, burn-in is 50% of total iterations.
For example,
```
./BA3 -i10000000 -b1000000 test.txt
```
will run the MCMC for 10 million iterations, discarding the first 1 million iterations. The option -n (--sampling) is
used to specify a positive integer that is the interval between samples. For example,
```
./BA3 -i10000000 -b1000000 -n1000 test.txt
```
will run the MCMC for 10 million iterations, discarding the first 1 million iterations and
sampling every 1000 iterations.

### 4.3 MCMC mixing parameters and auto-tuning

For continuous parameters such as migration rates, allele frequencies and inbreeding
coefficients, the size of the proposed change to the parameter value at each iteration of the
MCMC can be adjusted. There are 3 mixing parameter adjustments: -a
(--deltaA), -f (--deltaF) and -m (--deltaM) that adjust the proposal size for the allele
frequencies, inbreeding coefficients and migration rates, respectively.

**Auto-tuning (default):** By default, BA3 automatically adjusts the mixing parameters during
the burn-in phase to achieve optimal acceptance rates (target: 20-60%). The initial mixing
parameters (default 0.1) are shown in the pre-run summary, and the final auto-tuned values
are shown in the post-run summary. To disable auto-tuning, use the -N (--noautotune) option.

### 4.4 Options for printing output

By default, the output produced by BA3 is written to a file named BA3out.txt. An alternative name for the output file
can be specified using the option -o (--output). For example,
```
./BA3 -o myout.txt test.txt
```

**Allele frequencies file:** The option -F (--freqfile) specifies a separate file to write allele
frequencies in TSV (tab-separated values) format, which is convenient for further analysis:
```
./BA3 -F freqs.tsv -o results.txt mydata.txt
```
This creates a file with columns: Population, Locus, Allele, Frequency, SD

The option `-t (--trace)` specifies whether a trace output file is created that
lists all the parameter values at each iteration. The option `-g` causes detailed information
regarding the individual multilocus genotypes and posterior probabilities of migrant ancestries
to be written to a file named `BA3indiv.txt`.

### 4.5 VCF input options

For VCF file input, use the following options together:
- `-V filename.vcf` or `--vcf filename.vcf`: Specify the VCF input file
- `-M metadata.txt` or `--meta metadata.txt`: Specify the metadata file

Example:
```
./BA3 -V snps.vcf -M populations.txt -o results.txt -i1000000
```

## 5 Output file format

### 5.1 Main output file

The output file contains several sections:

#### Population Labels
```
 Population Labels:
  [0] pop0
  [1] pop1
  [2] pop2
```
Maps numeric indices to population names for compact display of the migration matrix.

#### Migration Rate Matrix
```
 Migration Rate Matrix m[source][dest]:
 Mean(SD)
                   [0]            [1]            [2]
  [0]  0.9718(0.0115) 0.0130(0.0100) 0.0152(0.0090)
  [1]  0.0878(0.0141) 0.7338(0.0396) 0.1784(0.0407)
  [2]  0.0870(0.0179) 0.2047(0.0326) 0.7083(0.0292)
```
Note that `m[i][j]` is the fraction of individuals in population i that
are migrants derived from population j (per generation).

#### Savage-Dickey Density Ratio Test

BA3 now includes a formal test for zero migration between populations:
```
 Savage-Dickey Density Ratio Test for Zero Migration:
 (Tests H0: m_ij = 0 vs H1: m_ij > 0)

 Source  Dest     Mean(SD)       BF_01  log10(BF)  KL(bits)  Interpretation
 ---------------------------------------------------------------------------
  [ 1]  [ 0]    0.01(0.01)     47.12    +1.67     2.80  Strong for H0
  [ 2]  [ 0]    0.02(0.02)     36.59    +1.56     2.46  Strong for H0
  ...
 ---------------------------------------------------------------------------
 Note: BF_01 > 1 supports H0 (no migration); BF_01 < 1 supports H1 (migration)
 KL(bits) = information gain from prior to posterior
```

The columns are:
- **Source/Dest**: Population indices for the migration rate being tested
- **Mean(SD)**: Posterior mean and standard deviation of the migration rate
- **BF_01**: Bayes factor comparing H0 (no migration) vs H1 (migration present)
- **log10(BF)**: Log10 of the Bayes factor
- **KL(bits)**: Kullback-Leibler information gain from prior to posterior
- **Interpretation**: Based on Jeffreys' scale (Weak, Substantial, Strong, Decisive)

#### Inbreeding Coefficients
```
 Inbreeding Coefficients:
 Index  Population                     F(SD)
 -----  ----------                     -----
 [ 0]  pop0                     0.2552(0.0367)
 [ 1]  pop1                     0.0810(0.0588)
 [ 2]  pop2                     0.2698(0.0891)
```

### 5.2 Allele frequencies file (-F option)

When using the -F option, allele frequencies are written to a separate TSV file:
```
Population	Locus	Allele	Frequency	SD
pop0	loc0	2	0.806	0.029
pop0	loc0	1	0.194	0.029
pop0	loc1	3	0.309	0.033
...
```

This format is convenient for import into R, Python, or spreadsheet programs.

## 6 Recommendations for running BA3

To generate correct results using BA3 it is important to use a sufficient number of iterations,
discard enough iterations as burn-in, and carry out several independent runs (started with different
random number seeds), examining the trace files for evidence of convergence and looking for
consistency of the estimates between independent runs.

### 6.1 Auto-tuning of mixing parameters

By default, BA3 automatically adjusts the mixing parameters during burn-in to achieve optimal
acceptance rates (20-60%). The pre-run summary shows the initial values with "(initial, autotune on)":
```
    Mixing:      dM=0.1 dA=0.1 dF=0.1 (initial, autotune on)
```

After the run completes, the final tuned values are displayed:
```
  Final mixing: dM=0.095 dA=0.312 dF=0.487
```

If you prefer to manually control the mixing parameters, use the -N option to disable auto-tuning
and specify your preferred values with -m, -a, and -f.

### 6.2 Diagnosing convergence

Two simple ways to examine convergence are:

- Conduct multiple runs initialized with different seeds and compare the posterior mean parameter estimates for concordance.
- Analyze the trace file for each run using the Tracer program.

### 6.3 Interpreting the individual ancestry output in BA3indiv.txt

If the command line option `-g` is used the file `BA3indiv.txt` is created containing an entry for each individual:
```
Individual: ind0 Source Popln: 0
Genotypes>>
loc0:2/2 loc1:3/3 loc2:2/2 loc3:2/2 loc4:?/? loc5:2/1 loc6:2/1
Migrant ancestry>>
[0,0]:0.931 [1,0]:0.000 [2,0]:0.000
[0,1]:0.000 [1,1]:0.002 [2,1]:0.001
[0,2]:0.000 [1,2]:0.035 [2,2]:0.031
```
The notation [i, j] indexes the population source i and generation j (0=nonmigrant, 1=1st
generation migrant, 2=second generation migrant) of migrant ancestry.

### 6.4 The Priors

The prior on allele frequencies is uniform Dirichlet. So, with two alleles the prior means
are $1/2$, with three alleles they are $1/3$ and so on. The prior distribution of the F statistic
is uniform on the interval (0, 1) with a mean of $1/2$. The prior on migration rates is
uniform with the constraint that $m_{ii} \geq 2/3$ and $\sum_i\sum_{j \neq i} m_{ij} \leq 1/3$. The prior means for $n$
populations are

$\overline{m}_{ii} = \frac{1}{n} + \frac{2}{3} \left(\frac{n-1}{n}\right),$

and

$\overline{m}_{ij} = \frac{1}{3n} \,\, \mathrm{for} \,\, i \neq j.$

The KL(bits) column in the Savage-Dickey test output indicates how much information the data
provides about each migration rate (higher values = more informative data).

## 7 Version History

- **3.4.0** (December 2025): Unified executable, VCF support, auto-tuning, Savage-Dickey test,
  improved output formatting, progress bar with ETA
- **3.0.5** (March 2023): Separate SNP/MSAT executables, bug fixes
- **3.0.0** (2007): Initial release of Edition 3

