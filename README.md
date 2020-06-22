# SHERPAS

**Screening Historical Events of Recombination in a Phylogeny via Ancestral Sequences**

A new, alignment-free genome recombination detection tool exploiting the idea of phylo-kmers (originally developed in [RAPPAS](https://github.com/phylo42/RAPPAS), Linard et al. 2019) to accelerate the process by several orders of magnitude while keeping comparable accuracy. 

__Reference:__

*Scholz GE, Linard B, Romashchenko N, Rivals E, Pardi F. Rapid screening and detection of inter-type viral recombinants using phylo-k-mers. (Submitted)*

 - [Overview](#overview)
 - [Usage](#usage)
   - [Compilation and rapid test](#compilation-and-rapid-test)
 - [SHERPAS Execution](#sherpas-execution)
   - [Mandatory options](#mandatory-options)
   - [Optional parameters](#optional-parameters)
   - [Advanced customization](#advanced-customization)
 - [Outputs](#outputs)


 # Overview

__Inputs:__
- A phylo-kmer database (usually having a __.rps__ extension) constructed from a reference alignment and tree (see below).
- A table associating each sequence in the reference alignment to a type (a.k.a. "strain"), in __csv__ format.
- A set of query sequences (e.g. whole genomes, long reads) for which recombination screening will be performed, in __fasta__ format.

__Outputs:__
- Table of recombination patterns detected in the query sequences.

### Construction of a phylo-kmer database (pkDB)

In the absence of a pre-computed pkDB, you will need the following to construct your own pkDB: 

- *Reference alignment*: a multiple alignment containing a number of reference sequences from each of the strains, in __fasta__ format. These sequences should be pure (i.e. non recombinant) with respect to their strain.
- *Reference tree*: a phylogenetic tree built from the reference alignment, in __newick__ format.

The construction of the pkDB should be performed with [RAPPAS2](https://github.com/phylo42/rappas2). Its code and documentation are available at https://github.com/phylo42/rappas2.

# Usage

## Compilation and rapid test

**Prerequisites:**

* Make sure Boost Librairies >=1.65 are present.
* Your GCC compiler must support c++17
* CMake >= 3.10 installed

In Debian, these can be installed with:
```
sudo apt install build-essential
sudo apt install cmake
sudo apt install libboost-all-dev
```

**Compilation:**

Clone the repository with recursive option:

```shell
git clone --recursive https://github.com/phylo42/sherpas.git        #do not forget the --recursive !!!
cd sherpas
mkdir release-build
cd release-build
cmake ..
make
```

**Rapid test:**

A rapid prediction of HBV (Hepathitis B Virus) recombinants can be performed.
From the `release-build` repository in which you compiled sources, execute the following commands:

```shell
# download a SHERPAS database already pre-built from HBV pure types: 
wget https://www.dropbox.com/s/m75hfo4mem4eb46/pkDB-HBV-full.zip
unzip pkDB-HBV-full.zip

# launch a prediction for 3000 HBV queries, using the pre-built HBV database:
sherpas/SHERPAS -d DB_k10_o1.5.rps -q ../examples/HBV_all/queries-3000.fasta -o prefix_ -g ref-groups.csv -c
```

These 3000 queries should be analyzed in less than 5 minutes (using a 3Ghz i7 CPU).
You should obtain the following file in the same directory :

``` shell
# the results themselves, e.g a the list of recombinant regions detected for each query:
prefix_res-queries-3000.txt 

# the queries in fasta format, matching the coordinates of prefix_res-queries-3000.txt 
prefix_queries-3000-circ300.fasta
```
More pre-built phylo-kmer databases (those used in the SHERPAS manuscript) can be downloaded from Dryad:
https://datadryad.org/stash/downloads/file_stream/373882.

# SHERPAS Execution

```shell
sherpas/SHERPAS [options] 
```

Command-line options are the following (see detailed description below):

Option | Description | Default value
--- | --- | ---
**-d** | path to the phylo-kmer database | None (mandatory field)
**-g** | path to the strain-assignment file | None (mandatory field)
**-q** | path to the queries file | None (mandatory field)
**-o** | prefix for the output file(s) | None (mandatory field)
**-w** | window size (>99)	 | 300
**-m** | method used (F or R) | F
**-t** | threshold for unassigned regions | 100 (F) or 0.99 (R)
**-c** | activates circularity options | None (none expected)
**-l** | changes output layout | None (none expected)
**-k** | disables N/A regions post treatment | None (none expected)
<!--- **-o** | path to the output directory | None (mandatory field) --->

## Mandatory options

**-d** : The path to the .rps file containing the phylo-kmer database. All the information on a reference alignment and the corresponding phylogeny used by SHERPAS are stored in this file, that needs to be built prior to using SHERPAS (with [RAPPAS2](https://github.com/phylo42/rappas2)). This building step only needs to be performed once for a given reference alignment. 
<!---
and a given kmer size k (note also that the database can independently be used with RAPPAS and SHERPAS).
--->

**-g** : The path to a .csv file that specifies the mapping between strains and sequences in the reference alignment.
This file consists of two columns, the first one containing the name of the sequences in the reference alignment (these names must match the names of in the alignment file used to build the database), the second one the corresponding strains.
Important note on typography: for technical reasons, neither the sequence names nor the strains names should contain a comma (‘,’). Moreover, strain names should not contain a star (‘*’).

*Example:*
```
ref_1,A
ref_2,A
ref_3,B
ref_4,C
```

**-q** : The path to a .fasta file of query sequences. The sequences in this dataset will be investigated individually by SHERPAS, using the information in the two files specified by the two options above.

**-o** : A prefix for the output file(s). This can also be used to specify an output directory (end the prefix with ‘/’).
See section [Outputs](#outputs) below for details on the output format.
<!--- The path to the output directory, where the result files will be created. --->

## Optional parameters

The parameters below influence the behavior of SHERPAS both in terms of accuracy and running times.
We refer the user to the manuscript for a discussion of how the parameters below should be adjusted, 
depending on the different use cases for SHERPAS. 

**-w** : The size, in number of *k*-mers, of the sliding window (default value 300). 

The window size controls how easily SHERPAS switches between different strain classifications along the query. 
Smaller windows tend to produce more fragmented partitions of the queries, thus increasing the detection of short recombinant segments, but also that of false positive recombinants.
The minimum possible value is 100, and the maximum is the length in k-mers* of the shortest sequence in the query file.
(*The number of kmers in a sequence of length L is L+1-k, where k is the k-mer size.) 
We advise against any value smaller than 200 or larger than 600. 

**-m** : The method used by SHERPAS (default value F). 

Currently two options are possible, “full” (called with F) and “reduced” (called with R). While “full” is usually more accurate, “reduced” is much faster.

**-t** : The threshold controlling unassigned regions (default value 100 for method F and 0.99 for method R). 

A section of a query can be returned as “unassigned” (N/A) if the evidence for any particular strain is too weak. This threshold governs the definition of “too weak”, with high thresholds generally producing larger unassigned regions, i.e. more conservative strain assignments. The precise meaning of this threshold depends on the chosen method: 
In method F, in order for a region to be assigned to a strain, the ratio of the best score over the second best score must be larger than the threshold. In method R, it is the ratio of the best score over the sum of all scores that must be larger than the threshold. 
<!---
For method F, the parameter controlled by the threshold is the ratio of the best score over the second best score, thus the threshold should be either zero or greater than one (any threshold between zero and one has the same effect as a threshold of zero). For method R, it is the ratio of the best score over the sum of all scores, so the threshold must be comprised between zero and one in that case (a threshold greater than one will return “unassigned” for the whole sequence, as the ratio computed is always smaller than one). For both methods, a threshold of zero means that no such control is operated.
--->

## Advanced customization

These are options that are disabled by default but may be useful in some circumstances.

**-c** : Use when queries are whole circular genomes.

**-l** : Changes the output mode to “linear” (see Output section below).

**-k** : Skips the post-processing part of unassigned segments.


#  Outputs

The results are written to a .txt file with the following format.
The header of the file summarizes information about the query file and the parameters used.
Then for each query, a line indicates the name of the query (its fasta header), 
and each subsequent line gives coordinates of a distinct sequence segment and the strain to which it was assigned, in the form:
```
[first_position] [last_position]	[strain]
```
The coordinates are 1-based and are relative to the query sequences with gaps (if any) removed. This is the same format as that used by [jpHMM](http://jphmm.gobics.de/), which facilitates comparisons.  

*Example:*
```
>query_1
1	4593	B
4594	7286	A
7287	8663	B
```

Alternatively (when option -l is activated), the output for one given query is reduced to a single comma separated line, in the form:
```
0,[strain_1],[breakpoint_1],[strain_2],...,[breakpoint_n-1],[strain_n],[end_of_sequence_position]
```

*Example:*
```
>query_1
0,B,4594,A,7287,B,8663
```

During the program execution, some information about its progress are written to the console. The running time employed to read all input files, and that for all subsequent operations are reported separately.
