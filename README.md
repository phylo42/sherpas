# SHERPAS

**Screening Historical Events of Recombination in a Phylogeny via Ancestral Sequences**

A new, alignment-free genome recombination detection tool exploiting the idea of phylo-kmers (originally developed in  from RAPPAS, Linard et al. 2019) to accelerate the process by several orders of magnitude while keeping comparable accuracy. 

__Reference:__

__*Scholz GE, Linard B, Romashchenko N, Rivals E, Pardi F. Rapid screening and detection of inter-(sub)type recombinants using phylo-k-mers. To be submitted to Bioinformatics (submitted)*__


 # Overview

__Inputs:__
- A reference phylogeny built with (ideally) non-recombinant genomes (or “pure types”) in __newick format__
- The multiple alignment from which was built this phylogeny, in __fasta format__.
- Table associating each tree leaf to a type (see examples below).
-  A set of query genomes for which recombination screening will be performed, in __fasta format__.


__Outputs:__
- Table of detected recombination patterns


## Construction of a phylo-kmer database

The construction of a **phylo-kmer database**, built from a reference phylogeny, is performed by RAPPAS (see https://www.ncbi.nlm.nih.gov/pubmed/30698645). The database builder and its documentation are available at https://github.com/phylo42/rappas2.

# Usage

## Compilation and rapid test

**Prerequisites:**

* Make sure Boost Librairies >=1.65 are present.
* Your GCC compiler must support c++17
* CMake >= 3.10 installed

In debian, these can be installed with:
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
Stil from the `release-build` repository in which you compiled sources, execute the following commands:

```shell
# download a SHERPAS database already pre-built from HBV pure types: 

wget https://www.dropbox.com/s/m75hfo4mem4eb46/pkDB-HBV-full.zip
unzip pkDB-HBV-full.zip

# launch a prediction for 3000 HBV queries, using the pre-built HBV database:

sherpas/SHERPAS -d DB_k10_o1.5.rps -q ../examples/HBV_all/queries-3000.fasta -o prefix_ -g ref-groups.csv -c
```

These 3000 queries should be analyze in less than 5 minutes (using a 3Ghz i7 CPU).
You should obtain the following file in the same directory :

``` shell
# the results itself, e.g a the list of recombinant regions detected for each query:
prefix_res-queries-3000.txt 

# the queries in fasta format, matching the coordinates of prefix_res-queries-3000.txt 
prefix_queries-3000-circ300.fasta
```
More pre-built database (those used in SHERPAS manuscript) can be downloaded from Dryad:
https://datadryad.org/stash/downloads/file_stream/373882.

## SHERPAS Execution

```shell
sherpas/SHERPAS [options] 
```

Command-line options are the following (see detailed description below):

Option | Description | Default value
--- | --- | ---
**-d** | path to the database | None (mandatory field)
**-q** | path to the queries file | None (mandatory field)
**-g** | path to the group-assignment file | None (mandatory field)
**-o** | path to the output directory | None (mandatory field)
**-w** |window size (>99)	 | 500
**-m** |method used (F or R) | F
**-t** | threshold for signal control | 10 (F) or 0.9 (R)
**-c** | activates circularity options | None (none expected)
**-l** | changes output layout | None (none expected)
**-k** | disables N/A regions post treatment | None (none expected)

## Mandatory options

These are the three necessary files used by SHERPAS to run, plus the desired output directory. Error messages will appear if one of them is missing, in which case the program will abort.

**-d** : The path to the .rps database. All the information on a reference alignment and the corresponding phylogeny used by SHERPAS are stored in a phylo-kmer database, that needs to be built prior to using SHERPAS (https://github.com/phylo42/rappas2). Once built, the database is serialized and stored as a .rps file. Thus, this building step only needs to be performed once for a given reference alignment and a given kmer size k (note also that the database can independently be used with RAPPAS and SHERPAS).

**-q** : The path to a .fasta file of query sequences. The sequences in this dataset will be investigated individually by SHERPAS, using the information contained in the database.

**-g** : The path to a .csv file defining the mapping between groups and genomes in the reference alignment used to build the database. This file consists of two columns, the first one containing the name of the sequences in the reference alignment (these names must match the names of in the alignment file used to build the database), the second one the corresponding groups.
Important note on typography: for technical reasons, neither the sequence names nor the groups names should contain a comma (‘,’). Moreover, group names should not contain a star (‘*’).

*Example:*
```
ref_1,A
ref_2,A
ref_3,B
ref_4,C
```

**-o** : The path to the output directory, where the result files will be created.

## Optional parameters

The following are numerical parameters with pre-set default values that produced good a good balance between sensitivity and precision for viral genomes such as HIV or HCV genomes. Different values may produce different levels of accuracy and with other references, it is advised to explore more values .

**-w** : The size, in k-mers, of the sliding window (default value 500). 

A larger window size will make the information more precise, but short recombinant segments might remain undetected. The minimum possible value is 100, and this value should not exceed the length in k-mers* of the shorter sequence in the query file.
*The number of kmers in a sequence of length L is L+1-k, where k is the k-mer size.

**-m** : The method used by SHERPAS (default value F). 

Currently two options are possible, “full” (called with F) and “reduced” (called with R). While “full” is usually more precise, “reduced” is much faster.

**-t** : The threshold for post-control (default value 10 for method F and 0.9 for method R). 

A section of a query can be returned as “unassigned” if the evidence for any particular group is too weak. This threshold governs the definition of “too weak”, in a way that depends on the method chosen. For method F, the parameter controlled by the threshold is the ratio of the best score over the second best score, thus the threshold should be either zero or greater than one (any threshold between zero and one has the same effect as a threshold of zero) . For method R, it is the ratio of the best score over the sum of all scores, so the threshold must be comprised between zero and one in that case (a threshold greater than one will return “unassigned” for the whole sequence, as the ratio computed is always smaller than one). For both methods, a threshold of zero means that no such control is operated.

## Advanced customization

These are options that are disabled by default but may be useful in some circonstances.

**-c** : Use when queries are full, circular genomes.

**-l** : Changes the output mode to “linear” (see “Output” section below).

**-k** : Skips the post-processing part of unassigned segments.


#  Output

A .txt file is written to the output directory given with -o.
The header of the file summarizes information about the query file and the parameters used.
Then for each query, a first line indicates the name of the query (e.g. its fasta header), and each subsequent line gives coordinates of  a distinct sequence segment and the group to which it was assigned, in the form:
```
[first_position]	[last_position]	[group]
```

*Example:*
```
>query_1
1	4593	B
4594	7286	A
7287	8663	B
```

Alternatively (when option -l is activated), the output for one given query is shrinked to a single comma separated line, in the form:
```
0,[group_1],[breakpoint_1],[group_2],...,[breakpoint_n-1],[group_n],[end_of_sequence_position]
```

*Example:*
```
>query_1
0,B,4594,A,7287,B,8663
```

-Once all the queries have been processed, the time taken by the program to run (not including the time taken to read all the necessary files or to build the database, but including the time taken to write the output file) is printed in the console.
