# SHERPAS

**Screening Historical Events of Recombination in a Phylogeny via Ancestral Sequences**
A new, alignment-free genome recombination detection tool exploiting the idea of phylo-kmers (originally developed in  from RAPPAS, Linard et al. 2019) to accelerate the process by several orders of magnitude while keeping comparable accuracy. 

__Reference:__
__Scholz GE, Linard B, Rivals E, Pardi F. Rapid screening and detection of inter-(sub)type recombinants using phylo-k-mers. To be submitted to Bioinformatics (submitted)__

→ si vous le mettez sur biorxiv après soumission, mettre le DOI associé ?


 # Overview

__Inputs:__
- A reference phylogeny built with (ideally) non-recombinant genomes (or “pure types”) in __newick format__
- The multiple alignment from which was built this phylogeny, in __fasta format__.
- Table associating each tree leaf to a type (see examples below).
-  A set of query genomes for which recombination screening will be performed, in __fasta format__.


__Outputs:__
- Table of detected recombination patterns


## Construction of a phylo-kmer database

The construction of a **phylo-kmer database**, built from a reference phylogeny, is performed by  RAPPAS ( see https://www.ncbi.nlm.nih.gov/pubmed/30698645 and . This approach is alignment-free, in the sense that future query sequences do not need to be aligned with the reference alignment. Instead, query k-mers will be searched against the phylo-kmer database, in which each k-mer is associated to phylogeny branches and  ta probability to belong to this branch. 
##  SHERPAS algorithm: general idea

In SHERPAS, the phylo-kmer database built by RAPPAS is used for the purpose of inter-subtype recombination detection. In essence, SHERPAS follows a "sliding-window" approach, that consists in investigating short subsequences (the so-called “windows”) of the query sequence in turn. For each window, SHERPAS analyzes k-mer matches between query sequence and phylo-kmer database to assign a best-matching clade among a set of user-defined clades (of the reference phylogeny). If, along the sequence, the best-matching clade changes significantly, this can be taken as a potential signal for recombination.





# Usage

##Compilation and execution

**Compilation:**

Clone the repository with recursive option:

```shell
git clone --recursive https://ADRESS_A_RAJOUTER
cdADRESS_A_RAJOUTER
```

In the cloned repository, successively run:

```shell
mkdir release-build
cmake
cd release-build
make
```

**Run SHERPAS:**

In the release-build repository:
```shell
recomb/rappas-recombn [options] 
```

Command-line options are the following (see detailed description below):

Option | Description | Default value
--- | --- | ---
-d | path to the database | None (mandatory field)
-q | path to the queries file | None (mandatory field)
-g | path to the group-assignment file | None (mandatory field)
-o | path to the output directory | None (mandatory field)
-w |window size (>99)	 | 500
-m |method used (F or R) | F
-t | threshold for signal control | 10 (F) or 0.9 (R)
-c | activates circularity options | None (none expected)
-l | changes output layout | None (none expected)
-k | disables N/A regions post treatment | None (none expected)

## Mandatory options

These are the three necessary files used by SHERPAS to run, plus the desired output directory. Error messages will appear if one of them is missing, in which case the program will abort.

-d : The path to the .rps database. All the information on a reference alignment and the corresponding phylogeny used by SHERPAS are stored in a phylo-kmer database, that needs to be built prior to using SHERPAS (https://gite.lirmm.fr/rappas-team/rappas-build, see [link to rappas-build documentation]). Once built, the database is serialized and stored as a .rps file. Thus, this building step only needs to be performed once for a given reference alignment and a given kmer size k (note also that the database can independently be used with RAPPAS and SHERPAS).

-q : The path to a .fasta file of query sequences. The sequences in this dataset will be investigated individually by SHERPAS, using the information contained in the database.

-g : The path to a .csv file defining the mapping between groups and genomes in the reference alignment used to build the database. This file consists of two columns, the first one containing the name of the sequences in the reference alignment (these names must match the names of in the alignment file used to build the database), the second one the corresponding groups.
Important note on typography: for technical reasons, neither the sequence names nor the groups names should contain a comma (‘,’). Moreover, group names should not contain a star (‘*’).

Example:
ref_1,A
ref_2,A
ref_3,B
ref_4,C

-o : The path to the output directory, where the result files will be created.

## Optional parameters

The following are numerical parameters with pre-set default values that produced good a good balance between sensitivity and precision for viral genomes such as HIV or HCV genomes. Different values may produce different levels of accuracy and with other references, it is advised to explore more values .

-w : The size, in k-mers, of the sliding window (default value 500). 

A larger window size will make the information more precise, but short recombinant segments might remain undetected. The minimum possible value is 100, and this value should not exceed the length in k-mers* of the shorter sequence in the query file.
*The number of kmers in a sequence of length L is L+1-k, where k is the k-mer size.

-m : The method used by SHERPAS (default value F). 

Currently two options are possible, “full” (called with F) and “reduced” (called with R). While “full” is usually more precise, “reduced” is much faster.

-t : The threshold for post-control (default value 10 for method F and 0.9 for method R). 

A section of a query can be returned as “unassigned” if the evidence for any particular group is too weak. This threshold governs the definition of “too weak”, in a way that depends on the method chosen. For method F, the parameter controlled by the threshold is the ratio of the best score over the second best score, thus the threshold should be either zero or greater than one (any threshold between zero and one has the same effect as a threshold of zero) . For method R, it is the ratio of the best score over the sum of all scores, so the threshold must be comprised between zero and one in that case (a threshold greater than one will return “unassigned” for the whole sequence, as the ratio computed is always smaller than one). For both methods, a threshold of zero means that no such control is operated.

## Advanced customization

These are options that are disabled by default but may be useful in some circonstances.

-c : Use when queries are full, circular genomes.

-l: Changes the output mode to “linear” (see “Output” section below).

-k : Skips the post-processing part of unassigned segments.


#  Output

A .txt file is written to the output directory given with -o.
The header of the file summarizes information about the query file and the parameters used.
Then for each query, a first line indicates the name of the query (e.g. its fasta header), and each subsequent line gives coordinates of  a distinct sequence segment and the group to which it was assigned, in the form:
[first_position]	[last_position]	[group]

Example:
>query_1
1	4593	B
4594	7286	A
7287	8663	B


Alternatively (when option -l is activated), the output for one given query is shrinked to a single comma separated line, in the form:
0,[group_1],[breakpoint_1],[group_2],...,[breakpoint_n-1],[group_n],[end_of_sequence_position]

Example:
>query_1
0,B,4594,A,7287,B,8663

-Once all the queries have been processed, the time taken by the program to run (not including the time taken to read all the necessary files or to build the database, but including the time taken to write the output file) is printed in the console.
