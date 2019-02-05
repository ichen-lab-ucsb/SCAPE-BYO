# Mapping a comprehensive ribozyme fitness landscape reveals a frustrated evolutionary network for self-aminoacylating RNA

Source files and python scripts used for the publication 
**Mapping a comprehensive ribozyme fitness landscape reveals a frustrated evolutionary network for self-aminoacylating RNA**
by Abe D. Pressman, Ziwei Liu, Evan Janzen, Celia Blanco, Ulrich F. Muller, Gerald F. Joyce, Robert Pascal and Irene A. Chen.

## Table of Contents

- [Prerequisites](#prerequisites)
- [k-seq Analysis](#installation)
- [Evolutionary Pathways](#features)
- [Correlation of Fitness Effect](#contributing)

## Prerequisites

What things you need to install the software and how to install them:

- Python 2.7  
- openpyxl package for python
- Levenshtein package for python 

## k-seq Analysis
A tool that calculates A and k times t, according to the equation A(1-Exp(-k[S]t)), for every sequence in a relevant k-seq data set. Requires a fairly large number of inputs, including sequence data/counts for a k-seq "start" round, a file for each tested k-seq round after selection has occurred, files describing the reaction conditions, and known or approximate normalization constants for each round based on the expected amount of DNA/RNA/protein present before and after the reaction in each sample.

### How to use the script to calculate k-seq results:

An example case to generate all k-seq data used in this publication are as follows:

```
python kseq_tools_v01.py R5c-counts.txt example-rounds.txt [output_file] example-normalization.txt example-subst-concs.txt example-rnds-to-avg.txt example-rnds-to-err.txt -v -p
```

It is NOT recommended to run the example case, as it will take 10-30 hours on a standard personal computer. Instead, consider one of the (much faster) following...

Calculate k-seq data for all sequences present at count ≥ 10 in round 5 of the selection:
```
python kseq_tools_v01.py R5c-counts.txt example-rounds.txt [output_file] example-normalization.txt example-subst-concs.txt example-rnds-to-avg.txt example-rnds-to-err.txt -v -m 10 -p
```

Calculate k-seq data for all sequences present in family 1A.1:
```
python kseq_tools_v01.py R5c-counts.txt example-rounds.txt [output_file] example-normalization.txt example-subst-concs.txt example-rnds-to-avg.txt example-rnds-to-err.txt -v -s center CTACTTCAAACAATCGGTCTG 3 -p
```

Current input file requirements for each rounds are "counts" file consisting of three lines of metadata followed by  one line per unique sequence in pool, of the format "sequence count" where count is an integer. Such files are produced by our Galaxy tools, currently available at http://galaxy-chen.cnsi.ucsb.edu:8080/. Future versions of this script will be included with tools to more easily process data from a number of other tools currently used to process high-throughput sequencing data; however, the version of the peak_pather script used in this publication (v0.1) will remain here for posterity.

For usage/argument details, run `python kseq_tools_v01.py -h`


## Evolutionary Pathways
The shortest pathway between two sequences, along a fitness landscape, can be found efficiently using the A* algorithm, and iterating over a single round's sequence population as a graph, with each sequence present as its own node and the edit distances between sequences as edges with distance. The script provided here finds the N best pathways between two sequences, ranked by 1) Shortest total path length, 2) Smallest maximum step size, 3) Smallest average step size, 4) Largest minimum sequence count, using a sequence counts file as the reference map (this could easily be replaced with highest minimum fitness, if the reference file is a list of sequences and their fitnesses instead).

Current input file requirements are a "counts" file consisting of three lines of metadata followed by  one line per unique sequence in pool, of the format "sequence count" where count is an integer. Such files are produced by our Galaxy tools, currently available at https://labs.chem.ucsb.edu/chen/irene/Chen_lab_at_UCSB/Galaxy_Tools.html. Future versions of this script will be included with tools to more easily process data from
a number of other tools currently used to process high-throughput sequencing data; however, the version of the peak_pather script used in this publication (v0.1) will remain here for posterity.

### How to use the script to calculate evolutionary pathways:

The pathways described in this publication are calculated as follows (with the flag -p to output progress data to the terminal, as the script can take a while to run):

```
python peak_pather_v01.py R5c-counts.txt [output_file] [start_sequence] [end_sequence] --min_count [variable] --max_step [variable] -n 5 --max_length 35 -p
```
Minimum sequence count was set at 3 for the highly-populated pathways between Motif 1A and 1B, and set at 2 for all other pathways. For all pairs of sequence endpoints investigated, maximum step size was set at 1, then incremented by 1 until 5 pathways were found; the script was then run again with min_step increased 1 further, to generate 5 additional pathways with larger step tolerance.

For usage/argument details, run `python peak_pather_v01.py -h`

## Correlation of Fitness Effects

The ruggedness of a ribozyme family can be measured usign the fitness correlation <a href="https://www.codecogs.com/eqnedit.php?latex=\gamma&space;_{d}" target="_blank"><img src="https://latex.codecogs.com/gif.latex?\gamma&space;_{d}" title="\gamma _{d}" /></a>, 
which is the average correlation of activity effects of single mutations in *d*-mutant neighbors
(where *d* is the Levenshtein edit distance, i.e., the number of substitutions, insertions or 
deletions between two related sequences).

The correlation of fitness effects, 
<a href="https://www.codecogs.com/eqnedit.php?latex=\gamma&space;_{d}" target="_blank"><img src="https://latex.codecogs.com/gif.latex?\gamma&space;_{d}" title="\gamma _{d}" /></a>, 
is calculated following [Ferretti et al., 2016](https://www.sciencedirect.com/science/article/pii/S0022519316000771?via%3Dihub).

The Python script (`ActivityCorrelationGamma.py`) and excel source file (`ActivityObservedData.xlsx`) can be found in the folder `Gamma_correlation`. Each spreadsheet in the excel file corresponds to a different family and follows the same format: column A) sequences, column B) calculated activity (including values under the baseline activity), column C) calculated activity if above the baseline activity (or baseline activity if below), and column D) logarithm of the values in column C).

### How to use the script to calculate correlation of fitness effects:

```
python ActivityCorrelationGamma.py ActivityObservedData.xlsx sheet_name correlation_distance
```

where `sheet_name` corresponds to either `Family_2.1`, `Family_1A.1`, `Family_1B.1`, `Family_1B.2` or `Family_1A.2`.

## Built With

* [Python](https://www.python.org/) - The Programming Language

## Authors

* **Abe Pressman** - *k-seq Analysis* - *Evolutionary Pathways* 
* **Celia Blanco** - *Correlation of Fitness Effects* 

