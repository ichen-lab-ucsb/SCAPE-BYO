# Systematic mapping of a ribozyme fitness landscape reveals a frustrated evolutionary network for self-aminoacylating RNA

Source files and python scripts used for the publication 
**Mapping a comprehensive ribozyme fitness landscape reveals a frustrated evolutionary network for self-aminoacylating RNA**
by Abe D. Pressman, Ziwei Liu, Evan Janzen, Celia Blanco, Ulrich F. Muller, Gerald F. Joyce, Robert Pascal and Irene A. Chen.

## Table of Contents

- [Prerequisites](#prerequisites)
- [k-seq Analysis](#installation)
- [Evolutionary Pathways](#features)
- [Correlation of Fitness Effect](#contributing)

## Prerequisites

- Python 2.7  
- openpyxl package for python
- Levenshtein package for python 
- numpy package for python 
- scipy package for python 
- argparse package for python 
- heapq package for python 

## *k*-Seq Analysis
This tool calculates *A* and *k&ast;t*, according to the equation *A(1-Exp(-k[S]t))*, for every sequence in a relevant *k*-Seq data set. It requires several input files: a file with sequence data/counts for a *k*-Seq "start" round, a file for each tested *k*-Seq round after selection has occurred, files describing the reaction conditions, and known or approximate normalization constants for each round based on the expected amount of DNA/RNA/protein present before and after the reaction in each sample.


### How to use the script to calculate *k*-Seq results:

To reproduce the numerical results reported in this publication, the python script `kseq_tools_v01.py` can be run as follows:

```
python kseq_tools_v01.py start_round kseq_rounds output_file normalization_list substrate_concs rounds_to_average rounds_to_error -v -p
```

where `start_round` is R5c-counts.txt, `kseq_rounds` is example-rounds.txt, `normalization_list` is example-normalization.txt, `substrate_concs` is example-subst-concs.txt, `rounds_to_average` is example-rnds-to-avg.txt and `rounds_to_error` is example-rnds-to-err.txt.

The code requires "counts" files as input, consisting of three lines of metadata followed by one line per unique sequence in the pool in the following format: sequences in the first column and counts (an integer number) in the second column. Such files are produced by our [Galaxy tools](https://labs.chem.ucsb.edu/chen/irene/Chen_lab_at_UCSB/Publications_files/Xulvi%20et%20al%20Methods%202016.pdf), currently available at the [Chen Lab website](https://labs.chem.ucsb.edu/chen/irene/Chen_lab_at_UCSB/Galaxy_Tools.html). The files corresponding to the start round (R5c-counts.txt) and each tested *k*-Seq round after selection (called in `kseq_rounds`) will be uploaded to a perpetual repository upon acceptance of the publication. A link to the repository will be added here. Every other input file can be found in the folder [kseq_tools](https://github.com/ichen-lab-ucsb/SCAPE-BYO/tree/master/kseq_tools). 

For more information on usage, see the [detailed readme file](https://github.com/ichen-lab-ucsb/SCAPE-BYO/blob/master/kseq_tools/README.md) or run `python kseq_tools_v01.py -h` in the terminal.


## Evolutionary Pathways

The shortest pathway between two sequences, along a fitness landscape, can be found efficiently using the [A* algorithm](https://en.wikipedia.org/wiki/A*_search_algorithm). The algorithm iterates over a single round's sequence population as a graph, with each sequence present as its own node and the edit distances between sequences as edges with distance. The python script `peak_pather_v01.py` finds the *N* best pathways between two sequences, ranked by: 1) shortest total path length, 2) smallest maximum step size, 3) smallest average step size, 4) largest minimum sequence count, using a sequence counts file as the reference map.

### How to use the script to calculate evolutionary pathways:

The pathways described in this publication are computed as follows:

```
python peak_pather_v01.py input_file output_file start_sequence end_sequence --min_count [variable] --max_step [variable] -n [variable] --max_length [variable] -p
```

where `input_file` is R5c-counts.txt. The code requires a "counts" file as input, consisting of three lines of metadata followed by  one line per unique sequence in the pool in the following format: sequences in the first column and counts (an integer number) in the second column. Such files are produced by our [Galaxy tools](https://labs.chem.ucsb.edu/chen/irene/Chen_lab_at_UCSB/Publications_files/Xulvi%20et%20al%20Methods%202016.pdf), currently available at the [Chen Lab website](https://labs.chem.ucsb.edu/chen/irene/Chen_lab_at_UCSB/Galaxy_Tools.html). 
The file corresponding to the round used to produce the pathways reported in the publication (R5c-counts.txt) will be uploaded to a perpetual repository upon acceptance of the publication. A link to the repository will be added here. 

Minimum sequence count (`min_count`) was set at 3 for the highly-populated pathways between Motif 1A and 1B, and set at 2 for all other pathways. For all pairs of sequence endpoints investigated, maximum step size (`max_step`) was set at 1, then incremented by 1 until 5 pathways were found. The script was then run again with min_step increased 1 further, to generate 5 additional pathways with larger step tolerance.

The number of best pathways to be generated with each run (`n`) was set to 5. The maximum length of any path searched before the path is abandoned (`max_length`) was set to 35.

For more information on usage, see the [detailed readme file](https://github.com/ichen-lab-ucsb/SCAPE-BYO/blob/master/peak_pather/README.md) or run `python peak_pather_v01.py -h` in the terminal.

## Correlation of Fitness Effects

The ruggedness of a ribozyme family can be measured using the fitness correlation <a href="https://www.codecogs.com/eqnedit.php?latex=\gamma&space;_{d}" target="_blank"><img src="https://latex.codecogs.com/gif.latex?\gamma&space;_{d}" title="\gamma _{d}" /></a>, 
which is the average correlation of activity effects of single mutations in *d*-mutant neighbors
(where *d* is the Levenshtein edit distance, i.e., the number of substitutions, insertions or 
deletions between two related sequences).

The correlation of fitness effects, 
<a href="https://www.codecogs.com/eqnedit.php?latex=\gamma&space;_{d}" target="_blank"><img src="https://latex.codecogs.com/gif.latex?\gamma&space;_{d}" title="\gamma _{d}" /></a>, 
is calculated following [Ferretti et al., 2016](https://www.sciencedirect.com/science/article/pii/S0022519316000771?via%3Dihub).

### How to use the script to calculate correlation of fitness effects:

To reproduce the numerical results reported in the publication, run:

```
python ActivityCorrelationGamma.py input_file sheet_name correlation_distance
```

where `input_file` corresponds to the excel file ActivityObservedData.xlsx and `sheet_name` corresponds to either Family_2.1, Family_1A.1, Family_1B.1, Family_1B.2 or Family_1A.2. `Correlation_distance` is an integer number {0..4}.

The Python script `ActivityCorrelationGamma.py` and excel source file `ActivityObservedData.xlsx` can be found in the folder [Gamma_correlation](https://github.com/ichen-lab-ucsb/SCAPE-BYO/tree/master/Gamma_correlation). Each spreadsheet in the excel file corresponds to a different family and follows the same format: column A) sequences, column B) calculated activity (including values under the baseline activity), column C) calculated activity if above the baseline activity (or baseline activity if below), and column D) logarithm of the values in column C).


## Built With

* [Python](https://www.python.org/) - The Programming Language

## Authors

* **Abe Pressman** - *k*-Seq Analysis - Evolutionary Pathways
* **Celia Blanco** - Correlation of Fitness Effects 

