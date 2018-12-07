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

Whatever goes here

## Evolutionary Pathways

Whatever goes here

## Correlation of Fitness Effects

The ruggedness of a ribozyme family can be measured usign the fitness correlation <a href="https://www.codecogs.com/eqnedit.php?latex=\gamma&space;_{d}" target="_blank"><img src="https://latex.codecogs.com/gif.latex?\gamma&space;_{d}" title="\gamma _{d}" /></a>, 
which is the average correlation of activity effects of single mutations in *d*-mutant neighbors
(where *d* is the Levenshtein edit distance, i.e., the number of substitutions, insertions or 
deletions between two related sequences).

The correlation of fitness effects, 
<a href="https://www.codecogs.com/eqnedit.php?latex=\gamma&space;_{d}" target="_blank"><img src="https://latex.codecogs.com/gif.latex?\gamma&space;_{d}" title="\gamma _{d}" /></a>, 
is calculated following [Ferretti et al., 2016](https://www.sciencedirect.com/science/article/pii/S0022519316000771?via%3Dihub).

The Python script (`ActivityCorrelationGamma.py`) and excel source file (`ActivityObservedData.xlsx`) can be found in the folder `Gamma_correlation`.

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

