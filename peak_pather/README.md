
## Tool for finding pathways between peak sequences in a fitness landscape

We expect future versions of this script to be updated/included with additional tools. As it is used in a publication, however, v0.1 of the peak_pather script will remain here for posterity.

### Goal:

The peak_pather script searches for the shortest pathway between two sequences, along a fitness landscape, using the A* algorithm (a decent description can be found on [Wikipedia](https://en.wikipedia.org/wiki/A*_search_algorithm). The script iterates over a single round's sequence population as a graph, with each sequence present as its own node and the edit distances between sequences as edges with distance. The script provided here finds the N best pathways between two sequences, ranked by 1) Shortest total path length, 2) Smallest maximum step size, 3) Smallest average step size, 4) Largest minimum sequence count, using a sequence counts file as the reference map (this could easily be replaced with highest minimum fitness, if the reference file is a list of sequences and their fitnesses instead).

### Input:

The script can be called as follows:

```
python peak_pather_v01 input output start_seq end_seq
```
#### Required arguments (positionally dependent):
`input`                 Location/name of file containing sequence counts for the round to be iterated over (e.g. `R5c-counts.txt`). Currently requires a "counts" file consisting of three lines of metadata followed by  one line per unique sequence in pool, of the format "sequence count" where count is an integer. Such files are produced by our Galaxy tools, currently available at https://labs.chem.ucsb.edu/chen/irene/Chen_lab_at_UCSB/Galaxy_Tools.html. 

`output`                Output file location/name (e.g. `paths-output.txt`)

`start_seq`             Path start sequence (e.g. 'CTACTTCAAACAATCGGTCTG')

`end_seq`               Path end sequence (e.g. 'CCACACTTCAAGCAATCGGTC')
  

#### Optional arguments:
These can be added to additionally configure pathfinding. Most require an additional argument. A comma indicates that an option can be called multiple ways.

 `-h`, `--help`            Show help message and exit

  `-i IN_TYPE`, `--in_type IN_TYPE`
                        Set input file type; default is 'counts', which
                        assumes three lines of header data followed by lines
                        of the format 'sequence count' where count is an
                        integer. Future versions will include additional options.
                        
  `-n NUM_PATHS`, `--num_paths NUM_PATHS`
                        Number of best pathways to be generated with each run (NUM_PATHS must be an integer). Default is 1.
                        
 `--max_length MAX_LENGTH`
                        Maximum length of path searched before giving up;
                        defaults to twice the length of the starting sequence. (MAX_LENGTH must be an integer)
 
 `--min_count MIN_COUNT`
                        Minimum count of sequences searched (the program
                        discards any sequences of lower count); default is 2
                        (setting lower than 2 may increase runtime
                        dramatically). (MIN_COUNT must be an integer)
  
  `--max_step MAX_STEP`   Maximum step size allowed by search paths; default is
                        1. (MAX_STEP must be an integer)
  
  `-d DIST_TYPE`, `--dist_type DIST_TYPE`
                        Distance metric used to find shortest pathway. `DIST_TYPE`
                        defaults to `edit` but `hamming` is also allowed.
  
  `-p`, `--track_progress`  If this flag is enabled, the terminal will output updates
                        on how much progress the code has made

### Output:
The output file from this script contains a list of the `NUM_PATHS` top pathways. Each is provided as a list of comma-separated values, with the data on each column (denoted by a header) corresponding to:

`step #` The number of steps from the sequence in this column to the start sequence

`sequence` The sequence of the node at this step along the pathway

`step size` The distance between this sequence along the pathway and the previous sequence

`total distance` The total distance of all steps up to this one

`sequence count` The count of this particular sequence in the counts file. May correspond to a different value, e.g. if the input file is a list of sequences and their fitness value instead of sequencing count.
