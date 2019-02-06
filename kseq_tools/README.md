
## Tool for calculating k-Seq data over a set of rounds whose corresponding DNA has been sequenced.

We expect future versions of this script to be updated/included with additional tools and a slightly improved algorithm. As it is used in a publication, however, v0.1 of the peak_pather script will remain here for posterity.

### Goal:

The kseq_tools calculates catalytic kinetics for a population of sequences, using the k-seq
methodology. Takes a k-seq 'start round' and a list of additional rounds (each corresponding to selection under known conditions);
gives predicted constants A and k*t for catalysis following [surviving
fraction] = A(1-Exp(-k*[S]*t)). If generally confused over use, see example
use case given on SCAPE-BYO readme [here](https://github.com/ichen-lab-ucsb/SCAPE-BYO/blob/master/README.md).


### Input:

The script can be called as follows:

```
python kseq_tools_v01 start_round kseq_rounds output normalization_list substrate_concs rounds_to_average rounds_to_error
```
#### Required arguments (positionally dependent):
`start_round`                 Location/name of file containing sequence counts for the pre-k-seq population (e.g. `R5c-counts.txt`). Currently requires a "counts" file consisting of three lines of metadata followed by  one line per unique sequence in pool, of the format "sequence count" where count is an integer. Such files are produced by our Galaxy tools, currently available at https://labs.chem.ucsb.edu/chen/irene/Chen_lab_at_UCSB/Galaxy_Tools.html. 

`kseq_rounds`                Location/name of file containing a list of additional
                        filenames, each of which contains sequence counts for
                        a post-k-seq population. (It is recommended to keep all counts files in the same location). Each of these files must share the same format as `start_round`.

`output`                Output file location/name (e.g. `kseq-data.csv`)

`normalization_list`             List of normalization factors for each round, starting
                        with the start round. One value per line. For the
                        start round, should typically be set equal to 1/[total
                        amount of DNA/RNA/protein] present at start of k-seq
                        selection rounds, and for all other rounds should be
                        1/[amount of DNA/RNA/protein] remaining after
                        selection step. The units here will correspond to the units used in the script's output.

`substrate_concs`             File containing list of substrate concentrations for
                        each set of kseq rounds. One row corresponds to each unique substrate concentration, not to each round/experimental sample (if duplicates are used). Must be the same number of rows as `rounds_to_averge`


`rounds_to_average`               File containing comma-separated lists of of sets of
                        rounds to average together for fits (e.g. row 1 as "1,2,3" and row 2 as "4,5,6"
                        will average the abundances of rounds 1,2, and 3, and
                        then also average 4,5, and 6). If only a single
                        replicate was carried out, then the file should only contain one number per round (e.g. row 1 as "1," row 2 as "2", etc.) If you just   want abundance without kseq calculations, leave blank. Must be the same number of rows as `substrate_concs`
                        See example use case (linked to at the top of README) if confused.

  

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
