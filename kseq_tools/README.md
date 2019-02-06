
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

`-s SEARCH_SET [SEARCH_SET ...]`, `--search_set SEARCH_SET [SEARCH_SET ...]` If used, this option will cause the script to only run k-seq over a subset of all sequences. Can be used to significantly speed up use. over. By default, search set
                        is set to `all`, generating k-seq data for all sequences in
                        the start round. If set to `center CENTER_SEQUENCE DISTANCE` (requires three arguments, the second
                        being a sequence and the third being an integer), will
                        only generate kseq data over all sequences within a
                        fixed distance of the center. If set to `list
                        SEQUENCE_LIST_LOCATION`, will only generate data over
                        the sequences listed in a file at
                        the given location (this file should contain one sequence per row).

  `-v`, `--verbose`         if this flag is used, output file will include data
                        on sequence concentration at every kseq round,
                        resulting in a larger output file. Default does not add this information.


  `-i IN_TYPE`, `--in_type IN_TYPE`
                        Set input file type; default is `counts`, which
                        assumes three lines of header data followed by lines
                        of the format 'sequence count' where count is an
                        integer. Future versions will include additional options.
                        
  `-o OUT_TYPE`, `--out_type OUT_TYPE`
                        Set output file type; default is `csv` or comma-
                        separated values, a format that can be used by a
                        variety of programs; other valid options are currently
                        `tsv` or tab-separated values
                        
 `--min_count MIN_COUNT`
                        Minimum count of sequences searched (the program
                        discards any sequences of lower count). Defaults to 1,
                        keeping all sequences present in the "k-seq start"
                        round, but can be set higher. (MIN_COUNT must be an integer)
    
  `-p`, `--track_progress`  If this flag is enabled, the terminal will output updates
                        on how much progress the code has made

### Output:
The output file from this script contains a csv (or other similarly-formatted) spreadsheet of all sequences for which k-Seq data is calculated. The first row is column headers, followed by a second row describing the number of unique sequences in each round, and a third row describing the total numbers of sequences. (We found keeping this data close by to be occasionally helpful in analysis).

For each sequence's k-Seq output (that is, each row), the values in each column correspond to the following:

`Seq. Name` The sequence identity. Same as sequence provided in input counts file.

`Sequence amount` (Only if `-v` is enabled) The total amount of this sequence present in the test tube at the start of k-seq. Units correspond to normalization constants (e.g. if the normalization units are 1/ng, this will be in ng of sequence).

`Surv. fraction ROUND_NAME` (Only if `-v` is enabled) The fraction of this sequence that survives selection to be sequenced in sample `ROUND_NAME`. Should ideally be < 1 for all sequences in all rounds.
