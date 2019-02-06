
## Tool for finding pathways between peak sequences in a fitness landscape

We expect future versions of this script to be updated/included with additional tools. As it is used in a publication, however, v0.1 of the peak_pather script will remain here for posterity.

### Goal:

The peak_pather script searches for the shortest pathway between two sequences, along a fitness landscape, using the A* algorithm (a decent description can be found on [Wikipedia](https://en.wikipedia.org/wiki/A*_search_algorithm). The script iterates over a single round's sequence population as a graph, with each sequence present as its own node and the edit distances between sequences as edges with distance. The script provided here finds the N best pathways between two sequences, ranked by 1) Shortest total path length, 2) Smallest maximum step size, 3) Smallest average step size, 4) Largest minimum sequence count, using a sequence counts file as the reference map (this could easily be replaced with highest minimum fitness, if the reference file is a list of sequences and their fitnesses instead).
