# VertexEliminationOrders

A package for computing vertex elimination orders with minimal treewidth for graphs.

The main function for computing elimination orders uses an anytime, branch-and-bound, 
depth-first search for computing the treewidth of a graph. If the search doesn't 
complete within the allocated time, an upper bound is returned representing the best 
solution the algorithm was able to find.

The basic algorithm is outlined by Gogate in the 2004 paper "A complete anytime 
algorithm for treewidth".

One of the heuristics mentioned in Yaun's 2011 paper "A Fast Parallel Branch and 
Bound Algorithm for Treewidth" is also used. Namely, the hash table for detecting
duplicates.
    
A discussion of relevant heuristics can be found in Bodlaender's "Treewidth 
Computations 1 & 2" papers but not yet used.

# Basic usage

```
using VertexEliminationOrders

# Load an example graph.
g = load_example_graph("sycamore_53_20_0")

# run the algorithm for 5 seconds.
treewidth, order, report = estimate_treewidth(g, 5.0);

# show the returned report containing the result and some metadata.
@show report
```
