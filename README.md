# Inference of Hierarchical Core—Periphery Structure in Temporal Networks

This code implements our Markov chain Monte Carlo approach for the identification of hierarchical core—periphery structure in temporal networks, as described in this paper.

# Acknowledgements

We thank [Christopher R. Anderson](http://www.math.ucla.edu/~anderson) for writing much of the MathVector and MathMatrix classes, as part of his Math 280 course.

# Usage

To run the code, one must first modify `main.cpp` by specifying:

`numGroups`: the number of groups in the initial group assignment

`numNodes`: the number of nodes in the network

`layersA`: the number of layers in the network

`dataset`: the start of the filenames for the imported network structure and for the outputted group structure

`mnMoveProb`: the probability that a multi-node move is proposed

`numMCSteps`: the number of Monte Carlo steps to perform

Additionally, for each layer of a network, one needs to add text files with filename `[dataset]layer[layer].txt` (where [dataset] was chosen previously and [layer] is the index of the given layer; for example, if `dataset ` was chosen to be `luke` and we are considering layer 3, the filename should be `lukelayer3.txt`) each of which contains the adjacency matrix for a given layer of the network. These text files should be such that each entry in a row of the adjacency matrix is seperated by a single space and such that each row is seperated by a newline. For example:
```
0 1 0
1 0 0
0 0 0
```

The code saves and outputs the group assignments for each node-layer every 10^4 steps. The group assignment of a given node-layer at a given step is represented as a number, namely \sum_{r=1}^{k-1} 2^{r-1} g^r_{(i,l)}. For example, if the number k of groups is 3, g^1_{(i,l)} = 0, and g^2_{(i,l)} = 1, we would denote the group assignment of node-layer (i,l) by 2.
The saved group assignments for a given layer are outputted to the files ```[dataset]layer[layer]groupAssigns.txt``` (where [dataset] was chosen previously and [layer] is the index of the given layer; for example, if `dataset ` was chosen to be `luke` and we are considering layer 3, the filename will be `lukelayer3groupAssigns.txt`).
