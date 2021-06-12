# How to create input files

The code uses an input file containing information about the swarm at each iteration. This input file could generated via a simulation of any swarm technique. In `./data/`, there are examples of four simulations of an implementation of the PSO algorihtm. 

The structure of this file is:

```
TAG ITERATION INFORMATION
```

For example, let's say that the tag is `it:#`, and the information represents the fitness, an input file should contain lines:

```
it:#1 102910.50
it:#2 92310.30
it:#3 82920.10
it:#4 72510.20
[...]
```

Note that you could have different types of information in the file: 

```
it:#1 102910.50
di:#1 60.3
it:#2 92310.30
di:#2 20.1
it:#3 82920.10
di:#3 50.1
it:#4 72510.20
di:#4 59.2
[...]
```

In the case of the interactions in the swarm, the line should contain a square matrix representing the interactions. 

For example, let's say that the tag is `ig:#`, the file of a swarm of size 5 should contain the following lines: 

```
ig:#1 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 1 1 1 0 1 0 0 0 1 0 
ig:#2 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 1 1 1 1 1 0 
ig:#3 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 1 1 1 1 1 0 
```

Note that each line contains 25 values, which enables the code to reconstruct a 5x5 matrix. This matrix represents the interactions in the swarm.

Each algorithm has its own definition of an interaction. More info here: 

- Uncovering the social interaction network in swarm intelligence algorithms. *Appl Netw Sci* 5, **24** (2020). https://doi.org/10.1007/s41109-020-00260-8

