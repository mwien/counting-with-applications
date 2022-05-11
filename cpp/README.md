# AAAI 2021 Experiments
A C++ implementation of the Clique-Picking algorithm (only counting, not sampling) is contained in this folder as `cliquepicking.cc`. It can be compiled with the following command:

```
g++ -O3 -march=native -lgmp -lgmpxx cliquepicking.cc
```

The solver can be run with `./a.out < <input-file>` and expects a chordal graph in the format defined below.

## Graph Format

Each graph on n vertices contains the vertices {1,2,...,n}. The first line of a graph file contains two space-separated integers: the number of vertices n and the number of edges m. It follows one blank line and m lines containing the edges as two space-separated integers.
