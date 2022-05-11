# Code for the paper "Polynomial-Time Algorithms for Counting and Sampling Markov Equivalent DAGs with Applications"
This repository contains all code used for the experiments in the paper.

## Usage
For replicating the experiments, the script "experiments.jl" has to be run:
```julia experiments.jl```. In order for this to work properly, you need to install the dct-policy package (https://github.com/csquires/dct-policy), add the files, which are under "/scripts", in the directory you have dct-policy, make the bash scripts executable and set appropriate paths in line 8 and 10 of experiments.jl.

Due to large file sizes this repository includes only one graph per choice of parameters for the sampling and active learning experiments. If you are interested in all the graphs or if you have any questions, do not hesitate to contact us.

## Functionality
In case you want to use the provided functions (in particular counting and sampling) apart from replicating these experiments, we encourage you to take a look at the repository CliquePicking: https://github.com/mwien/CliquePicking.
It includes some more explanation on using the code directly and will be kept up to date.