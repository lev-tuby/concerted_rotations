# Concerted rotations of protein backbones using Wolfram Mathematica

This repository contains the source code to generate concerted rotations of a protein backbone using the approach introduced by [Zamuner et al](http://journals.plos.org/plosone/article?id=10.1371/journal.pone.0118342).
In our implementation, the C code is precomputed using Wolfram Mathematica, making it significantly fast. While the library in this repository implements a move of 7 consecutive dihedrals on a protein backbone, the 
user can modify the Mathematica notebook and Python wrapper to use a different set of free variables, or the C source code to reproduce a different backbone. The geometry of the backbone is known to the concerted rotation 
algorithm only through the function which convert the backbone into a series of Denavit-Hartenberg parameters.

The general scheme of our library is reported in the following diagram:

<div style="text-align:center">
  <img src="Images/from_mathematica_to_C.png" width="600">
</div>

## Generating the code

To generate the precompiled code for a generic concerted rotation, open and run the Mathematica notebook 
in folder `Code/generic_move`.

To compile the test program, in folder `Code` run:

`make concerted_rot.x`

Our implementation of the concerted rotations is based on the [Gnu Scientific Library](https://www.gnu.org/software/gsl/), which needs to be installed before compiling the code.

To run the test program:

`./concerted_rot.x aapot`

The test code simulates a phantom protein backbone. Nonetheless, the code can be easily modified to 
sample according to the [Caterpillar model](http://journals.plos.org/plosone/article?id=10.1371/journal.pone.0112852), an implementation of which is contained in directory `lib`. The file 
aapot contains the Caterpillar interaction matrix, as published in [Coluzza, PLoS 2014](http://journals.plos.org/plosone/article?id=10.1371/journal.pone.0112852).

In the main function of `concerted_rot.c` the value of the following parameters can be easily modified 
to test their influence on the concerted rotation sampling.

* Number of iterations: `maxIter`
* Variance of normal distribution used to obtain step sizes in tangent space: `sigma`

Sequence and sequence length are hardcoded in the code at beginning of `concerted_rot.c`.

## Modifying the concerted rotation

* Instructions on how to change the 7 free variables are included in the Mathematica notebooke in `Code/generic_move`. 
* The library produced by Mathematica+Python works using the DH hartenberg convention. In order to simulate a different backbone one has to write a function which constructs the DH bases starting from that. Check `Code/minimal.c` for a minimal example, and `Code/lib/CAT_moves.c` for the functions mapping a protein backbone into a set of DH parameters.


