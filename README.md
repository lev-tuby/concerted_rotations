# Concerted rotations of protein backbones using Wolfram Mathematica

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


