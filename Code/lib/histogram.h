#ifndef __HISTOGRAM__
#define __HISTOGRAM__

#include <stddef.h>
#include <stdlib.h>
#include <string.h>
#include <stdio.h>
#include "my_memory.h"
#include "messages.h"

#define MAX_HISTOGRAM_DIMENSION 8
#define BYTES_TO_MB 0.000000954

//TODO: manage error codes ... for "methods"

typedef struct {
    // main data structures
    int *frequency;
    double *free_energy;

    // properties of histogram/free energy (dimensionality)
    int     dimension; // Determine dimension of histogram
    // these arreys are defined explicitly to stay close in memory ... and assuming maximal dimension to be lets say 8 is pretty generous ...
    int     number_of_bins   [ MAX_HISTOGRAM_DIMENSION ]; // number of bins in all dimension
    double  min              [ MAX_HISTOGRAM_DIMENSION ]; // minimal value in all direction
    double  max              [ MAX_HISTOGRAM_DIMENSION ]; // maximal value in all direction
    double  bin_size_inv     [ MAX_HISTOGRAM_DIMENSION ]; // invers of bin size in all dimension

    // aditional stuff for free energy calculation
    double N; //number of samples taken
    double max_frequency; // maximal hight of histogram
    double min_frequency; // minimal height of histogram
    double alpha; // current alpha
}histogram;


histogram * histogram_init  (double min, double max, double number_of_bins, double alpha    );
void histogram_add_dimension(double min, double max, double number_of_bins, histogram *histo);
void histogram_free(histogram *histo);

void histogram_add     (double *point, histogram *histo);
void histogram_add_safe(double *point, histogram *histo);

// function return value of free energy at given point in histogram
double histogram_free_energy(double *point, histogram *histo);
double histogram_free_energy_safe(double *point, histogram *histo);
// function return value of frequency at given point in histogram
double histogram_frequency(double *point, histogram *histo);
double histogram_frequency_safe(double *point, histogram *histo);



void histogram_print(char *filename, char *format, histogram *histo);
void histogram_rescale(histogram *histo); // function decide if it is suitable to 

#endif
