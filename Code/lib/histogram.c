#include "histogram.h"


histogram * histogram_init(double min, double max, double number_of_bins, double alpha) 
{
/* Errors and Warning Handling */
    if (min >= max)
    {
        printf("Maximal value in histogram have to be larger then minimal!!!\n");
        return NULL;
    }
    if (number_of_bins <= 0)
    {
        printf("Number of bins can not be lower then 0!!\n");
        return NULL;
    }
    printf("INFO: Initial size of histogram and free energy is 2x %g\n", sizeof(double)*number_of_bins);
/*-----------------------------*/

    histogram 
        *histo              = (histogram *)malloc(sizeof(histogram));

    int
        *frequency          = i1t(number_of_bins);                      // actual histogram 

    double
        *free_energy        = d1t(number_of_bins);                      // free energy surface

    for (int i = 0; i < number_of_bins; i++)
    {
        frequency[i]        = 0;
        free_energy[i]      = 0.0;
    }

    // main data structures
    histo->frequency        = frequency;
    histo->free_energy      = free_energy;

    // properties of histogram/free energy (dimensionality)
    histo->dimension        = 1;
    histo->number_of_bins[0]= number_of_bins;
    histo->min[0]           = min;
    histo->max[0]           = max;
    histo->bin_size_inv[0]  = number_of_bins/(max-min);

    // aditional stuff for free energy calculation
    histo->N                = 0;
    histo->max_frequency    = 0;                                        // frequncy of maximaly populated bin (used to rescale alpha)
    histo->min_frequency    = 0;                                        // frequncy of minimaly populated bin (used to rescale alpha)
    histo->alpha            = alpha;                                    // alpha that is used in wang landau sampling

    return histo;
}

void histogram_add_dimension(double min, double max, double number_of_bins, histogram *histo)
{
    // resize 1D array
    // adding dimension to histogram with values will rise warning since it arase all date ... since we cant asigne old date to any new value ... might be changed but do we need it ?
/* Errors and Warning Handling */
    if ( histo->N != 0 ) // WARRNING
    {
        printf("We are adding dimension to non-empty histogram ... do not do that ... it is something u probably do not want to do!\n");
    }
    if ( histo->dimension >= MAX_HISTOGRAM_DIMENSION ) //ERROR
    {
        printf("blabla\n");
        return;
    }
/*-----------------------------*/

    int
        *new_frequency;

    double
        number_of_all_bins = number_of_bins,
        *new_free_energy;

    // calculate current size of histogram and freeEnergy
    for(int i = 0; i < histo->dimension; i++){
        number_of_all_bins *= histo->number_of_bins[i];
    }

    //now expand 1D arrays
    new_frequency = (int *) realloc ( histo->frequency, (size_t) number_of_all_bins * sizeof(int));
    if ( new_frequency == NULL) //ERROR
    {
        printf("Memory: was not possible to realloc histogram in histogram_add_dimension(...)!\nSize of new histogram was: %lf MB\n", number_of_all_bins*sizeof(int)*2*BYTES_TO_MB);
        return;
    }
    histo->frequency   = new_frequency;

    new_free_energy = (double *) realloc ( histo->free_energy, (size_t) number_of_all_bins * sizeof(double));
    if ( new_free_energy == NULL) //ERROR
    {
        printf("Memory: was not possible to realloc histogram in histogram_add_dimension(...)!\nSize of new histogram was: %lf MB\n", number_of_all_bins*sizeof(double)*2*BYTES_TO_MB);
        return;
    }
    histo->free_energy = new_free_energy;

    for (int i = 0; i < number_of_all_bins; i++)
    {
        histo->frequency[i]     = 0;
        histo->free_energy[i]   = 0.0;
    }

    // properties of histogram/free energy (dimensionality)
    // here histo->dimension is index +1 since array starts from 0
    histo->number_of_bins[histo->dimension] = number_of_bins;
    histo->min[histo->dimension]            = min;
    histo->max[histo->dimension]            = max;
    histo->bin_size_inv[histo->dimension]   = number_of_bins/(max-min);
    histo->dimension++;
}

void histogram_free(histogram *histo)
{
    free_i1t(histo->frequency);
    free_d1t(histo->free_energy);
    free(histo);
}

/*
histogram * histogram_load_text(char *filename)
{

}

histogram * histogram_load_bin(char *filename)
{

}

void histogram_save_text(char *filename, histogram *histo)
{

}

void histogram_save_bin(char *filename, histogram *histo)
{

}
*/

void histogram_add(double *point, histogram *histo)
{
// function take value in N-dimensional space and transfer and add 1 to particular bin in histogram
    int
        index   = 0;// index in 1D array coresponding to coordinates in N-dimensional array

    double
        multi = 1.0; // multiplication prefactor coresponding to size of preceding dimensions 

    for (int i = 0; i < histo->dimension; i++)
    {
        index += (int)( ((point[i] - histo->min[i]) * histo->bin_size_inv[i])) * multi;
        multi *= histo->number_of_bins[i];
    }
    histo->frequency[index]++;
    histo->free_energy[index] -= histo->alpha;
    histo->N++;
}

void histogram_add_safe(double *point, histogram *histo)
{
// function take value in N-dimensional space and transfer and add 1 to particular bin in histogram
// check if value is in range of histogram 
    int
        index   = 0,// index in 1D array coresponding to coordinates in N-dimensional array
        intermediate;// index of value[i] in i-th dimension (used for checking if value belongs to interval histo->min[i] to histo->max[i]
    double
        multi =1.0;     // multiplication prefactor coresponding to size of preceding dimensions 

    for (int i = 0; i < histo->dimension; i++)
    {
        intermediate = (int)((point[i] - histo->min[i]) * histo->bin_size_inv[i]);
        if( intermediate >= histo->number_of_bins[i] )
        {
            printf("EPIC TERRORR we are out of histogram ... \n");
        }
        index += intermediate * multi;
        multi *= histo->number_of_bins[i];
    }
    histo->frequency[index]++;
    histo->free_energy[index] -= histo->alpha;
    histo->N++;
}

int histogram_add_soft(double *point, histogram *histo) //function return 0 if point is inside histogram boundaries ... if not it returns -1
{
// function take value in N-dimensional space and transfer and add 1 to particular bin in histogram
// check if value is in range of histogram 
    int
        index = 0,// index in 1D array coresponding to coordinates in N-dimensional array
        intermediate;// index of value[i] in i-th dimension (used for checking if value belongs to interval histo->min[i] to histo->max[i]
    double
        multi =1.0;     // multiplication prefactor coresponding to size of preceding dimensions 

    for (int i = 0; i < histo->dimension; i++)
    {
        intermediate = (int)((point[i] - histo->min[i]) * histo->bin_size_inv[i]);
        if( intermediate >= histo->number_of_bins[i] )
        {
            return -1;
        }
        index += intermediate * multi;
        multi *= histo->number_of_bins[i];
    }
    histo->frequency[index]++;
    histo->free_energy[index] -= histo->alpha;
    histo->N++;
    return 0;
}


double histogram_free_energy(double *point, histogram *histo)
{
    int
        index=0;// index in 1D array coresponding to coordinates in N-dimensional array

    double
        multi = 1.0; // multiplication prefactor coresponding to size of preceding dimensions 
    for (int i = 0; i < histo->dimension; i++)
    {
        index += (int)((point[i] - histo->min[i]) * histo->bin_size_inv[i]) * multi;
        multi *= histo->number_of_bins[i];
    }
    return histo->free_energy[index];
}

double histogram_free_energy_safe(double *point, histogram *histo)
{
    int
        index=0,// index in 1D array coresponding to coordinates in N-dimensional array
        intermediate;// index of value[i] in i-th dimension (used for checking if value belongs to interval histo->min[i] to histo->max[i]
    double
        multi =1.0;     // multiplication prefactor coresponding to size of preceding dimensions 


    for (int i = 0; i < histo->dimension; i++)
    {
        intermediate = (int)((point[i] - histo->min[i]) * histo->bin_size_inv[i]);
        if( intermediate >= histo->number_of_bins[i] )
        {
            printf("EPIC TERRORR we are out of histogram ... \n");
        }
        index += intermediate * multi;
        multi *= histo->number_of_bins[i];
    }
    return histo->free_energy[index];
}

double histogram_frequency(double *point, histogram *histo)
{
    int
        index = 0;  // index in 1D array coresponding to coordinates in N-dimensional array
    double
        multi = 1.0; // multiplication prefactor coresponding to size of preceding dimensions 
    for (int i = 0; i < histo->dimension; i++)
    {
        index += (int)((point[i] - histo->min[i]) * histo->bin_size_inv[i]) * multi;
        multi *= histo->number_of_bins[i];
    }
    return histo->frequency[index];
}

double histogram_frequency_safe(double *point, histogram *histo)
{
    int
        index=0,// index in 1D array coresponding to coordinates in N-dimensional array
        intermediate;   // index of value[i] in i-th dimension (used for checking if value belongs to interval histo->min[i] to histo->max[i]

    double
        multi =1.0;     // multiplication prefactor coresponding to size of preceding dimensions 

    for (int i = 0; i < histo->dimension; i++)
    {
        intermediate = (int)((point[i] - histo->min[i]) * histo->bin_size_inv[i]);
        if( intermediate >= histo->number_of_bins[i] )
        {
            printf("EPIC TERRORR we are out of histogram ... \n");
        }
        index += intermediate * multi;
        multi *= histo->number_of_bins[i];
    }
    return histo->frequency[index];
}


void histogram_print(char *filename, char *format, histogram *histo)//would be nice to specifie format of data ... that can be easyli printed via gnuplot etc.
{
    char
        error[1024];

    FILE
        *histogram_out;

    if ((histogram_out=fopen(filename,"w"))==NULL)
    {
        sprintf(error,"Could not open file %s (%s).\n",filename,"w");
        failed(error);
    }
    // printing is done in format for gnuplot matrix ... so it is just reprint of matrix in memory
    if (strcmp(format, "matrix2d") == 0)
    {
        for(int i=0; i < histo->number_of_bins[0]; i++)
        {
            for(int ii=0; ii< histo->number_of_bins[1]; ii++)
            {
                fprintf(histogram_out, "%+08i\t", histo->frequency[(i*histo->number_of_bins[0])+ii]);
            }
            fprintf(histogram_out, "\n");
        }
    }else{
        printf("No other methods implemented sofar\n");
    }

    fclose(histogram_out);
}


