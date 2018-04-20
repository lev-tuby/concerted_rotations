/**
 * @file
 * @brief Source file for all functions for work with histogram structure
 * Well not much else to say ...
 * @todo Add methods to save and load histogram in binary and text format
 */
#include "histogram.h"


/**
 * @brief Initialization of 1D histogram
 *
 * Function return histogram structure of 1D histogram.
 * To add more dimensions to existing histogram see histogram_add_dimension().
 *
 * @param[in]  min              Minimal value of histogram in 1D histogram.
 * @param[in]  max              Maximal value of histogram in 1D histogram.
 * @param[in]  number_of_bins   Number of bins in 1D histogram.
 * @param[in]  alpha            Initial value of microcanonical entropy change \f$\alpha\f$, used in Wang and Landau algorithm.
 *
 * @return histogram            Initialized 1D histogram structure.
 */
histogram * histogram_init(double min, double max, double number_of_bins, double alpha) 
{
/* Errors and Warning Handling */
    char
        msg[1024];

    if (min >= max)
    {
        snprintf(msg, sizeof(msg), "Maximal value in histogram have to be larger then minimal!\nminimal: %g | maximal %g", min, max);
        error(msg, __FILE__, __LINE__);
    }
    if (number_of_bins <= 0)
    {
        error("Number of bins can not be lower then 0!\n", __FILE__, __LINE__);
    }
    snprintf(msg, sizeof(msg), "Initial size of histogram and free energy is 2x %g.\n", sizeof(double)*number_of_bins);
    info(msg, __FILE__, __LINE__);
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

/**
 * @brief Function add one dimension to existing histogram
 *
 * Function modifie histogram structure by adding one aditional dimension.
 * <b>Function set all values in histogram and free energy to 0.0!</b>
 *
 * @param[in]       min              Minimal value of histogram in additional dimension.
 * @param[in]       max              Maximal value of histogram in additional dimension.
 * @param[in]       number_of_bins   Number of bins in additional dimension.
 * @param[in,out]   *histo           Histogram strucure to whcih dimension is added.
 *
 * @return \c void
 */
void histogram_add_dimension(double min, double max, double number_of_bins, histogram *histo)
{
/* Errors and Warning Handling */
    char
        msg[1024];

    if (min >= max)
    {
        snprintf(msg, sizeof(msg), "Maximal value in histogram have to be larger then minimal!\nminimal: %g | maximal %g", min, max);
        error(msg, __FILE__, __LINE__);
    }
    if (number_of_bins <= 0)
    {
        error("Number of bins can not be lower then 0!\n", __FILE__, __LINE__);
    }
    if ( histo->N != 0 )
    {
        warning("We are adding dimension to non-empty histogram ... do not do that ... it is something you probably do not want to do!", __FILE__, __LINE__);
    }
    if ( histo->dimension >= MAX_HISTOGRAM_DIMENSION )
    {
        snprintf(msg, sizeof(msg), "You reach maximal histogram dimension (%d)!\nMaximal histogram dimension could be changed in histogram.h", MAX_HISTOGRAM_DIMENSION);
        error(msg, __FILE__, __LINE__);
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
    if ( new_frequency == NULL)
    {
        snprintf(msg, sizeof(msg), "Memory: was not possible to realloc histogram in histogram_add_dimension(...)!\nSize of new histogram was: %lf MB\n", number_of_all_bins*sizeof(int)*2*0.001);
        error(msg, __FILE__, __LINE__);
    }
    histo->frequency   = new_frequency;

    new_free_energy = (double *) realloc ( histo->free_energy, (size_t) number_of_all_bins * sizeof(double));
    if ( new_free_energy == NULL)
    {
        snprintf(msg, sizeof(msg), "Memory: was not possible to realloc histogram in histogram_add_dimension(...)!\nSize of new histogram was: %lf MB\n", number_of_all_bins*sizeof(int)*2*0.001);
        error(msg, __FILE__, __LINE__);
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

/**
 * @brief Free memory ocupied by histogram structure
 *
 * @param[in,out]   *histo           Pointer to histogram structure to be deallocated.
 *
 * @return \c void
 */
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

/**
 * @brief Function add sample to histogram
 *
 * Function add 1 to bin coresponding to coordinates in point in N-dimensional histogram.
 * Function also update <c>free_energy</c> by adding alpha to particular bin.
 * Function also increase counter of samples in histogram <c>N</c>.
 * <b>Point should have same dimesnion as histogram! Function do not check if point has same dimesnion as histogram!</b>
 * <b>Function do not check if point belongs to histogram (migh cause acces to wrong memory)!</b>
 *
 * @param[in]       *point      N-dimensional value to be added in particular bin in histogram.
 * @param[in,out]   *histo      Histogram strucure to whcih sample is added.
 *
 * @return \c void
 */
void histogram_add(double *point, histogram *histo)
{
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

/**
 * @brief Function add sample to histogram
 *
 * Function add 1 to bin coresponding to coordinates in point in N-dimensional histogram if point belong to histogram.
 * Function also update <c>free_energy</c> by adding alpha to particular bin if point belong to histogram.
 * Function also increase counter of samples in histogram <c>N</c> if point belong to histogram.
 * <b>Function do check if point belongs to histogram.</b> If point do not belong to histogram then rise ERROR.
 *
 * @param[in]       *point      N-dimensional value to be added in particular bin in histogram.
 * @param[in,out]   *histo      Histogram strucure to whcih sample is added.
 *
 * @return \c void
 */
void histogram_add_safe(double *point, histogram *histo)
{
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
            error("Sample can not be added outside of histogram!", __FILE__, __LINE__);
        }
        index += intermediate * multi;
        multi *= histo->number_of_bins[i];
    }
    histo->frequency[index]++;
    histo->free_energy[index] -= histo->alpha;
    histo->N++;
}

/**
 * @brief Function add sample to histogram
 *
 * Function add 1 to bin coresponding to coordinates in point in N-dimensional histogram if point belong to histogram.
 * Function also update <c>free_energy</c> by adding alpha to particular bin if point belong to histogram.
 * Function also increase counter of samples in histogram <c>N</c> if point belong to histogram.
 * <b>Function do check if point belongs to histogram.</b>
 * If point do not belong to histogram then function return <c>OUT_HISTOGRAM</c> if point belong to histogram then return <c>IN_HISTOGRAM</c>.
 *
 * @param[in]       *point      N-dimensional value to be added in particular bin in histogram.
 * @param[in,out]   *histo      Histogram strucure to whcih sample is added.
 *
 * @return <c>OUT_HISTOGRAM</c> or <c>IN_HISTOGRAM</c>
 */
int histogram_add_soft(double *point, histogram *histo)
{
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
            return OUT_HISTOGRAM;
        }
        index += intermediate * multi;
        multi *= histo->number_of_bins[i];
    }
    histo->frequency[index]++;
    histo->free_energy[index] -= histo->alpha;
    histo->N++;
    return IN_HISTOGRAM;
}

/**
 * @brief Function return free energy for given bin
 *
 * Function read free energy at given point of histogram.
 * <b>Function do not check if point belongs to histogram.</b> So that nonsense value can be returned if point is outside of histogram.
 *
 * @param[in]       *point      N-dimensional value at which free energy is returned.
 * @param[in]       *histo      Histogram strucure from where free energy is obtained.
 *
 * @return Free energy at <c>point</c>
 */
double histogram_free_energy(double *point, const histogram *histo)
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

/**
 * @brief Function return free energy for given bin
 *
 * Function read free energy at given point of histogram.
 * <b>Function do check if point belongs to histogram.</b>
 * If point do not belong to histogram then function return <c>OUT_HISTOGRAM</c> if point belong to histogram then return <c>IN_HISTOGRAM</c>.
 *
 * @param[in]       *point      N-dimensional value at which free energy is returned.
 * @param[out]      *value      Free energy at <c>point</c>.
 * @param[in]       *histo      Histogram strucure from where free energy is obtained.
 *
 * @return <c>OUT_HISTOGRAM</c> or <c>IN_HISTOGRAM</c>
 */
int histogram_free_energy_safe(double *point, double *value, const histogram *histo)
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
            error("Sample can not be read from outside of histogram!", __FILE__, __LINE__);
        }
        index += intermediate * multi;
        multi *= histo->number_of_bins[i];
    }
    return histo->free_energy[index];
}

/**
 * @brief Function return frequency for given bin
 *
 * Function read frequency at given point of histogram.
 * <b>Function do not check if point belongs to histogram!</b> So that nonsense value can be returned if point is outside of histogram.
 *
 * @param[in]       *point      N-dimensional value at which free energy is returned.
 * @param[in]       *histo      Histogram strucure from where free energy is obtained.
 *
 * @return Frequency at <c>point</c>
 */
double histogram_frequency(double *point, const histogram *histo)
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

/**
 * @brief Function return frequency for given bin
 *
 * Function read frequency at given point of histogram.
 * <b>Function do check if point belongs to histogram.</b>
 * If point do not belong to histogram then function return <c>OUT_HISTOGRAM</c> if point belong to histogram then return <c>IN_HISTOGRAM</c>.
 *
 * @param[in]       *point      N-dimensional value at which free energy is returned.
 * @param[out]      *value      Frequency at <c>point</c>.
 * @param[in]       *histo      Histogram strucure from where free energy is obtained.
 *
 * @return <c>OUT_HISTOGRAM</c> or <c>IN_HISTOGRAM</c>
 */
double histogram_frequency_safe(double *point, double *value, const histogram *histo)
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
            error("Sample can not be read form outside of histogram!", __FILE__, __LINE__);
        }
        index += intermediate * multi;
        multi *= histo->number_of_bins[i];
    }
    return histo->frequency[index];
}

/**
 * @brief Function print histogram in selected format to file
 *
 * Function so far support only ploting to GNUPLOT matrix format <c>matrix2d</c>/
 *
 * @param[in]       *filename       Name of file used for printing histogram.
 * @param[in]       *format         Constant string that specifies format in which output data are printed.
 * @param[in]       *histo          Histogram strucure from where data are obtained.
 *
 * @return \c void
 */
void histogram_print(const char *filename, const char *format, const histogram *histo)//would be nice to specifie format of data ... that can be easyli printed via gnuplot etc.
{
    char
        msg[1024];

    FILE
        *histogram_out;

    if ((histogram_out=fopen(filename,"w"))==NULL)
    {
        snprintf(msg, sizeof(msg), "Could not open file %s (%s).\n",filename,"w");
        error(msg, __FILE__, __LINE__);
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
        warning("No other methods implemented so far!", __FILE__, __LINE__);
    }

    fclose(histogram_out);
}




