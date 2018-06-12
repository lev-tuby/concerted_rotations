#ifndef __HISTOGRAM__
#define __HISTOGRAM__
#include <stdio.h>
/**
 * @file
 * @brief Header file for all functions for work with histogram structure
 * Well not much else to say ...
 */


/** @brief Maximal value of histogram dimensions used in initialization of histogram structure. */
#define MAX_HISTOGRAM_DIMENSION 8

/** @brief Return values of histogram add/read if value end up inside histogram range.*/
#define IN_HISTOGRAM  0
/** @brief Return values of histogram add/read if value end up outside histogram range.*/
#define OUT_HISTOGRAM 1

/**
 * @brief Main structure to store data of histograms.
 *
 * Except for normal histogram also microcanonical entropy based on WL algorithm is stored for each bin.
 */
typedef struct {
    // main data structures
    int *frequency;                                         /**< Histogram data stored in 1D array. */
    double *free_energy;                                    /**< Free energy data stored in 1D array. */

    // properties of histogram/free energy (dimensionality)
    int     dimension;                                      /**< Determine dimension of histogram. */

    // these arreys are defined explicitly to stay close in memory ... and assuming maximal dimension to be lets say 8 is pretty generous ...
    int     number_of_bins   [ MAX_HISTOGRAM_DIMENSION ];   /**< Number of bins in all dimension. */
    double  min              [ MAX_HISTOGRAM_DIMENSION ];   /**< Minimal value in all dimension. */
    double  max              [ MAX_HISTOGRAM_DIMENSION ];   /**< Maximal value in all direction. */
    double  bin_size_inv     [ MAX_HISTOGRAM_DIMENSION ];   /**< Invers of bin size in all dimension. */

    // aditional stuff for free energy calculation
    double N;                                               /**< Number of samples taken. Used resclae of alpha in WL free energy calculation. */
    double max_frequency;                                   /**< Maximal hight of histogram. Used resclae of alpha in WL free energy calculation. */
    double min_frequency;                                   /**< Minimal height of histogram. Used resclae of alpha in WL free energy calculation. */
    double alpha;                                           /**< Current value of microcanonical entropy change. Used for WL free energy calculation. */
}histogram;


/**
 * @brief Initialization of 1D histogram
 *
 * Function return histogram structure of 1D histogram.
 * To add more dimensions to existing histogram see histogram_add_dimension().
 */
histogram * histogram_init  (double min, double max, double number_of_bins, double alpha    );

/**
 * @brief Function add one dimension to existing histogram
 *
 * Function modifie histogram structure by adding one aditional dimension.
 * <b>Function set all values in histogram and free energy to 0.0!</b>
 */
void histogram_add_dimension(double min, double max, double number_of_bins, histogram *histo);

/**
 * @brief Free memory ocupied by histogram structure
 */
void histogram_free(histogram *histo);

/**
 * @brief Function add sample to histogram
 *
 * Function add 1 to bin coresponding to coordinates in point in N-dimensional histogram.
 * Function also update <c>free_energy</c> by adding alpha to particular bin.
 * Function also increase counter of samples in histogram <c>N</c>.
 * <b>Point should have same dimesnion as histogram! Function do not check if point has same dimesnion as histogram!</b>
 * <b>Function do not check if point belongs to histogram (migh cause acces to wrong memory)!</b>
 */
void histogram_add     (double *point, histogram *histo);
/**
 * @brief Function add sample to histogram
 *
 * Function add 1 to bin coresponding to coordinates in point in N-dimensional histogram if point belong to histogram.
 * Function also update <c>free_energy</c> by adding alpha to particular bin if point belong to histogram.
 * Function also increase counter of samples in histogram <c>N</c> if point belong to histogram.
 * <b>Function do check if point belongs to histogram.</b> If point do not belong to histogram then rise ERROR.
 */
void histogram_add_safe(double *point, histogram *histo);

/**
 * @brief Function add sample to histogram
 *
 * Function add 1 to bin coresponding to coordinates in point in N-dimensional histogram if point belong to histogram.
 * Function also update <c>free_energy</c> by adding alpha to particular bin if point belong to histogram.
 * Function also increase counter of samples in histogram <c>N</c> if point belong to histogram.
 * <b>Function do check if point belongs to histogram.</b>
 * If point do not belong to histogram then function return <c>OUT_HISTOGRAM</c> if point belong to histogram then return <c>IN_HISTOGRAM</c>.
 */
int histogram_add_soft(double *point, histogram *histo);

/**
 * @brief Function return free energy for given bin
 *
 * Function read free energy at given point of histogram.
 * <b>Function do not check if point belongs to histogram.</b> So that nonsense value can be returned if point is outside of histogram.
 */
double histogram_free_energy(double *point, const histogram *histo);

/**
 * @brief Function return free energy for given bin
 *
 * Function read free energy at given point of histogram.
 * <b>Function do check if point belongs to histogram.</b>
 * If point do not belong to histogram then function return <c>OUT_HISTOGRAM</c> if point belong to histogram then return <c>IN_HISTOGRAM</c>.
 */
int histogram_free_energy_safe(double *point, double *value, const histogram *histo);

/**
 * @brief Function return frequency for given bin
 *
 * Function read frequency at given point of histogram.
 * <b>Function do not check if point belongs to histogram!</b> So that nonsense value can be returned if point is outside of histogram.
 */
double histogram_frequency(double *point, const histogram *histo);

/**
 * @brief Function return frequency for given bin
 *
 * Function read frequency at given point of histogram.
 * <b>Function do check if point belongs to histogram.</b>
 * If point do not belong to histogram then function return <c>OUT_HISTOGRAM</c> if point belong to histogram then return <c>IN_HISTOGRAM</c>.
 */
double histogram_frequency_safe(double *point, double *value, const histogram *histo);


/**
 * @brief Function print histogram in selected format to file
 *
 * Function so far support only ploting to GNUPLOT matrix format <c>matrix2d</c>/
 */
void histogram_print(const char *filename, const char *format, const histogram *histo);


//void histogram_rescale(const histogram *histo);

#endif //__HISTOGRAM__
