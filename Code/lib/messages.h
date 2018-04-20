#ifndef __MSG_HDR__
#define __MSG_HDR__
/**
 * @file
 * @brief Header file for function for typical STDOUT
 * Contain function for reading and writing PDB files
 */

/** @brief Global counter of warnings raised.*/
extern int global_warn;

/** @brief Maximal number of warnings that are tolerated.*/
extern int global_max_warn;

/** @brief Fail message*/
void failed (char message[]);

/** @brief Checkpoint message*/
void checkpoint(char message[]);

/**
 * @brief Standard ERROR message
 *
 * Call of function results in program termination
 */
void error(char message[], char *path, int line);

/**
 * @brief Standard WARNING message
 *
 * Call of function results in rising warrning message if number of already rised warning is smaller then maximal number of allowed warnings.
 * If number of warnings reach maximal number <c>global_max_warn</c> program is terminated and place from where warning was generated is shown.
 */
void warning(char message[], char *path, int line);

/**
 * @brief Standard INFO message
 *
 * Call of this function print out INFO message ... should contain inportant informations or hints for simulation.
 */
void info(char message[], char *path, int line);

#endif //__MSG_HDR__

