#ifndef __MSG_HDR__
#define __MSG_HDR__
/**
 * @file
 * @brief Header file for function for typical STDERR
 * Contain function for reading and writing PDB files
 */

/** @brief Global counter of warnings raised.*/
extern int global_warn;

/** @brief Maximum number of warnings that are tolerated.*/
extern int global_max_warn;

/** @brief Aborts program and prints a fail message*/
void failed (char message[]);

/** @brief Prints checkpoint message*/
void checkpoint(char message[]);

/** @brief Aborts with a detailed Fatal ERROR message */
void error(char message[], char *path, int line);

/** @brief Prints a  standard WARNING message */
void warning(char message[], char *path, int line);

/** @brief Prings a Standard INFO message */
void info(char message[], char *path, int line);

#endif //__MSG_HDR__

