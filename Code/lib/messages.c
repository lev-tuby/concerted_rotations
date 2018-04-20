#include <stdlib.h>
#include <stdio.h>
#include "messages.h"

/** @brief Global counter of warnings raised.*/
int global_warn;

/** @brief Maximal number of warnings that are tolerated.*/
int global_max_warn;

/**
 * @brief Fail message
 *
 * Wrapper around program termination
 *
 * @param[in]  message         Message describing reason for program failure.
 *
 * @return \c void
 */
void failed (char message[])
{
    fprintf (stderr,"Error. \n FAILED *** %s ***\n \n", message);
    abort();
}

/**
 * @brief Checkpoint message
 *
 * Wrapper around checkpoint output
 *
 * @param[in]  message         Message describing state of program.
 *
 * @return \c void
 */
void checkpoint(char message[])
{
    printf ("\n CHECKPOINT *** %s ***\n \n", message);
}

/**
 * @brief Standard ERROR message
 *
 * Call of function results in program termination
 *
 * @param[in]  msg              Text of error message printed in comand line.
 * @param[in]  *path            Path to file where error was raised.
 * @param[in]  line             Line number from where error was rised.
 *
 * @return \c void
 */
void error(char msg[], char *path, int line)
{
    fprintf(stderr,"ERROR:\t%s\nIn file:%s line:%d\n", msg, path, line);
    fflush(stderr);
    abort();
}

/**
 * @brief Standard WARNING message
 *
 * Call of function results in rising warrning message if number of already rised warning is smaller then maximal number of allowed warnings.
 * If number of warnings reach maximal number <c>global_max_warn</c> program is terminated and place from where warning was generated is shown.
 *
 * @param[in]  msg              Text of error message printed in comand line.
 * @param[in]  *path            Path to file where error was raised.
 * @param[in]  line             Line number from where error was rised.
 *
 * @return \c void
 */
void warning(char msg[], char *path, int line)
{
    if (global_warn < global_max_warn)
    {
        fprintf(stderr, "WARNING:\t%s\n", msg);
        fflush(stderr);
        global_warn++;
    }
    else
    {
        fprintf(stderr, "WARNING:\t%s\nIn file:%s line:%d\n", msg, path, line);
        fflush(stderr);
        abort();
    }
}

/**
 * @brief Standard INFO message
 *
 * Call of this function print out INFO message ... should contain inportant informations or hints for simulation.
 *
 * @param[in]  msg              Text of error message printed in comand line.
 * @param[in]  *path            Path to file where error was raised.
 * @param[in]  line             Line number from where error was rised.
 *
 * @return \c void
 */
void info(char msg[], char *path, int line)
{
    fprintf(stdout, "INFO:\t%s\n", msg);
    fflush(stdout);
}

