#include <stdlib.h>
#include <stdio.h>
#include "messages.h"

/** @brief Global counter of warnings raised.*/
int global_warn;

/** @brief Maximum number of warnings that are tolerated.*/
int global_max_warn;

/**
 *
 * Wrapper around program termination
 *
 * @param[in]  message         Message describing the reason for program failure.
 *
 * @return \c void
 */
void failed (char message[])
{
    fprintf (stderr,"Error. \n FAILED *** %s ***\n \n", message);
    abort();
}

/**
 *
 * Wrapper around checkpoint output
 *
 * @param[in]  message         Message describing program state.
 *
 * @return \c void
 */
void checkpoint(char message[])
{
    fprintf (stderr,"\n CHECKPOINT *** %s ***\n \n", message);
}

/**
 *
 * Terminates the program
 *
 * @param[in]  msg              Error message to be printed on stderr.
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
 *
 * Rises a warning message if the number of already rised warning is smaller then maximal number of allowed warnings.
 * Otherwise, if the number of warnings has reached the maximum <c>global_max_warn</c>, the program is terminated and the origin of the last warning is shown.
 *
 * @param[in]  msg              Error message to be printed on stderr.
 * @param[in]  *path            Path to the file where the error was raised.
 * @param[in]  line             Line number from where the error was rised.
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
 *
 * Prints out an INFO message ... should contain inportant informations or hints for simulation.
 *
 * @param[in]  msg              Error message to be  printed on stderr.
 * @param[in]  *path            Path to the file where the error was raised.
 * @param[in]  line             Line number from where the error was raised.
 *
 * @return \c void
 */
void info(char msg[], char *path, int line)
{
    fprintf(stderr, "INFO:\t%s\n", msg);
    fflush(stderr);
}

