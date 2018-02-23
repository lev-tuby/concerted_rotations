#ifndef __MSG_HDR__
#define __MSG_HDR__

extern int global_warn;
extern int global_max_warn;

void failed (char message[]);
void checkpoint(char message[]);

void error(char message[], char *path, int line);
void warning(char message[], char *path, int line);
void info(char message[], char *path, int line);

#endif //__MSG_HDR

