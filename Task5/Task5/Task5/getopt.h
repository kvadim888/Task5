#ifndef GETOPT_H
#define GETOPT_H

#include <stdint.h>

extern char	*optarg;
extern int	optind;

int getopt(int argc, char *const argv[], const char *optstring);

#endif
