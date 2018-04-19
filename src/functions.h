#ifndef _FUNCTIONS_
#define _FUNCTIONS_

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>
#include <math.h>
#include <limits.h>

#include "structs.c"

#define PATH_LENGTH 32
#define LINE_LENGTH 26000
#define DECIMAL_ACCURACY 1000000
#define KVEC 3
#define W 4
#define DELIM ":"

#define FRECHET 1
#define DTW 2

#define THRESHOLD 0.2


// Declare functions of general use
double** find_grid_curve(curve *, double, double);

node *create_hashtable_node(double **, curve *, int);

double *grid_curve_to_1D(double **, int, int);

int is_centroid(node **, node *, int);

#endif