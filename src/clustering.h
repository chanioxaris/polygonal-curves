#ifndef _CLUSTERING_
#define _CLUSTERING_

// Declare functions of clustering.c file
void clustering(char *, hashtable **, int, int, int, int, int, int, int, int);

double centroids_transposition(node **, node **, int, int);

double *silhouette(node **, cluster **, int, int);

#endif