#ifndef _CLUSTERING_UPDATE_
#define _CLUSTERING_UPDATE_

// Declare functions of clustering_update.c file
void PAM(node **, cluster **, int, int);

void Mean_Discrete_Frechet(node **, cluster **, int);

node *frechet_traversal_node(double **, double **, int, int, int);

int most_similar_curve_ID(cluster *, node *);

int min_3_index(double, double, double);

double *mean_coordinates(double *, double *, int);

#endif