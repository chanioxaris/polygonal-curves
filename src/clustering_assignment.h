#ifndef _CLUSTERING_ASSIGNMENT_
#define _CLUSTERING_ASSIGNMENT_

// Declare functions of clustering_assignment.c file
cluster **Assignment_1(hashtable *, int , node **, int , int );

cluster **Assignment_2(hashtable **, int, int, node **, int, int);

double initial_min_range(node **, int , int );

int **centroids_index_table(hashtable **, node **, int, int);

int is_centroid_LSH(node **, node *, int);

void true_curve_assignments(hashtable *, cluster **, node **, int);

void LSH_assignments_lloyd(hashtable *, cluster **, node **, int, int);

#endif