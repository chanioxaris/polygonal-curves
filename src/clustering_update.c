#include "functions.h"
#include "clustering_update.h"
#include "metric_functions.h"
#include "binary_tree.h"

// Implementation of Partitioning Around Medoids (PAM) – Improved update
void PAM(node **centroids, cluster **clusters, int k_cluster, int metric_function)
{
	int i, j, z, min_index;	
	double min;
	double *distances;
	
	node *tmp;	
	
	// For every cluster..
	for (i = 0 ; i < k_cluster ; i++)
	{
		distances = (double*) malloc(clusters[i]->cluster_size * sizeof(double));
		
		if (distances == NULL)
		{	
			printf("Malloc: memory allocation error!\n");
			exit(3);
		}
		
		// For every object in this cluster..
		for (j = 0 ; j < clusters[i]->cluster_size ; j++)
		{
			distances[j] = 0;
			
			// Find the sum distance with every other object, and update the array distances..
			for (z = 0 ; z < clusters[i]->cluster_size ; z++)
			{
				if (metric_function == FRECHET)
					distances[j] += frechet_distance(clusters[i]->curves_in_cluster[j]->real_curve->coordinates, clusters[i]->curves_in_cluster[z]->real_curve->coordinates,
														clusters[i]->curves_in_cluster[j]->real_curve->points, clusters[i]->curves_in_cluster[z]->real_curve->points, clusters[i]->curves_in_cluster[j]->real_curve->dimension);
				else
					distances[j] += DTW_distance(clusters[i]->curves_in_cluster[j]->real_curve->coordinates, clusters[i]->curves_in_cluster[z]->real_curve->coordinates,
														clusters[i]->curves_in_cluster[j]->real_curve->points, clusters[i]->curves_in_cluster[z]->real_curve->points, clusters[i]->curves_in_cluster[j]->real_curve->dimension);	
			}				
		}	

		min = INFINITY;
		
		// Find out the min sum distance and define this object as the new centroid
		for (j = 0 ; j < clusters[i]->cluster_size ; j++)
		{
			if (distances[j] < min)
			{
				min = distances[j];
				min_index = j;
			}
		}
		
		// Update the array with centroids
		centroids[i] = clusters[i]->curves_in_cluster[min_index];
		
		free(distances);
	}
	return;
}



void Mean_Discrete_Frechet(node **centroids, cluster **clusters, int k_cluster)
{
	int i, height;
	
	node **tree;
	
	for (i = 0 ; i < k_cluster ; i++)
	{	
		height = floor(log((double)(clusters[i]->cluster_size)) / log(2));
		
		if( pow(2.0, height) < (double)(clusters[i]->cluster_size) )   // 2^h <= n < 2^h+1.
			height++;
		
		// Binary Tree's creation
		tree = binary_tree_init(height);	

		// Insert cluster's curve in trees leaves
		binary_tree_insert_leaf(tree, clusters[i], height);
			
		// Assign the new centroid into 'centroids' and put the old centroid into 'clusters'
		centroids[i] = postOrderTraversal(tree, 0, height);	
	
		centroids[i]->real_curve->ID_curve = most_similar_curve_ID(clusters[i], centroids[i]);
		
		binary_tree_destroy(tree, height);
	}
	return;
}


node *frechet_traversal_node(double **curve1, double **curve2, int m1, int m2, int dimension)
{
  	int i, j, min_index, index = 0, P = m1, Q = m2;
	double distance;
	double **array, **traversal, **traversal_reversed;
	
	traversal_reversed = (double**) malloc((m1+m2) * sizeof(double*));

	if (traversal_reversed == NULL)
	{	
		printf("Malloc: memory allocation error!\n");
		exit(3);
	}

	for (i = 0 ; i < m1+m2 ; i++)
	{
		traversal_reversed[i] = (double*) malloc(dimension * sizeof(double));

		if (traversal_reversed[i] == NULL)
		{
			printf("Malloc: memory allocation error!\n");
			exit(3);
		}
	}

  	array = (double**) malloc(m1 * sizeof(double*));

	if (array == NULL)
	{
		printf("Malloc: memory allocation error!\n");
		exit(3);
	}

    for (i = 0 ; i < m1 ; i++)
	{
    	array[i] = (double*) malloc(m2 * sizeof(double));

		if (array[i] == NULL)
		{
			printf("Malloc: memory allocation error!\n");
			exit(3);
		}
	}

	// Frechet distance
    array[0][0] = euclidean(curve1[0], curve2[0], dimension);

    for (i = 1 ; i < m1 ; i++)
    	array[i][0] = max_2(array[i-1][0], euclidean(curve1[i], curve2[0], dimension));

    for (j = 1 ; j < m2 ; j++)
    	array[0][j] = max_2(array[0][j-1], euclidean(curve1[0], curve2[j], dimension));

    for (i = 1 ; i < m1 ; i++)
      	for (j = 1 ; j < m2 ; j++)
        	array[i][j] = max_2(min_3(array[i-1][j], array[i][j-1], array[i-1][j-1]), euclidean(curve1[i], curve2[j], dimension));

	traversal_reversed[index++] = mean_coordinates(curve1[--P], curve2[--Q], dimension);

	// Get traversal out of frechet array
	while (P != 0 || Q != 0)
	{
		if (!Q)
			min_index = 0;
		else if (!P)
			min_index = 1;
		else
			min_index = min_3_index(array[P-1][Q], array[P][Q-1], array[P-1][Q-1]);
		
		if (min_index == 0)
			traversal_reversed[index++] = mean_coordinates(curve1[--P], curve2[Q], dimension);
		else if (min_index == 1)
			traversal_reversed[index++] = mean_coordinates(curve1[P], curve2[--Q], dimension);
		else
			traversal_reversed[index++] = mean_coordinates(curve1[--P], curve2[--Q], dimension);	
	}

	traversal = (double**) malloc(index * sizeof(double*));

	if (traversal == NULL)
	{
		printf("Malloc: memory allocation error!\n");
		exit(3);
	}

	for (i = 0 ; i < index-1 ; i++)
	{
    	traversal[i] = (double*) malloc(dimension * sizeof(double));
		
		if (traversal[i] == NULL)
		{
			printf("Malloc: memory allocation error!\n");
			exit(3);
		}
	}

	for (i = 0 ; i < index ; i++)
		traversal[i] = traversal_reversed[index-i-1];
	
	
	// Create node
	node *new_node;

	// Allocate memory for new_node
	new_node = (node*) malloc(sizeof(node));

	if (new_node == NULL)
	{
		printf("Malloc: memory allocation error!\n");
		exit(3);
	}

	new_node->real_curve = (curve*) malloc(sizeof(curve));

	if (new_node->real_curve == NULL)
	{
		printf("Malloc: memory allocation error!\n");
		exit(3);
	}
	
	// Fill struct with data	
	new_node->real_curve->ID_curve = -1;
	new_node->real_curve->dimension = dimension;
	new_node->real_curve->points = index;
	new_node->real_curve->coordinates = traversal;
	

	for (i = 0 ; i < m1 ; i++)
		free(array[i]);
	free(array);
		
	for (i = 0 ; i < m1+m2 ; i++)
		traversal_reversed[i] = NULL;
	
	for (i = 0 ; i < m1+m2 ; i++)
		free(traversal_reversed[i]);
	free(traversal_reversed);
	
	return new_node;
}


// Auxilary function to find nearest dataset curve on the same cluster of mean discrete travesal 
int most_similar_curve_ID(cluster *clusters, node *centroid)
{
	int i, min_index;
	double distance, min_distance = INFINITY;
	
	for (i = 0 ; i < clusters->cluster_size ; i++)
	{
		distance = frechet_distance(clusters->curves_in_cluster[i]->real_curve->coordinates, centroid->real_curve->coordinates,
									clusters->curves_in_cluster[i]->real_curve->points, centroid->real_curve->points, 
									clusters->curves_in_cluster[i]->real_curve->dimension);
				
		if (distance < min_distance)
		{
			min_distance = distance;
			min_index = clusters->curves_in_cluster[i]->real_curve->ID_curve;
		}
	}
	return min_index;
}


// Auxilary function to get minimum value out of three
int min_3_index(double up, double left, double diagonal)
{
	if (up <= left && up <= diagonal)
      	return 0;
  	else if (left <= up && left <= diagonal)
      	return 1;
    else
      	return 2;
}


// Auxilary function to calculate mean coordinates out of two curves
double *mean_coordinates(double *coord1, double *coord2, int dimension)
{
	int i;
	double *mean;
	
	mean = (double*) malloc(dimension * sizeof(double));
	
	if (mean == NULL)
	{	
		printf("Malloc: memory allocation error!\n");
		exit(3);
	}
	
	for (i = 0 ; i < dimension ; i++)
		mean[i] = (double) ((coord1[i] + coord2[i])/2);
		
	return mean;
}