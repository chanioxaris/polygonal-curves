#include "functions.h"
#include "clustering_init.h"
#include "hashtable.h"
#include "preprocessing.h"
#include "metric_functions.h"

double frechet_distance(double **, double **, int, int, int);
double DTW_distance(double **, double **, int, int, int);

// Implementation of K-means++
node **Init_1(hashtable *ht, int k_cluster, int dataset_size, int metric_function)
{
	int i, j, z, index, max_P;	
	double distance, x, max;
	double *D, *P;	
	node *tmp;
	node **centroids;

	// Array with min distances
	D = (double*) malloc(dataset_size * sizeof(double));
	
	if (D == NULL)
	{	
		printf("Malloc: memory allocation error!\n");
		exit(3);
	}
	
	// Array with probabilities
	P = (double*) malloc(dataset_size * sizeof(double));

	if (P == NULL)
	{	
		printf("Malloc: memory allocation error!\n");
		exit(3);
	}	
	
	centroids = (node**) malloc(k_cluster * sizeof(node*));
	
	if (centroids == NULL)
	{	
		printf("Malloc: memory allocation error!\n");
		exit(3);
	}
	
	//Choose first centroid uniformly at random
	index = (int) (rand_gaussian() * dataset_size);
	//Assign it into array with centroids
	centroids[0] = hashtable_crossing(ht, index+1);	

	// Initialization of arrays
	for(i = 1; i < k_cluster ; i++)
	{ 
		for(j = 0; j < dataset_size; j++)
		{
			D[j] = INFINITY;		
			P[j] = 0.0;			
		}

		// For every single curve
		for(j = 0; j < dataset_size; j++)
		{
			// Get from hashtable every single object
			tmp = hashtable_crossing(ht, j+1);
				
			// If its not already a centroid
			if (!is_centroid(centroids, tmp, i))
			{
				// Calculate the distances from every centroid 
				for(z = 0 ; z < i ; z++)
				{
					if (metric_function == FRECHET)
						distance = frechet_distance(centroids[z]->real_curve->coordinates, tmp->real_curve->coordinates,
													centroids[z]->real_curve->points, tmp->real_curve->points, tmp->real_curve->dimension);
					else
						distance = DTW_distance(centroids[z]->real_curve->coordinates, tmp->real_curve->coordinates,
													centroids[z]->real_curve->points, tmp->real_curve->points, tmp->real_curve->dimension);
				}
				// and keep the minimum distance
				if (distance < D[j])
					D[j] = distance;
			}
			
		}

		max = -INFINITY;

		// Find the max value of minimum distances in D
		for(j = 0; j < dataset_size; j++)
			if (D[j] != INFINITY && D[j] > max)
				max = D[j];
	
		// Creation of P(r) array
		for(j = 0; j < dataset_size; j++)
		{
			if (D[j] != INFINITY)
			{
				for(z = 0; z < j; z++){
					if (D[z] != INFINITY)
						P[j] += pow(D[z]/max, 2);
				}					
			}				
		}

		// Keep in max_P variable the max value of P array
		for(j = dataset_size-1 ; j >= 0; j--)
			if (P[j] != 0)
			{
				max_P = P[j];
				break;
			}
		// Choose randomly a real number in space [0, max_P]
		x = ((float)rand()/(float)(RAND_MAX)) * max_P;
		
		for(j = 0; j < dataset_size; j++)
		{				
			if (x <= P[j]){
				//this is the curve that we'll choose for centroid
				centroids[i] = hashtable_crossing(ht, j+1);
				break;
			}	
		}
	}
	
	free(D);
	free(P);
	
	return centroids;
}

// Implementation of Random selection of k-points (simplest)
node **Init_2(hashtable *ht, int k_cluster, int dataset_size)
{
	int i, index;
	node **centroids;
	
	centroids = (node**) malloc(k_cluster * sizeof(node*));
	
	if (centroids == NULL)
	{	
		printf("Malloc: memory allocation error!\n");
		exit(3);
	}

	for(i = 0; i < k_cluster; i++)
	{
		do{
			index = (int) (rand_gaussian() * dataset_size) + 1;
			centroids[i] = hashtable_crossing(ht, index);
			
		}while(same_centroids(centroids, i));  // Keep choosing centroids till they are all different
	}
	return centroids;
}


// Auxiliary function that returns True if current curve is a centroid
int same_centroids(node **centroids, int index)
{
	int i;
	
	for(i = 0; i < index; i++)
		if(centroids[i]->real_curve->ID_curve == centroids[index]->real_curve->ID_curve)
			return 1;

	return 0;
}