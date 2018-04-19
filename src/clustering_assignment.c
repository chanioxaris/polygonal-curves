#include "functions.h"
#include "clustering_assignment.h"
#include "hashtable.h"
#include "metric_functions.h"

double frechet_distance(double **, double **, int, int, int);
double DTW_distance(double **, double **, int, int, int);

// Implementation of Lloyd’s assignment (simplest approach)
cluster **Assignment_1(hashtable *ht, int k_cluster, node **centroids, int dataset_size, int metric_function)
{
	int i, j, ID;
	double min_D, distance;
	
	node *tmp;
	cluster **clusters;
		
	// Initilization of k clusters
	clusters = (cluster**) malloc(k_cluster * sizeof(cluster*));
	
	if (clusters == NULL)
	{
		printf("Malloc: memory allocation error!\n");
		exit(3);
	}
	
	for(i = 0; i < k_cluster; i++)
	{
		clusters[i] = (cluster*) malloc(sizeof(cluster));
		
		if (clusters[i] == NULL)
		{
			printf("Malloc: memory allocation error!\n");
			exit(3);
		}
		
		clusters[i]->cluster_size = 1;
		clusters[i]->curves_in_cluster = (node**) malloc(sizeof(node*));
		
		if (clusters[i]->curves_in_cluster == NULL)
		{
			printf("Malloc: memory allocation error!\n");
			exit(3);
		}
		
		clusters[i]->curves_in_cluster[0] = centroids[i];
	}
	
	// For every curve of dataset
	for(j = 0; j < dataset_size; j++)
	{		
		// Find each curve in hashtable
		tmp = hashtable_crossing(ht, j+1);
	
		// If its not a centroid
		if (!is_centroid(centroids, tmp, i))
		{	
			min_D = INFINITY;
	
			// For every centroid
			for(i = 0; i < k_cluster; i++)
			{	
				// Calculate the distances from every centroid 
				if (metric_function == FRECHET)
					distance = frechet_distance(centroids[i]->real_curve->coordinates, tmp->real_curve->coordinates,
												centroids[i]->real_curve->points, tmp->real_curve->points, tmp->real_curve->dimension);
				else
					distance = DTW_distance(centroids[i]->real_curve->coordinates, tmp->real_curve->coordinates,
											centroids[i]->real_curve->points, tmp->real_curve->points, tmp->real_curve->dimension);

				// And keep the minimum distance
				if (distance < min_D){
					min_D = distance;
					ID = i;
				}
			}
						
			clusters[ID]->cluster_size++;
			clusters[ID]->curves_in_cluster = (node**) realloc(clusters[ID]->curves_in_cluster, clusters[ID]->cluster_size * sizeof(node));
			
			if (clusters[ID]->curves_in_cluster == NULL)
			{
				printf("Malloc: memory allocation error!\n");
				exit(3);
			}
			
			clusters[ID]->curves_in_cluster[clusters[ID]->cluster_size-1] = tmp;
		}
	}	
	return clusters;
}


// Implementation of LSH assignment
cluster **Assignment_2(hashtable **ht, int L, int k_cluster, node **centroids, int dataset_size, int metric_function)
{
	int i, j, assigned_curves_counter;
	int **centroids_index;
	double range, distance;

	node *tmp;
	cluster **clusters;

	// Calculate the min distance between two centroids
	range = initial_min_range(centroids, k_cluster, metric_function)/2;
	
	// Initilization of k clusters
	clusters = (cluster**) malloc(k_cluster * sizeof(cluster*));

	if (clusters == NULL)
	{
		printf("Malloc: memory allocation error!\n");
		exit(3);
	}

	for(i = 0; i < k_cluster; i++)
	{
		clusters[i] = (cluster*) malloc(sizeof(cluster));
		
		if (clusters[i] == NULL)
		{
			printf("Malloc: memory allocation error!\n");
			exit(3);
		}
		
		clusters[i]->cluster_size = 1;
		clusters[i]->curves_in_cluster = (node**) malloc(sizeof(node*));
		
		if (clusters[i]->curves_in_cluster == NULL)
		{
			printf("Malloc: memory allocation error!\n");
			exit(3);
		}
		
		clusters[i]->curves_in_cluster[0] = centroids[i];
	}

	// Every single curve..
	for(i = 0; i < dataset_size; i++)
	{
		tmp = hashtable_crossing(ht[0], i+1);
	
		tmp->real_curve->assigned = 0;		// ..initially is not assigned to any centroid
		tmp->real_curve->assigned_number = 0;
		tmp->real_curve->assigned_ID = (int *) malloc(tmp->real_curve->assigned_number * sizeof(int));
		
		if (tmp->real_curve->assigned_ID == NULL)
		{
			printf("Malloc: memory allocation error!\n");
			exit(3);
		}
		
	}

	for(i = 0; i < k_cluster; i++)
		centroids[i]->real_curve->assigned = 1;		
	
	centroids_index = centroids_index_table(ht, centroids, k_cluster, L);
	
	do
	{
		assigned_curves_counter = 0;

		for(i = 0; i < k_cluster; i++)
		{
			for(j = 0; j < L; j++)
			{
				tmp = ht[j]->table[centroids_index[i][j]];

				while (tmp != NULL)
				{
					if (tmp->real_curve->assigned == 1)
					{
						tmp = (tmp->next);
						continue;
					}
					
					// Calculate distance
					if (metric_function == FRECHET)
						distance = frechet_distance(centroids[i]->real_curve->coordinates, tmp->real_curve->coordinates, 
									centroids[i]->real_curve->points, tmp->real_curve->points, tmp->real_curve->dimension);
					else
						distance = DTW_distance(centroids[i]->real_curve->coordinates, tmp->real_curve->coordinates, 
									centroids[i]->real_curve->points, tmp->real_curve->points, tmp->real_curve->dimension);	
									
					
					// R diameter neighbors						
					if (distance <= range)
					{
						tmp->real_curve->assigned_ID = (int*) realloc(tmp->real_curve->assigned_ID, ++(tmp->real_curve->assigned_number) * sizeof(int));
						
						if (tmp->real_curve->assigned_ID == NULL)
						{
							printf("Realloc: memory allocation error1!\n");
							exit(3);
						}
						
						tmp->real_curve->assigned_ID[tmp->real_curve->assigned_number-1] = i;
						assigned_curves_counter++;
					}	
					
					// Step to the next node in bucket
					tmp = (tmp->next);
				}				
			}
		}
		// Choose the real assignment of a curve among many centroids ( if this objects belongs in the range of many centroids )
		true_curve_assignments(ht[0], clusters, centroids, metric_function);
		range *= 2;
		
	}while(assigned_curves_counter != 0);

	// For those curves that havent assigned yet.. find the closest centroid and assign them
	LSH_assignments_lloyd(ht[0], clusters, centroids, k_cluster, metric_function);
	
	for (i = 0 ; i < k_cluster ; i++)
		free(centroids_index[i]);
	free(centroids_index);
	
	return clusters;
}


// Auxiliary function that calculates the minimum distance among all centroids
double initial_min_range(node **centroids, int k_cluster, int metric_function)
{
	int i, j;
	double distance, minD[k_cluster], min = INFINITY;

	// For every single centroid
	for(i = 0; i < k_cluster; i++)
	{	
		minD[i] = INFINITY;
		
		// For every other centroid
		for(j = 0; j < k_cluster; j++)
		{	
			if ( i == j ) // If the centroid is the same skip it
				continue;
			
			// Calculate the distances
			if (metric_function == FRECHET)
				distance = frechet_distance(centroids[i]->real_curve->coordinates, centroids[j]->real_curve->coordinates,
													centroids[i]->real_curve->points, centroids[j]->real_curve->points, centroids[j]->real_curve->dimension);
			else
				distance = DTW_distance(centroids[i]->real_curve->coordinates, centroids[j]->real_curve->coordinates,
													centroids[i]->real_curve->points, centroids[j]->real_curve->points, centroids[j]->real_curve->dimension);
			// and keep the minimum
			if (distance < minD[i])
				minD[i] = distance;
		}
	}
	
	// Finally find the minimum distance among all centroids
	for(i = 0; i < k_cluster; i++)
	{
		if (minD[i] < min)
			min = minD[i];	
	}	
	return min;
}


// Auxilary function that finds the index in all hashtables, that centroids , belong
int **centroids_index_table(hashtable **ht, node **centroids, int k_cluster, int L)
{
	int i, j, index;
	int **centroids_index;
	node *tmp_node;

	// Centroids_index will be an array k_cluster x L that coints the index, for every centroid , in every hashtable
	centroids_index = (int**) malloc(k_cluster * sizeof(int*));
	
	if (centroids_index == NULL)
	{
		printf("Malloc: memory allocation error!\n");
		exit(3);
	}

	for (i = 0 ; i < k_cluster ; i++)
	{
		centroids_index[i] = (int*) malloc(L * sizeof(int));
		
		if (centroids_index[i] == NULL)
		{
			printf("Malloc: memory allocation error!\n");
			exit(3);
		}
	}
	
	// For every hashtable
	for (i = 0 ; i < L ; i++)
	{	
		// For every bucket
		for (j = 0 ; j < ht[i]->size ; j++)
		{	
			// For every node
			tmp_node = ht[i]->table[j];		

			while (tmp_node != NULL)
			{
				// Check if node is a centroid
				index = is_centroid_LSH(centroids, tmp_node, k_cluster);
				
				if (index != -1)
					centroids_index[index][i] = j;

				tmp_node = tmp_node->next;
			}
		}	
	}	
	return centroids_index;
}


// Auxilary function that checks if a curve is centroid , and returns the index of its position in array centroids
int is_centroid_LSH(node **centroids, node *new, int k_cluster)
{
	int i;
	
	for(i = 0; i < k_cluster ; i++)
		if (centroids[i]->real_curve->ID_curve == new->real_curve->ID_curve)
			return i;
		
	return -1;
}


// Auxilary function that finds the real assignment of a curve, to a centroid ( in the case that it is assigned to many centroids )
void true_curve_assignments(hashtable *ht, cluster **clusters, node **centroids, int metric_function)
{
	int i, j, index, min_index;
	double distance, min_distance;
	
	node *tmp;
	
	// In every bucket	
	for (i = 0 ; i < ht->size ; i++)
	{
		tmp = ht->table[i];		
		
		while (tmp != NULL)
		{
			// If a curve isnt assigned yet..
			if (tmp->real_curve->assigned == 0)
			{	
				// and it belongs in more than one centroid's range 
				if (tmp->real_curve->assigned_number > 1)
				{
					min_distance = INFINITY;
					
					// find out which centroid is closer
					for (j = 0 ; j < tmp->real_curve->assigned_number ; j++) 
					{
						index = tmp->real_curve->assigned_ID[j]; 
						
						if (metric_function == FRECHET)
							distance = frechet_distance(centroids[index]->real_curve->coordinates, tmp->real_curve->coordinates,
														centroids[index]->real_curve->points, tmp->real_curve->points, tmp->real_curve->dimension);
						else
							distance = DTW_distance(centroids[index]->real_curve->coordinates, tmp->real_curve->coordinates,
													centroids[index]->real_curve->points, tmp->real_curve->points, tmp->real_curve->dimension);													
													
						if (distance < min_distance)
						{
							min_distance = distance;
							min_index = index;
						}																				
					}

					// From now on this curve will be assigned to only one centroid ( the closest )
					tmp->real_curve->assigned = 1;
					
					// Update clusters' info
					clusters[min_index]->cluster_size++;
					clusters[min_index]->curves_in_cluster = (node**) realloc(clusters[min_index]->curves_in_cluster, clusters[min_index]->cluster_size * sizeof(node));
					
					if (clusters[min_index]->curves_in_cluster == NULL)
					{
						printf("Malloc: memory allocation error!\n");
						exit(3);
					}
					
					
					clusters[min_index]->curves_in_cluster[clusters[min_index]->cluster_size-1] = tmp;		
				}
			}
			tmp = tmp->next;
		}
	}	
	return;
}


// Auxilary function that assigns all unassigned curves after LSH assignments
void LSH_assignments_lloyd(hashtable *ht, cluster **clusters, node **centroids, int k_cluster, int metric_function)
{	
	int i, j, min_index;
	double distance, min_distance;
	
	node *tmp;
	
	// For the whole dataset..
	for (i = 0 ; i < ht->size ; i++)
	{
		tmp = ht->table[i];		

		while (tmp != NULL)
		{	
			// If a curve isnt assigned yet..
			if (tmp->real_curve->assigned == 0)
			{				
				min_distance = INFINITY;
				
				for (j = 0 ; j < k_cluster ; j++) 
				{
					// Compare the distances between every other centroid..
					if (metric_function == FRECHET)
						distance = frechet_distance(centroids[j]->real_curve->coordinates, tmp->real_curve->coordinates,
													centroids[j]->real_curve->points, tmp->real_curve->points, tmp->real_curve->dimension);
					else
						distance = DTW_distance(centroids[j]->real_curve->coordinates, tmp->real_curve->coordinates,
												centroids[j]->real_curve->points, tmp->real_curve->points, tmp->real_curve->dimension);													
					
					// and choose the closest
					if (distance < min_distance)
					{
						min_distance = distance;
						min_index = j;
					}																				
				}

				// From now on this curve will be assigned to only one centroid ( the closest )
				tmp->real_curve->assigned = 1;
				
				// Update clusters' info
				clusters[min_index]->cluster_size++;
				clusters[min_index]->curves_in_cluster = (node**) realloc(clusters[min_index]->curves_in_cluster, clusters[min_index]->cluster_size * sizeof(node));
				
				if (clusters[min_index]->curves_in_cluster == NULL)
				{
					printf("Malloc: memory allocation error!\n");
					exit(3);
				}
				
				clusters[min_index]->curves_in_cluster[clusters[min_index]->cluster_size-1] = tmp;		
			}
			tmp = tmp->next;
		}
	}	
	return;
}