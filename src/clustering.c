#include "functions.h"
#include "clustering.h"
#include "clustering_assignment.h"
#include "clustering_init.h"
#include "clustering_update.h"
#include "hashtable.h"
#include "metric_functions.h"
#include "output_functions.h"

double frechet_distance(double **, double **, int, int, int);
double DTW_distance(double **, double **, int, int, int);

// Basic function for clustering
void clustering(char *output_file, hashtable **ht, int L, int k_cluster, int dataset_size, int metric_function, int init, int assign, int update, int complete) {
	int i, iterations = 0;
	double difference, difference_old = -1;
	double *sil;

	time_t start, end;
	node **centroids, **centroids_old;
	cluster **clusters;

	centroids_old = (node**) malloc(k_cluster * sizeof(node*));

	if (centroids_old == NULL) {
		printf("Malloc: memory allocation error!\n");
		exit(3);
	}

	start = time(NULL);

	if (init == 1)
		centroids = Init_1(ht[0], k_cluster, dataset_size, metric_function);
	else
		centroids = Init_2(ht[0], k_cluster, dataset_size);

	if (assign == 1)
		clusters = Assignment_1(ht[0], k_cluster, centroids, dataset_size, metric_function);
	else
		clusters = Assignment_2(ht, L, k_cluster, centroids, dataset_size, metric_function);

	if (update == 1)
		Mean_Discrete_Frechet(centroids, clusters, k_cluster);
	else
		PAM(centroids, clusters, k_cluster, metric_function);

	for (i = 0 ; i < k_cluster ; i++)
		free(clusters[i]);
	free(clusters);

	do {
		for (i = 0 ; i < k_cluster ; i++)
			centroids_old[i] = centroids[i];

		if (assign == 1)
			clusters = Assignment_1(ht[0], k_cluster, centroids, dataset_size, metric_function);
		else
			clusters = Assignment_2(ht, L, k_cluster, centroids, dataset_size, metric_function);

		if (update == 1)
			Mean_Discrete_Frechet(centroids, clusters, k_cluster);
		else
			PAM(centroids, clusters, k_cluster, metric_function);

		difference = centroids_transposition(centroids, centroids_old, k_cluster, metric_function);

		if (difference == difference_old)
			iterations++;
		else
			iterations = 0;

		difference_old = difference;

	} while(difference > THRESHOLD && iterations < 3); //no centroids transposition

	end = time(NULL);

	// Calculate silhouette for the final clusters
	sil = silhouette(centroids, clusters, k_cluster, metric_function);

	// Write to output file
	output(output_file, init, assign, update, metric_function, k_cluster, end-start, centroids, clusters, sil, complete);

	printf("\n");

	free(centroids);
	free(centroids_old);
	free(sil);

	for (i = 0 ; i < k_cluster ; i++)
		free(clusters[i]);
	free(clusters);

	return;
}

// Calculate transposition of old and new centroids
double centroids_transposition(node **centroids_new, node **centroids_old, int k_cluster, int metric_function) {
	int i;
	double transposition = 0;

	for (i = 0 ; i < k_cluster ; i++) {
		// If transposition is zero, centroid didn't change
		if (metric_function == FRECHET)
			transposition += frechet_distance(centroids_old[i]->real_curve->coordinates, centroids_new[i]->real_curve->coordinates,
										centroids_old[i]->real_curve->points, centroids_new[i]->real_curve->points, centroids_old[i]->real_curve->dimension);
		else
			transposition += DTW_distance(centroids_old[i]->real_curve->coordinates, centroids_new[i]->real_curve->coordinates,
										centroids_old[i]->real_curve->points, centroids_new[i]->real_curve->points, centroids_old[i]->real_curve->dimension);
	}

	printf("DIFFERENCE: %f\n", transposition);

	return transposition;
}


// Function that calculates the silhouette for each cluster and the total silhouette of the whole clustering procedure
double *silhouette(node **centroids, cluster **clusters, int k_cluster, int metric_function) {
	int i, j, k, index;
	double s, sum, avg_a, avg_b, s_total, min_distance;

	double *total_s = (double*) malloc((k_cluster+1) * sizeof(double));

	if (total_s == NULL) {
		printf("Malloc: memory allocation error!\n");
		exit(3);
	}

	// For every cluster..
	for (i = 0 ; i < k_cluster ; i++) {
		s_total = 0;

		// For every object in this cluster..
		for (j = 0 ; j < clusters[i]->cluster_size ; j++) {
			// Calculate the distance between this object and its centroid
			if (metric_function == FRECHET)
				sum = frechet_distance(centroids[i]->real_curve->coordinates, clusters[i]->curves_in_cluster[j]->real_curve->coordinates,
													centroids[i]->real_curve->points, clusters[i]->curves_in_cluster[j]->real_curve->points, centroids[i]->real_curve->dimension);
			else
				sum = DTW_distance(centroids[i]->real_curve->coordinates, clusters[i]->curves_in_cluster[j]->real_curve->coordinates,
													centroids[i]->real_curve->points, clusters[i]->curves_in_cluster[j]->real_curve->points, centroids[i]->real_curve->dimension);

			// Calculate the sum distance among this object and every other object in this cluster
			for (k = 0 ; k < clusters[i]->cluster_size ; k++) {
				if (metric_function == FRECHET)
					sum += frechet_distance(clusters[i]->curves_in_cluster[j]->real_curve->coordinates, clusters[i]->curves_in_cluster[k]->real_curve->coordinates,
														clusters[i]->curves_in_cluster[j]->real_curve->points, clusters[i]->curves_in_cluster[k]->real_curve->points, clusters[i]->curves_in_cluster[j]->real_curve->dimension);
				else
					sum += DTW_distance(clusters[i]->curves_in_cluster[j]->real_curve->coordinates, clusters[i]->curves_in_cluster[k]->real_curve->coordinates,
														clusters[i]->curves_in_cluster[j]->real_curve->points, clusters[i]->curves_in_cluster[k]->real_curve->points, clusters[i]->curves_in_cluster[j]->real_curve->dimension);
			}

			// Calculate the average distance for this specific object
			avg_a = (double) sum / clusters[i]->cluster_size;

			min_distance = INFINITY;

			// Calculate the second nearest centroid from this sepcific object
			for (k = 0 ; k < k_cluster ; k++)
				if (k != i) {
					if (metric_function == FRECHET)
						sum = frechet_distance(clusters[i]->curves_in_cluster[j]->real_curve->coordinates, centroids[k]->real_curve->coordinates,
															clusters[i]->curves_in_cluster[j]->real_curve->points, centroids[k]->real_curve->points, centroids[k]->real_curve->dimension);
					else
						sum = DTW_distance(clusters[i]->curves_in_cluster[j]->real_curve->coordinates, centroids[k]->real_curve->coordinates,
															clusters[i]->curves_in_cluster[j]->real_curve->points, centroids[k]->real_curve->points, centroids[k]->real_curve->dimension);

					if (sum < min_distance) {
						min_distance = sum;
						index = k;
					}
				}

			sum = min_distance;

			// And for this specific object calculate its shilhoutte for the second nearest centroid
			for (k = 0 ; k < clusters[index]->cluster_size ; k++) {
				if (metric_function == FRECHET)
					sum += frechet_distance(clusters[i]->curves_in_cluster[j]->real_curve->coordinates, clusters[index]->curves_in_cluster[k]->real_curve->coordinates,
														clusters[i]->curves_in_cluster[j]->real_curve->points, clusters[index]->curves_in_cluster[k]->real_curve->points, clusters[i]->curves_in_cluster[j]->real_curve->dimension);
				else
					sum += DTW_distance(clusters[i]->curves_in_cluster[j]->real_curve->coordinates, clusters[index]->curves_in_cluster[k]->real_curve->coordinates,
														clusters[i]->curves_in_cluster[j]->real_curve->points, clusters[index]->curves_in_cluster[k]->real_curve->points, clusters[i]->curves_in_cluster[j]->real_curve->dimension);
			}

			// Calculate the average distance for this specific object
			avg_b = (double) sum / (clusters[index]->cluster_size);

			if (avg_a >= avg_b)
				s = (double) (avg_b - avg_a) / avg_a;
			else
				s = (double) (avg_b - avg_a) / avg_b;

			s_total += s;
		}

		// Update total_s array with #k_clusters shilhouttes ( for every cluster respectively )
		total_s[i] = s_total/clusters[i]->cluster_size;
	}

	total_s[k_cluster] = 0;
	for (i = 0 ; i < k_cluster ; i++)
		total_s[k_cluster] += total_s[i];

	// Update total_s array's last position that is the mean shilhoutte of all clustering procedure
	total_s[k_cluster] = (double) (total_s[k_cluster] / k_cluster);

	return total_s;
}
