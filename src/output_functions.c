#include "functions.h"
#include "output_functions.h"

// Function that appends to output file
void output(char *output_file, int init, int assign, int update, int metric_function, int k_cluster, double time,
			node **centroids, cluster **clusters, double *silhouette, int complete) {
  	int i, j, z;

  	// Open the output_file to append new information
	FILE *output = fopen(output_file, "a");

  	// Check if the file opened successfully
  	if (output == NULL) {
		printf("Fopen: error opening output file!\n");
      	exit(1);
    }


  	// Let's start writing in file..
	fprintf(output, "Algorithm: I%dA%dU%d\n", init, assign, update);

  	if (metric_function == FRECHET)
    	fprintf(output, "Metric: Frechet\n");
    else
        fprintf(output, "Metric: DTW\n");


	if (update == 1) {
		for (i = 0 ; i < k_cluster ; i++) {
			fprintf(output, "CLUSTER-%d {size: %d}\n", (i+1), clusters[i]->cluster_size);
			fprintf(output, "{centroid}\n");
			for (j = 0 ; j < centroids[i]->real_curve->points ; j++)
				for (z = 0 ; z < centroids[i]->real_curve->dimension; z++)
					fprintf(output, "%.15f ", centroids[i]->real_curve->coordinates[j][z]);

		fprintf(output, "\n");
		}
	}
	else {
		for (i = 0 ; i < k_cluster ; i++)
			fprintf(output, "CLUSTER-%d {size: %d, centroid: %d}\n", (i+1), clusters[i]->cluster_size, centroids[i]->real_curve->ID_curve);
	}

	fprintf(output, "clustering_time: %.1f seconds\n", time);

	fprintf(output, "Silhouette: [");
	for (i = 0 ; i < k_cluster ; i++)
      	fprintf(output, "%.5f, ", silhouette[i]);
	fprintf(output, "%.5f]\n", silhouette[k_cluster]);

	if (complete) {
		for (i = 0 ; i < k_cluster ; i++) {
			fprintf(output, "CLUSTER-%d { ", (i+1));
			for (j = 0 ; j < clusters[i]->cluster_size-1 ; j++)
				fprintf(output, "%d, ", clusters[i]->curves_in_cluster[j]->real_curve->ID_curve);
			fprintf(output, "%d }\n", clusters[i]->curves_in_cluster[clusters[i]->cluster_size-1]->real_curve->ID_curve);
		}
	}

  	fprintf(output, "\n\n");

  	// When the writing is finished, we close the file!
	fclose(output);

	return;
}
