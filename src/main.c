#include "functions.h"
#include "preprocessing.h"
#include "clustering.h"
#include "clustering_init.h"
#include "clustering_assignment.h"
#include "clustering_update.h"
#include "hashtable.h"

int main(int argc, char* argv[]) {
	int i, j, K, L, k_cluster, complete = 0, metric_function = 0;
	char input_file[PATH_LENGTH] = "NULL", conf_file[128] = "NULL", output_file[PATH_LENGTH] = "NULL", user_input[128] = "NULL";

	config *config_info;
	dataset_info *dataset_information;
	hashtable **ht;

    // Parse arguments from command line
    if (argc >= 1 && argc <= 10) {
        for (i = 1 ; i <= (argc-1) ; i++) {
            if (!strcmp(argv[i], "-i")) { // If flag is -i ..
                stpcpy(input_file, argv[++i]); // Copy of argument's value to var
                continue;
            }
            if (!strcmp(argv[i], "-c")) { // If flag is -c ..
                stpcpy(conf_file, argv[++i]); // Copy of argument's value to var
                continue;
            }
            if (!strcmp(argv[i], "-o")) { // If flag is -o ..
                stpcpy(output_file, argv[++i]); // Copy of argument's value to var
                continue;
            }
            if (!strcmp(argv[i], "-complete")) { // If flag is -stats ..
                complete = 1;	// Set this var True
                continue;
            }
            if (!strcmp(argv[i], "-d")) { // If flag is -d ..
                if (!strcmp(argv[++i], "Frechet")) // Check if metric is "Frechet"
                    metric_function = FRECHET;
                else if (!strcmp(argv[i], "DTW")) // or if it's "DTW"
                    metric_function = DTW;
                else {	// Otherwise .. the input is wrong!
                    printf("Wrong metric function argument! \n");
                    return -1;
                }

                continue;
            }
            // If flag of any argument is wrong .. exit!
            printf("Wrong input arguments! \n");
            return -2;
        }
    }
    else {	// If number of arguments is wrong .. exit!
        printf("Wrong number of arguments! \n");
        return -3;
    }

    // Parse input from user
    if (!strcmp(input_file, "NULL")) {
        printf("Please insert the path to input file: ");
        scanf("%s", input_file);
    }

	if (!strcmp(conf_file, "NULL")) {
        printf("Please insert the path to configuration file: ");
        scanf("%s", conf_file);
    }

	if (!strcmp(output_file, "NULL")) {
        printf("Please insert the path to output file: ");
        scanf("%s", output_file);
    }

    if (!metric_function) {
        printf("Please insert the metric function: ");
        scanf("%s", user_input);

        if (!strcmp(user_input, "Frechet")) // Check if metric is "DFT"
            metric_function = FRECHET;
        else if (!strcmp(user_input, "DTW")) // or if it's "DTW"
            metric_function = DTW;
        else {	//Otherwise .. the input is wrong!
            printf("Wrong metric function argument! \n");
            return -1;
        }
    }

	// Read config file to retrieve information	(k_cluster, K, L)
	config_info = get_config_information(conf_file);

	// Seed random number generator
	srand(time(NULL));

	// Parse of dataset and creation of L hashtables
	ht = preprocessing(input_file, config_info->L, config_info->K);

	dataset_information = get_dataset_information(input_file);

	for (i = 1 ; i < 3 ; i++) {
		for (j = 1 ; j < 3 ; j++) {
			printf("i: %d ,j: %d\n",i,j);

			if (metric_function == DTW)
				clustering(output_file, ht, config_info->L, config_info->k_cluster, dataset_information->number_of_curves, metric_function, i, j, 2, complete);
			else {
				clustering(output_file, ht, config_info->L, config_info->k_cluster, dataset_information->number_of_curves, metric_function, i, j, 1, complete);
				clustering(output_file, ht, config_info->L, config_info->k_cluster, dataset_information->number_of_curves, metric_function, i, j, 2, complete);
			}
		}
	}

	hashtable_destroy(ht, config_info->L);
	free(config_info);
	free(dataset_information);

	return 0;
}
