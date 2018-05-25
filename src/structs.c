typedef struct _curve_ {
	int ID_curve;
	int dimension;
	int points;
	double **coordinates;

	// LSH Assignment
	int assigned;
	int assigned_number;
	int *assigned_ID;
}curve;

// Nodes for hashtable
typedef struct _node_ {
	curve *real_curve;
	int grid_points;
	double *grid_curve_1D;
	struct _node_ *next;
}node;


typedef struct _hashtable_ {
	int size;
	struct _node_ **table;
}hashtable;


typedef struct _dataset_info_ {
	int min_points;
	int max_points;
	int number_of_curves;
}dataset_info;


typedef struct _config_ {
	int k_cluster;
	int K;
	int L;
}config;


typedef struct _cluster_ {
	int cluster_size;
	node **curves_in_cluster;
}cluster;
