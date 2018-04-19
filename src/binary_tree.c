#include "functions.h"
#include "binary_tree.h"
#include "clustering_update.h"

// Function that initializes the binary tree ( implementation with array )
node **binary_tree_init(int h)
{
	// Tree's size is 2^(h+1) - 1
	int i, size = pow(2, h+1) - 1;
	
	node **tree = (node**) malloc(size * sizeof(node*));
	
	if (tree == NULL)
	{
		printf("Malloc: memory allocation error!\n");
		exit(3);
	}

	// Initialize every tree node to NULL
	for (i = 0 ; i < size ; i++)
		tree[i] = NULL;

	return tree;
}


// Function that inserts objects in leafs
void binary_tree_insert_leaf(node **tree, cluster *clusters, int h)
{
	// Offset is actually the first node of the array that is a leaf
	int i, offset = pow(2, h) - 1;
	
	// Update the array's indexes that are leaves
	for (i = 0 ; i < clusters->cluster_size ; i++)
		tree[i+offset] = clusters->curves_in_cluster[i];

	return;
}


// Function to free the memory allocated for binary tree
void binary_tree_destroy(node **tree, int h)
{
	int i, j, size = pow(2, h) - 1;
		
	for (i = 1 ; i < size ; i++)
	{	
		if (tree[i] != NULL)
		{
			for (j = 0 ; j < tree[i]->real_curve->points ; j++)
				free(tree[i]->real_curve->coordinates[j]);
			free(tree[i]->real_curve->coordinates);
			free(tree[i]->real_curve);
		}
				
		free(tree[i]);
	}
	free(tree);
	
	return;
}


// Function that recursively finds the mean (frechet) curve among many curves
node *postOrderTraversal(node **tree, int index, int h)
{
	node *leftCurve, *rightCurve;

	// Return if node is leaf
	if(index >= (pow(2, h) - 1))
		return tree[index];
	
	else
	{	// Firstly we check the left child until reach leaf node
		leftCurve = postOrderTraversal(tree, 2*index+1, h);
				
		// Secondly we check the rigth child until reach leaf node
		rightCurve = postOrderTraversal(tree, 2*index+2, h);
		
		if (leftCurve == NULL)
			return NULL;
		
		if (rightCurve == NULL)
			return leftCurve;
	
		// Get the mean optimal traversal of both childs
		tree[index] = frechet_traversal_node(leftCurve->real_curve->coordinates, rightCurve->real_curve->coordinates, 
											leftCurve->real_curve->points, rightCurve->real_curve->points, 
											leftCurve->real_curve->dimension);
										
		return tree[index];
	}
}