#ifndef _BINARY_TREE_
#define _BINARY_TREE_

// Declare functions of binary_tree.c file
node **binary_tree_init(int);

void binary_tree_insert_leaf(node **, cluster *, int);

void binary_tree_destroy(node **, int);

node *postOrderTraversal(node **, int, int);

#endif