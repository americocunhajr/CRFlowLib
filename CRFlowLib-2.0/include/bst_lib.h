
/*
* -----------------------------------------------------------------
*  bst_lib.h
*  Binary Search Tree Library
*  Version: 2.0
*  Last Update: Oct 9, 2019
*  
*  Programmer: Americo Barbosa da Cunha Junior
*              americo.cunhajr@gmail.com
* ----------------------------------------------------------------- 
*  Copyright (c) 2010-2019, Americo Barbosa da Cunha Junior
*  All rights reserved.
* -----------------------------------------------------------------
*  This is the header file for BST_LIB module, a computational 
*  library with routines to manipulate a binary search tree.
* -----------------------------------------------------------------
*/




#ifndef __BST_LIB_H__
#define __BST_LIB_H__


#include <gsl/gsl_vector.h>
#include <gsl/gsl_matrix.h>


#undef  BST_RIGHT
#define BST_RIGHT 1

#undef  BST_LEFT
#define BST_LEFT 0




/*
*------------------------------------------------------------
*   structure of binary search tree leaf
*
*   last update: Feb 2, 2009
*------------------------------------------------------------
*/

typedef struct leaf
{
    gsl_vector *phi;   /* composition */
    gsl_vector *Rphi;  /* reaction mapping */
    gsl_matrix *A;     /* mapping gradient matrix */
    gsl_matrix *L;     /* ellipsoid Cholesky matrix */
} bst_leaf;
/*------------------------------------------------------------*/





/*
*------------------------------------------------------------
*   structure of binary search tree node
*
*   last update: Feb 2, 2009
*------------------------------------------------------------
*/

typedef struct node
{
    gsl_vector *v;          /* cutting plane normal vector */
    double a;               /* scalar a */
    struct node *r_node;    /* right node */
    struct node *l_node;    /* left node */
    struct leaf *r_leaf;    /* right leaf */
    struct leaf *l_leaf;    /* left leaf */
} bst_node;
/*------------------------------------------------------------*/





/*
*------------------------------------------------------------
* function prototypes
*------------------------------------------------------------
*/

bst_leaf *bst_leaf_alloc();

bst_node* bst_node_alloc();

void bst_leaf_free(void **leaf_bl);

void bst_node_free(void **node_bl);

int bst_leaf_set(gsl_vector *phi,
                 gsl_vector *Rphi,
                 gsl_matrix *A,
                 gsl_matrix *L,
                 bst_leaf *leaf);

double bst_cutplane(gsl_vector *phi_l,
                    gsl_vector *phi_r,
                    gsl_vector *v);

int bst_node_set(bst_leaf *l_leaf,
                 bst_leaf *r_leaf,
                 bst_node *node);

int bst_node_add(unsigned int side,
                 bst_node *old_node,
                 bst_node *new_node);

int bst_search(bst_node *root,
               gsl_vector *phi,
               bst_node **end_node,
               bst_leaf **end_leaf);

int bst_height(bst_node *root);

#endif /* __BST_LIB_H__ */
