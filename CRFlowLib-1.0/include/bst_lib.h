
/*
* -----------------------------------------------------------------
*  Binary Search Tree Library --- bst_lib.h
*  Version: 1.6180
*  Date: Feb 19, 2010
* ----------------------------------------------------------------- 
*  Programmer: Americo Barbosa da Cunha Junior
*              americo.cunhajr@gmail.com
* -----------------------------------------------------------------
*  Copyright (c) 2010 by Americo Barbosa da Cunha Junior
*
*  This program is free software: you can redistribute it and/or
*  modify it under the terms of the GNU General Public License as
*  published by the Free Software Foundation, either version 3 of
*  the License, or (at your option) any later version.
*
*  This program is distributed in the hope that it will be useful,
*  but WITHOUT ANY WARRANTY; without even the implied warranty of
*  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
*  GNU General Public License for more details.
*
*  A copy of the GNU General Public License is available in
*  LICENSE.txt or http://www.gnu.org/licenses/.
* -----------------------------------------------------------------
*  This is the header file for a library with
*  to operates a binary search tree.
* -----------------------------------------------------------------
*/




#ifndef __BST_LIB_H__
#define __BST_LIB_H__


#include <gsl/gsl_vector.h>
#include <gsl/gsl_matrix.h>


#undef  RIGHT
#define RIGHT 1

#undef  LEFT
#define LEFT 0




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

int bst_node_add(int side,
                 bst_node *old_node,
                 bst_node *new_node);

int bst_search(bst_node *root,
               gsl_vector *phi,
               bst_node **end_node,
               bst_leaf **end_leaf);

int bst_height(bst_node *root);

#endif /* __BST_LIB_H__ */
