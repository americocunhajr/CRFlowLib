
/*
* -----------------------------------------------------------------
*  Binary Search Tree Library --- bst_lib.c
*  Version: 1.6180
*  Date: Mar 12, 2010
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
*  This is the implementation file of a library
*  to work with binary search trees.
* -----------------------------------------------------------------
*/




#include <stdlib.h>
#include <gsl/gsl_blas.h>
#include <gsl/gsl_errno.h>

#include "../include/bst_lib.h"




/*
*------------------------------------------------------------
*   bst_leaf_alloc
*
*   This function allocates a struct leaf.
*
*   Input:
*   void
*
*   Output:
*   leaf - pointer to a struct leaf
*
*   last update: May 10, 2009
*------------------------------------------------------------
*/

bst_leaf *bst_leaf_alloc()
{
    bst_leaf *leaf = NULL;
    
    /* memory allocation for bst_leaf */
    leaf = (bst_leaf *) malloc(sizeof(bst_leaf));
    if ( leaf == NULL )
        return NULL;
    
    /* setting bst_leaf elements equal to NULL */
    leaf->phi  = NULL;
    leaf->Rphi = NULL;
    leaf->A    = NULL;
    leaf->L    = NULL;
    
    return leaf;
}
/*------------------------------------------------------------*/





/*
*------------------------------------------------------------
*   bst_node_alloc
*
*   This function allocates a struct node.
*
*   Input:
*   void
*
*   Output:
*   node - pointer to a struct node
*
*   last update: May 10, 2009
*------------------------------------------------------------
*/

bst_node* bst_node_alloc()
{
    bst_node *node = NULL;
    
    /* memory allocation for bst_node */
    node = (bst_node *) malloc(sizeof(bst_node));
    if ( node == NULL )
        return NULL;
    
    /* setting bst_node elements equal to NULL */
    node->v      = NULL;
    node->r_node = NULL;
    node->l_node = NULL;
    node->r_leaf = NULL;
    node->l_leaf = NULL;

    return node;
}
/*------------------------------------------------------------*/





/*
*------------------------------------------------------------
*   bst_leaf_free
*
*   This function frees the memory used by a struct leaf.
*
*   Input:
*   leaf - pointer to a struct leaf
*
*   Output:
*   void
*
*   last update: May 9, 2009
*------------------------------------------------------------
*/

void bst_leaf_free(void **leaf_bl)
{
    bst_leaf *leaf = NULL;
    
    /* checking if leaf_bl is NULL */
    if ( *leaf_bl == NULL )
        return;
    
    leaf = (bst_leaf *) (*leaf_bl);
    
    /* releasing bst_leaf elements */
    if( leaf->phi != NULL )
    {
        gsl_vector_free(leaf->phi);
        leaf->phi = NULL;
    }
    
    if( leaf->Rphi != NULL )
    {
        gsl_vector_free(leaf->Rphi);
        leaf->Rphi = NULL;
    }
    
    if( leaf->A != NULL )
    {
        gsl_matrix_free(leaf->A);
        leaf->A = NULL;
    }
    
    if( leaf->L != NULL )
    {
        gsl_matrix_free(leaf->L);
        leaf->L = NULL;
    }
    
    /* releasing allocated memory by bst_leaf */
    free(*leaf_bl);
    *leaf_bl = NULL;
    
    return;
}
/*------------------------------------------------------------*/





/*
*------------------------------------------------------------
*   bst_node_free
*
*   This function frees the memory used by a struct node.
*
*   Input:
*   node - pointer to a struct node
*
*   Output:
*   void
*
*   last update: Feb 19, 2010
*------------------------------------------------------------
*/

void bst_node_free(void **node_bl)
{
    bst_node *node = NULL;
    
    /* checking if node_bl is NULL */
    if ( *node_bl == NULL )
        return;
    
    node = (bst_node *) (*node_bl);
    
    /* releasing allocated memory by bst_node elements */
    if( node->l_leaf != NULL )
    {
        bst_leaf_free((void **)&(node->l_leaf));
        node->l_leaf = NULL;
    }
    
    if( node->r_leaf != NULL )
    {
        bst_leaf_free((void **)&(node->r_leaf));
        node->r_leaf = NULL;
    }
    
    if( node->l_node != NULL )
    {
        bst_node_free((void **)&(node->l_node));
        node->l_node = NULL;
    }
    
    if( node->r_node != NULL )
    {
        bst_node_free((void **)&(node->r_node));
        node->r_node = NULL;
    }
    
    if( node->v != NULL )
    {
        gsl_vector_free(node->v);
        node->v = NULL;
    }
    
    /* releasing allocated memory by bst_node */
    free(*node_bl);
    *node_bl = NULL;

    return;
}
/*------------------------------------------------------------*/





/*
*------------------------------------------------------------
*   bst_leaf_set
*
*   This function sets the bst_leaf elements equal to
*   the input elements and returns GSL_SUCCESS if
*   there is no problem during its execution.
*
*   Input:
*   phi  - composition
*   Rphi - reaction mapping
*   A    - mapping gradient matrix
*   L    - ellipsoid Cholesky matrix
*
*   Output:
*   leaf - pointer to a struct leaf
*   success or error
*
*   last update: Feb 22, 2009
*------------------------------------------------------------
*/

int bst_leaf_set(gsl_vector *phi,gsl_vector *Rphi,
                    gsl_matrix *A,gsl_matrix *L,bst_leaf *leaf)
{
    /* checking if leaf is NULL */
    if ( leaf == NULL )
        return GSL_EFAILED;
    
    /* setting phi equal to the input vector phi */
    if ( leaf->phi == NULL )
        leaf->phi = gsl_vector_calloc(phi->size);
    
    gsl_vector_memcpy(leaf->phi,phi);
    
    /* setting Rphi equal to the input vector Rphi */
    if ( leaf->Rphi == NULL )
        leaf->Rphi = gsl_vector_calloc(Rphi->size);
    
    gsl_vector_memcpy(leaf->Rphi,Rphi);
    
    /* setting A equal to the input matrix A */
    if ( leaf->A == NULL )
        leaf->A = gsl_matrix_calloc(A->size1,A->size2);
    
    gsl_matrix_memcpy(leaf->A,A);
    
    /* setting L equal to the input matrix L */
    if ( leaf->L == NULL )
        leaf->L = gsl_matrix_calloc(L->size1,L->size2);
    
    gsl_matrix_memcpy(leaf->L,L);
    
    return GSL_SUCCESS;
}
/*------------------------------------------------------------*/





/*
*------------------------------------------------------------
*   bst_cutplane
*
*   This function compute the vector ( v = phi_r - phi_l )
*   and the scalar ( a = (|phi_r|^2-|phi_l|^2)/2 )
*   which define a cutting plane used to search any
*   information in the binary search tree.
*   
*   Input:
*   phi_l - left  composition
*   phi_r - right compositon
*
*   Output:
*   v    - cutting plane normal vector
*   a    - cutting plane scalar
*
*   last update: Feb 27, 2009
*------------------------------------------------------------
*/

double bst_cutplane(gsl_vector *phi_l,gsl_vector *phi_r,gsl_vector *v)
{
    double a1, a2;

    /* v := phi_r */
    gsl_vector_memcpy(v,phi_r);
    
    /* v := phi_r - phi_l */
    gsl_vector_sub(v,phi_l);
    
    /* a1 := |phi_r|^2 */
    gsl_blas_ddot(phi_r,phi_r,&a1);

    /* a2 := |phi_l|^2 */
    gsl_blas_ddot(phi_l,phi_l,&a2);
    
    return 0.5*(a1-a2);
}
/*------------------------------------------------------------*/





/*
*------------------------------------------------------------
*   bst_node_set
*
*   This function set the bst_node elements equal to
*   the input elements and define the vector and
*   scalar to do the search for a information.
*
*   Input:
*   l_leaf - left  leaf
*   r_leaf - right leaf
*
*   Output:
*   node - pointer to a struct node
*   success or error
*
*   last update: Feb 27, 2009
*------------------------------------------------------------
*/

int bst_node_set(bst_leaf *l_leaf,bst_leaf *r_leaf,bst_node *node)
{
    /* checking if node is NULL */
    if ( node == NULL )
        return GSL_EFAILED;
    
    /* memory allocation for v */
    if ( node->v == NULL )
    	node->v = gsl_vector_calloc(r_leaf->phi->size);
    
    /* setting v and a */
    node->a = bst_cutplane(l_leaf->phi,r_leaf->phi,node->v);

    /* setting leaf elements equal to the input leaves */
    node->l_leaf = l_leaf;
    node->r_leaf = r_leaf;
    
    /* setting node elements equal to NULL */
    node->l_node = NULL;
    node->r_node = NULL;
    
    return GSL_SUCCESS;
}
/*------------------------------------------------------------*/





/*
*------------------------------------------------------------
*   bst_node_add
*
*   This function adds a new node to a given side
*   of a old node and returns GSL_SUCCESS if
*   there is no problem during its execution.
*
*   Input:
*   side     - new node side
*   old_node - pointer to the old node
*   new_node - pointer to the new node
*
*   Output:
*   success or error
*
*   last update: Feb 22, 2009
*------------------------------------------------------------
*/

int bst_node_add(int side,bst_node *old_node,bst_node *new_node)
{
    if( side == RIGHT )
    {
        if( old_node->r_node != NULL )
            return GSL_EFAILED;
        else
        {
            /* adding a new right node */
            old_node->r_node = new_node;
            old_node->r_leaf = NULL;
            
            return GSL_SUCCESS;
        }
    }
    else if( side == LEFT )
    {
        if( old_node->l_node != NULL )
            return GSL_EFAILED;
        else
        {
            /* adding a new left node */
            old_node->l_node = new_node;
            old_node->l_leaf = NULL;
            
            return GSL_SUCCESS;
        }
    }
    else
        return GSL_EFAILED;
}
/*------------------------------------------------------------*/





/*
*------------------------------------------------------------
*   bst_search
*
*   This function searchs for a near composition phi0 in the
*   binary search tree.
*
*   Input:
*   root     - binary tree root
*   phi      - composition
*
*   Output:
*   end_leaf - leaf with the near composition
*   end_node - node with the near composition
*   side     - near compositon leaf side
*
*   last update: Feb 3, 2009
*------------------------------------------------------------
*/

int bst_search(bst_node *root,gsl_vector *phi,
                        bst_node **end_node,bst_leaf **end_leaf)
{
    double dot;
    
    if( root == NULL )
    {
        *end_leaf = NULL;
        *end_node = NULL;
        
        return GSL_EFAILED;
    }
    
    /* dot := v^T*phi */
    gsl_blas_ddot(root->v,phi,&dot);
    
    /* if dot > a then right node else left node */
    if( dot > root->a )
    {
        /* extern node */
        if( root->r_node == NULL )
        {
            *end_leaf = root->r_leaf;
            *end_node = root;
            
            return RIGHT;
        }
        else
            return bst_search(root->r_node,phi,end_node,end_leaf);
    }
    else
    {
        /* extern node */
        if( root->l_node == NULL )
        {
            *end_leaf = root->l_leaf;
            *end_node = root;
            
            return LEFT;
        }
        else
            return bst_search(root->l_node,phi,end_node,end_leaf);
    }
}
/*------------------------------------------------------------*/




/*
*------------------------------------------------------------
*   bst_height
*
*   This function computes the binary search tree height.
*
*   Input:
*   root - binary tree root
*
*   Output:
*   height - binary search tree height
*
*   last update: Nov 2, 2009
*------------------------------------------------------------
*/

int bst_height(bst_node *root)
{
    int l_tree_height = 0;
    int r_tree_height = 0;
    
    if ( root != NULL )
    {
        l_tree_height = bst_height(root->l_node);
        r_tree_height = bst_height(root->r_node);
        
        if ( l_tree_height > r_tree_height )
            return l_tree_height + 1;
        else
            return r_tree_height + 1;
    }
    else
        return 0;
    
}
/*------------------------------------------------------------*/

