
/*
* -----------------------------------------------------------------
*  isat_lib.c
*  In Situ Adaptive Tabulation Library
*  Version: 2.0
*  Last Update: Nov 3, 2019
*  
*  Programmer: Americo Barbosa da Cunha Junior
*              americo.cunhajr@gmail.com
* -----------------------------------------------------------------
*  Copyright (c) 2010-2019, Americo Barbosa da Cunha Junior
*  All rights reserved.
* -----------------------------------------------------------------
*  This is the implementation file for ISAT_LIB module, a 
*  computational library with In Situ Adaptive Tabulation (ISAT) 
*  algorithm routines.
* -----------------------------------------------------------------
*/




#include <stdlib.h>
#include <math.h>
#include <time.h>
#include <gsl/gsl_math.h>
#include <gsl/gsl_errno.h>

#include "../include/thrm_lib.h"
#include "../include/ell_lib.h"
#include "../include/ode_lib.h"
#include "../include/isat_lib.h"




/*
*------------------------------------------------------------
*   isat_alloc
*
*   This function alocates memory for an ISAT workspace
*   structure and initialize its elements.
*
*   Output:
*   isat_mem - pointer to ISAT workspace
*
*   last update: Oct 9, 2019
*------------------------------------------------------------
*/

isat_wrk *isat_alloc()
{
    /* create ISAT workspace */
    isat_wrk *isat_mem = NULL;

    /* memory allocation for ISAT workspace */
    isat_mem = (isat_wrk *) malloc(sizeof(isat_wrk));
    if ( isat_mem == NULL )
        return NULL;
    
    /* initialize ISAT workspace elements */
    isat_mem->root      = NULL;
    isat_mem->lf        = 0;
    isat_mem->nd        = 0;
    isat_mem->add       = 0;
    isat_mem->grw       = 0;
    isat_mem->rtv       = 0;
    isat_mem->dev       = 0;
    isat_mem->hgt       = 0;
    isat_mem->maxleaves = 0;
    isat_mem->time_add  = 0.0;
    isat_mem->time_grw  = 0.0;
    isat_mem->time_rtv  = 0.0;
    isat_mem->time_dev  = 0.0;
    
    return isat_mem;
}
/*------------------------------------------------------------*/




/*
*------------------------------------------------------------
*   isat_free
*
*   This function release the memory used by ISAT workspace.
*
*   Input:
*   isat_bl - pointer to ISAT workspace
*
*   last update: Oct 9, 2019
*------------------------------------------------------------
*/

void isat_free(void **isat_bl)
{
    /* create ISAT workspace */
    isat_wrk *isat_mem = NULL;
    
    /* check if ISAT workspace memory block is NULL */
    if ( *isat_bl == NULL )
        return;
    
    isat_mem = (isat_wrk *) (*isat_bl);
    
    /* release the memory allocated for ISAT workspace elements */
    if( isat_mem->root != NULL )
    {
        bst_node_free((void **)&(isat_mem->root));
        isat_mem->root = NULL;
    }
    
    /* release the memory allocated for ISAT workspace */
    free(*isat_bl);
    *isat_bl = NULL;
    
    return;
}
/*------------------------------------------------------------*/




/*
*------------------------------------------------------------
*   isat_bst_init
*
*   This function initiates ISAT binary search tree.
*
*   Input:
*   isat_mem - pointer to ISAT workspace
*
*   last update: Oct 9, 2019
*------------------------------------------------------------
*/

int isat_bst_init(isat_wrk *isat_mem)
{
    /* check if ISAT workspece is allocated */
    if ( isat_mem == NULL )
        return GSL_EINVAL;
    
    /* memory allocation for ISAT binary search tree root */
    if ( isat_mem->root == NULL )
    {
        isat_mem->root = bst_node_alloc();
        if ( isat_mem->root == NULL )
        {
            free(isat_mem);
            isat_mem = NULL;
            return GSL_EINVAL;
        }
    }
    
    return GSL_SUCCESS;
}
/*------------------------------------------------------------*/




/*
*------------------------------------------------------------
*   isat_statistics
*
*   This function prints on screen ISAT workspace
*   statistics of usage.
*
*   Input:
*   isat_mem - pointer to ISAT workspace
*
*   last update: Oct 1, 2019
*------------------------------------------------------------
*/

void isat_statistics(isat_wrk *isat_mem)
{
    double time_add;
    double time_grw;
    double time_rtv;
    double time_dev;
    
    time_add = (double) (isat_mem->time_add / CLOCKS_PER_SEC) / isat_mem->add;
    time_grw = (double) (isat_mem->time_grw / CLOCKS_PER_SEC) / isat_mem->grw;
    time_rtv = (double) (isat_mem->time_rtv / CLOCKS_PER_SEC) / isat_mem->rtv;
    time_dev = (double) (isat_mem->time_dev / CLOCKS_PER_SEC) / isat_mem->dev;
    
    printf("\n ISAT statistics:");
    printf("\n # of adds         = %d",   isat_mem->add);
    printf("\n # of grows        = %d",   isat_mem->grw);
    printf("\n # of retrieves    = %d",   isat_mem->rtv);
    printf("\n # of dir. eval.   = %d",   isat_mem->dev);
    printf("\n # of leaves       = %d",   isat_mem->lf);
    printf("\n # of nodes        = %d",   isat_mem->nd);
    printf("\n tree height       = %d\n", isat_mem->hgt);
    
    printf("\n average values for CPU time (s):");
    printf("\n add: %+.6e",time_add);
    printf("\n grw: %+.6e",time_grw);
    printf("\n rtv: %+.6e",time_rtv);
    printf("\n dev: %+.6e",time_dev);
    
    return;
}
/*------------------------------------------------------------*/




/*
* -----------------------------------------------------------------
*   isat_input
*
*   This function receives ISAT input parameters from the user.
*
*   Input:
*   maxleaves - maximum number of ISAT tree leaves
*   etol      - ISAT error tolerance
*
*   Output:
*   success or error
*
*   last update: Nov 3, 2019
* -----------------------------------------------------------------
*/

int isat_input(unsigned int *maxleaves,double *etol)
{
    printf("\n Input ISAT parameters:\n");
    
    printf("\n maximum of tree leaves:");
    scanf("%d", maxleaves);
    printf("\n %d\n", *maxleaves);
    if( *maxleaves <= 0 )
	GSL_ERROR(" maxleaves must be a positive integer",GSL_EINVAL);

    printf("\n ISAT error tolerance:");
    scanf("%lf", etol);
    printf("\n %+.1e\n", *etol);
    if( *etol < 0.0 )
	GSL_ERROR(" etol must be grather than zero",GSL_EINVAL);
    
    return GSL_SUCCESS;
}
/*----------------------------------------------------------------*/




/*
*------------------------------------------------------------
*   isat_eoa_matrix
*
*   This function computes the EOA matrix in Cholesky form.
*   
*   Input:
*   A    - gradient matrix
*   etol - error tolerance
*
*   Output:
*   L    - EOA Cholesky matrix
*
*   last update: Oct 9, 2019
*------------------------------------------------------------
*/

void isat_eoa_mtrx(gsl_matrix *A,
                   double etol,
                   gsl_matrix *L)
{
    unsigned int i;
    double  eps_max   = 0.5;
    gsl_matrix *Aetol = NULL;
    gsl_matrix *V     = NULL;
    gsl_vector *sigma = NULL;

    /* memory allocation */
    Aetol = gsl_matrix_calloc(A->size2,A->size2);
    V     = gsl_matrix_calloc(A->size2,A->size2);
    sigma = gsl_vector_calloc(A->size2);
    
    /* Aetol := A */
    gsl_matrix_memcpy(Aetol,A);
    
    /* Aetol := (1/etol).A */
    gsl_matrix_scale (Aetol,1.0/etol);

    /* Aetol = U*sigma*V^T */
    ell_psd2eig(Aetol,V,sigma);

    /* eliminate small and large singular values */
    for ( i = 0; i < sigma->size; i++ )
        sigma->data[i] = GSL_MAX(sigma->data[i],eps_max);
    
    /* V*sigma^2*V^T = L*L^T */
    ell_eig2chol(V,sigma,L);
    
    /* release allocated memory */
    gsl_vector_free(sigma);
    gsl_matrix_free(V);
    gsl_matrix_free(Aetol);
    sigma = NULL;
    V     = NULL;
    Aetol = NULL;
    
    return;
}
/*------------------------------------------------------------*/




/*
*------------------------------------------------------------
*   isat_lerror
*
*   This function computes the local error defined as
*
*   eps = ||Rl(phi)- R(phi)||_2 where
*
*   Rl(phi) = R(phi0) + A*(phi-phi0).
*
*   Input:
*   Rphi  - reaction mapping of phi
*   Rphi0 - reaction mapping of phi0
*   A     - mapping gradient matrix
*   phi   - query composition
*   phi0  - initial composition
*
*   Output:
*   eps   - local error
*
*   last update: Feb 19, 2010
*------------------------------------------------------------
*/

double isat_lerror(gsl_vector *Rphi,
                   gsl_vector *Rphi0,
                   gsl_matrix *A,
                   gsl_vector *phi,
                   gsl_vector *phi0)
{
    double eps;
    gsl_vector *Rlphi = NULL;
    
    /* memory allocation for Rlphi */
    Rlphi = gsl_vector_calloc(phi->size);
    
    /* Rlphi := R(phi0) + A*(phi-phi0) */
    linear_approx(phi,phi0,Rphi0,A,Rlphi);
    
    /* Rlphi := Rl(phi) - R(phi) */
    gsl_vector_sub(Rlphi,Rphi);
    
    /* eps := 2-norm(Rl(phi) - R(phi)) */
    eps = gsl_blas_dnrm2(Rlphi);
    
    /* releasing allocated memory */
    gsl_vector_free(Rlphi);
    Rlphi = NULL;
    
    return eps;
}
/*------------------------------------------------------------*/




/*
*------------------------------------------------------------
*   isat4
*
*   This function executes the 4-th version of the
*   In Situ Adaptive Tabulation (ISAT) algorithm.
*
*   
*   Input:
*   isat_mem  - ISAT workspace
*   cvode_mem - ODE solver workspace
*   etol      - error tolerance
*   t0        - initial time
*   delta_t   - time step
*   phi       - query composition
*   A         - mapping gradient matrix
*   L         - EOA Cholesky matrix
*
*   Output:
*   Rphi      - reaction mapping
*   success or error
*
*   last update: Nov 3, 2019
*------------------------------------------------------------
*/

int isat4(isat_wrk *isat_mem,
          void *cvode_mem,
          double etol,
          double t0,
          double delta_t,
          gsl_vector *phi,
          gsl_matrix *A,
          gsl_matrix *L,
          gsl_vector *Rphi)
{

    /* CPU clock start counter */
    clock_t cpu_start = clock();

    int flag;
    unsigned int bst_side;
    bst_leaf *end_leaf     = NULL;
    bst_node *end_node     = NULL;
    
    /****** First call for ISAT algorithm ******/

    /* check if there is no leaf in the binary search tree */
    if( isat_mem->lf == 0 )
    {
        int flag;
        bst_leaf *first_leaf = NULL;
	
	/* memory allocation for BST first leaf */
        first_leaf = bst_leaf_alloc();
	if ( first_leaf == NULL )
                return GSL_ENOMEM;
        
        /* direct integration */
        flag = odesolver_reinit(t0,phi,cvode_mem);
        if ( flag != GSL_SUCCESS )
            return flag;
        
        flag = odesolver(cvode_mem,
                         delta_t,
                         Rphi);
        if ( flag != GSL_SUCCESS )
            return flag;
        
        /* compute the mapping gradient matrix */
        flag = gradient(cvode_mem,t0,delta_t,phi,Rphi,A);
        if ( flag != GSL_SUCCESS )
            return flag;
        
        /* compute EOA Cholesky matrix */
        isat_eoa_mtrx(A,etol,L);
        
        /* define first leaf elements */
        bst_leaf_set(phi,Rphi,A,L,first_leaf);
        isat_mem->root->r_leaf = first_leaf;
        
        /* update ISAT workspace counters */
        isat_mem->lf++;
        isat_mem->hgt = bst_height(isat_mem->root);
        
        return GSL_SUCCESS;
    }
    
    
    /****** Further calls for ISAT algorithm ******/  
    
    /* search for the near composition in binary search tree */
    if( isat_mem->lf > 1 )
	bst_side = bst_search(isat_mem->root,
                              phi,
                              &end_node,
                              &end_leaf);
    else
    {
        end_node = isat_mem->root;
        end_leaf = isat_mem->root->r_leaf;
        bst_side = BST_RIGHT;
    }
    
    /* check if the near composition is inside the ellipsoid */
    flag = ell_pt_in(phi,end_leaf->phi,end_leaf->L);
    if( flag == ELL_TRUE )
    {
        /* compute the linear approximation */
        linear_approx(phi,end_leaf->phi,end_leaf->Rphi,end_leaf->A,Rphi);
        
        /* update ISAT workspace counters */
        isat_mem->rtv++;
	    isat_mem->time_rtv += clock() - cpu_start;
        
        return GSL_SUCCESS;
    }
    else
    {
        double lerror = 0.0;
	
        /* performe direct integration */
        flag = odesolver_reinit(t0,phi,cvode_mem);
        if ( flag != GSL_SUCCESS )
            return flag;
        
        flag = odesolver(cvode_mem,delta_t,Rphi);
        if ( flag != GSL_SUCCESS )
            return flag;
        
        /* compute ISAT local error */
        lerror = isat_lerror(Rphi,end_leaf->Rphi,
                             end_leaf->A,phi,end_leaf->phi);
        
        /* check if the local error is greater than etol */
        if( lerror < etol )
        {
            /* grow the ellipsoid */
            ell_pt_modify(phi,end_leaf->phi,end_leaf->L);
            
            /* update ISAT workspace counters */
            isat_mem->grw++;
	        isat_mem->time_grw += clock() - cpu_start;
            
            return GSL_SUCCESS;
        }
        else
        {
	    /* check if the maximum number of leaves is excced */
	    if( isat_mem->lf > isat_mem->maxleaves )
	    {
	        /*  update ISAT workspace counters */
	        isat_mem->dev++;
		    isat_mem->time_dev +=  clock() - cpu_start;
		    return GSL_SUCCESS;
	    }
	    
            bst_leaf *new_leaf = NULL;
            
	    /* memory allocation */
            new_leaf = bst_leaf_alloc();
            if ( new_leaf == NULL )
                return GSL_ENOMEM;
            
            /* compute the mapping gradient matrix */
            flag = gradient(cvode_mem,t0,delta_t,phi,Rphi,A);
            if ( flag != GSL_SUCCESS )
                return flag;
            
            /* compute ellipsoid matrix */
            isat_eoa_mtrx(A,etol,L);
            
            /* define the new leaf elements */
            bst_leaf_set(phi,Rphi,A,L,new_leaf);
            
            /* check if the binary search tree has more than one leaf */
            if( isat_mem->lf > 1 )
            {
                bst_node *new_node = NULL;
		
		        /* memory allocation for the new node */
                new_node = bst_node_alloc();
                if ( new_node == NULL )
                    return GSL_ENOMEM;
                
		        /* define the new node leaves */
                bst_node_set(end_leaf,new_leaf,new_node);
		
		        /* add the new node to the binary search tree */
                bst_node_add(bst_side,end_node,new_node);
            }
	        else
                bst_node_set(end_leaf,new_leaf,end_node);
            
            /* update ISAT workspace counters */
            isat_mem->add++;
            isat_mem->lf++;
            isat_mem->nd++;
            isat_mem->hgt = bst_height(isat_mem->root);
	        isat_mem->time_add += clock() - cpu_start;
            
            return GSL_SUCCESS;
        }
    }
}
/*------------------------------------------------------------*/
