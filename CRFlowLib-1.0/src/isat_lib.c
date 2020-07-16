
/*
* -----------------------------------------------------------------
*  In Situ Adaptive Tabulation Library --- isat_lib.c
*  Version: 1.6180
*  Date: Nov 15, 2010
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
*  to work with ISAT algorithm.
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
*   This function alocates memory for a struct isat.
*
*   Output:
*   isat - pointer to a struct isat
*
*   last update: Jul 13, 2010
*------------------------------------------------------------
*/

isat_wrk *isat_alloc()
{
    /* memory allocation for isat_wrk */
    isat_wrk *isat = NULL;
    isat = (isat_wrk *) malloc(sizeof(isat_wrk));
    if ( isat == NULL )
        return NULL;
    
    /* setting isat_wrk elements equal NULL and 0.0 */
    isat->root     = NULL;
    isat->lf       = 0;
    isat->nd       = 0;
    isat->add      = 0;
    isat->grw      = 0;
    isat->rtv      = 0;
    isat->dev      = 0;
    isat->hgt      = 0;
    isat->max_lf   = 0;
    isat->time_add = 0.0;
    isat->time_grw = 0.0;
    isat->time_rtv = 0.0;
    isat->time_dev = 0.0;
    
    return isat;
}
/*------------------------------------------------------------*/




/*
*------------------------------------------------------------
*   isat_free
*
*   This function frees the memory used by a struct isat.
*
*   Input:
*   isat - pointer to a struct isat
*
*   last update: May 10, 2009
*------------------------------------------------------------
*/

void isat_free(void **isat_bl)
{
    isat_wrk *isat = NULL;
    
    /* checking if isat_bl is NULL */
    if ( *isat_bl == NULL )
        return;
    
    isat = (isat_wrk *) (*isat_bl);
    
    /* releasing allocated memory by isat_wrk elements */
    if( isat->root != NULL )
    {
        bst_node_free((void **)&(isat->root));
        isat->root = NULL;
    }
    
    /* releasing allocated memory by isat_wrk */
    free(*isat_bl);
    *isat_bl = NULL;
    
    return;
}
/*------------------------------------------------------------*/




/*
*------------------------------------------------------------
*   isat_set
*
*   This function initiates ISAT binary search tree.
*
*   Input:
*   isat - pointer to a struct isat
*
*   last update: May 13, 2009
*------------------------------------------------------------
*/

int isat_set(isat_wrk *isat)
{
    /* checking if isat is NULL */
    if ( isat == NULL )
        return GSL_EINVAL;
    
    /* memory allocation for binary search tree root */
    if ( isat->root == NULL )
    {
        isat->root = bst_node_alloc();
        if ( isat->root == NULL )
        {
            free(isat);
            isat = NULL;
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
*   This function prints on the screen informations
*   about the isat workspace.
*
*   Input:
*   isat - pointer to a struct isat
*
*   last update: Jul 13, 2010
*------------------------------------------------------------
*/

void isat_statistics(isat_wrk *isat)
{
    double time_add;
    double time_grw;
    double time_rtv;
    double time_dev;
    
    time_add = (double) (isat->time_add / CLOCKS_PER_SEC) / isat->add;
    time_grw = (double) (isat->time_grw / CLOCKS_PER_SEC) / isat->grw;
    time_rtv = (double) (isat->time_rtv / CLOCKS_PER_SEC) / isat->rtv;
    time_dev = (double) (isat->time_dev / CLOCKS_PER_SEC) / isat->dev;
    
    printf("\n ISAT statistics:");
    printf("\n # of adds         = %d",   isat->add);
    printf("\n # of grows        = %d",   isat->grw);
    printf("\n # of retrieves    = %d",   isat->rtv);
    printf("\n # of dir. eval.   = %d",   isat->dev);
    printf("\n # of leaves       = %d",   isat->lf);
    printf("\n # of nodes        = %d",   isat->nd);
    printf("\n tree height       = %d\n", isat->hgt);
    
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
*   This function receives ISAT parameters.
*
*   Input:
*   max_lf - maximum value of tree leaves
*   etol   - ISAT error tolerance
*   n0     - factor to multiply by unit roundoff
*
*   Output:
*   success or error
*
*   last update: Feb 25, 2010
* -----------------------------------------------------------------
*/

int isat_input(unsigned int *max_lf,double *etol,double *n0)
{
    printf("\n Input ISAT parameters:\n");
    
    printf("\n maximum of tree leaves:");
    scanf("%d", max_lf);
    printf("\n %d\n", *max_lf);
    if( *max_lf <= 0 )
	GSL_ERROR(" max_lf must be a positive integer (max_lf > 0)",GSL_EINVAL);

    printf("\n ISAT error tolerance:");
    scanf("%lf", etol);
    printf("\n %+.1e\n", *etol);
    if( *etol < 0.0 )
	GSL_ERROR(" etol must be grather than zero (etol > 0.0)",GSL_EINVAL);
    
    printf("\n Unit roundoff multiple:");
    scanf("%lf", n0);
    printf("\n %+.1e\n", *n0);
    if( *n0 < 0.0 )
	GSL_ERROR(" n0 must be grather than zero (n0 > 0.0)",GSL_EINVAL);
    
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
*   n0   - factor to multiply by unit roundoff
*   etol - error tolerance
*
*   Output:
*   L    - EOA Cholesky matrix
*
*   last update: Feb 19, 2010
*------------------------------------------------------------
*/

void isat_eoa_mtrx(gsl_matrix *A,double etol,double n0,
						gsl_matrix *L)
{
    unsigned int i;
    double  eps_max   = 0.5;
    double  eps_min   = etol/(n0*DBL_EPSILON);
    gsl_matrix *Aetol = NULL;
    gsl_matrix *V     = NULL;
    gsl_vector *sig   = NULL;

    /* memory allocation */
    Aetol = gsl_matrix_calloc(A->size2,A->size2);
    V     = gsl_matrix_calloc(A->size2,A->size2);
    sig   = gsl_vector_calloc(A->size2);
    
    /* Aetol := A */
    gsl_matrix_memcpy(Aetol,A);
    
    /* Aetol := (1/etol).A */
    gsl_matrix_scale (Aetol,1.0/etol);

    /* Aetol = U*sig*V^T */
    ell_psd2eig(Aetol,V,sig);

    /* eliminating small and large singular values */
    for ( i = 0; i < sig->size; i++ )
        sig->data[i] = GSL_MIN( GSL_MAX(sig->data[i],eps_max), eps_min);
    
    /* V*sig^2*V^T = L*L^T */
    ell_eig2chol(V,sig,L);
    
    /* releasing allocated memory */
    gsl_vector_free(sig);
    gsl_matrix_free(V);
    gsl_matrix_free(Aetol);
    sig   = NULL;
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
*   eps = 2-norm(Rl(phi)- R(phi)) where
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

double isat_lerror(gsl_vector *Rphi,gsl_vector *Rphi0,
		gsl_matrix *A,gsl_vector *phi,gsl_vector *phi0)
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
*   in situ adaptive tabulation algorithm.
*
*   
*   Input:
*   isat      - isat workspace
*   thrm_data - thermochemistry workspace
*   etol      - error tolerance
*   n0        - factor to multiply by unit roundoff
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
*   last update: Nov 15, 2010
*------------------------------------------------------------
*/

int isat4(isat_wrk *isat,void *thrm_data,void *cvode_mem,
	    double etol,double n0,double t0,double delta_t,
    gsl_vector *phi,gsl_matrix *A,gsl_matrix *L,gsl_vector *Rphi)
{
    clock_t cpu_start = clock();
    
    /* ISAT first step */
    if( isat->lf == 0 )
    {
        int flag;
        bst_leaf *first_leaf = NULL;
	
	/* memory allocation for first_leaf */
        first_leaf = bst_leaf_alloc();
	if ( first_leaf == NULL )
                return GSL_ENOMEM;
        
        /* performing direct integration */
        flag = odesolver_reinit(thrm_eqs,thrm_data,t0,
                                        ATOL,RTOL,phi,cvode_mem);
        if ( flag != GSL_SUCCESS )
            return flag;
        
        flag = odesolver(cvode_mem,delta_t,Rphi);
        if ( flag != GSL_SUCCESS )
            return flag;
        
        /* computing the mapping gradient matrix */
        flag = gradient(thrm_eqs,thrm_data,cvode_mem,t0,delta_t,
                                                ATOL,RTOL,phi,Rphi,A);
        if ( flag != GSL_SUCCESS )
            return flag;
        
        /* computing EOA Cholesky matrix */
        isat_eoa_mtrx(A,etol,n0,L);
        
        /* setting first leaf elements */
        bst_leaf_set(phi,Rphi,A,L,first_leaf);
        isat->root->r_leaf = first_leaf;
        
        /* updating leaves and height counters */
        isat->lf++;
        isat->hgt = bst_height(isat->root);
        
        return GSL_SUCCESS;
    }
    
    
    /* ISAT second and following steps */
    int side, flag;
    double dot         = 0.0;
    bst_leaf *end_leaf = NULL;
    bst_node *end_node = NULL;
    
    
    /* searching for the near composition in binary search tree */
    if( isat->lf > 1 )
	side = bst_search(isat->root,phi,&end_node,&end_leaf);
    else
    {
        end_node = isat->root;
        end_leaf = isat->root->r_leaf;
        side     = RIGHT;
    }
    
    /* checking if the near composition is inside EOA */
    dot = ell_pt_in(phi,end_leaf->phi,end_leaf->L);
    if( gsl_fcmp(dot,1.0,ATOL) < 1 )
    {
        /* computing linear approximation */
        linear_approx(phi,end_leaf->phi,end_leaf->Rphi,end_leaf->A,Rphi);
        
        /* updating counters */
        isat->rtv++;
	isat->time_rtv += clock() - cpu_start;
        
        return GSL_SUCCESS;
    }
    else
    {
        double lerror = 0.0;
	
        /* performing direct integration */
        flag = odesolver_reinit(thrm_eqs,thrm_data,t0,
                                        ATOL,RTOL,phi,cvode_mem);
        if ( flag != GSL_SUCCESS )
            return flag;
        
        flag = odesolver(cvode_mem,delta_t,Rphi);
        if ( flag != GSL_SUCCESS )
            return flag;
        
        /* computing ISAT local error */
        lerror = isat_lerror(Rphi,end_leaf->Rphi,
                            end_leaf->A,phi,end_leaf->phi);
        
        /* checking if lerror is greater than etol */
        if( gsl_fcmp(lerror,etol,ATOL) < 1 )
        {
            /* growing the EOA */
            ell_pt_modify(phi,end_leaf->phi,end_leaf->L);
            
            /* updating counters */
            isat->grw++;
	    isat->time_grw += clock() - cpu_start;
            
            return GSL_SUCCESS;
        }
        else
        {
	    /* checking if the maximum # of leaves was excced */
	    if( isat->lf > isat->max_lf )
	    {
	        /*  updating counter */
	        isat->dev++;
		isat->time_dev +=  clock() - cpu_start;
		
		return GSL_SUCCESS;
	    }
	    
            bst_leaf *new_leaf = NULL;
            
	    /* memory allocation */
            new_leaf = bst_leaf_alloc();
            if ( new_leaf == NULL )
                return GSL_ENOMEM;
            
            /* computing the mapping gradient matrix */
            flag = gradient(thrm_eqs,thrm_data,cvode_mem,t0,delta_t,
                                                ATOL,RTOL,phi,Rphi,A);
            if ( flag != GSL_SUCCESS )
		return flag;
            
            /* computing EOA Cholesky matrix */
            isat_eoa_mtrx(A,etol,n0,L);
            
            /* setting new_leaf elements */
            bst_leaf_set(phi,Rphi,A,L,new_leaf);
            
            /* cheking if binary search tree has more than one leaf */
            if( isat->lf > 1 )
            {
                bst_node *new_node = NULL;
		
		/* memory allocation for new node */
                new_node = bst_node_alloc();
                if ( new_node == NULL )
                    return GSL_ENOMEM;
                
		/* setting the leaves of the new node */
                bst_node_set(end_leaf,new_leaf,new_node);
		
		/* adding the new node to the tree */
                bst_node_add(side,end_node,new_node);
            }
	    else
                bst_node_set(end_leaf,new_leaf,end_node);
            
            /* updating counters */
            isat->add++;
            isat->lf++;
            isat->nd++;
            isat->hgt = bst_height(isat->root);
	    isat->time_add += clock() - cpu_start;
            
            return GSL_SUCCESS;
        }
    }
}
/*------------------------------------------------------------*/
