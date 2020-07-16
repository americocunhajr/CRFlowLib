
/*
* -----------------------------------------------------------------
*  Ellipsoid Library --- ell_lib.c
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
*  to work with n-dimensional ellipsoids.
* -----------------------------------------------------------------
*/




#include <stdlib.h>
#include <gsl/gsl_linalg.h>
#include <gsl/gsl_blas.h>

#include "../include/ell_lib.h"




/*
*------------------------------------------------------------
*   ell_psd2eig
*
*  The n x n positive symetric defined matrix A is given
*  by A = B^T*B. This function computes the eigendecomposition
*  of A given the matrix B.
*
*  Method:
*          The SVD of B is: B = U*S*V^T, so that
*          B^T*B = V*S^2*V^T.
*
*   Input:
*   B   - ellipsoid matrix
*
*   Output:
*   V   - orthogonal matrix
*   sig - diagonal elments of S
*
*   last update: Feb 26, 2009
*------------------------------------------------------------
*/

void ell_psd2eig(gsl_matrix *B,gsl_matrix *V,gsl_vector *sig)
{    
    gsl_matrix *A    = NULL;
    gsl_vector *work = NULL;
    
    /* memory allocation */
    A    = gsl_matrix_calloc(B->size1,B->size2);
    work = gsl_vector_calloc(sig->size);
    
    /* A := B */
    gsl_matrix_memcpy(A,B);
    
    /* B = U*S*V^T */
    gsl_linalg_SV_decomp(A,V,sig,work);

    /* releasing allocated memory */
    gsl_vector_free(work);
    gsl_matrix_free(A);
    A    = NULL;
    work = NULL;
    
    return;
}
/*------------------------------------------------------------*/




/*
*------------------------------------------------------------
*   ell_psd2chol
*
*  The n x n positive symetric defined matrix A is given
*  by A = B^T*B. This function computes the Cholesky
*  decomposition of A given the matrix B.
*  
*  Method:
*          The matrix B has QR decomposition Q*R,
*          where Q is orthogonal and R is upper-triangular,
*          so that the Cholesky decomposition of A is given by
*          L*L^T where L := R^T.
*
*   Input:
*   B - ellipsoid matrix
*
*   Output:
*   L - lower triangular matrix
*
*   last update: Feb 22, 2009
*------------------------------------------------------------
*/

void ell_psd2chol(gsl_matrix *B,gsl_matrix *L)
{
    /* L := B^T*B */
    gsl_blas_dgemm(CblasTrans,CblasNoTrans,1.0,B,B,0.0,L);
    
    /* B^T*B = L*L^T */
    gsl_linalg_cholesky_decomp(L);
    
    return;
}
/*------------------------------------------------------------*/




/*
*------------------------------------------------------------
*   ell_eig2chol
*
*  The n x n matrix A has the eigendecomposition
*  A = V*S^2*V^T and Cholesky decomposition A = L*L^T.
*  This function computes the Cholesky decomposition of A
*  given V and sig = diag(S).
*
*  Method:
*          The eigendecomposition of is:
*          A = V*S^2*V^T = B^T*B where B = S*V, so that
*          B^T*B = L*L^T where L is a Cholesky matrix.
*
*   Input:
*   V   - orthogonal matrix
*   sig - singular values vector
*
*   Output:
*   L - lower triangular matrix
*
*   last update: Mar 1, 2010
*------------------------------------------------------------
*/

void ell_eig2chol(gsl_matrix *V,gsl_vector *sig,gsl_matrix *L)
{
    unsigned int i;
    gsl_matrix *S = NULL;
    gsl_matrix *B = NULL;
    
    /* memory allocation */
    S = gsl_matrix_calloc(sig->size,sig->size);
    B = gsl_matrix_calloc(L->size1,L->size2);
    
    /* S := diag(S) */
    for ( i = 0; i < sig->size; i++ )
        S->data[i*S->tda+i] = sig->data[i];
    
    /* B := S*V */
    gsl_blas_dgemm(CblasNoTrans,CblasNoTrans,1.0,S,V,0.0,B);
    
    /* B^T*B = L*L^T */
    ell_psd2chol(B,L);
    
    /* releasing allocated memory */
    gsl_matrix_free(S);
    gsl_matrix_free(B);
    S = NULL;
    B = NULL;
    
    return;
}
/*------------------------------------------------------------*/



/*
*------------------------------------------------------------
*   ell_pt_in
*  
*  This function checks if a vector x is covered by an
*  n-ellipsoid centered at the vector c and returns
*  TRUE if the answer is afirmative or FALSE otherwise.
*
*  2-norm(L^T*(x-c)) <= 1.0
*
*   Input:
*   x - point to be tested
*   c - ellipsoid center
*   L - ellipsoid Cholesky matrix
*
*   Output:
*   dot - 2-norm of ellipsoid equation
*
*   last update: Oct 6, 2009
*------------------------------------------------------------
*/

double ell_pt_in(gsl_vector *x,gsl_vector *c,gsl_matrix *L)
{
    double dot;
    gsl_vector *aux = NULL;
    
    /* memory allocation */
    aux = gsl_vector_calloc(x->size);
    
    /* aux := x */
    gsl_vector_memcpy(aux,x);

    /* aux := x - c */
    gsl_vector_sub(aux,c);
    
    /* aux := L^T*(x-c) */
    gsl_blas_dtrmv(CblasLower,CblasTrans,CblasNonUnit,L,aux);

    /* dot := 2-norm(L^T*(x-c)) */
    dot = gsl_blas_dnrm2(aux);
    
    /* releasing allocated memory */
    gsl_vector_free(aux);
    aux = NULL;

    return dot;
}
/*------------------------------------------------------------*/




/*
*------------------------------------------------------------
*   ell_map2uhs
*
*   This function maps a vector x into other vector y in the
*   transformed space where the n-ellipsoid is a unity
*   n-sphere centered at origin.
*   
*   y = L^T*(x-c)
*   
*   Input:
*   x - vector to be mapped
*   c - ellipsoid center
*   L - ellipsoid Cholesky matrix
*
*   Output:
*   y - mapped vector
*
*   last update: Feb 22, 2009
*------------------------------------------------------------
*/

void ell_map2uhs(gsl_vector *x,gsl_vector *c,gsl_matrix *L,gsl_vector *y)
{
    /* y := x */
    gsl_vector_memcpy(y,x);
    
    /* y := x - c */
    gsl_vector_sub(y,c);
    
    /* y := L^T*(x-c) */
    gsl_blas_dtrmv(CblasLower,CblasTrans,CblasNonUnit,L,y);

    return;
}
/*------------------------------------------------------------*/




/*
*------------------------------------------------------------
*   ell_pt_modify
*
*  This function changes an n-ellipsoid centered at point c
*  to cover a given point p and the old n-ellipsoid.
*  Given p, c and the ellipsoid Cholesky matrix L,
*  this function chnages L such that it becomes to be
*  associated with the new n-ellipsoid.
*
*   Input:
*   p - point to be covered
*   c - ellipsoid center
*
*   Output:
*   L - new ellipsoid Cholesky matrix
*
*   last update: Feb 23, 2009
*------------------------------------------------------------
*/

void ell_pt_modify(gsl_vector *p,gsl_vector *c,gsl_matrix *L)
{
    double ypnorm, gamma;
    gsl_vector *yp = NULL;
    gsl_matrix *G  = NULL;
    
    /* memory allocation */
    yp = gsl_vector_calloc(p->size);
    G  = gsl_matrix_calloc(p->size,p->size);
    
    /* mapping p to the y-space, the unit n-sphere space */
    ell_map2uhs(p,c,L,yp);
    ypnorm = gsl_blas_dnrm2(yp);
    
    /* checking the degenerate case where p = c */
    if( ypnorm < TOL )
    {
        gsl_matrix_set_identity(L);
        
        /* releasing allocated memory */
        gsl_matrix_free(G);
        gsl_vector_free(yp);
        
        return;
    }
    
    /* rank-one modification algorithm */
    gamma = (1.0/ypnorm - 1.0)/(ypnorm*ypnorm);
    
    /* G := I */
    gsl_matrix_set_identity(G);
    
    /* G := I + gamma.yp*yp^T */
    gsl_blas_dsyr(CblasUpper,gamma,yp,G);
    
    /* G := G*L^T */
    gsl_blas_dtrmm(CblasRight,CblasLower,CblasTrans,CblasNonUnit,1.0,L,G);
    
    /* L := L*G^T*G*L^T = L*L^T */
    ell_psd2chol(G,L);
    
    /* releasing allocated memory */
    gsl_matrix_free(G);
    gsl_vector_free(yp);
    G  = NULL;
    yp = NULL;
    
    return;
}
/*------------------------------------------------------------*/

