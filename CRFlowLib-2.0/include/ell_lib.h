
/*
* -----------------------------------------------------------------
*  ell_lib.h
*  Ellipsoid Library
*  Version: 2.0
*  Last Update: July 22, 2020
*  
*  Programmer: Americo Barbosa da Cunha Junior
*              americo.cunhajr@gmail.com
* -----------------------------------------------------------------
*  Copyright (c) 2010-2019, Americo Barbosa da Cunha Junior
*  All rights reserved.
* -----------------------------------------------------------------
*  This is the header file for BST_LIB module, a computational
*  library to work with n-dimensional ellipsoids.
* -----------------------------------------------------------------
*/



#ifndef __ELL_LIB_H__
#define __ELL_LIB_H__

#include <stdlib.h>
#include <gsl/gsl_linalg.h>


#undef  ELL_TOL
#define ELL_TOL 1.0e-7

#undef  ELL_TRUE
#define ELL_TRUE 1

#undef  ELL_FALSE
#define ELL_FALSE 0



/*
*------------------------------------------------------------
* function prototypes
*------------------------------------------------------------
*/

void ell_psd2eig(gsl_matrix *B,
                 gsl_matrix *V,
                 gsl_vector *lam);

void ell_psd2chol(gsl_matrix *B,
                  gsl_matrix *L);

void ell_eig2chol(gsl_matrix *V,
                  gsl_vector *sig,
                  gsl_matrix *L);

unsigned int ell_pt_in(gsl_vector *x,
                       gsl_vector *c,
                       gsl_matrix *L);

void ell_map2uhs(gsl_vector *x,
                 gsl_vector *c,
                 gsl_matrix *L,
                 gsl_vector *y);

void ell_pt_modify(gsl_vector *p,
                   gsl_vector *c,
                   gsl_matrix *L);


#endif /* __ELL_LIB_H__ */
