
/*
* -----------------------------------------------------------------
*  Ellipsoid Library --- ell_lib.h
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
*  This is the header file for a library that operates an
*  n-dimensional ellipsoid.
* -----------------------------------------------------------------
*/




#ifndef __ELL_LIB_H__
#define __ELL_LIB_H__


#include <gsl/gsl_vector.h>
#include <gsl/gsl_matrix.h>


#undef  TOL
#define TOL 1.0e-7

#undef  TRUE
#define TRUE 1

#undef  FALSE
#define FALSE 0




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

double ell_pt_in(gsl_vector *x,
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
