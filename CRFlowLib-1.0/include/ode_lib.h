
/*
* -----------------------------------------------------------------
*  ODE Solver Library --- ode_lib.h
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
*  This is the header file for a library
*  with the ODE solution tools.
* -----------------------------------------------------------------
*/



#ifndef __ODESOLVER_LIB_H__
#define __ODESOLVER_LIB_H__


#include <gsl/gsl_vector.h>
#include <gsl/gsl_matrix.h>
#include <cvode/cvode.h>


#undef  RTOL
#define RTOL 1.0e-6

#undef  ATOL
#define ATOL 1.0e-15




/*
*------------------------------------------------------------
* function prototypes
*------------------------------------------------------------
*/

double uround(void);

double wnorm(int n,
             double *v1,
             double *v2);

void ewtset(int n,
            double *v,
            double atol,
            double rtol,
            double *ewt);

int jacobian(CVRhsFn f,
             void *f_data,
             gsl_vector *Fx,
             gsl_vector *x,
             double t,
             double atol,
             double rtol,
             gsl_matrix *J);

int gradient(CVRhsFn f,
             void *f_data,
             void *cvode_mem,
             double t0,
             double delta_t,
             double atol,
             double rtol,
             gsl_vector *phi,
             gsl_vector *Rphi,
             gsl_matrix *A);

void linear_approx(gsl_vector *x,
                   gsl_vector *x0,
                   gsl_vector *Fx0,
                   gsl_matrix *DFx0,
                   gsl_vector *Fx);

int odesolver_init(CVRhsFn f,
                   void *data,
                   double t0,
                   int mxsteps,
		           double atol,
                   double rtol,
                   gsl_vector *x,
                   void *cvode_mem);

int odesolver_reinit(CVRhsFn f,
                     void *f_data,
                     double t0,
		             double atol,
                     double rtol,
                     gsl_vector *x,
                     void *cvode_mem);

int odesolver(void *cvode_mem,
              double tf,
              gsl_vector *Fx);


#endif /* __ODESOLVER_LIB_H__ */
