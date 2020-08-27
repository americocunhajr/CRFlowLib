
/*
* -----------------------------------------------------------------
*  ode_lib.h
*  ODE Solver Library
*  Version: 2.0
*  Last Update: Oct 3, 2019
*  
*  This is the header file for computational library with
*  routines to deal with Ordinary Differential Equations (ODE).
* ----------------------------------------------------------------- 
*  Programmer: Americo Barbosa da Cunha Junior
*              americo.cunhajr@gmail.com
* -----------------------------------------------------------------
*  Copyright (c) 2019 by Americo Barbosa da Cunha Junior
* -----------------------------------------------------------------
*/



#ifndef __ODESOLVER_LIB_H__
#define __ODESOLVER_LIB_H__


#include <gsl/gsl_vector.h>
#include <gsl/gsl_matrix.h>
#include <cvode/cvode.h>




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

int gradient(void *cvode_mem,
             double t0,
             double delta_t,
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

int odesolver_reinit(double t0,
                     gsl_vector *x,
                     void *cvode_mem);

int odesolver(void *cvode_mem,
              double tf,
              gsl_vector *Fx);


#endif /* __ODESOLVER_LIB_H__ */
