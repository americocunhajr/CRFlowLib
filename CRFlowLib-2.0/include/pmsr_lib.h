
/*
* -----------------------------------------------------------------
*  pmsr_lib.h
*  Pairwise Mixing Stirred Reactor Library
*  Version: 2.0
*  Last Update: Oct 1, 2019
*  
*  This is the header file of a computational library with
*  routines to implement a Pairwise Mixing Stirred Reactor
*  (PMSR) model.
* ----------------------------------------------------------------- 
*  Programmer: Americo Barbosa da Cunha Junior
*              americo.cunhajr@gmail.com
* -----------------------------------------------------------------
*  Copyright (c) 2019 by Americo Barbosa da Cunha Junior
* -----------------------------------------------------------------
*/



#ifndef __PMSR_LIB_H__
#define __PMSR_LIB_H__


#include <gsl/gsl_vector.h>
#include <gsl/gsl_rng.h>




/*
* ------------------------------------------------------------
*  Types: struct pmsr, pmsr_wrk
* ------------------------------------------------------------
* This structure contains fields of PMSR model.
* 
* last update: Nov 5, 2009
* ------------------------------------------------------------
*/

typedef struct pmsr
{
    int *pmsr_idx;         /* particles index vector */
    gsl_vector **pmsr_ps;  /* particle system        */
} pmsr_wrk;
/*------------------------------------------------------------*/





/*
*------------------------------------------------------------
* function prototypes
*------------------------------------------------------------
*/

void pmsr_title();

int pmsr_input(unsigned long int *seed,
               unsigned int *Np,
               unsigned int *Ndt,
	       double *t0,
               double *delta_t,
               double *tau_res,
               double *tau_mix,
	       double *tau_pair,
               double *Tmax,
               double *Tmin,
               double *atol,
               double *rtol);

pmsr_wrk *pmsr_alloc();

int pmsr_init(int Np,
              int Neq,
              pmsr_wrk *pmsr);

void pmsr_free(int Np,
               void **pmsr_bl);

void pmsr_set_all(int Np,
                  gsl_vector *phi,
                  pmsr_wrk *pmsr);

int pmsr_Nexc(int Np,
              double delta_t,
              double tau_res);

int pmsr_Npair(int Np,
               double delta_t,
               double tau_pair);

void pmsr_meanvar(int Np,
                  int k,
                  gsl_vector **ps,
                  double *mean,
                  double *var);

void pmsr_mixture(int Np,
                  int Neq,
                  double delta_t,
                  double tau_mix,
                  int *idx,
                  gsl_vector **ps);

void pmsr_iops(int Np,
               int Nexc,
               int Npair,
               gsl_rng *r,
               gsl_vector *phi,
               pmsr_wrk *pmsr);


#endif /* __PMSR_LIB_H__ */
