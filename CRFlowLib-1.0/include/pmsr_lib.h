
/*
* -----------------------------------------------------------------
*  Pairwise Mixing Stirred Reaction Library --- pmsr_lib.h
*  Version: 1.6180
*  Date: Nov 16, 2010
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
*  This is the header file for the PMSR model.
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
               unsigned int *K,
               double *t0,
               double *delta_t,
               double *tau_res,
               double *tau_mix,
	           double *tau_pair,
               double *Tmax,
               double *Tmin);

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
