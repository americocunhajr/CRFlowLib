
/*
* -----------------------------------------------------------------
*  Thermochemistry Library --- thrm_lib.h
*  Version: 1.6180
*  Date: Dec 29, 2010
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
*  This is the header file for a library with
*  thermochemistry routines.
* -----------------------------------------------------------------
*/



#ifndef __THRM_LIB_H__
#define __THRM_LIB_H__


#include <gsl/gsl_vector.h>
#include <nvector/nvector_serial.h>




/*
* ------------------------------------------------------------
*  Types: struct thrm_struct, thrm_wrk
* ------------------------------------------------------------
* This structure contains fields of the Chemkin workspace.
* 
* last update: Mar 5, 2009
* ------------------------------------------------------------
*/

typedef struct thrm_struct
{
    int    *iwk;   /* integer work vector */
    double *rwk;   /* real work vector */
    char   *cwk;   /* character work vector */
    int leniwk;    /* iwk dimension */
    int lenrwk;    /* rwk dimension */
    int lencwk;    /* cwk dimension */
    int n_e;       /* # of chemical elements */
    int n_s;       /* # of elementary species */
    int n_r;       /* # of elementary reactions */
} thrm_wrk;
/*------------------------------------------------------------*/




/*
*------------------------------------------------------------
* function prototypes
*------------------------------------------------------------
*/

void thrm_title();

void thrm_mech(int n_e,
               int n_s,
               int n_r,
               char *mech);

int thrm_input(int *steps,
               double *t0,
               double *delta_t,
               double *Tmax,
               double *Tmin,
               double *tol);

thrm_wrk* thrm_alloc();

void thrm_free(void **thrm_bl);

int thrm_init(thrm_wrk *thrm);

void thrm_composition(thrm_wrk *thrm,
                      gsl_vector *phi);

void thrm_page_heading0(FILE *file);

double thrm_h2T(thrm_wrk *thrm,
                double *phi,
                double Tmax,
                double Tmin,
                double tol);

void thrm_temp_meanvar(int Np,
                       gsl_vector **ps,
                       thrm_wrk *thrm,
                       double Tmax,
                       double Tmin,
                       double tol,
                       double *mean,
                       double *var);

double thrm_rrsum(gsl_vector *phi);

double thrm_mfsum(gsl_vector *phi);

int thrm_mfsign(gsl_vector *phi);

int thrm_eqs(double t,
             N_Vector phi,
             N_Vector Rphi,
             void *thrm_data);

int thrm_eqs2(double t,
              N_Vector phi,
              N_Vector Sphi,
              void *thrm_data);

#endif /* __THRM_LIB_H__ */
