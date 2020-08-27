
/*
* -----------------------------------------------------------------
*  thrm_lib.h
*  Thermochemistry Library
*  Version: 2.0
*  Last Update: Nov 3, 2019
*
*  Programmer: Americo Barbosa da Cunha Junior
*              americo.cunhajr@gmail.com
* -----------------------------------------------------------------
*  Copyright (c) 2010-2019, Americo Barbosa da Cunha Junior
*  All rights reserved.
* -----------------------------------------------------------------
*  This is the header file for THRM_LIB module, a computational
*  library with thermochemistry routines.
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
