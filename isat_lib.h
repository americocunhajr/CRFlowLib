
/*
* -----------------------------------------------------------------
*  In Situ Adaptive Tabulation Library --- isat_lib.h
*  Version: 1.6180
*  Date: Feb 25, 2010
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
*  with the ISAT algorithm.
* -----------------------------------------------------------------
*/



#ifndef __ISAT_LIB_H__
#define __ISAT_LIB_H__


#include <gsl/gsl_blas.h>
#include <gsl/gsl_math.h>

#include "../include/bst_lib.h"


#undef  RTOL
#define RTOL 1.0e-6

#undef  ATOL
#define ATOL 1.0e-9




/*
*------------------------------------------------------------
*   structure of ISAT workspace
*
*   last update: Jul 13, 2010
*------------------------------------------------------------
*/

typedef struct isat_struct
{
    bst_node *root;         /* binary search tree root      */
    unsigned int lf;        /* # of leaves                  */
    unsigned int nd;        /* # of nodes                   */
    unsigned int add;       /* # of adds                    */
    unsigned int grw;       /* # of grows                   */
    unsigned int rtv;       /* # of retrieves               */
    unsigned int dev;       /* # of discarted evaluations   */
    unsigned int hgt;       /* tree height                  */
    unsigned int max_lf;    /* maximum value of tree leaves */
    clock_t time_add;       /* cpu clock spent by additions */
    clock_t time_grw;       /* cpu clock spent by growths   */
    clock_t time_rtv;       /* cpu clock spent by retrieves */
    clock_t time_dev;       /* cpu clock spent by dev       */
    
} isat_wrk;
/*------------------------------------------------------------*/




/*
*------------------------------------------------------------
* function prototypes
*------------------------------------------------------------
*/

isat_wrk *isat_alloc();

void isat_free(void **isat_bl);

int isat_set(isat_wrk *isat);

void isat_statistics(isat_wrk *isat);

int isat_input(unsigned int *max_lf,
               double *etol,
               double *n0);

void isat_eoa_mtrx(gsl_matrix *A,
                   double etol,
                   double n0,
                   gsl_matrix *L);

double isat_lerror(gsl_vector *Rphi,
                   gsl_vector *Rphi0,
		           gsl_matrix *A,
                   gsl_vector *phi,
                   gsl_vector *phi0);

int isat4(isat_wrk *isat,
          void *thrm_data,
          void *cvode_mem,
          double etol,
          double n0,
          double t0,
          double delta_t,
          gsl_vector *phi,
          gsl_matrix *A,
          gsl_matrix *L,
          gsl_vector *Rphi);


#endif /* __ISAT_LIB_H__ */
