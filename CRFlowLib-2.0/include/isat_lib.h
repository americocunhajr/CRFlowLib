
/*
* -----------------------------------------------------------------
*  isat_lib.h
*  In Situ Adaptive Tabulation Library
*  Version: 2.0
*  Last Update: Nov 3, 2019
*  
*  Programmer: Americo Barbosa da Cunha Junior
*              americo.cunhajr@gmail.com
* ----------------------------------------------------------------- 
*  Copyright (c) 2010-2019, Americo Barbosa da Cunha Junior
*  All rights reserved.
* -----------------------------------------------------------------
*  This is the header file for ISAT_LIB module, a computational 
*  library with In Situ Adaptive Tabulation (ISAT) algorithm 
*  routines.
* -----------------------------------------------------------------
*/



#ifndef __ISAT_LIB_H__
#define __ISAT_LIB_H__


#include <gsl/gsl_blas.h>
#include <gsl/gsl_math.h>

#include "../include/bst_lib.h"




/*
*------------------------------------------------------------
*   structure for ISAT workspace
*
*   last update: Oct 9, 2019
*------------------------------------------------------------
*/

typedef struct isat_struct
{
    bst_node *root;         /* binary search tree root       */
    unsigned int lf;        /* # of leaves                   */
    unsigned int nd;        /* # of nodes                    */
    unsigned int add;       /* # of adds                     */
    unsigned int grw;       /* # of grows                    */
    unsigned int rtv;       /* # of retrieves                */
    unsigned int dev;       /* # of direct evaluations       */
    unsigned int hgt;       /* tree height                   */
    unsigned int maxleaves; /* max. # of leaves in the tree  */
    clock_t time_add;       /* CPU clock spent by additions  */
    clock_t time_grw;       /* CPU clock spent by growths    */
    clock_t time_rtv;       /* CPU clock spent by retrieves  */
    clock_t time_dev;       /* CPU clock spent by dir. eval. */
    
} isat_wrk;
/*------------------------------------------------------------*/




/*
*------------------------------------------------------------
* function prototypes
*------------------------------------------------------------
*/

isat_wrk *isat_alloc();

void isat_free(void **isat_bl);

int isat_bst_init(isat_wrk *isat);

void isat_statistics(isat_wrk *isat);

int isat_input(unsigned int *max_lf,
               double *etol);

void isat_eoa_mtrx(gsl_matrix *A,
                   double etol,
                   gsl_matrix *L);

double isat_lerror(gsl_vector *Rphi,
                   gsl_vector *Rphi0,
		           gsl_matrix *A,
                   gsl_vector *phi,
                   gsl_vector *phi0);

int isat4(isat_wrk *isat_mem,
          void *cvode_mem,
	      double etol,
          double t0,
          double delta_t,
          gsl_vector *phi,
          gsl_matrix *A,
          gsl_matrix *L,
          gsl_vector *Rphi);


#endif /* __ISAT_LIB_H__ */
