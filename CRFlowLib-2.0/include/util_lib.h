
/*
* -----------------------------------------------------------------
*  util_lib.h
*  Utilities Library
*  Version: 2.0
*  Last Update: Nov 3, 2019
*
*  Programmer: Americo Barbosa da Cunha Junior
*              americo.cunhajr@gmail.com
* -----------------------------------------------------------------
*  Copyright (c) 2010-2019, Americo Barbosa da Cunha Junior
*  All rights reserved.
* -----------------------------------------------------------------
*  This is the header file for UTIL_LIB module, a computational
*  library with routines for miscelaneous applications.
* -----------------------------------------------------------------
*/




#ifndef __UTIL_LIB_H__
#define __UTIL_LIB_H__

#include <time.h>

#include <gsl/gsl_vector.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_rng.h>




/*
*------------------------------------------------------------
* function prototypes
*------------------------------------------------------------
*/

int gsl_fprint_vector(FILE* file,
                      gsl_vector *v);

int gsl_print_vector(gsl_vector *v);

int gsl_print_matrix_line(unsigned int i,
                          gsl_matrix *M);

int gsl_print_matrix_column(unsigned int j,
                            gsl_matrix *M);

int gsl_print_matrix(gsl_matrix *M);

int gsl_rand_vector(gsl_rng *r,
                    gsl_vector *v);

int gsl_rand_matrix(gsl_rng *r,
                    gsl_matrix *M);

int gsl_diagonal_matrix(gsl_vector *v,
                        gsl_matrix *D);

void util_time_used(clock_t cpu_start,
                    clock_t cpu_end,
                    time_t wall_start,
                    time_t wall_end);

#endif /* __UTIL_LIB_H__ */
