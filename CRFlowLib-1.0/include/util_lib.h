
/*
* -----------------------------------------------------------------
*  Utilities Library --- util_lib.h
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
*  This is the header file for a library with miscelaneous
*  functions witch have a lot of applications.
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
