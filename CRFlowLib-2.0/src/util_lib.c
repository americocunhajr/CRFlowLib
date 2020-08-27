
/*
* -----------------------------------------------------------------
*  Utilities Library --- util_lib.c
*  Version: 1.6180
*  Date: Mar 12, 2010
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
*  This is the implementation file of a library
*  with miscelaneous functions.
* -----------------------------------------------------------------
*/




#include <stdio.h>
#include <gsl/gsl_errno.h>

#include "../include/util_lib.h"




/*
*------------------------------------------------------------
*   gsl_fprint_vector
*
*   This function prints a vector in a file.
*
*   Input:
*   file - output file
*   v    - GSL vector
*
*   Output:
*
*   Return:
*   success or error
*
*   last update: Nov 6, 2008
*------------------------------------------------------------
*/

int gsl_fprint_vector(FILE* file,gsl_vector *v)
{
    if( v == NULL )
	return GSL_EINVAL;
    
    unsigned int i;
    
    for( i = 0; i < v->size; i++ )
	fprintf(file," %+.3e", v->data[i]);

    return GSL_SUCCESS;
}
/*------------------------------------------------------------*/




/*
*------------------------------------------------------------
*   gsl_print_vector
*
*   This function prints a vector on the screen
*
*   Input:
*   v  - GSL vector
*
*   Output:
*
*   Return:
*   success or error
*
*   last update: Dec 26, 2008
*------------------------------------------------------------
*/

int gsl_print_vector(gsl_vector *v)
{
    if( v == NULL )
	return GSL_EINVAL;
    
    unsigned int i;
    
    for( i = 0; i < v->size; i++ )
	printf(" %+.3e", v->data[i]);

    return GSL_SUCCESS;
}
/*------------------------------------------------------------*/




/*
*------------------------------------------------------------
*   gsl_print_matrix_line
*
*   This function prints a matrix line on the screen.
*
*   Input:
*   i  - matrix line to be printed
*   M  - GSL matrix
*
*   Output:
*
*   Return:
*   success or error
*
*   last update: Dec 26, 2008
*------------------------------------------------------------
*/

int gsl_print_matrix_line(unsigned int i,gsl_matrix *M)
{
    if( M == NULL )
	return GSL_EINVAL;
    
    unsigned int j;
    
    printf("\n");
    for( j = 0; j < M->size2; j++ )
	printf(" %+.3e ", M->data[i*M->tda+j]);
    
    return GSL_SUCCESS;
}
/*------------------------------------------------------------*/




/*
*------------------------------------------------------------
*   gsl_print_matrix_column
*
*   This function prints a matrix column on the screen
*
*   Input:
*   M  - GSL matrix
*   i  - matrix column to be printed
*
*   Output:
*
*   Return:
*   success or error
*
*   last update: Dec 26, 2008
*------------------------------------------------------------
*/

int gsl_print_matrix_column(unsigned int j,gsl_matrix *M)
{
    if( M == NULL )
    	return GSL_EINVAL;
    
    unsigned int i;
    
    printf("\n");
    for( i = 0; i < M->size1; i++ )
	printf(" %+.3e\n", M->data[i*M->tda+j]);
    
    return GSL_SUCCESS;
}
/*------------------------------------------------------------*/




/*
*------------------------------------------------------------
*   gsl_print_matrix
*
*   This function prints a matrix on the screen
*
*   Input:
*   M  - GSL matrix
*
*   Outpur:
*
*   Return:
*   success or error
*
*   last update: Dec 26, 2008
*------------------------------------------------------------
*/

int gsl_print_matrix(gsl_matrix *M)
{
    if( M == NULL )
    	return GSL_EINVAL;
    
    unsigned int i;
    
    for( i = 0; i < M->size1; i++ )
	gsl_print_matrix_line(i,M);
    
    return GSL_SUCCESS;
}
/*------------------------------------------------------------*/




/*
*------------------------------------------------------------
*   gsl_rand_vector
*
*   This function creates a random vector.
*
*   Input:
*   seed - random number generator seed
*
*   Output:
*   v  - GSL vector
*
*   Return:
*   success or error
*
*   last update: Feb 19, 2010
*------------------------------------------------------------
*/

int gsl_rand_vector(gsl_rng *r,gsl_vector *v)
{
    if( v == NULL )
    	return GSL_EINVAL;
    
    unsigned int i;
    
    for( i = 0; i < v->size; i++ )
	v->data[i] = gsl_rng_uniform(r) + 1.0;
    
    return GSL_SUCCESS;
}
/*------------------------------------------------------------*/




/*
*------------------------------------------------------------
*   gsl_rand_matrix
*
*   This function creates a radom matrix.
*
*   Input:
*   seed - random number generator seed
*
*   Output:
*   M    - GSL matrix
*
*   Return:
*   success or error
*
*   last update: Feb 27, 2009
*------------------------------------------------------------
*/

int gsl_rand_matrix(gsl_rng *r,gsl_matrix *M)
{
    if( M == NULL )
    	return GSL_EINVAL;
    
    unsigned int i, j;
    
    for( i = 0; i < M->size1; i++ )
        for( j = 0; j < M->size2; j++ )
            M->data[i*M->tda+j] = gsl_rng_uniform(r) + 1.0;
    
    return GSL_SUCCESS;
}
/*------------------------------------------------------------*/




/*
*------------------------------------------------------------
*   gsl_diagonal_matrix
*
*   This function creates a diagonal matrix with the elements
*   of the vector v.
*
*   Input:
*   v - GSL vector
*
*   Output:
*   D  - GSL matrix
*
*   Return:
*   success or error
*
*   last update: Feb 26, 2009
*------------------------------------------------------------
*/

int gsl_diagonal_matrix(gsl_vector *v,gsl_matrix *D)
{
    if( v == NULL || D == NULL )
    	return GSL_EINVAL;
    
    unsigned int i, j;
    
    for( i = 0; i < D->size1; i++ )
        for( j = 0; j < D->size2; j++ )
            if( i == j )
                D->data[i*D->tda+j] = v->data[i];
            else
                D->data[i*D->tda+j] = 0.0;
    
    return GSL_SUCCESS;
}
/*------------------------------------------------------------*/




/*
*------------------------------------------------------------
*   util_time_used
*
*   This function computes the CPU time and Wall time
*   used in a simulation.
*
*   Input:
*   cpu_start  - cpu clock start flag
*   cpu_end    - cpu clock end flag
*   wall_start - wall clock start flag
*   wall_end   - wall clock end flag
*
*   last update: Nov 10, 2009
*------------------------------------------------------------
*/

void util_time_used(clock_t cpu_start,clock_t cpu_end,
			    time_t wall_start,time_t wall_end)
{
    double cpu_used;
    double wall_used;
    
    cpu_used  = ((double) (cpu_end - cpu_start)) / CLOCKS_PER_SEC;
    wall_used = difftime(wall_end,wall_start);
    
    printf("\n");
    printf("\n CPU  time (s): %+.6e",cpu_used);
    printf("\n Wall time (s): %+.6e",wall_used);
    printf("\n");
    
    return;
}
/*------------------------------------------------------------*/
