
/*
* -----------------------------------------------------------------
*  Pope case equivalent composition --- main__pope-case.c
*  Version: 1.6180
*  Date: Oct 13, 2010
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
*  This is the implementation file of a pairwise mixing
*  stirred reactor using direct integration.
* -----------------------------------------------------------------
*/




#include <math.h>
#include <time.h>

#include "../include/util_lib.h"
#include "../include/thrm_lib.h"



/*
* -----------------------------------------------------------------
*   main
* -----------------------------------------------------------------
*/

int main(void)
{
    unsigned int flag;       /* flag              */
    unsigned int Neq;        /* # of equations    */
    double T;                /* temperature                  */
    double Tmax;             /* temperature upper bound      */
    double Tmin;             /* temperature lower bound      */
    double tol;              /* tolerance                    */
    gsl_vector *phi;         /* composition       */
    thrm_wrk *thrm;          /* chemkin workspace */
    
    /* initiating variables */
    flag       = GSL_SUCCESS;
    Neq        = 0;
    T          = 0.0;
    Tmax       = 5000.0;
    Tmin       = 100.0;
    tol        = 1.0e-4;
    phi        = NULL;
    thrm       = NULL;

    
    /* memory allocation for chemkin workspace */
    thrm = thrm_alloc();
    if ( thrm == NULL )
        GSL_ERROR(" can't alloc thrm_wrk ",flag);    
    
    /* initiating chemkin workspace */
    flag = thrm_init(thrm);
    if ( flag != GSL_SUCCESS )
        GSL_ERROR(" can't initiate thrm_wrk ",flag);
    
    /* setting the number of equations */
    Neq = thrm->n_s + 2;

    /* memory allocation for composition vectors */
    phi = gsl_vector_calloc(Neq);
    
    /* inputing inflow parameters */
    printf("\n input composition:\n");
    thrm_composition(thrm,phi);
    
    phi->data[0] = -2.5459e9;
    
    T = thrm_h2T(thrm,phi->data,Tmax,Tmin,tol);
    
    printf("\n mixture temperature: %f \n",T);
    
    /* releasing allocated memory */
    gsl_vector_free(phi);
    thrm_free((void **)&thrm);
    
    phi  = NULL;
    thrm = NULL;

    return GSL_SUCCESS;
}
