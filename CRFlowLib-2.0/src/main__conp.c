
/*
* -----------------------------------------------------------------
*  Adiabatic Constant Pressure Problem --- main_conp.c
*  Version: 3.14
*  Date: Dec 28, 2010
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
#include "../include/cklib.h"
#include "../include/ckthrm.h"
#include "../include/thrm_lib.h"
#include "../include/ode_lib.h"





/*
* -----------------------------------------------------------------
*   pmsr_title
*
*   This function prints the program title on the screen.
*
*   last update: Dec 29, 2010
* -----------------------------------------------------------------
*/

void conp_title()
{
    time_t curtime = time (NULL);
    struct tm *loctime = localtime (&curtime);
    
    printf("\n");
    printf("==================================================");
    printf("\n Adiabatic Constant Pressure Problem --- conp \n");
    printf("\n Author:");
    printf("\n Americo Barbosa da Cunha Junior --- PUC-Rio");
    printf("\n americo.cunhajr@gmail.com\n");
    printf("\n Date:");
    printf("\n %s", asctime(loctime) );
    printf("==================================================");
    printf("\n");
    
    return;
}
/*----------------------------------------------------------------*/




/*
* -----------------------------------------------------------------
*   conp_input
*
*   This function receives conp parameters
*   and returns GSL_SUCCESS if there is no error.
*
*   Input:
*   K            - # of time steps
*   initial time - time step
*   delta_t      - time step
*
*   Output:
*   success or error
*
*   last update: Dec 29, 2010
* -----------------------------------------------------------------
*/

int conp_input(unsigned int *K,double *t0,double *delta_t)
{
    printf("\n Input CONP parameters:\n");
    
    printf("\n # of time steps:");
    scanf("%d", K);
    printf("\n %d\n", *K);
    if( *K <= 0 )
	GSL_ERROR(" K must be a positive integer (K > 0)",GSL_EINVAL);

    printf("\n initial time:");
    scanf("%lf", t0);
    printf("\n %+.1e\n", *t0);
    if( *t0 < 0.0 )
	GSL_ERROR(" t0 must be a grather than zero (t0 > 0.0)",GSL_EINVAL);
    
    printf("\n time step:");
    scanf("%lf", delta_t);
    printf("\n %+.1e\n", *delta_t);
    if( *delta_t <= 0.0 )
	GSL_ERROR(" delta_t must be a grather than zero (delta_t > 0.0)",GSL_EINVAL);
    
    return GSL_SUCCESS;
}
/*----------------------------------------------------------------*/



/*
* -----------------------------------------------------------------
*   main
*
*   This is the main function of conp program.
* -----------------------------------------------------------------
*/

int main(void)
{
    unsigned int k;          /* events counter               */
    unsigned int flag;       /* flag                         */
    unsigned int Neq;        /* # of equations               */
    unsigned int K;          /* # of time steps              */
    unsigned int mxsteps;    /* max of solver iterations     */
    clock_t cpu_start;       /* cpu clock start flag         */
    clock_t cpu_end;         /* cpu clock end flag           */
    time_t wall_start;       /* wall clock start flag        */
    time_t wall_end;         /* wall clock end flag          */
    double atol;             /* absolute tolerance           */
    double rtol;             /* relative tolerance           */
    double t0;               /* initial time                 */
    double delta_t;          /* time step                    */
    double T;                /* temperature                  */
    double Tmax;             /* temperature upper bound      */
    double Tmin;             /* temperature lower bound      */
    double tol;              /* tolerance                    */
    gsl_vector *phi;         /* state vector                 */
    thrm_wrk *thrm;          /* chemkin workspace            */
    void *cvode_mem;         /* CVODE workspace              */
    FILE *file1;
    
    
    /* initiating variables */
    k          = 0;
    flag       = GSL_SUCCESS;
    Neq        = 0;
    K          = 0;
    mxsteps    = 5000;
    atol       = ATOL;
    rtol       = RTOL;
    t0         = 0.0;
    delta_t    = 0.0;
    T          = 0.0;
    Tmax       = 5000.0;
    Tmin       = 100.0;
    tol        = 1.0e-4;
    phi        = NULL;
    thrm       = NULL;
    cvode_mem  = NULL;
    file1      = NULL;
    
    
    /* CPU and Wall clocks flags in the program begining */
    cpu_start  = clock();
    wall_start = time(NULL);
    
    /* program title */
    conp_title();
    
    /* memory allocation for chemkin workspace */
    thrm = thrm_alloc();
    if ( thrm == NULL )
        GSL_ERROR(" can't alloc thrm_wrk ",flag);
    
    /* memory allocation for ODE solver workspace */
    cvode_mem = CVodeCreate(CV_BDF,CV_NEWTON);
    if ( cvode_mem == NULL )
        GSL_ERROR(" can't alloc cvode_mem ",flag);
    
    /* inputing conp parameters */
    conp_input(&K,&t0,&delta_t);
    
    /* initiating chemkin workspace */
    flag = thrm_init(thrm);
    if ( flag != GSL_SUCCESS )
        GSL_ERROR(" can't initiate thrm_wrk ",flag);
    
    /* setting the number of equations */
    Neq = thrm->n_s + 2;

    /* memory allocation for state vector */
    phi = gsl_vector_calloc(Neq);
    
    /* printing mechanism information */
    thrm_mech(thrm->n_e,thrm->n_s,thrm->n_r,"---");
    
    /* inputing reagents parameters */
    printf("\n Input reagents parameters:\n");
    thrm_composition(thrm,phi);
        
    /* initiating ODE solver workspace */
    flag = odesolver_init(thrm_eqs,thrm,t0,mxsteps,
                            atol,rtol,phi,cvode_mem);
    if ( flag != GSL_SUCCESS )
        GSL_ERROR(" can't initiate cvode_mem ",flag);
    
    /* opening file to write the properties evolution */
    file1 = fopen("conp.dat","w");
    if ( file1 == NULL )
    {
        printf("can't open an output file (*.dat)\n");
        exit(1);
    }
    
    /* printing initial time */
    fprintf(file1,"\n %+.3e",t0);
    
    /* printing initial properties */
    T = thrm_h2T(thrm,phi->data,Tmax,Tmin,tol);
    fprintf(file1," %+.3e",T);
    gsl_fprint_vector(file1,phi);
    
    /* solving governing equations for K events */
    while ( k < K && flag == GSL_SUCCESS )
    {
        fprintf(file1,"\n %+.3e",(k+1)*delta_t);
        
        /* ode integrattion */
        t0 += delta_t;
        gsl_vector_set_zero(phi);
        flag = odesolver(cvode_mem,t0,phi);
        T = thrm_h2T(thrm,phi->data,Tmax,Tmin,tol);
        
        /* printing reactive flow properties */
        fprintf(file1," %+.3e",T);
        gsl_fprint_vector(file1,phi);
        
        /* updating events counter */
        k++;
    }
    
    
    /* closing file1 */
    fprintf(file1,"\n");
    fclose(file1);
    
    /* CPU and Wall clocks at program end */
    cpu_end  = clock();
    wall_end = time(NULL);
    
    /* printing CPU and Wall time spent at program execution */
    util_time_used(cpu_start,cpu_end,wall_start,wall_end);
    
    /* checking for errors in simalation */
    if ( flag != GSL_SUCCESS )
        printf("\n\n Simualtion finish with error\n");
    else
        printf("\n\n Simulation finish without error\n");
    
    /* releasing allocated memory */
    gsl_vector_free(phi);
    thrm_free((void **)&thrm);
    CVodeFree(&cvode_mem);
    
    phi        = NULL;
    thrm       = NULL;
    cvode_mem  = NULL;
    file1      = NULL;
    
    return GSL_SUCCESS;
}
