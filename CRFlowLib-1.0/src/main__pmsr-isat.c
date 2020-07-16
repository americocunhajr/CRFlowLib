
/*
* -----------------------------------------------------------------
*  Pairwise Mixing Stirred Reaction ISAT --- main__pmsr-isat.c
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
*  This is the implementation file of a pairwise mixing
*  stirred reactor using ISAT approach.
* -----------------------------------------------------------------
*/




#include <math.h>
#include <time.h>

#include "../include/util_lib.h"
#include "../include/thrm_lib.h"
#include "../include/pmsr_lib.h"
#include "../include/ode_lib.h"
#include "../include/isat_lib.h"




/*
* -----------------------------------------------------------------
*   main
*
*   This is the main function of PMSR/ISAT program.
* -----------------------------------------------------------------
*/

int main(void)
{
    unsigned int i;
    unsigned int k;          /* events counter               */
    unsigned int flag;       /* flag                         */
    unsigned int Np;         /* # of particles               */
    unsigned int Neq;        /* # of equations               */
    unsigned int K;          /* # of time steps              */
    unsigned int Nexc;       /* particles to exchange        */
    unsigned int Npair;      /* particles to pairwize        */
    unsigned int mxsteps;    /* max of solver iterations     */
    unsigned int max_lf;     /* maximun value of tree leaves */
    unsigned long int seed;  /* rng seed                     */
    clock_t cpu_start;       /* cpu clock start flag         */
    clock_t cpu_end;         /* cpu clock end flag           */
    time_t wall_start;       /* wall clock start flag        */
    time_t wall_end;         /* wall clock end flag          */
    double atol;             /* absolute tolerance           */
    double rtol;             /* relative tolerance           */
    double t0;               /* initial time                 */
    double delta_t;          /* time step                    */
    double tau_res;          /* residence time               */
    double tau_mix;          /* mixture time                 */
    double tau_pair;         /* pairwise time                */
    double T;                /* temperature                  */
    double Tmax;             /* temperature upper bound      */
    double Tmin;             /* temperature lower bound      */
    double tol;              /* tolerance                    */
    double etol;             /* error tolerance              */
    double kappa;            /* ISAT lower bound             */
    double mean;             /* properties mean value        */
    double var;              /* properties variance          */
    gsl_vector *phi_inflow;  /* fresh gases composition      */
    gsl_vector *phi_inside;  /* burned gases composition     */
    gsl_vector *Rphi;        /* reaction mapping             */
    gsl_matrix *A;           /* mapping gradient matrix      */
    gsl_matrix *L;           /* EOA Cholesky matrix          */
    gsl_rng *rand;           /* rng workspace                */
    pmsr_wrk *pmsr;          /* pmsr workspace               */
    thrm_wrk *thrm;          /* chemkin workspace            */
    isat_wrk *isat;          /* isat workspace               */
    void *cvode_mem;         /* CVODE workspace              */
    FILE *file1;
    FILE *file2;
    FILE *file3;
    FILE *file4;
    
    
    /* initiating variables */
    k          = 0;
    flag       = GSL_SUCCESS;
    Np         = 0;
    Neq        = 0;
    K          = 0;
    Nexc       = 0;
    Npair      = 0;
    mxsteps    = 5000;
    max_lf     = 0;
    seed       = 0;
    atol       = ATOL;
    rtol       = RTOL;
    t0         = 0.0;
    delta_t    = 0.0;
    tau_res    = 0.0;
    tau_mix    = 0.0;
    tau_pair   = 0.0;
    T          = 0.0;
    Tmax       = 0.0;
    Tmin       = 0.0;
    tol        = 0.0;
    etol       = 0.0;
    kappa      = 0.0;
    mean       = 0.0;
    var        = 0.0;
    phi_inflow = NULL;
    phi_inside = NULL;
    Rphi       = NULL;
    A          = NULL;
    L          = NULL;
    rand       = NULL;
    pmsr       = NULL;
    thrm       = NULL;
    isat       = NULL;
    cvode_mem  = NULL;
    file1      = NULL;
    file2      = NULL;
    file3      = NULL;
    file4      = NULL;
    
    
    /* CPU and Wall clocks flags in the program begining */
    cpu_start  = clock();
    wall_start = time(NULL);
    
    /* program title */
    pmsr_title();
    
    /* memory allocation for chemkin workspace */
    thrm = thrm_alloc();
    if ( thrm == NULL )
        GSL_ERROR(" can't alloc thrm_wrk ",flag);
    
    /* memory allocation for pmsr workspace */
    pmsr = pmsr_alloc();
    if ( pmsr == NULL )
        GSL_ERROR(" can't alloc pmsr_wrk ",flag);
    
    /* memory allocation for ODE solver workspace */
    cvode_mem = CVodeCreate(CV_BDF,CV_NEWTON);
    if ( cvode_mem == NULL )
        GSL_ERROR(" can't alloc cvode_mem ",flag);
    
    /* memory allocation for isat workspace */
    isat = isat_alloc();
    if ( isat == NULL )
        GSL_ERROR(" can't alloc isat_wrk ",flag);
    
    /* memory allocation for rng workspace */
    rand = gsl_rng_alloc(gsl_rng_taus);
    
    /* inputing isat parameters */
    isat_input(&max_lf,&etol,&kappa);
    
    /* inputing pmsr parameters */
    pmsr_input(&seed,&Np,&K,&t0,&delta_t,
        &tau_res,&tau_mix,&tau_pair,&Tmin,&Tmax);
    
    /* initiating chemkin workspace */
    flag = thrm_init(thrm);
    if ( flag != GSL_SUCCESS )
        GSL_ERROR(" can't initiate thrm_wrk ",flag);
    
    /* setting the number of equations */
    Neq = thrm->n_s + 2;

    /* initiating pmsr workspace */
    flag = pmsr_init(Np,Neq,pmsr);
    if ( flag != GSL_SUCCESS )
        GSL_ERROR(" can't initiate pmsr_wrk ",flag);
    
    /* setting rng seed */
    gsl_rng_set(rand,seed);

    /* initiating isat workspace */
    flag = isat_set(isat);
    if ( flag != GSL_SUCCESS )
        GSL_ERROR(" can't initiate isat_wrk ",flag);
    
    /* setting maximum value of tree leaves */
    isat->max_lf = max_lf;
    
    /* memory allocation for composition vectors */
    phi_inflow = gsl_vector_calloc(Neq);
    phi_inside = gsl_vector_calloc(Neq);
    Rphi       = gsl_vector_calloc(Neq);
    
    /* memory allocation for ISAT matrices */
    A = gsl_matrix_calloc(Neq,Neq);
    L = gsl_matrix_calloc(Neq,Neq);
    
    /* computing Nexc and Npair */
    Nexc  = pmsr_Nexc(Np,delta_t,tau_res);
    Npair = pmsr_Npair(Np,delta_t,tau_pair);
    
    /* printing mechanism information */
    thrm_mech(thrm->n_e,thrm->n_s,thrm->n_r,"---");
    
    /* inputing inflow parameters */
    printf("\n Input inflow parameters:\n");
    thrm_composition(thrm,phi_inflow);
    
    /* inputing initial parametersfr */
    printf("\n Input initial parameters:\n");
    thrm_composition(thrm,phi_inside);

    /* setting the system of particles with the burned gases composition */
    pmsr_set_all(Np,phi_inside,pmsr);
    
    /* initiating ODE solver workspace */
    flag = odesolver_init(thrm_eqs,thrm,t0,mxsteps,
                            atol,rtol,phi_inflow,cvode_mem);
    if ( flag != GSL_SUCCESS )
        GSL_ERROR(" can't initiate cvode_mem ",flag);
    
    /* opening file to write the properties evolution */
    file1 = fopen("pmsr_isat_mean_vs_t.dat","w");
    file2 = fopen("pmsr_isat_var_vs_t.dat","w");
    file3 = fopen("pmsr_isat_output_vs_t.dat","w");
    file4 = fopen("pmsr_isat_phi_vs_index.dat","w");
    if ( file1 == NULL || file2 == NULL || file3 == NULL || file4 == NULL )
    {
        printf("can't open an output file (*.dat)\n");
        exit(1);
    }
    
    /* printing initial time */
    fprintf(file1,"\n %+.3e",t0);
    fprintf(file2,"\n %+.3e",t0);
    fprintf(file3,"\n %+.3e",t0);
    
    /* printing initial mean and variance */
    thrm_temp_meanvar(Np,pmsr->pmsr_ps,thrm,Tmax,Tmin,tol,&mean,&var);
    fprintf(file1," %+.3e",mean);
    fprintf(file2," %+.3e",var);
    for ( i = 0; i < Neq; i++ )
    {
        pmsr_meanvar(Np,i,pmsr->pmsr_ps,&mean,&var);
        fprintf(file1," %+.3e",mean);
        fprintf(file2," %+.3e",var);
    }
    
    /* printing initial ISAT output values */
    fprintf(file3," %7d",isat->add);
    fprintf(file3," %7d",isat->grw);
    fprintf(file3," %7d",isat->rtv);
    fprintf(file3," %7d",isat->dev);
    fprintf(file3," %7d",isat->hgt);
    
    
    /* solve governing equations for Nev events */
    while ( k < K && flag == GSL_SUCCESS )
    {
        fprintf(file1,"\n %+.3e",(k+1)*delta_t);
        fprintf(file2,"\n %+.3e",(k+1)*delta_t);
        fprintf(file3,"\n %+.3e",(k+1)*delta_t);
        
        /* inflow, outflow and pairing */
        pmsr_iops(Np,Nexc,Npair,rand,phi_inflow,pmsr);
        
        /* mixture step */
        pmsr_mixture(Np,Neq,delta_t,tau_mix,pmsr->pmsr_idx,pmsr->pmsr_ps);
        
        /* reaction step */
        for ( i = 0; i < Np && flag == GSL_SUCCESS; i++ )
        {            
            /* in situ adaptive tabulation */
            flag = isat4(isat,thrm,cvode_mem,etol,kappa,
                            t0,delta_t,pmsr->pmsr_ps[i],A,L,Rphi);
            gsl_vector_memcpy(pmsr->pmsr_ps[i],Rphi);
            gsl_vector_set_zero(Rphi);
            gsl_matrix_set_zero(A);
            gsl_matrix_set_zero(L);
        }
        
        /* checking if ISAT exceed the memory limit */
        if ( flag == GSL_ENOMEM )
        {   
            printf("\n");
            isat_statistics(isat);
            GSL_ERROR(" ISAT exceed the limit of tree leaves",flag);
        }
        
        /* printing temperature mean and variance */
        thrm_temp_meanvar(Np,pmsr->pmsr_ps,thrm,Tmax,Tmin,tol,&mean,&var);
        fprintf(file1," %+.3e",mean);
        fprintf(file2," %+.3e",var);
        
        /* printing other properties mean and variance */
        for ( i = 0; i < Neq; i++ )
        {
            pmsr_meanvar(Np,i,pmsr->pmsr_ps,&mean,&var);
            fprintf(file1," %+.3e",mean);
            fprintf(file2," %+.3e",var);
        }
        
        /* printing ISAT output values */
        fprintf(file3," %7d",isat->add);
        fprintf(file3," %7d",isat->grw);
        fprintf(file3," %7d",isat->rtv);
        fprintf(file3," %7d",isat->dev);
        fprintf(file3," %7d",isat->hgt);
        
        /* printing the properties of all particles */
        if ( k > K - ceil(50*tau_res/delta_t) )
            for ( i = 0; i < Np; i++ )
            {
                fprintf(file4,"\n %6d",pmsr->pmsr_idx[i]);
                T = thrm_h2T(thrm,pmsr->pmsr_ps[i]->data,Tmax,Tmin,tol);
                fprintf(file4," %+.3e",T);
                gsl_fprint_vector(file4,pmsr->pmsr_ps[i]);
            }
        
        /* updating events counter */
        k++;
    }
    
    
    /* closing file1, file2, file3, and file4 */
    fprintf(file1,"\n");
    fprintf(file2,"\n");
    fprintf(file3,"\n");
    fprintf(file4,"\n");
    fclose(file1);
    fclose(file2);
    fclose(file3);
    fclose(file4);
    
    /* printing isat statistics */
    printf("\n");
    isat_statistics(isat);
    
    /* CPU and Wall clocks at program end */
    cpu_end  = clock();
    wall_end = time(NULL);
    
    /* printing cpu and wall time spent at program execution */
    util_time_used(cpu_start,cpu_end,wall_start,wall_end);
    
    /* checking for errors in simalation */
    if ( flag != GSL_SUCCESS )
        printf("\n\n Simualtion finish with error\n");
    else
        printf("\n\n Simulation finish without error\n");
    
    
    /* releasing allocated memory */
    gsl_vector_free(phi_inflow);
    gsl_vector_free(phi_inside);
    gsl_vector_free(Rphi);
    gsl_matrix_free(A);
    gsl_matrix_free(L);
    gsl_rng_free(rand);
    pmsr_free(Np,(void **)&pmsr);
    thrm_free((void **)&thrm);
    isat_free((void **)&isat);
    CVodeFree(&cvode_mem);
    
    phi_inflow = NULL;
    phi_inside = NULL;
    Rphi       = NULL;
    A          = NULL;
    L          = NULL;
    rand       = NULL;
    pmsr       = NULL;
    thrm       = NULL;
    isat       = NULL;
    cvode_mem  = NULL;
    file1      = NULL;
    file2      = NULL;
    file3      = NULL;
    file4      = NULL;
    
    return GSL_SUCCESS;
}
