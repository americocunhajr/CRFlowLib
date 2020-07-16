
/*
* -----------------------------------------------------------------
*  Pairwise Mixing Stirred Reaction ISAT-DI --- main__pmsr-error.c
*  Version: 1.6180
*  Date: Dec 15, 2010
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
*  stirred reactor using ISAT and DI approaches, which 
*  estimate the error in ISAT calculations.
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
    double epsG;             /* ISAT global error            */
    double epsk;             /* auxiliar global error        */
    gsl_vector *phi_inflow;  /* fresh gases composition      */
    gsl_vector *phi_inside;  /* burned gases composition     */
    gsl_vector *Rphi;        /* reaction mapping             */
    gsl_vector *Rphi_DI;     /* reaction mapping for DI      */
    gsl_vector *b;           /* scaling matrix elements      */
    gsl_matrix *A;           /* mapping gradient matrix      */
    gsl_matrix *L;           /* EOA Cholesky matrix          */
    gsl_rng *rand_di;        /* rng workspace                */
    gsl_rng *rand_isat;      /* rng workspace                */
    pmsr_wrk *pmsr_di;       /* pmsr workspace               */
    pmsr_wrk *pmsr_isat;     /* pmsr workspace               */
    thrm_wrk *thrm;          /* chemkin workspace            */
    isat_wrk *isat;          /* isat workspace               */
    void *cvode_mem_di;      /* CVODE workspace              */
    void *cvode_mem_isat;    /* CVODE workspace              */
    FILE *file1;
    FILE *file2;
    FILE *file3;
    FILE *file4;
    
    
    /* initiating variables */
    k              = 0;
    flag           = GSL_SUCCESS;
    Np             = 0;
    Neq            = 0;
    K              = 0;
    Nexc           = 0;
    Npair          = 0;
    mxsteps        = 5000;
    max_lf         = 0;
    seed           = 0;
    atol           = ATOL;
    rtol           = RTOL;
    t0             = 0.0;
    delta_t        = 0.0;
    tau_res        = 0.0;
    tau_mix        = 0.0;
    tau_pair       = 0.0;
    T              = 0.0;
    Tmax           = 0.0;
    Tmin           = 0.0;
    tol            = 0.0;
    etol           = 0.0;
    kappa          = 0.0;
    mean           = 0.0;
    var            = 0.0;
    epsG           = 0.0;
    epsk           = 0.0;
    phi_inflow     = NULL;
    phi_inside     = NULL;
    Rphi           = NULL;
    Rphi_DI        = NULL;
    b              = NULL;
    A              = NULL;
    L              = NULL;
    rand_di        = NULL;
    rand_isat      = NULL;
    pmsr_di        = NULL;
    pmsr_isat      = NULL;
    thrm           = NULL;
    isat           = NULL;
    cvode_mem_di   = NULL;
    cvode_mem_isat = NULL;
    file1          = NULL;
    file2          = NULL;
    file3          = NULL;
    file4          = NULL;
    
    
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
    pmsr_di = pmsr_alloc();
    if ( pmsr_di == NULL )
        GSL_ERROR(" can't alloc pmsr_wrk ",flag);
    
    pmsr_isat = pmsr_alloc();
    if ( pmsr_isat == NULL )
        GSL_ERROR(" can't alloc pmsr_wrk ",flag);
    
    /* memory allocation for ODE solver workspace */
    cvode_mem_di = CVodeCreate(CV_BDF,CV_NEWTON);
    if ( cvode_mem_di == NULL )
        GSL_ERROR(" can't alloc cvode_mem_di ",flag);
    
    cvode_mem_isat = CVodeCreate(CV_BDF,CV_NEWTON);
    if ( cvode_mem_isat == NULL )
        GSL_ERROR(" can't alloc cvode_mem_isat ",flag);

    
    /* memory allocation for isat workspace */
    isat = isat_alloc();
    if ( isat == NULL )
        GSL_ERROR(" can't alloc isat_wrk ",flag);
    
    /* memory allocation for rng workspace */
    rand_di   = gsl_rng_alloc(gsl_rng_taus);
    rand_isat = gsl_rng_alloc(gsl_rng_taus);
    
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
    flag = pmsr_init(Np,Neq,pmsr_di);
    if ( flag != GSL_SUCCESS )
        GSL_ERROR(" can't initiate pmsr_wrk ",flag);
    
    flag = pmsr_init(Np,Neq,pmsr_isat);
    if ( flag != GSL_SUCCESS )
        GSL_ERROR(" can't initiate pmsr_wrk ",flag);
    
    /* setting rng sead */
    gsl_rng_set(rand_di,seed);
    gsl_rng_set(rand_isat,seed);

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
    Rphi_DI    = gsl_vector_calloc(Neq);
    b          = gsl_vector_calloc(Neq);
    
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
    
    /* inputing initial parameters */
    printf("\n Input initial parameters:\n");
    thrm_composition(thrm,phi_inside);
    
    /* defining scalling elements */
    gsl_vector_set_all(b,1.0);
    b->data[0] = 1.0/phi_inside->data[0];

    /* setting the system of particles with the burned gases composition */
    pmsr_set_all(Np,phi_inside,pmsr_di);
    pmsr_set_all(Np,phi_inside,pmsr_isat);
    
    /* initiating ODE solver workspace */
    flag = odesolver_init(thrm_eqs,thrm,t0,mxsteps,
                            atol,rtol,phi_inflow,cvode_mem_di);
    if ( flag != GSL_SUCCESS )
        GSL_ERROR(" can't initiate cvode_mem_di ",flag);
    
    flag = odesolver_init(thrm_eqs,thrm,t0,mxsteps,
                            atol,rtol,phi_inflow,cvode_mem_isat);
    if ( flag != GSL_SUCCESS )
        GSL_ERROR(" can't initiate cvode_mem_isat ",flag);

    
    /* opening file to write the properties evolution */
    file1 = fopen("pmsr_isat-di_mean_vs_t.dat","w");
    file2 = fopen("pmsr_isat-di_var_vs_t.dat","w");
    file3 = fopen("pmsr_isat-di_output_vs_t.dat","w");
    file4 = fopen("pmsr_isat-di_phi_vs_index.dat","w");
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
    thrm_temp_meanvar(Np,pmsr_isat->pmsr_ps,thrm,Tmax,Tmin,tol,&mean,&var);
    fprintf(file1," %+.3e",mean);
    fprintf(file2," %+.3e",var);
    for ( i = 0; i < Neq; i++ )
    {
        pmsr_meanvar(Np,i,pmsr_isat->pmsr_ps,&mean,&var);
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
        pmsr_iops(Np,Nexc,Npair,rand_di,phi_inflow,pmsr_di);
        pmsr_iops(Np,Nexc,Npair,rand_isat,phi_inflow,pmsr_isat);
        
        /* mixture step */
        pmsr_mixture(Np,Neq,delta_t,tau_mix,pmsr_di->pmsr_idx,pmsr_di->pmsr_ps);
        pmsr_mixture(Np,Neq,delta_t,tau_mix,pmsr_isat->pmsr_idx,pmsr_isat->pmsr_ps);
        
        
        /* reaction step */
        for ( i = 0; i < Np && flag == GSL_SUCCESS; i++ )
        {
            /* initiating ODE solver */
            flag = odesolver_reinit(thrm_eqs,thrm,t0,
                            ATOL,RTOL,pmsr_di->pmsr_ps[i],cvode_mem_di);
            
            /* direct integration */
            if ( flag == GSL_SUCCESS )
                flag = odesolver(cvode_mem_di,delta_t,Rphi_DI);
            
            /* in situ adaptive tabulation */
            flag = isat4(isat,thrm,cvode_mem_isat,etol,kappa,
                            t0,delta_t,pmsr_isat->pmsr_ps[i],A,L,Rphi);
            gsl_vector_memcpy(pmsr_isat->pmsr_ps[i],Rphi);
            
            /* Rphi_DI := R(phi)_DI - R(phi)_ISAT */
            gsl_vector_sub(Rphi_DI,Rphi);
            /* Rphi_DI := B*(R(phi)_DI - R(phi)_ISAT) */
            gsl_vector_mul(b,Rphi_DI);
            /* computing k-th global error */
            epsk += gsl_blas_dnrm2(Rphi_DI);
            
            gsl_vector_set_zero(Rphi);
            gsl_vector_set_zero(Rphi_DI);
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
        thrm_temp_meanvar(Np,pmsr_isat->pmsr_ps,thrm,Tmax,Tmin,tol,&mean,&var);
        fprintf(file1," %+.3e",mean);
        fprintf(file2," %+.3e",var);
        
        /* printing other properties mean and variance */
        for ( i = 0; i < Neq; i++ )
        {
            pmsr_meanvar(Np,i,pmsr_isat->pmsr_ps,&mean,&var);
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
                fprintf(file4,"\n %6d",pmsr_isat->pmsr_idx[i]);
                T = thrm_h2T(thrm,pmsr_isat->pmsr_ps[i]->data,Tmax,Tmin,tol);
                fprintf(file4," %+.3e",T);
                gsl_fprint_vector(file4,pmsr_isat->pmsr_ps[i]);
            }
        
        /* computing ISAT global error */
        epsG += epsk;
        epsk  = 0.0;
        
        /* updating events counter */
        k++;
    }
    
    /* computing ISAT global error */
    epsG = epsG/(K*Np);
    
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
    
    printf("\n\n ISAT global error = %+.6e",epsG);
    
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
    gsl_vector_free(Rphi_DI);
    gsl_matrix_free(A);
    gsl_matrix_free(L);
    gsl_rng_free(rand_di);
    gsl_rng_free(rand_isat);
    pmsr_free(Np,(void **)&pmsr_di);
    pmsr_free(Np,(void **)&pmsr_isat);
    thrm_free((void **)&thrm);
    isat_free((void **)&isat);
    CVodeFree(&cvode_mem_di);
    CVodeFree(&cvode_mem_isat);
    
    phi_inflow     = NULL;
    phi_inside     = NULL;
    Rphi           = NULL;
    Rphi_DI        = NULL;
    A              = NULL;
    L              = NULL;
    rand_di        = NULL;
    rand_isat      = NULL;
    pmsr_di        = NULL;
    pmsr_isat      = NULL;
    thrm           = NULL;
    isat           = NULL;
    cvode_mem_di   = NULL;
    cvode_mem_isat = NULL;
    file1          = NULL;
    file2          = NULL;
    file3          = NULL;
    file4          = NULL;
    
    return GSL_SUCCESS;
}
