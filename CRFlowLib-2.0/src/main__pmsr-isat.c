
/*
* -----------------------------------------------------------------
*  main__pmsr-isat.c
*  Pairwise Mixing Stirred Reactor via ISAT
*  Version: 2.0
*  Last Update: Nov 3, 2019
* 
*  This is the main program for an implementation of
*  a Pairwise Mixing Stirred Reactor (PMSR) employing 
*  In Situ Adaptive Tabulation (ISAT) algorithm for 
*  integration of the governing equations.
* ----------------------------------------------------------------- 
*  Programmer: Americo Barbosa da Cunha Junior
*              americo.cunhajr@gmail.com
* -----------------------------------------------------------------
*  Copyright (c) 2019 by Americo Barbosa da Cunha Junior
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
*   This is the main function for PMSR-ISAT program.
* -----------------------------------------------------------------
*/

int main(void)
{
    unsigned int i;
    unsigned int events;     /* events counter               */
    unsigned int flag;       /* return value flag            */
    unsigned int Np;         /* # of particles               */
    unsigned int Neq;        /* # of equations               */
    unsigned int Ndt;        /* # of time steps              */
    unsigned int Nexc;       /* # of particles to exchange   */
    unsigned int Npair;      /* # of particles to pairwise   */
    unsigned int maxsteps;   /* max of solver iterations     */
    unsigned int maxleaves;  /* maximun value of tree leaves */
    unsigned long int seed;  /* RNG seed                     */
    clock_t cpu_start;       /* CPU clock start flag         */
    clock_t cpu_end;         /* CPU clock end flag           */
    time_t wall_start;       /* Wall clock start flag        */
    time_t wall_end;         /* Wall clock end flag          */
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
    double mean;             /* properties mean value        */
    double var;              /* properties variance          */
    gsl_vector *phi_inflow;  /* fresh gases composition      */
    gsl_vector *phi_inside;  /* burned gases composition     */
    gsl_vector *Rphi;        /* reaction mapping             */
    gsl_matrix *A;           /* mapping gradient matrix      */
    gsl_matrix *L;           /* EOA Cholesky matrix          */
    gsl_rng *rand;           /* RNG workspace                */
    pmsr_wrk *pmsr;          /* PMSR workspace               */
    thrm_wrk *thrm;          /* Chemkin workspace            */
    isat_wrk *isat;          /* ISAT workspace               */
    void *cvode_mem;         /* CVODE workspace              */
    FILE *file1;
    FILE *file2;
    FILE *file3;
    FILE *file4;
    
    
    /* initiating variables */
    events     = 0;
    flag       = GSL_SUCCESS;
    Np         = 0;
    Neq        = 0;
    Ndt        = 0;
    Nexc       = 0;
    Npair      = 0;
    maxsteps   = 5000;
    maxleaves  = 0;
    seed       = 0;
    atol       = 0.0;
    rtol       = 0.0;
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
    
    
    /* program title */
    pmsr_title();
    
    /* CPU clock and Wall clock at the program begining */
    cpu_start  = clock();
    wall_start = time(NULL);
    
    /* allocate memory for Chemkin workspace */
    thrm = thrm_alloc();
    if ( thrm == NULL )
        GSL_ERROR(" thrm_wrk cannot be allocated",GSL_ENOMEM);
    
    /* allocate memory for PMSR workspace */
    pmsr = pmsr_alloc();
    if ( pmsr == NULL )
        GSL_ERROR(" pmsr_wrk cannot be allocated",GSL_ENOMEM);
    
   /* creates CVODE solver object and defines solution method */
    cvode_mem = CVodeCreate(CV_BDF);
    if ( cvode_mem == NULL )
        GSL_ERROR(" cvode_mem cannot be allocated",GSL_ENOMEM);
    
    /* allocate memory for ISAT workspace */
    isat = isat_alloc();
    if ( isat == NULL )
        GSL_ERROR(" isat_wrk cannot be allocated",flag);
    
    /* allocate memory for RNG workspace */
    rand = gsl_rng_alloc(gsl_rng_taus);
    
    /* input ISAT parameters */
    isat_input(&maxleaves,&etol);
    
    /* input PMSR parameters */
    pmsr_input(&seed,&Np,&Ndt,&t0,&delta_t,
                  &tau_res,&tau_mix,&tau_pair,
                              &Tmin,&Tmax,&atol,&rtol);
    
    /* initiate Chemkin workspace */
    flag = thrm_init(thrm);
    if ( flag != GSL_SUCCESS )
        GSL_ERROR(" thrm_wrk cannot be initiated",flag);
    
    /* set the number of equations */
    Neq = thrm->n_s + 2;

    /* initiate PMSR workspace */
    flag = pmsr_init(Np,Neq,pmsr);
    if ( flag != GSL_SUCCESS )
        GSL_ERROR(" pmsr_wrk cannot be initiated",flag);
    
    /* set RND seed */
    gsl_rng_set(rand,seed);

    /* initiate ISAT workspace */
    flag = isat_bst_init(isat);
    if ( flag != GSL_SUCCESS )
        GSL_ERROR(" isat_wrk cannot be initiated",flag);
    
    /* set the maximum value of tree leaves */
    isat->maxleaves = maxleaves;
    
    /* allocate memory for composition vectors */
    phi_inflow = gsl_vector_calloc(Neq);
    phi_inside = gsl_vector_calloc(Neq);
    Rphi       = gsl_vector_calloc(Neq);
    
    /* allocate memory for ISAT matrices */
    A = gsl_matrix_calloc(Neq,Neq);
    L = gsl_matrix_calloc(Neq,Neq);
    
    /* compute Nexc and Npair */
    Nexc  = pmsr_Nexc(Np,delta_t,tau_res);
    Npair = pmsr_Npair(Np,delta_t,tau_pair);
    
    /* print mechanism information */
    thrm_mech(thrm->n_e,thrm->n_s,thrm->n_r,"---");
    
    /* input inflow parameters */
    printf("\n Input inflow parameters:\n");
    thrm_composition(thrm,phi_inflow);
    
    /* input initial parameters fr */
    printf("\n Input initial parameters:\n");
    thrm_composition(thrm,phi_inside);

    /* set the system of particles with the burned gases composition */
    pmsr_set_all(Np,phi_inside,pmsr);
    
    /* initiate ODE solver workspace */
    flag = odesolver_init(thrm_eqs,thrm,t0,maxsteps,
                            atol,rtol,phi_inflow,cvode_mem);
    if ( flag != GSL_SUCCESS )
        GSL_ERROR(" cvode_mem cannot be initiated",flag);
    
    /* open files to write the properties evolution */
    file1 = fopen("results/pmsr-isat_mean_vs_t.dat","w");
    file2 = fopen("results/pmsr-isat_var_vs_t.dat","w");
    file3 = fopen("results/pmsr-isat_output_vs_t.dat","w");
    file4 = fopen("results/pmsr-isat_phi_vs_index.dat","w");
    if ( file1 == NULL || file2 == NULL || file3 == NULL || file4 == NULL )
    {
        printf("can't open an output file (*.dat)\n");
        exit(1);
    }
    
    /* print initial time */
    fprintf(file1,"\n %+.3e",t0);
    fprintf(file2,"\n %+.3e",t0);
    fprintf(file3,"\n %+.3e",t0);
    
    /* print initial mean and variance */
    thrm_temp_meanvar(Np,pmsr->pmsr_ps,thrm,Tmax,Tmin,tol,&mean,&var);
    fprintf(file1," %+.3e",mean);
    fprintf(file2," %+.3e",var);
    for ( i = 0; i < Neq; i++ )
    {
        pmsr_meanvar(Np,i,pmsr->pmsr_ps,&mean,&var);
        fprintf(file1," %+.3e",mean);
        fprintf(file2," %+.3e",var);
    }
    
    /* print initial ISAT output values */
    fprintf(file3," %7d",isat->add);
    fprintf(file3," %7d",isat->grw);
    fprintf(file3," %7d",isat->rtv);
    fprintf(file3," %7d",isat->dev);
    fprintf(file3," %7d",isat->hgt);
    
    
    /* solve governing equations for Ndt events */
    while ( events < Ndt && flag == GSL_SUCCESS )
    {
        fprintf(file1,"\n %+.3e",(events+1)*delta_t);
        fprintf(file2,"\n %+.3e",(events+1)*delta_t);
        fprintf(file3,"\n %+.3e",(events+1)*delta_t);
        
        /* inflow, outflow and pairing */
        pmsr_iops(Np,Nexc,Npair,rand,phi_inflow,pmsr);
        
        /* mixture step */
        pmsr_mixture(Np,Neq,delta_t,tau_mix,pmsr->pmsr_idx,pmsr->pmsr_ps);
        
        /* reaction step */
        for ( i = 0; i < Np && flag == GSL_SUCCESS; i++ )
        {            
            /* in situ adaptive tabulation */
            flag = isat4(isat,cvode_mem,etol,
                            t0,delta_t,pmsr->pmsr_ps[i],A,L,Rphi);
            gsl_vector_memcpy(pmsr->pmsr_ps[i],Rphi);
            gsl_vector_set_zero(Rphi);
            gsl_matrix_set_zero(A);
            gsl_matrix_set_zero(L);
        }
        
        /* check if ISAT exceed the memory limit */
        if ( flag == GSL_ENOMEM )
        {   
            printf("\n");
            isat_statistics(isat);
            GSL_ERROR(" ISAT exceed the limit of tree leaves",flag);
        }
        
        /* print temperature mean and variance */
        thrm_temp_meanvar(Np,pmsr->pmsr_ps,thrm,Tmax,Tmin,tol,&mean,&var);
        fprintf(file1," %+.3e",mean);
        fprintf(file2," %+.3e",var);
        
        /* print other properties mean and variance */
        for ( i = 0; i < Neq; i++ )
        {
            pmsr_meanvar(Np,i,pmsr->pmsr_ps,&mean,&var);
            fprintf(file1," %+.3e",mean);
            fprintf(file2," %+.3e",var);
        }
        
        /* print ISAT output values */
        fprintf(file3," %7d",isat->add);
        fprintf(file3," %7d",isat->grw);
        fprintf(file3," %7d",isat->rtv);
        fprintf(file3," %7d",isat->dev);
        fprintf(file3," %7d",isat->hgt);
        
        /* print the properties of all particles */
        if ( events > Ndt - ceil(50*tau_res/delta_t) )
            for ( i = 0; i < Np; i++ )
            {
                fprintf(file4,"\n %6d",pmsr->pmsr_idx[i]);
                T = thrm_h2T(thrm,pmsr->pmsr_ps[i]->data,Tmax,Tmin,tol);
                fprintf(file4," %+.3e",T);
                gsl_fprint_vector(file4,pmsr->pmsr_ps[i]);
            }
        
        /* update events counter */
        events++;
    }
    
    
    /* close file1, file2, file3, and file4 */
    fprintf(file1,"\n");
    fprintf(file2,"\n");
    fprintf(file3,"\n");
    fprintf(file4,"\n");
    fclose(file1);
    fclose(file2);
    fclose(file3);
    fclose(file4);
    
    /* print isat statistics */
    printf("\n");
    isat_statistics(isat);
    
    /* CPU clock and Wall clock at the program end */
    cpu_end  = clock();
    wall_end = time(NULL);
    
    /* print CPU time and Wall time spent at program execution */
    util_time_used(cpu_start,cpu_end,wall_start,wall_end);
    
    /* check for errors in simulation */
    if ( flag != GSL_SUCCESS )
        printf("\n\n Simualtion finish with error\n");
    else
        printf("\n\n Simulation finish without error\n");
    
    
    /* release allocated memory */
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
