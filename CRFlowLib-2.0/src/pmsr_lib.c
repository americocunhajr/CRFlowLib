
/*
* -----------------------------------------------------------------
*  pmsr_lib.c
*  Pairwise Mixing Stirred Reactor Library
*  Version: 2.0
*  Last Update: Oct 1, 2019
* 
*  This is a computational library with routines to implement
*  a Pairwise Mixing Stirred Reactor (PMSR) model.
* ----------------------------------------------------------------- 
*  Programmer: Americo Barbosa da Cunha Junior
*              americo.cunhajr@gmail.com
* -----------------------------------------------------------------
*  Copyright (c) 2019 by Americo Barbosa da Cunha Junior
* -----------------------------------------------------------------
*/




#include <math.h>
#include <time.h>
#include <gsl/gsl_randist.h>
#include <gsl/gsl_errno.h>

#include "../include/pmsr_lib.h"




/*
* -----------------------------------------------------------------
*   pmsr_title
*
*   This function prints the program title on the screen.
*
*   last update: Oct 20, 2009
* -----------------------------------------------------------------
*/

void pmsr_title()
{
    time_t curtime = time (NULL);
    struct tm *loctime = localtime (&curtime);
    
    printf("\n");
    printf("=============================================");
    printf("\n Pairwise Mixing Stirred Reactor --- PMSR\n");
    printf("\n by");
    printf("\n Americo Barbosa da Cunha Junior");
    printf("\n americo.cunhajr@gmail.com\n");
    printf("\n Date:");
    printf("\n %s", asctime(loctime) );
    printf("=============================================");
    printf("\n");
    
    return;
}
/*----------------------------------------------------------------*/




/*
* -----------------------------------------------------------------
*   pmsr_input
*
*   This function receives PMSR parameters from user
*   and returns GSL_SUCCESS if there is no error.
*
*   Input:
*   seed         - RNG seed
*   Np           - # of particles
*   Ndt          - # of time steps
*   t0           - initial time
*   delta_t      - time step
*   tau_res      - residence time
*   tau_mix      - mixture time
*   tau_pair     - pairwise time
*
*   Output:
*   success or error
*
*   last update: Oct 1, 2019
* -----------------------------------------------------------------
*/

int pmsr_input(unsigned long int *seed,unsigned int *Np,unsigned int *Ndt,
	       double *t0,double *delta_t,double *tau_res,double *tau_mix,
	       double *tau_pair,double *Tmax,double *Tmin,
               double *atol,double *rtol)
{
    printf("\n Input PMSR parameters:\n");

    printf("\n RNG seed:");
    scanf("%lu", seed);
    printf("\n %lu\n", *seed);
    if( *seed <= 0 )
	*seed = time(NULL);
    
    printf("\n # of particles:");
    scanf("%d", Np);
    printf("\n %d\n", *Np);
    if( *Np <= 0 || *Np % 2 != 0 )
	GSL_ERROR(" Np must be a positive even integer",GSL_EINVAL);
    
    printf("\n # of time steps:");
    scanf("%d", Ndt);
    printf("\n %d\n", *Ndt);
    if( *Ndt <= 0 )
	GSL_ERROR(" Ndt must be a positive integer",GSL_EINVAL);

    printf("\n initial time:");
    scanf("%lf", t0);
    printf("\n %+.1e\n", *t0);
    if( *t0 < 0.0 )
	GSL_ERROR(" t0 must be a grather than zero",GSL_EINVAL);
    
    printf("\n time step:");
    scanf("%lf", delta_t);
    printf("\n %+.1e\n", *delta_t);
    if( *delta_t <= 0.0 )
	GSL_ERROR(" delta_t must be a grather than zero",GSL_EINVAL);

    printf("\n residence time:");
    scanf("%lf", tau_res);
    printf("\n %+.1e\n", *tau_res);
    if( *tau_res <= 0.0 )
	GSL_ERROR(" tau_res must be a grather than zero",GSL_EINVAL);
    
    printf("\n mixture time:");
    scanf("%lf", tau_mix);
    printf("\n %+.1e\n", *tau_mix);
    if( *tau_mix <= 0.0 )
	GSL_ERROR(" tau_mix must be a grather than zero",GSL_EINVAL);
    
    printf("\n pairwise time:");
    scanf("%lf", tau_pair);
    printf("\n %+.1e\n", *tau_pair);
    if( *tau_pair <= 0.0 )
	GSL_ERROR(" tau_pair must be a grather than zero",GSL_EINVAL);
    
    printf("\n temperature upper bound:");
    scanf("%lf", Tmax);
    printf("\n %+.1e\n", *Tmax);
    if( *Tmax <= 0.0 )
	GSL_ERROR(" Tmax must be a grather than zero",GSL_EINVAL);
    
    printf("\n temperature lower bound:");
    scanf("%lf", Tmin);
    printf("\n %+.1e\n", *Tmin);
    if( *Tmin < 0.0 || *Tmin > *Tmax )
	GSL_ERROR(" Tmin must be less than Tmax and a grather than zero",GSL_EINVAL);

    printf("\n absolute tolerance:");
    scanf("%lf", atol);
    printf("\n %+.1e\n", *atol);
    if( *atol <= 0.0 )
	GSL_ERROR(" atol must be a grather than zero",GSL_EINVAL);

    printf("\n relative tolerance:");
    scanf("%lf", rtol);
    printf("\n %+.1e\n", *rtol);
    if( *rtol <= 0.0 )
	GSL_ERROR(" rtol must be a grather than zero",GSL_EINVAL);
    
    return GSL_SUCCESS;
}
/*----------------------------------------------------------------*/




/*
* -----------------------------------------------------------------
*  pmsr_alloc
*
*  This function allocates an internal memory block for
*  pmsr workspace struct. If successful, returns a pointer to
*  pmsr_wrk. If a startup error occurs returns NULL.
*
*   last update: Nov 5, 2009
* -----------------------------------------------------------------
*/

pmsr_wrk *pmsr_alloc()
{
    pmsr_wrk *pmsr = NULL;
    
    /* memory allocation for pmsr struct */
    pmsr = (pmsr_wrk *) malloc(sizeof(pmsr_wrk));
    if ( pmsr == NULL )
        return NULL;
    
    /* setting pmsr_wrk elements */
    pmsr->pmsr_idx = NULL;
    pmsr->pmsr_ps  = NULL;
    
    return pmsr;
}
/*----------------------------------------------------------------*/




/*
* -----------------------------------------------------------------
*  pmsr_init
*
*  This function initiates the system of particles
*  with an initial composition phi.
*
*   Input:
*   Np   - # of particles
*   Neq  - # of equations
*   pmsr - pointer to a struct pmsr
*
*   Output:
*   success or error
*
*  last update: Sep 22, 2009
* -----------------------------------------------------------------
*/

int pmsr_init(int Np,int Neq,pmsr_wrk *pmsr)
{
    int i;
    
    /* checking if pmsr is NULL */
    if ( pmsr == NULL )
        return GSL_EINVAL;
    
    /* memory allocation for the system of particles */
    pmsr->pmsr_ps = (gsl_vector **) malloc(Np*sizeof(gsl_vector *));
    if ( pmsr->pmsr_ps == NULL )
	return GSL_ENOMEM;
    
    for ( i = 0; i < Np; i++ )
        pmsr->pmsr_ps[i] = gsl_vector_calloc(Neq);
    
    /* memory allocation for the vector of indices */
    pmsr->pmsr_idx = (int *) malloc(Np*sizeof(int));
    if ( pmsr->pmsr_idx == NULL )
    {
	for ( i = 0; i < Np; i++ )
	{
	    gsl_vector_free(pmsr->pmsr_ps[i]);
	    pmsr->pmsr_ps[i] = NULL;
	}
	
        free(pmsr->pmsr_ps);
	pmsr->pmsr_ps = NULL;
        
	return GSL_ENOMEM;
    }
    
    /* initial configuration of particles pairwise */
    for ( i = 0; i < Np; i++ )
        pmsr->pmsr_idx[i] = i;
    
    return GSL_SUCCESS;
}
/*----------------------------------------------------------------*/




/*
* -----------------------------------------------------------------
*   pmsr_free
*
*   This routine frees the memory allocated by pmsr_alloc.
*
*   Input:
*   pmsr - pointer to a struct pmsr
*
*   Output:
*   void
*
*   last update: May 10, 2009
* -----------------------------------------------------------------
*/

void pmsr_free(int Np,void **pmsr_bl)
{
    int i;
    pmsr_wrk *pmsr = NULL;
    
    /* checking if pmsr_bl is NULL */
    if ( *pmsr_bl == NULL )
        return;
    
    /* checking if Np is positive */
    if ( Np <= 0 )
        return;
    
    pmsr = (pmsr_wrk *) (*pmsr_bl);
    
    /* releasing pmsr elements */
    if ( pmsr->pmsr_ps != NULL )
    {
	for ( i = 0; i < Np; i++ )
	    if ( pmsr->pmsr_ps[i] != NULL )
	    {
	        gsl_vector_free(pmsr->pmsr_ps[i]);
	        pmsr->pmsr_ps[i] = NULL;
	    }
	/* releasing pmsr struct */
	free(pmsr->pmsr_ps);
	pmsr->pmsr_ps = NULL;
    }
    
    if ( pmsr->pmsr_idx != NULL )
    {
	free(pmsr->pmsr_idx);
	pmsr->pmsr_idx = NULL;
    }
    
    free(*pmsr_bl);
    *pmsr_bl = NULL;
    
    return;
}
/*----------------------------------------------------------------*/




/*
* -----------------------------------------------------------------
*  pmsr_set_all
*
*  This function sets the system of particles with a given
*  composition phi.
*
*   Input:
*   pmsr - pointer to a struct pmsr
*   phi  - composition
*
*   Output:
*   success or error
*
*  last update: Mar 2, 2010
* -----------------------------------------------------------------
*/

void pmsr_set_all(int Np,gsl_vector *phi,pmsr_wrk *pmsr)
{
    int i;
    
    /* setting all particles with an initial phi */
    for ( i = 0; i < Np; i++ )
        gsl_vector_memcpy(pmsr->pmsr_ps[i],phi);
    
    return;
}
/*----------------------------------------------------------------*/




/*
* -----------------------------------------------------------------
*  pmsr_Nexc
*
*  This function computes the number of particles to be exchanged,
*  Nexc, corresponding to outﬂow/inﬂow particles.
*
*  Nexc  = ceil( 1/2 Np delta_t/taus_res )
*
*   Input:
*   Np      - # of particles
*   delta_t - time step
*   tau_res - residence time
*
*   Output:
*   Nexc    - # of particles to be exchanged
*
*   last update: Feb 12, 2009
* -----------------------------------------------------------------
*/

int pmsr_Nexc(int Np,double delta_t,double tau_res)
{   
    return ceil( (Np/2)*(delta_t/tau_res) );
}
/*----------------------------------------------------------------*/




/*
* -----------------------------------------------------------------
*  pmsr_Npair
*
*  This function computes the number of particles to be pairwised,
*  Npair, corresponding to pairing particles.
*
*  Npair = ceil( 1/2 Np delta_t/taus_pair )
*
*   Input:
*   Np       - # of particles
*   delta_t  - time step
*   tau_pair - pairwise time
*
*   Output:
*   Npair    - # of particles to be pairwised
*
*   last update: Feb 12, 2009
* -----------------------------------------------------------------
*/

int pmsr_Npair(int Np,double delta_t,double tau_pair)
{   
    return ceil( (Np/2)*(delta_t/tau_pair) );
}
/*----------------------------------------------------------------*/




/*
* -----------------------------------------------------------------
*   pmsr_meanvar
*
*   This function computes the mean and the variance
*   of a particular property in a system of particles
*   using the algorithm of Knuth/Welford.
*
*   Ref.:
*   D. E. Knuth (1998)
*   The Art of Computer Programming,
*   volume 2: Seminumerical Algorithms
*   3rd ed., Addison-Wesley.
*   
*   B. P. Welford (1962).
*   "Note on a method for calculating corrected sums
*   of squares and products". Technometrics 4(3):419–420
*
*   Input:
*   Np - # of particles
*   k  - property index
*   ps - system of particle
*
*   Output:
*   mean - propertie mean value
*   var  - propertie variance
*
*   last update: Mar 4, 2010
* -----------------------------------------------------------------
*/

void pmsr_meanvar(int Np,int k,gsl_vector **ps,double *mean,double *var)
{
    int i;
    double M     = 0.0;
    double S     = 0.0;
    double delta = 0.0;

    for ( i = 0; i < Np; i++ )
    {
        delta = ps[i]->data[k] - M;
	M    += delta/(double) (i+1);
        S    += delta*(ps[i]->data[k] - M);
    }

    /* mean */
    *mean  = M;
    
    /* variance */
    *var   = S/(double) (i-1);
   
    return;
}
/*----------------------------------------------------------------*/




/*
* -----------------------------------------------------------------
*   pmsr_mixture
*
*   This function computes the mixture step in PMSR model.
*
*   Input:
*   Np      - # of particles
*   Neq     - # of equations
*   delta_t - time step
*   tau_mix - mixture time
*   idx     - particles index vector
*
*   Output:
*   ps      - system of particles
*
*   last update: Mar 9, 2009
* -----------------------------------------------------------------
*/

void pmsr_mixture(int Np,int Neq,double delta_t,double tau_mix,
					    int *idx,gsl_vector **ps)
{
    int k, i;
    double sum, dif;
    
    /* particle index */
    k = 0;
    
    /* mixture of the Np particles */
    while( k < Np/2 )
    {
        /* mixture of two particles */
        for ( i = 0; i < Neq; i++ )
        {
            /* sum = 1/2(phi_p0 + phi_q0) */
            sum = 0.5*( ps[idx[2*k]]->data[i] + ps[idx[2*k+1]]->data[i] );
            
            /* dif = 1/2(phi_p0 - phi_q0) */
            dif = 0.5*( ps[idx[2*k]]->data[i] - ps[idx[2*k+1]]->data[i] );
            
            /* phi_p =  dif.exp(-2/tau delta_t) + sum */
            ps[idx[2*k]]->data[i] =
                        dif*exp(-(2.0/tau_mix)*delta_t) + sum;
            
            /* phi_q = -dif.exp(-2/tau delta_t) + sum */
            ps[idx[2*k+1]]->data[i] =
                        -dif*exp(-(2.0/tau_mix)*delta_t) + sum;
        }
        k++;
    }
    
    return;
}
/*----------------------------------------------------------------*/




/*
* -----------------------------------------------------------------
*   pmsr_iops
*
*   This function executes the inflow/outflow,
*   pairwise and shuflle steps.
*
*   Input:
*   Np      - # of particles
*   Nexc    - # of particles to be exchanged
*   Npair   - # of particles to be pairwised
*   rand    - rng workspace
*   phi     - inflow composition
*
*   Output:
*   p       - system of particles
*   pairs   - pairs of particles index
*
*   last update: Mar 1, 2010
* -----------------------------------------------------------------
*/

void pmsr_iops(int Np,int Nexc,int Npair,gsl_rng *rand,
				    gsl_vector *phi,pmsr_wrk *pmsr)
{
    int i;
    int idx1[Nexc+Npair];
    int idx2[2*(Nexc+Npair)];
    int *aux = NULL;
    
    /* memory allocation */
    aux = (int *) malloc((Np/2)*sizeof(int));
    
    /* setting aux vector with integers from 0 to Np/2-1 */
    for ( i = 0; i < Np/2; i++ )
        aux[i] = i;
    
    /* selecting Nexc + Npair pairs at random */
    gsl_ran_choose(rand,idx1,Nexc+Npair,aux,Np/2,sizeof(int));
    for ( i = 0; i < Nexc + Npair; i++ )
    {
        idx2[2*i  ] = pmsr->pmsr_idx[2*idx1[i]  ];
        idx2[2*i+1] = pmsr->pmsr_idx[2*idx1[i]+1];
    }
    
    
    /* changing composition of Nexc pairs (inflow / outflow) */
    for ( i = 0; i < Nexc; i++ )
    {
        gsl_vector_memcpy(pmsr->pmsr_ps[idx2[2*i]  ],phi);
        gsl_vector_memcpy(pmsr->pmsr_ps[idx2[2*i+1]],phi);
    }
    
    /* shuffling the selected pairs */
    gsl_ran_shuffle(rand,idx2,2*(Nexc+Npair),sizeof(int));
    
    /* setting the new partners in the system */
    for ( i = 0; i < Nexc + Npair; i++ )
    {
        pmsr->pmsr_idx[2*idx1[i]  ] = idx2[2*i  ];
        pmsr->pmsr_idx[2*idx1[i]+1] = idx2[2*i+1];
    }
    
    /* releasing allocated memory */
    free(aux);
    aux = NULL;
    
    return;
}
/*----------------------------------------------------------------*/

