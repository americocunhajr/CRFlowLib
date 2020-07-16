
/*
* -----------------------------------------------------------------
*  Partially Stirred Reaction Library --- pasr_lib.c
*  Version: 1.6180
*  Date: Mar 15, 2010
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
*  This is the implementation file of PaSR model.
* -----------------------------------------------------------------
*/




#include <math.h>
#include <time.h>
#include <gsl/gsl_randist.h>
#include <gsl/gsl_errno.h>

#include "../include/pasr_lib.h"




/*
* -----------------------------------------------------------------
*   pasr_title
*
*   This function prints the program title on the screen.
*
*   last update: Nov 4, 2009
* -----------------------------------------------------------------
*/

void pmsr_title()
{
    time_t curtime = time (NULL);
    struct tm *loctime = localtime (&curtime);
    
    printf("\n");
    printf("=============================================");
    printf("\n Partially Stirred Reactor --- PaSR\n");
    printf("\n Author:");
    printf("\n Americo Barbosa da Cunha Junior --- PUC-Rio");
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
*   pasr_input
*
*   This function receives pasr parameters
*   and returns GSL_SUCCESS if there is no error.
*
*   Input:
*   seed         - rng seed
*   Np           - # of particles
*   K            - # of time steps
*   initial time - time step
*   delta_t      - time step
*   tau_res      - residence time
*   tau_mix      - mixture time
*   tau_pair     - pairwize time
*
*   Output:
*   success or error
*
*   last update: Nov 16, 2010
* -----------------------------------------------------------------
*/

int pasr_input(unsigned long int *seed,unsigned int *Np,unsigned int *K,
			double *t0,double *delta_t,double *tau_res,
				double *tau_mix,double *Tmax,double *Tmin)
{
    printf("\n Input PaSR parameters:\n");
    
    printf("\n rng seed:");
    scanf("%lu", seed);
    printf("\n %lu\n", *seed);
    if( *seed <= 0 )
	*seed = time(NULL);
    
    printf("\n # of particles:");
    scanf("%d", Np);
    printf("\n %d\n", *Np);
    if( *Np <= 0 || *Np % 2 != 0 )
	return GSL_EINVAL;
    
    printf("\n # of time steps to advance in time:");
    scanf("%d", K);
    printf("\n %d\n", *K);
    if( *K <= 0 )
	return GSL_EINVAL;

    printf("\n initial time:");
    scanf("%lf", t0);
    printf("\n %+.1e\n", *t0);
    if( *t0 < 0.0 )
	return GSL_EINVAL;
    
    printf("\n time step:");
    scanf("%lf", delta_t);
    printf("\n %+.1e\n", *delta_t);
    if( *delta_t <= 0.0 )
	return GSL_EINVAL;

    printf("\n residence time:");
    scanf("%lf", tau_res);
    printf("\n %+.1e\n", *tau_res);
    if( *tau_res <= 0.0 )
	return GSL_EINVAL;
    
    printf("\n mixture time:");
    scanf("%lf", tau_mix);
    printf("\n %+.1e\n", *tau_mix);
    if( *tau_mix <= 0.0 )
	return GSL_EINVAL;
    
    printf("\n temperature upper bound:");
    scanf("%lf", Tmax);
    printf("\n %+.1e\n", *Tmax);
    if( *Tmax < 0.0 )
	return GSL_EINVAL;
    
    printf("\n temperature lower bound:");
    scanf("%lf", Tmin);
    printf("\n %+.1e\n", *Tmin);
    if( *Tmin < 0.0 || *Tmin > *Tmax )
	return GSL_EINVAL;
    
    return GSL_SUCCESS;
}
/*----------------------------------------------------------------*/




/*
* -----------------------------------------------------------------
*  pasr_alloc
*
*  This function allocates an internal memory block for
*  pasr workspace struct. If successful, returns a pointer to
*  pasr_wrk. If a startup error occurs returns NULL.
*  
*  last update: Nov 5, 2009
* -----------------------------------------------------------------
*/

pasr_wrk *pasr_alloc()
{
    pasr_wrk *pasr = NULL;
    
    /* memory allocation for pasr struct */
    pasr = (pasr_wrk *) malloc(sizeof(pasr_wrk));
    if ( pasr == NULL )
        return NULL;
    
    /* setting pasr_wrk elements */
    pasr->pasr_idx = NULL;
    pasr->pasr_ps  = NULL;
    
    return pasr;
}
/*----------------------------------------------------------------*/




/*
* -----------------------------------------------------------------
*  pasr_init
*
*  This function initiates the system of particles
*  with an initial composition phi.
*
*   Input:
*   Np   - # of particles
*   Neq  - # of equations
*   pasr - pointer to a struct pasr
*
*   Output:
*   success or error
*
*  last update: Nov 4, 2009
* -----------------------------------------------------------------
*/

int pasr_init(int Np,int Neq,pasr_wrk *pasr)
{
    int i;
    
    /* checking if pasr is NULL */
    if ( pasr == NULL )
        return GSL_EINVAL;
    
    /* memory allocation for the system of particles */
    pasr->pasr_ps = (gsl_vector **) malloc(Np*sizeof(gsl_vector *));
    if ( pasr->pasr_ps == NULL )
	return GSL_ENOMEM;
    
    for( i = 0; i < Np; i++ )
        pasr->pasr_ps[i] = gsl_vector_calloc(Neq);
    
    /* memory allocation for the vector of indices */
    pasr->pasr_idx = (int *) malloc(Np*sizeof(int));
    if ( pasr->pasr_idx == NULL )
    {
	for( i = 0; i < Np; i++ )
	{
	    gsl_vector_free(pasr->pasr_ps[i]);
	    pasr->pasr_ps[i] = NULL;
	}
	
        free(pasr->pasr_ps);
	pasr->pasr_ps = NULL;
        
	return GSL_ENOMEM;
    }
    
    /* initial configuration of particles pairwize */
    for( i = 0; i < Np; i++ )
        pasr->pasr_idx[i] = i;
    
    return GSL_SUCCESS;
}
/*----------------------------------------------------------------*/




/*
* -----------------------------------------------------------------
*   pasr_free
*
*   This routine frees the memory allocated by pasr_alloc.
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

void pasr_free(int Np,void **pasr_bl)
{
    int i;
    pasr_wrk *pasr = NULL;
    
    /* checking if pasr_bl is NULL */
    if ( *pasr_bl == NULL )
        return;
    
    /* checking if Np is positive */
    if ( Np <= 0 )
        return;
    
    pasr = (pasr_wrk *) (*pasr_bl);
    
    /* releasing pasr elements */
    if ( pasr->pasr_ps != NULL )
    {
	for( i = 0; i < Np; i++ )
	    if ( pasr->pasr_ps[i] != NULL )
	    {
	        gsl_vector_free(pasr->pasr_ps[i]);
	        pasr->pasr_ps[i] = NULL;
	    }
	/* releasing pasr struct */
	free(pasr->pasr_ps);
	pasr->pasr_ps = NULL;
    }
    
    if ( pasr->pasr_idx != NULL )
    {
	free(pasr->pasr_idx);
	pasr->pasr_idx = NULL;
    }
    
    free(*pasr_bl);
    *pasr_bl = NULL;
    
    return;
}
/*----------------------------------------------------------------*/




/*
* -----------------------------------------------------------------
*  pasr_set_all
*
*  This function sets the system of particles with a given
*  composition phi.
*
*   Input:
*   pasr - pointer to a struct pasr
*   phi  - composition
*
*   Output:
*   success or error
*
*  last update: Nov 4, 2009
* -----------------------------------------------------------------
*/

int pasr_set_all(int Np,gsl_vector *phi,pasr_wrk *pasr)
{
    int i;
    
    /* check if pasr is NULL */
    if ( pasr == NULL )
        return GSL_EINVAL;
    
    /* set all particles to initial phi */
    for( i = 0; i < Np; i++ )
        gsl_vector_memcpy(pasr->pasr_ps[i],phi);
    
    return GSL_SUCCESS;
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
*   last update: Nov 4, 2009
* -----------------------------------------------------------------
*/

int pasr_Nexc(int Np,double delta_t,double tau_res)
{   
    return ceil( (Np/2)*(delta_t/tau_res) );
}
/*----------------------------------------------------------------*/




/*
* -----------------------------------------------------------------
*   pasr_meanvar
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
*   last update: Jan 5, 2010
* -----------------------------------------------------------------
*/

void pasr_meanvar(int Np,int k,gsl_vector **ps,double *mean,double *var)
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
*   pasr_iem
*
*   This function computes the IEM mixture step in PaSR model.
*
*   Input:
*   Np      - # of particles
*   Neq     - # of equations
*   delta_t - time step
*   tau_mix - mixture time
*   phi_m   - mean value of phi
*   idx     - particles index vector
*
*   Output:
*   ps      - system of particles
*
*   last update: Nov 4, 2009
* -----------------------------------------------------------------
*/

void pasr_iem(int Np,int Neq,double delta_t,double tau_mix,
			    double phi_m,int *idx,gsl_vector **ps)
{
    int k, i;
    double dif;
    
    /* particle index */
    k = 0;
    
    /* mixture of the Np particles */
    while( k < Np )
    {
        for( i = 0; i < Neq; i++ )
        {
            /* dif = phi_0 - phi_m */
            dif = ps[idx[k]]->data[i] - phi_m;
            
            /* phi =  phi_m + dif.exp(- tau / delta_t) */
            ps[idx[k]]->data[i] = phi_m + dif*exp(- delta_t / tau_mix);
        }
        k++;
    }
    
    return;
}
/*----------------------------------------------------------------*/




/*
* -----------------------------------------------------------------
*   pasr_iops
*
*   This function executes the inflow/outflow,
*   and shuflle steps.
*
*   Input:
*   Np      - # of particles
*   Nexc    - # of particles to be exchanged
*   rand    - rng workspace
*   phi     - inflow composition
*
*   Output:
*   p       - system of particles
*   pairs   - pairs of particles index
*
*   last update: Feb 19, 2010
* -----------------------------------------------------------------
*/

void pasr_iops(int Np,int Nexc,gsl_rng *rand,
				    gsl_vector *phi,pasr_wrk *pasr)
{
    int i;
    int idx1[Nexc];
    int idx2[2*Nexc];
    int *aux = NULL;
    
    /* memory allocation */
    aux = (int *) malloc((Np/2)*sizeof(int));
    
    /* setting aux vector with integers from 0 to Np/2-1 */
    for( i = 0; i < Np/2; i++ )
        aux[i] = i;
    
    /* selecting Nexc pairs at random */
    gsl_ran_choose(rand,idx1,Nexc,aux,Np/2,sizeof(int));
    for( i = 0; i < Nexc; i++ )
    {
        idx2[2*i  ] = pasr->pasr_idx[2*idx1[i]  ];
        idx2[2*i+1] = pasr->pasr_idx[2*idx1[i]+1];
    }
    
    
    /* changing composition of Nexc pairs (inflow / outflow) */
    for( i = 0; i < Nexc; i++ )
    {
        gsl_vector_memcpy(pasr->pasr_ps[idx2[2*i]  ],phi);
        gsl_vector_memcpy(pasr->pasr_ps[idx2[2*i+1]],phi);
    }
    
    /* shuffling the selected pairs */
    gsl_ran_shuffle(rand,idx2,2*Nexc,sizeof(int));
    
    /* setting the new partners in the system */
    for( i = 0; i < Nexc; i++ )
    {
        pasr->pasr_idx[2*idx1[i]  ] = idx2[2*i  ];
        pasr->pasr_idx[2*idx1[i]+1] = idx2[2*i+1];
    }
    
    /* releasing allocated memory */
    free(aux);
    
    return;
}
/*----------------------------------------------------------------*/
