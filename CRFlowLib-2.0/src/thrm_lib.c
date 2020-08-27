
/*
* -----------------------------------------------------------------
*  Thermochemistry Library --- thrm_lib.c
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
*  This is the implementation file of a library
*  with thermochemistry routines.
* -----------------------------------------------------------------
*/



#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>
#include <gsl/gsl_errno.h>
#include <gsl/gsl_math.h> 

#include "../include/cklib.h"
#include "../include/ckthrm.h"
#include "../include/thrm_lib.h"




/*
*------------------------------------------------------------
*   thrm_title
*
*   This function prints the program title on the screen
*
*   last update: Oct 20, 2009
* -----------------------------------------------------------------
*/

void thrm_title()
{
    time_t curtime = time (NULL);
    struct tm *loctime = localtime (&curtime);
    
    printf("\n");
    printf("=============================================");
    printf("\n Themochemistry Direct Simulation\n");
    printf("\n Author:");
    printf("\n Americo Barbosa da Cunha Junior --- PUC-Rio");
    printf("\n americo.cunhajr@gmail.com\n");
    printf("\n Date:");
    printf("\n %s", asctime(loctime) );
    printf("=============================================");
    printf("\n");
    
    return;
}




/*
*------------------------------------------------------------
*   thrm_mech
*
*   This function prints on the screen informations about the
*   chemical kinetics mechanism
*
*   Input:
*   n_e  - # of elements
*   n_s  - # of species
*   n_r  - # of reactions
*   mech - mechanism name
*
*   last update: Jul 22, 2009
*------------------------------------------------------------
*/

void thrm_mech(int n_e,int n_s,int n_r,char *mech)
{    
    printf("\n");
    printf("\n Thermochemistry Mechanism: %s", mech);
    printf("\n # of elements  = %d", n_e);
    printf("\n # of species   = %d", n_s);
    printf("\n # of reactions = %d", n_r);
    printf("\n");
    printf("\n");
    
    return;
}
/*------------------------------------------------------------*/




/*
* -----------------------------------------------------------------
*   thrm_input
*
*   This function receives the pmsr parameters
*   and returns GSL_SUCCESS if there is no error.
*
*   Input:
*   Np       - # of particles
*   k        - # of discrete time steps
*   delta_t  - time step
*   tau_res  - residence time
*   tau_mix  - mixture time
*   tau_pair - pairwize time
*
*   Output:
*   success or error
*
*   last update: Feb 5, 2009
* -----------------------------------------------------------------
*/

int thrm_input(int *steps,double *t0,double *delta_t,
                            double *Tmax,double *Tmin,double *tol)
{
    printf("\n Input simulation parameters:\n");
    
    printf("\n # of integration steps:");
    scanf("%d", steps);
    printf("\n %d\n", *steps);
    if( *steps <= 0 )
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
    
    printf("\n tolerance:");
    scanf("%lf", tol);
    printf("\n %+.1e\n", *tol);
    if( *tol <= 0.0 )
	return GSL_EINVAL;
    
    return GSL_SUCCESS;
}
/*----------------------------------------------------------------*/




/*
*------------------------------------------------------------
*   thrm_alloc
*
*   This function allocates a struct thermochemistry.
*
*   Output:
*   thrm - pointer to a struct thermochemistry
*
*   last update: May 10, 2009
*------------------------------------------------------------
*/

thrm_wrk* thrm_alloc()
{
    thrm_wrk *thrm = NULL;
    
    /* memory allocation for ck_wrk struct */
    thrm = (thrm_wrk *) malloc(sizeof(thrm_wrk));
    if ( thrm == NULL )
        return NULL;
    
    /* setting ck_wrk elements equal NULL and 0.0 */
    thrm->iwk    = NULL;
    thrm->rwk    = NULL;
    thrm->cwk    = NULL;
    thrm->leniwk = 0;
    thrm->lenrwk = 0;
    thrm->lencwk = 0;
    thrm->n_e    = 0;
    thrm->n_s    = 0;
    thrm->n_r    = 0;
    
    return thrm;
}
/*------------------------------------------------------------*/




/*
*------------------------------------------------------------
*   thrm_free
*
*   This function frees the memory used by CK work vectors
*
*   Input:
*   thrm - pointer to a struct chemkin
*
*   Output:
*   void
*
*   last update: Jul 21, 2009
*------------------------------------------------------------
*/

void thrm_free(void **thrm_bl)
{
    thrm_wrk* thrm = NULL;
    
    /* checking if thrm_bl is NULL */
    if ( *thrm_bl == NULL )
        return;
    
    thrm = (thrm_wrk *) (*thrm_bl);
    
    /* releasing ck work vectors */
    if( thrm->iwk != NULL )
    {
        free(thrm->iwk);
        thrm->iwk = NULL;
    }
    
    if( thrm->rwk != NULL )
    {
        free(thrm->rwk);
        thrm->rwk = NULL;
    }
    
    if( thrm->cwk != NULL )
    {
        free(thrm->cwk);
        thrm->cwk = NULL;
    }
    
    /* releasing thrm struct */
    free(*thrm_bl);
    *thrm_bl = NULL;
    
    return;
}
/*------------------------------------------------------------*/




/*
*------------------------------------------------------------
*   thrm_init
*
*   This function initiates chemkin work vectors and the
*   mechanism parameters.
*
*   Input:
*   thrm - pointer to a struct thermochemistry
*
*   Output:
*   success or error
*
*   last update: Feb 19, 2010
*------------------------------------------------------------
*/

int thrm_init(thrm_wrk *thrm)
{
    if ( thrm == NULL )
        return GSL_EINVAL;
    
    /* calling ck_wrk_len subroutine */
    ck_wrk_len(&(thrm->leniwk),&(thrm->lenrwk),&(thrm->lencwk));
    
    /* memory allocation for ck work vectors */
    thrm->iwk =    (int *) malloc(   (thrm->leniwk)*sizeof(int));
    thrm->rwk = (double *) malloc(   (thrm->lenrwk)*sizeof(double));
    thrm->cwk =   (char *) malloc(16*(thrm->lencwk)*sizeof(char));
    if ( thrm->iwk == NULL || thrm->rwk == NULL || thrm->cwk == NULL )
    {
        thrm_free((void **)&thrm);
        return GSL_ENOMEM;
    }
    
    /* calling ck_wrk_init subroutine */
    ck_wrk_init(&(thrm->leniwk),&(thrm->lenrwk),&(thrm->lencwk),
        thrm->iwk,thrm->rwk,thrm->cwk);
    
    /* calling ck_wrk_indx subroutine */
    ck_wrk_indx(thrm->iwk,thrm->rwk,&(thrm->n_e),&(thrm->n_s),&(thrm->n_r));
    
    return GSL_SUCCESS;
}
/*------------------------------------------------------------*/




/*
*------------------------------------------------------------
*   thrm_composition
*
*   This function sets composition vector with values.
*
*   Input:
*   thrm - pointer to a struct thermochemistry
*
*   Output:
*   phi - composition vector
*   success or error
*
*   last update: Mar 2, 2010
*------------------------------------------------------------
*/

void thrm_composition(thrm_wrk *thrm,gsl_vector *phi)
{   
    /* calling ck_wrk_comp subroutine */
    ck_composition(&(thrm->n_s),thrm->iwk,thrm->rwk,thrm->cwk,phi->data);
    
    return;
}
/*------------------------------------------------------------*/




/*
*------------------------------------------------------------
*   thrm_page_heading0
*
*   This function prints the properties names on the screen.
*
*   Input:
*   thrm - pointer to a struct thermochemistry
*
*   last update: Nov 6, 2009
*------------------------------------------------------------
*/

void thrm_page_heading0(FILE *file)
{    
    fprintf(file," t(s)      ");
    fprintf(file," T(K)      ");
    fprintf(file," h         ");
    fprintf(file," p         ");
    fprintf(file," CO2       ");
    fprintf(file," O2        ");
    fprintf(file," CO        ");
    fprintf(file," O         ");
    
    return;
}
/*------------------------------------------------------------*/




/*
*------------------------------------------------------------
*   thrm_h2T
*
*   This function computes the mixture temperature for a 
*   given enthalpy, pressure and species mass fraction.
*
*   Input:
*   ck  - pointer to a struct chemkin
*   phi - composition vector
*
*   Output:
*   T   - mixture temperature
*
*   last update: Mar 2, 2010
*------------------------------------------------------------
*/

double thrm_h2T(thrm_wrk *thrm,double *phi,double Tmax,double Tmin,double tol)
{
    double T = 0.0;
 
    ck_h2t_bisec(&(thrm->n_s),phi,&Tmax,&Tmin,&tol,thrm->iwk,thrm->rwk,&T);
    /*ck_h2t_newt(&(thrm->n_s),phi,&tol,thrm->iwk,thrm->rwk,&T,&T);*/
    
    return T;
}
/*------------------------------------------------------------*/





/*
* -----------------------------------------------------------------
*   thrm_temp_meanvar
*
*   This function computes the temperature mean and variance
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
*   Np   - # of particles
*   ps   - system of particle
*   thrm - pointer to a struct thermochemistry
*   Tmax - temperature upper bound
*   Tmin - temperature lower bound
*   tol  - tolerance
*   mean - temperature mean
*   var  - temperature variance
*
*   Output:
*   success or error
*
*   last update: Mar 4, 2010
* -----------------------------------------------------------------
*/

void thrm_temp_meanvar(int Np,gsl_vector **ps,thrm_wrk *thrm,
        double Tmax,double Tmin,double tol,double *mean,double *var)
{
    int i;
    double T     = 0.0;
    double S     = 0.0;
    double M     = 0.0;
    double delta = 0.0;
    
    for ( i = 0; i < Np; i++ )
    {
        T     = thrm_h2T(thrm,ps[i]->data,Tmax,Tmin,tol);
	delta = T - M;
	M    += delta/(double) (i+1);
        S    += delta*(T - M);
    }
    
    /* mean */
    *mean  = M;
    
    /* variance */
    *var   = S/(double) (i-1);
    
    return;
}
/*----------------------------------------------------------------*/




/*
*------------------------------------------------------------
*   thrm_rrsum
*
*   This function computes the reaction rates sum given the
*   reaction rate vector.
*
*   Input:
*   phi - composition vector
*
*   Output:
*   sum - reaction rates sum
*
*   last update: Sep 28, 2009
*------------------------------------------------------------
*/

double thrm_rrsum(gsl_vector *phi)
{
    unsigned int i;
    double sum = 0.0;
 
    for (i = 2; i < phi->size; i++)
        sum += phi->data[i];
    
    return sum;
}
/*------------------------------------------------------------*/




/*
*------------------------------------------------------------
*   thrm_mfsum
*
*   This function computes the mass fraction sum given a
*   composition vector.
*
*   Input:
*   phi - composition vector
*
*   Output:
*   sum - mass fraction sum
*
*   last update: Sep 24, 2009
*------------------------------------------------------------
*/

double thrm_mfsum(gsl_vector *phi)
{
    unsigned int i;
    double sum = 0.0;
 
    for (i = 2; i < phi->size; i++)
        sum += phi->data[i];
    
    return sum;
}
/*------------------------------------------------------------*/




/*
*------------------------------------------------------------
*   thrm_mfsign
*
*   This function checks if the mass fraction
*   components are all positive in a composition vector.
*
*   Input:
*   phi  - composition vector
*
*   Output:
*   flag - success or fail
*
*   last update: Sep 28, 2009
*------------------------------------------------------------
*/

int thrm_mfsign(gsl_vector *phi)
{
    unsigned int i;
    int flag = GSL_SUCCESS;
 
    for (i = 2; i < phi->size; i++)
        if ( phi->data[i] < 0.0 )
            flag = GSL_EFAILED;
    
    return flag;
}
/*------------------------------------------------------------*/




/*
*------------------------------------------------------------
*   thrm_eqs
*
*   This function defines the set of governing diferential
*   equations, one for enthalpy, one for pressure and ns
*   for species mass fractions.
*
*   Input:
*   t         - instant time
*   phi       - composition vector
*   Sphi      - reaction rate vector
*   thrm_data - pointer to external data
*
*   Output:
*   success or error
*
*   last update: Feb 19, 2010
*------------------------------------------------------------
*/

int thrm_eqs(double t,N_Vector phi,N_Vector Sphi,void *thrm_data)
{
    int i, Neq;
    double tol = 1.0e-4;
    double Tmax = 5000.0, Tmin = 100.0;
    double p, T, rho;
    double *phi_data  = NULL;
    double *Sphi_data = NULL;
    double *wdot      = NULL;
    double *wt        = NULL;
    thrm_wrk *thrm    = NULL;
    
    /* initiating ck workspace */
    thrm = (thrm_wrk *)thrm_data;
    if ( thrm == NULL )
        return GSL_EINVAL;
    
    /* obtaining phi and Rphi length and components */
    Neq       = NV_LENGTH_S(phi);
    phi_data  = NV_DATA_S(phi);
    Sphi_data = NV_DATA_S(Sphi);
    
    /* memory allocation for species molar mass and reaction rate vectors */
    wt   = (double *) malloc((Neq-2)*sizeof(double));
    wdot = (double *) malloc((Neq-2)*sizeof(double));
    if( wt == NULL || wdot == NULL )
    {
        free(wt);
        free(wdot);
        return GSL_ENOMEM;
    }
    
    /* pressure (dynes/cm^2) */
    p = phi_data[1];
    
    /* mixture temperature (K) */
    T = thrm_h2T(thrm,phi_data,Tmax,Tmin,tol);
    if( T < tol )
        return GSL_EFAILED;
    
    /* mixture density (g/cm^3)*/
    ckrhoy(&p,&T,&(phi_data[2]),thrm->iwk,thrm->rwk,&rho);
    
    /* species molar weight (g/mol) */
    ckwt(thrm->iwk,thrm->rwk,wt);
    
    /* species molar production rates (mol/(cm^3.s)) */
    ckwyp(&p,&T,&(phi_data[2]),thrm->iwk,thrm->rwk,wdot);
    
    /* dh/dt = 0 (isoenthalpic) */
    Sphi_data[0] = 0.0;
    
    /* dp/dt = 0 (isobaric) */
    Sphi_data[1] = 0.0;
    
    /* dY[i]/dt = wdot[i-2].wt[i-2]/rho */
    for ( i = 2; i < Neq; i++ )
        Sphi_data[i] = wdot[i-2]*wt[i-2]/rho;
    
    /* releasing allocated memory */
    free(wdot);
    free(wt);
    phi_data  = NULL;
    Sphi_data = NULL;
    wdot      = NULL;
    wt        = NULL;
    thrm      = NULL;
    
    return GSL_SUCCESS;
}
/*------------------------------------------------------------*/



/*
*------------------------------------------------------------
*   thrm_eqs2
*
*   This function defines the set of governing diferential
*   equations, one for temperature and ns
*   for species mass fractions.
*
*   Input:
*   t         - instant time
*   phi       - state vector
*   Sphi      - reaction rate vector
*   thrm_data - pointer to external data
*
*   Output:
*   success or error
*
*   last update: Dec 28, 2010
*------------------------------------------------------------
*/

int thrm_eqs2(double t,N_Vector phi,N_Vector Sphi,void *thrm_data)
{
    int i, Neq;
    double T, p, rho, cpb, sum = 0.0;
    double *phi_data  = NULL;
    double *Sphi_data = NULL;
    double *wdot      = NULL;
    double *wt        = NULL;
    double *h         = NULL;
    thrm_wrk *thrm    = NULL;
    
    /* initiating ck workspace */
    thrm = (thrm_wrk *)thrm_data;
    if ( thrm == NULL )
        return GSL_EINVAL;
    
    /* obtaining phi and Rphi length and components */
    Neq       = NV_LENGTH_S(phi);
    phi_data  = NV_DATA_S(phi);
    Sphi_data = NV_DATA_S(Sphi);
    
    /* memory allocation for species molar mass and reaction rate vectors */
    wt   = (double *) malloc((Neq-2)*sizeof(double));
    wdot = (double *) malloc((Neq-2)*sizeof(double));
    h    = (double *) malloc((Neq-2)*sizeof(double));
    if( wt == NULL || wdot == NULL || h == NULL)
    {
        free(wt);
        free(wdot);
	free(h);
        return GSL_ENOMEM;
    }

    /* mixture temperature (K) */
    T = phi_data[0];
    
    /* pressure (dynes/cm^2) */
    p = phi_data[1];
    
    /* mixture density (g/cm^3)*/
    ckrhoy(&p,&T,&(phi_data[2]),thrm->iwk,thrm->rwk,&rho);

    /* mean specific heat at constant pressure (ergs/g.K) */
    ckcpbs(&T,&(phi_data[2]),thrm->iwk,thrm->rwk,&cpb);
    
    /* species molar weight (g/mol) */
    ckwt(thrm->iwk,thrm->rwk,wt);
    
    /* species molar production rates (mol/(cm^3.s)) */
    ckwyp(&p,&T,&(phi_data[2]),thrm->iwk,thrm->rwk,wdot);
    
    /* species specific enthalpy (ergs/g) */
    ckhms(&T,thrm->iwk,thrm->rwk,h);
    
    /* dY[i+2]/dt = wdot[i].wt[i]/rho */
    for ( i = 0; i < thrm->n_s; i++ )
    {
        Sphi_data[i+2] = wdot[i]*wt[i]/rho;
	sum += h[i]*wdot[i]*wt[i];
    }

    /* dT/dt =  -1/(rho.cpb) . sum_i=1^ns h[i].wdot[i].wt[i] */
    Sphi_data[0] = -sum/(rho*cpb);
    
    /* dp/dt = 0 (isobaric) */
    Sphi_data[1] = 0.0;
    
    
    /* releasing allocated memory */
    free(wdot);
    free(wt);
    free(h);
    phi_data  = NULL;
    Sphi_data = NULL;
    wdot      = NULL;
    wt        = NULL;
    h         = NULL;
    thrm      = NULL;
    
    return GSL_SUCCESS;
}
/*------------------------------------------------------------*/
