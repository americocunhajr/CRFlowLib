
/*
* -----------------------------------------------------------------
*  ode_lib.c
*  ODE Solver Library
*  Version: 2.0
*  Last Update: Oct 1, 2019
* 
*  This is the implementation of a computational library with
*  numerical tools for Ordinary Differential Equations (ODE).
* ----------------------------------------------------------------- 
*  Programmer: Americo Barbosa da Cunha Junior
*              americo.cunhajr@gmail.com
* -----------------------------------------------------------------
*  Copyright (c) 2019 by Americo Barbosa da Cunha Junior
* -----------------------------------------------------------------
*/




#include <stdlib.h>
#include <math.h>
#include <float.h>

#include <cvode/cvode.h>               /* prototypes for CVODE functions and const */
#include <nvector/nvector_serial.h>    /* access to serial NVector                 */
#include <sunmatrix/sunmatrix_dense.h> /* access to dense SUNMatrix                */
#include <sunlinsol/sunlinsol_dense.h> /* access to dense SUNLinearSolver          */
#include <sundials/sundials_types.h>   /* defs. of realtype, sunindextype          */

#include <gsl/gsl_blas.h>
#include <gsl/gsl_math.h>
#include <gsl/gsl_errno.h>

#include "../include/ode_lib.h"



/*
*------------------------------------------------------------
*   uround
*
*   This function computes the machine unit roundoff
*   defined as the smallest u such that 1.0 + u > 1.0
*
*   Output:
*   uround - machine unit roundoff
*
*   last update: May 20, 2009
*------------------------------------------------------------
*/

double uround(void)
{
    double u = 1.0;
    double comp = 0.0;

    while(comp != 1.0)
    {
	u *= 0.5;
	comp = 1.0 + u;
    }
    
    return 2.0*u;
}
/*------------------------------------------------------------*/




/*
*------------------------------------------------------------
*   wnorm
*
*   This function computes v1 and v2 vector weight norm
*   defined as
*   wnorm = sqrt( (1/n)*(sum (v1[i]*v2[i])^2) )
*
*   Input:
*   n     - vectors dimension
*   v     - vector 1
*   v     - vector 2
*
*   Output:
*   wnorm - weight norm
*
*   last update: May 20, 2009
*------------------------------------------------------------
*/

double wnorm(int n,double *v1,double *v2)
{
    int i;
    double wnorm = 0.0;
    
    for ( i = 0; i < n; i++ )
	wnorm += v1[i]*v1[i]*v2[i]*v2[i];
    
    return sqrt( wnorm / (double)n );
}
/*------------------------------------------------------------*/




/*
*------------------------------------------------------------
*   ewset
*
*   This function creates a error weight vector
*   defined as ewv[i] = rtol*abs(v[i]) + atol
*
*   Input:
*   n    - vector dimension
*   v    - vector
*   atol - absolute tolerance
*   rtol - relative tolerance
*
*   Output:
*   ewv  - error weight vector
*
*   last update: May 22, 2009
*------------------------------------------------------------
*/

void ewtset(int n,double *v,double atol,double rtol,double *ewt)
{
    int i;
    
    for ( i = 0; i < n; i++ )
        ewt[i] =  rtol*fabs(v[i]) + atol;
    
    return;
}
/*------------------------------------------------------------*/




/*------------------------------------------------------------
*   jacobian
*
*   This function computes the jacobian matrix of ydot = f(y,t)
*   defiend as  Jij = dfi/dyj
*
*   Input:
*   n      - vector dimension
*   f_data - pointer to external data
*   Fy     - reaction mapping
*   y      - composition vector
*   t      - time step
*   atol   - absolute tolerance
*   rtol   - relative tolerance
*
*   Output:
*   J       - jacobian matrix
*   success or error
*
*   last update: Oct 7, 2009
------------------------------------------------------------*/

int jacobian(CVRhsFn f,void *f_data,gsl_vector *Fy,gsl_vector *y,
                    double t,double atol,double rtol,gsl_matrix *J)
{
    unsigned int i, j, N1, N2;
    double roundoff, min_inc_mult;
    double fnorm, minInc, inc, inc_inv, yjsaved, srur;
    double *yd_data  = NULL;
    double *Fyd_data = NULL;
    N_Vector yd      = NULL;
    N_Vector Fyd     = NULL;
    gsl_vector *ewt  = NULL;
    
    
    min_inc_mult = 1.0e3;
    
    /* vectors dimensions */
    N1 = y->size;
    N2 = Fy->size;
    
    /* error weight vector */
    ewt = gsl_vector_calloc(N1);
    
    /* computing the unit roundoff */
    roundoff = uround();
    
    /* memory allocation and startup of yd and Fyd */
    yd  = N_VNew_Serial(N1);
    Fyd = N_VNew_Serial(N2);
    if ( yd == NULL || Fyd == NULL )
        return GSL_ENOMEM;
    
    /* obtaining yd and Fyd compoments  */
    yd_data  = NV_DATA_S(yd);
    Fyd_data = NV_DATA_S(Fyd);
    
    /* setting yd equal y vector */
    for ( i = 0; i < N1; i++ )
        yd_data[i] = y->data[i];
    
    /* setting Fyd equal the null vector */
    for ( i = 0; i < N2; i++ )
        Fyd_data[i] = 0.0;
    
    /* defing error weight vector */
    ewtset(N1,y->data,atol,rtol,ewt->data);
    
    /* computing weight norm */
    fnorm = wnorm(N2,Fy->data,ewt->data);
    
    /* square root of the machine roundoff */
    srur = sqrt(roundoff);
    
    /* computing disturbance parameter */
    minInc = (fnorm != 0.0) ?
	    (min_inc_mult*fabs(t)*roundoff*((double)N1)*fnorm) : 1.0;
    
    for ( j = 0; j < N1; j++ )
    {
        /* saving y[j] value */
        yjsaved = yd_data[j];
        
        /* disturbance */
        inc = GSL_MAX(srur*fabs(yjsaved),minInc/ewt->data[j]);
    	
    	/* disturbing y */
        yd_data[j] += inc;
	
        /* computing Fy disturbed */
        f(t,yd,Fyd,f_data);
	
	/* restoring yd[j] original value */
        yd_data[j] = yjsaved;
	
	/* computing the step */
        inc_inv = 1.0/inc;
        
        /* computing jacobian matrix column j */
    	for ( i = 0; i < N2; i++ )
            J->data[i*J->tda+j] = inc_inv*(Fyd_data[i] - Fy->data[i]);
    }
    
    /* releasing allocated memory */
    N_VDestroy_Serial(yd);
    N_VDestroy_Serial(Fyd);
    gsl_vector_free(ewt);
    yd       = NULL;
    Fyd      = NULL;
    yd_data  = NULL;
    Fyd_data = NULL;
    ewt      = NULL;
   
    return GSL_SUCCESS;
}
/*------------------------------------------------------------*/




/*------------------------------------------------------------
*   gradient
*
*   This function computes the gradient matrix of R(phi)
*   defiend as: Aij = d R_i/d phi_j
*
*   Input:
*   cvode_mem - ODE solver workspace
*   t0        - initial time
*   delta_t   - time step
*   phi     - composition vector
*   Rphi    - reac tion mapping
*
*   Output:
*   A       - gradient matrix
*   success or error
*
*   last update: Oct 1, 2019
------------------------------------------------------------*/

int gradient(void *cvode_mem,double t0,double delta_t,
             gsl_vector *phi,gsl_vector *Rphi,gsl_matrix *A)
{
    unsigned int i, j, flag;
    double inc, inc_inv, phijsaved, srur;
    gsl_vector *Rphid = NULL;
    
    /* memory allocation for Rphid */
    Rphid = gsl_vector_calloc(Rphi->size);
    
    /* square root of 1.0e6 times the machine roundoff */
    srur = sqrt(1.0e6*DBL_EPSILON);
    
    for ( j = 0; j < phi->size; j++ )
    {
        /* save phi_j value */
        phijsaved = phi->data[j];
        
        /* disturbance */
        inc = phijsaved*srur + srur;
    	
    	/* disturbe phi_j */
        phi->data[j] += inc;
	
        /* compute R(phi) disturbed */
	flag = odesolver_reinit(t0,phi,cvode_mem);
        if ( flag != GSL_SUCCESS )
            return flag;
        
	flag = odesolver(cvode_mem,delta_t,Rphid);
        if ( flag != GSL_SUCCESS )
            return flag;
	
	/* restore phi_j original value */
        phi->data[j] = phijsaved;
	
	/* compute the step */
        inc_inv = 1.0/inc;
        
        /* compute gradient matrix column j */
    	for ( i = 0; i < Rphi->size; i++ )
            A->data[i*A->tda+j] = inc_inv*(Rphid->data[i] - Rphi->data[i]);
    }
    
    
    /* release allocated memory */
    gsl_vector_free(Rphid);
    Rphid = NULL;
   
    return GSL_SUCCESS;
}
/*------------------------------------------------------------*/




/*
*------------------------------------------------------------
*   linear_approx
*
*   This function computes the linear approximation
*   for a vector function F(x) near the point x0.
*   
*   F(x) = F(x0) + DF(x0)*(x-x0)
*   
*   Input:
*   x    - vector x
*   x0   - vector x0
*   Fx0  - function at x0
*   DFx0 - jacobian matrix at x0
*
*   Output:
*   Fx   - linear approximation for F(x)
*
*   last update: Jun 9, 2009
*------------------------------------------------------------
*/

void linear_approx(gsl_vector *x,gsl_vector *x0,
                gsl_vector *Fx0,gsl_matrix *DFx0,gsl_vector *Fx)
{
    gsl_vector *aux = NULL;
    
    /* memory allocation */
    aux = gsl_vector_alloc(x0->size);

    /* aux := x */
    gsl_vector_memcpy(aux,x);
    
    /* aux := x - x0 */
    gsl_vector_sub(aux,x0);
    
    /* Fx := F(x0) */
    gsl_vector_memcpy(Fx,Fx0);
    
    /* Fx := F(x0) + DF(x0)*(x-x0) */
    gsl_blas_dgemv(CblasNoTrans,1.0,DFx0,aux,1.0,Fx);
    
    /* releasing allocated memory */
    gsl_vector_free(aux);
    aux = NULL;
    
    return;
}
/*------------------------------------------------------------*/




/*
*------------------------------------------------------------
*   odesolver_init
*
*   This function initiates the ODE solver workspace.
*
*   Input:
*   f         - right hand side function
*   f_data    - right hand side function data
*   t0        - initial time
*   mxsteps   - max # of solver iterations
*   atol      - absolute tolerance
*   rtol      - relative tolerance
*   x         - initial condition
*   cvode_mem - ODE solver workspace
*
*   Output:
*   success or error
*
*   last update: Oct 1, 2019
*------------------------------------------------------------
*/

int odesolver_init(CVRhsFn f,void *f_data,double t0,int mxsteps,
		   double atol,double rtol,gsl_vector *x,void *cvode_mem)
{
    int         maxnef = 100;           /* max. # of error test failures  */
    int           flag = CV_SUCCESS;    /* return value flag              */
    N_Vector        y0 = NULL;          /* initial condition vector       */
    SUNMatrix        A = NULL;          /* Matrix for linear solver use   */
    SUNLinearSolver LS = NULL;          /* dense linear solver object     */

    
    /* memory allocation for initial condition y0 */
    y0 = N_VMake_Serial(x->size,x->data);
    if ( y0 == NULL )
	return GSL_ENOMEM;

    /* memory allocation for CVODE workspace */
    flag = CVodeInit(cvode_mem,f,t0,y0);
    if( flag != CV_SUCCESS )
    	return GSL_ENOMEM;

    /* specifies scalar absolute and relative tolerances */
    flag = CVodeSStolerances(cvode_mem,rtol,atol);
    if( flag != CV_SUCCESS )
    	return flag;

    /* set user data for right hand side function */ 
    flag = CVodeSetUserData(cvode_mem,f_data);
    if( flag != CV_SUCCESS )
        return flag;
    
    /* Create dense SUNMatrix for use in linear solver */
    A = SUNDenseMatrix(x->size,x->size);
    if ( A == NULL )
	return GSL_ENOMEM;

    /* Create dense SUNLinearSolver object for use by CVode */
    LS = SUNLinSol_Dense(y0, A);
    if ( LS == NULL )
	return GSL_ENOMEM;

    /* Call CVodeSetLinearSolver to attach the matrix and linear solver to CVode */
    flag = CVodeSetLinearSolver(cvode_mem, LS, A);
    if( flag != CVLS_SUCCESS )
    	return flag;

    /* set solver max # of steps */
    flag = CVodeSetMaxNumSteps(cvode_mem,mxsteps);
    if( flag != CV_SUCCESS )
        return flag;    

    /* set max # of error test failures permitted */
    flag = CVodeSetMaxErrTestFails(cvode_mem,maxnef);
    if( flag != CV_SUCCESS )
        return flag;
    
    /* release allocated memory for y0 */
    N_VDestroy_Serial(y0);
    y0 = NULL;

    return GSL_SUCCESS;
}
/*------------------------------------------------------------*/




/*
*------------------------------------------------------------
*   odesolver_reinit
*
*   This function reinitiates the ODE solver workspace.
*
*   Input:
*   f         - right hand side function
*   f_data    - right hand side function data
*   t0        - initial time
*   x         - initial condition
*   cvode_mem - ODE solver workspace
*
*   Output:
*   success or error
*
*   last update: Oct 1, 2019
*------------------------------------------------------------
*/

int odesolver_reinit(double t0,gsl_vector *x,void *cvode_mem)
{
    int flag    = CV_SUCCESS;  /* return value flag              */
    N_Vector y0 = NULL;        /* initial condition vector       */
    
    /* memory allocation for y0 */
    y0 = N_VMake_Serial(x->size,x->data);
    if ( y0 == NULL )
	return GSL_ENOMEM;
    
    /* restart CVODE workspace */
    flag = CVodeReInit(cvode_mem,t0,y0);
    if( flag != CV_SUCCESS )
    	return flag;
    
    /* release allocated memory */
    N_VDestroy_Serial(y0);
    y0 = NULL;

    return GSL_SUCCESS;
}
/*------------------------------------------------------------*/




/*
*------------------------------------------------------------
*   odesolver
*
*   This function performs the direct integration of the
*   governing equations.
*
*   Input:
*   cvode_mem - ODE solver workspace
*   tf   - final time
*   Fx   - solution vector
*
*   Output:
*   success or error
*
*   last update: Oct 31, 2009
*------------------------------------------------------------
*/

int odesolver(void *cvode_mem,double tf,gsl_vector *Fx)
{
    double t;
    int flag   = CV_SUCCESS;
    N_Vector y = NULL;
    
    /* memory allocation and startup of y */
    y = N_VMake_Serial(Fx->size,Fx->data);
    if ( y == NULL )
	return GSL_ENOMEM;
    
    /* calling CVode solver */
    flag = CVode(cvode_mem,tf,y,&t,CV_NORMAL);
    if( flag != CV_SUCCESS )
        return flag;
    
    /* releasing allocated memory */
    N_VDestroy_Serial(y);
    y = NULL;

    return GSL_SUCCESS;
}
/*------------------------------------------------------------*/
