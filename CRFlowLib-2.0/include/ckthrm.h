
/*
* -----------------------------------------------------------------
*  ckthrm.h
*  Chemkin Thermochemistry Library
*  Version: 2.0
*  Last Update: Nov 3, 2019
*
*  Programmer: Americo Barbosa da Cunha Junior
*              americo.cunhajr@gmail.com
* -----------------------------------------------------------------
*  Copyright (c) 2010-2019, Americo Barbosa da Cunha Junior
*  All rights reserved.
* -----------------------------------------------------------------
*  This is the header file for CKTHRM module, a computational
*  library with routines to do the interface with Chemkin-II code.
* -----------------------------------------------------------------
*/




#ifndef __CKTHRM_H__
#define __CKTHRM_H__

/*
#define ck_wrk_len         ck_wrk_len__
#define ck_wrk_init        ck_wrk_init__
#define ck_wrk_indx        ck_wrk_indx__
#define ck_composition     ck_composition__
#define ck_species_heading ck_species_heading__
#define ck_normalize       ck_normalize__
#define ck_h2t_newt3       ck_h2t_newt3__
#define ck_h2t_newt        ck_h2t_newt__
#define ck_h2t_bisec       ck_h2t_bisec__ */


#define ck_wrk_len         ck_wrk_len_
#define ck_wrk_init        ck_wrk_init_
#define ck_wrk_indx        ck_wrk_indx_
#define ck_composition     ck_composition_
#define ck_species_heading ck_species_heading_
#define ck_normalize       ck_normalize_
#define ck_h2t_newt3       ck_h2t_newt3_
#define ck_h2t_newt        ck_h2t_newt_
#define ck_h2t_bisec       ck_h2t_bisec_


/*
*------------------------------------------------------------
* function prototypes
*------------------------------------------------------------
*/

void ck_wrk_len(int *leniwk,
                int *lenrwk,
                int *lencwk);

void ck_wrk_init(int *leniwk,
                 int *lenrwk,
                 int *lencwk,
                 int *iwrk,
                 double *rwrk,
                 char *cwrk);

void ck_wrk_indx(int *iwrk,
                 double *rwrk,
                 int *mm,
                 int *kk,
                 int *ii);

void ck_composition(int *kk,
                    int *iwrk,
                    double *rwrk,
                    char *cwrk,
                    double *phi);

void ck_species_heading(int *kk,
                        char *cwrk,
                        char *ksym);

void ck_normalize(int *nn,
                  double *vv);

void ck_h2t_newt(int *kk,
                 double *phi,
                 double *tol,
                 int *iwrk,
                 double *rwrk,
                 double *Tguest,
                 double *Tsafe);

void ck_h2t_bisec(int *kk,
                  double *phi,
                  double *Tmax,
                  double *Tmin,
                  double *tol,
                  int *iwrk,
                  double *rwrk,
                  double *Tsafe);


#endif /* __CKTHRM_H__ */
