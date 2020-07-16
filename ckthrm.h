
/*
* -----------------------------------------------------------------
*  Chemkin Thermochemistry Library --- ckthrm.h
*  Version: 1.6180
*  Date: Mar 4, 2010
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
*  This is the header file for a library to initialize
*  Chemkin workspace.
* -----------------------------------------------------------------
*/




#ifndef __CKTHRM_H__
#define __CKTHRM_H__


#define ck_wrk_len         ck_wrk_len__
#define ck_wrk_init        ck_wrk_init__
#define ck_wrk_indx        ck_wrk_indx__
#define ck_composition     ck_composition__
#define ck_species_heading ck_species_heading__
#define ck_normalize       ck_normalize__
#define ck_h2t_newt3       ck_h2t_newt3__
#define ck_h2t_newt        ck_h2t_newt__
#define ck_h2t_bisec       ck_h2t_bisec__



/*
*------------------------------------------------------------
* function prototypes
*------------------------------------------------------------
*/

void ck_wrk_len(int *leniwk,int *lenrwk,int *lencwk);

void ck_wrk_init(int *leniwk,int *lenrwk,int *lencwk,
                    int *iwrk,double *rwrk,char *cwrk);

void ck_wrk_indx(int *iwrk,double *rwrk,int *mm,int *kk,int *ii);

void ck_composition(int *kk,int *iwrk,double *rwrk,
                                    char *cwrk, double *phi);

void ck_species_heading(int *kk,char *cwrk,char *ksym);

void ck_normalize(int *nn,double *vv);

void ck_h2t_newt(int *kk, double *phi, double *tol,
                int *iwrk, double *rwrk, double *Tguest, double *Tsafe);

void ck_h2t_bisec(int *kk,double *phi,double *Tmax,double *Tmin,
                double *tol,int *iwrk,double *rwrk, double *Tsafe);


#endif /* __CKTHRM_H__ */
