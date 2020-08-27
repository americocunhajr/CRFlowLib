
/*
* -----------------------------------------------------------------
*  Chemkin Library --- cklib.h
*  Version: 1.6180
*  Date: Dec 29, 2010
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
*  This is the header file with the prototypes of Chemkin-II
*  package routines writte in Fortran 77.
* -----------------------------------------------------------------
*/




#ifndef __CK_LIB_H__
#define __CK_LIB_H__


#define ckindx ckindx_
#define ckinit ckinit_
#define cklen  cklen_
#define cksyms cksyms_
#define ckwt   ckwt_
#define ckrp   ckrp_
#define ckrhoy ckrhoy_
#define ckxty  ckxty_
#define ckytx  ckytx_
#define ckhbms ckhbms_
#define ckubms ckubms_
#define ckhms  ckhms_
#define ckums  ckums_
#define ckcpbs ckcpbs_
#define ckcvbs ckcvbs_
#define ckwyp  ckwyp_
#define cksnum cksnum_




/*
*------------------------------------------------------------
* function prototypes
*------------------------------------------------------------
*/

/* Initialization */
void ckindx(int *ICKWRK,double *RCKWRK,int *MM,int *KK,int *II,int *NFIT);

void ckinit(int *LENICK,int *LENRCK,int *LENCCK,int *LINC,int *LOUT,
            int *ICKWRK,double *RCKWRK,char *CCKWRK,int *IFLAG,int *LENCHAR);

void cklen(int* LINC,int* LOUT,int* LENI,int *LENR,int *LENC,int *IFLAG);

/* Information about Species */
void cksyms(char *CCKWRK,int *LOUT,char *KNAME,int *KERR,
                                            int *LENCHAR1,int *LENCHAR2);

void ckwt(int *ICKWRK,double *RCKWRK,double *WT);


/* Gas Constants and Units */
void ckrp(int *ICKWRK,double *RCKWRK,double *RU,double *RUC,double *PA);


/* Equation of State */
void ckrhoy(double *P,double *T,double *Y,int *ICKWRK,
                                            double *RCKWRK,double *RHO);


/* Mole-Mass Convertion */
void ckxty(double *X,int *ICKWRK,double *RCKWRK,double *Y);

void ckytx(double *Y,int *ICKWRK,double *RCKWRK,double *X);


/* Thermodynamic Properties (mass units) */
void ckhms(double *T,int *ICKWRK,double *RCKWRK,double *HMS);

void ckums(double *T,int *ICKWRK,double *RCKWRK,double *UMS);


/* Mean Thermodynamic Properties (mass units) */
void ckhbms(double *T,double *Y,int *ICKWRK,double *RCKWRK,double *HBMS);

void ckubms(double *T,double *Y,int *ICKWRK,double *RCKWRK,double *UBMS);
  
void ckhms(double *T,int *ICKWRK,double *RCKWRK,double *HMS);

void ckums(double *T,int *ICKWRK,double *RCKWRK,double *UMS);

void ckcpbs(double *T,double *Y,int *ICKWRK,double *RCKWRK,double *CPBMS);
  
void ckcvbs(double *T,double *Y,int *ICKWRK,double *RCKWRK,double *CVBMS);

/* Chemical Production Rates */
void ckwyp(double *P,double *T,double *Y,int *ICKWRK,
                                            double *RCKWRK,double *WDOT);


/* Utilities */
void cksnum(char *LINE,int *NEXP,int *LOUT,char *KRAY,int *NN,
    int *KNUM,int *NVAL,double *RVAL,int *KERR,int *LENCHAR1,int *LENCHAR2);


#endif /* __CK_LIB_H__ */
