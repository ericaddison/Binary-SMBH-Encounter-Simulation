/********************************************************************
*********************************************************************
 
	File: ODE.h
	Purpose: Main header file for ODE solvers
 
    Author: Eric Addison
    Initial Date: 12 May 2011
    Affiliation: USU Astrophysics


*********************************************************************/


#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <sys/time.h>
#include <errno.h>
#include <ctype.h>
#include "matrix.h"

/* Function Declarations */
double mid(double y[], int nvar, double xs, double h,double yout[],void (*f)(double, double[], double[]));
double mmid(double y[],int nvar, double xs, double htot, int nstep,double yout[],void (*f)(double, double[], double[]));
double RK4(double y[], double dy[], int nvar, double xs, double h,double yout[],void (*f)(double, double[], double[]));
double RK4_adaptive(double y[], int nvar, double xs, double *hh, double yout[],double tol,double hmin,void (*f)(double, double[], double[])); 
void RK4A_const_step(double y[], int nvar, double *xs, double htot, double *htry,double eps, double hmin, void (*f)(double, double[], double[]));
double bsstep(double y[], int nvar, double *xs, double *htry,double eps, double hmin, double yscal[], int method, void (*f)(double, double[], double[]));
void bs_const_step(double y[], int nvar, double *xs, double htot, double *htry,double eps,double hmin, double yscal[], int method, void (*f)(double, double[], double[]));
void pzextr(int iest, double xest, double yest[], double yz[], double dy[], int nvar);
void rzextr(int iest, double xest, double yest[], double yz[], double dy[], int nvar);
double mleap(double y[], int nvar, double xs, double htot, double nstep, double yout[],void (*F)(double, double[], double[]));
void scale_vars(double s, double h, int nvar, double y[], double yscal[], void (*F)(double, double[], double[]));

/* macros */
#define SIGN(X) ( X >= 0 ? X : -X )
#define SQR(X) (X*X)
#define CUBE(X) (X*X*X)
#define FMAX(X,Y) ( X >= Y ? X : Y )
#define FMIN(X,Y) ( X <= Y ? X : Y )

/* constants used by Burlisch-Stoer method */
#define KMAXX 8					//Maximum row number used in the extrapolation.
#define IMAXX (KMAXX+1)
#define SAFE1 0.25				//Safety factors.
#define SAFE2 0.7
#define REDMAX 1.0e-5			//Maximum factor for stepsize reduction.
#define REDMIN 0.7				//Minimum factor for stepsize reduction.
#define TINY 1.0e-30			//Prevents division by zero.
#define SCALMX 0.1				// 1/SCALMX is the maximum factor by which a stepsize can be increased.
double **dBS, *xBS;				//pointers used by extrapolation routines for Burilsch-Stoer

/* constants used by Adaptive RK4 */
#define S1 0.9
#define S2 4.0

