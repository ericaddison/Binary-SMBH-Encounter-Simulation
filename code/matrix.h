/********************************************************************
*********************************************************************
 
	File: matrix.h
	Purpose: Header file for Matrix utilities
 
	Author: Eric Addison
	Initial Date: Summer 2009
	Affiliation: Utah State University

*********************************************************************/

#include <stdio.h>
#include <stdlib.h>
#include <math.h>

//vector routines
void vect_mult_element(double * x, double * y, double * prod, int size);
void vect_add(double * x, double * y, double * sum, int size);
void vect_mult_scalar(double * x, double a, double * result, int size);
void vect_add2(double * x, double * y, double * sum, int size);
void outer_prod(double * x, double * y, int size, double prod[][size]);
double dot_prod(double * x, double * y, int size);
void cross_prod(double * A, double * B, double * X);
double mag(double * x, int size);
void vec_rot(double x[][3], int axis, double theta);
double vec_min(int s, double *A, int *min1);
void Vcopy(int size, double A[size], double B[size]);
void Vcopy_int(int size, int A[size], int B[size]);
double min(double x, double y);

//matrix routines
void mat_mult(int rowsA, int colsA, int rowsB, int colsB, double A[][colsA], double B[][colsB], double prod[][colsB]);
void transpose(int rows, int cols, double A[][cols], double At[][rows]);
int LUdcmp(int size, double A[][size], double LU[][size], double * indx);
void LUbksub(int size, double LU[][size], double * b, double * indx, double * res);
void MInv(int size, double A[][size], double Y[][size]);
void Mcopy(int rows, int cols, double A[][cols], double B[][rows]);
double trace(int n, double A[][n]);
void mat_vec_mult(int rows, int cols, double A[][cols], double * v, double * prod);
void vec_mat_mult(int rows, int cols, double A[][cols], double * v, double * prod);
double mat_min(int rows, int cols, double A[rows][cols], int *min1, int*min2);
void mat_add(int rows, int cols, double A[][cols], double B[][cols], double sum[][cols]);
void mat_scalar_mult(int rows, int cols, double A[][cols], double a, double prod[][cols]);


//display routines
void Vdisplay(double * x, int size);
void Mdisplay(int rows, int cols, double x[][cols]);
void Vdisplay_int(int * x, int size);
