/********************************************************************
*********************************************************************

    CHAIN algorithm attempt

    Author: Eric Addison
    Initial Date: 20Feb2011
    Affiliation: USU Astrophysics

    File: chain.h
    Purpose: Primary header file for CHAIN algorithm

*********************************************************************
*********************************************************************/


#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <sys/time.h>
#include <errno.h>
#include <ctype.h>
#include "matrix.h"
#include "ODE.h"
#include "cosmo.h"

#define N_MAX 100


typedef struct{
	//to use this structure, first fill in N,mass, r, and q with intial values
	//then run init_chain to initialize the chain stucture
	int N;						//number of ptcls
	double M;					//total mass of the chain
	double mass[N_MAX];			//masses of ptcls
	double r[N_MAX][3];			//positions of ptcls
	double v[N_MAX][3];			//velocities of ptcls
	double q[N_MAX][3];			//center of mass positions
	double p[N_MAX][3];			//center of mass momenta
	double r0[3];				//center of mass position
	double v0[3];				//center of mass velocity
	double p0[3];				//center of mass momentum
	int chain[N_MAX];			//chain vector
	int old_chain[N_MAX];			//previous chain vector
	double R[N_MAX][3];			//chain position vectors R
	double W[N_MAX][3];			//chain momentum vectors W
	double T0;					//initial kinetic energy term
	double U0;					//initial potential energy term
	double E;					//initial numerical energy term = T+U
	int init;					//flag to see if it's the initial chain construction
}CHAIN;



/* chain.c functions */
void init_chain(CHAIN * C);
void set_N(CHAIN * C);
void build_chain(CHAIN * C);
void set_infs(int N,double (*DM)[N],int min1, int min2);
void set_more_infs(int N,double (*DM)[N],int old);
void append_start(int *chain,int smin,int chained);
void distance_matrix(CHAIN * C, double (*DM)[C->N]);
double distance(double x[3], double y[3]);
void print_masses(CHAIN * C);
void print_q(CHAIN * C);
void print_p(CHAIN * C);
int check_chain( CHAIN *C);
int check_chain1(CHAIN *C);
int check_chain2(CHAIN *C);
void init_R_W(CHAIN *C);
void set_p_q(CHAIN * C);
void set_energy(CHAIN * C);
int get_Rij(CHAIN * C,int i, int j, double *Rij);
void output_positions_C(CHAIN * C, FILE * fp, double t);
void update_positions(CHAIN * C);
void update_momenta(CHAIN * C);
void update_R_W(CHAIN * C);
double chain_avg(CHAIN * C);

/* leapfrog.c functions */
void calc_dU1(CHAIN * C, double du1[][3]);
void calc_dU2(CHAIN * C, double du2[][3]);
void calc_dT(CHAIN * C, double dT[][3]);
double calc_T(CHAIN * C);
double calc_U(CHAIN * C);
double leap(CHAIN * C, double h);
