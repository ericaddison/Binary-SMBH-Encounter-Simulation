/*******************************************************************
********************************************************************

    BEMRI Tidal Disruption Simulation

    Author: Eric Addison
    Initial Date: 8 Feb 2010
    Affiliation: USU Astrophysics

    File: BEMRI.h
    Purpose: Main project header file. Contains all global function
             declarations, #includes, etc.

********************************************************************
********************************************************************/

#ifndef _BEMRI_H_
#define _BEMRI_H_

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <time.h>
#include <sys/time.h>
#include <errno.h>
#include <ctype.h>
#include "kb.h"
#include "matrix.h"
#include "ODE.h"
#include "chain.h"

int fcount;

//macros for ease of coding
#define m1 b->binary_b.mass1
#define m2 b->binary_b.mass2
#define m3 b->binary_h.mass1
#define eb b->binary_b.ecc
#define Pb b->binary_b.P_orb
#define ab b->binary_b.a
#define ia b->binary_b.i_a
#define aop b->binary_b.ao_p
#define lan b->binary_b.la_n
#define rpb b->binary_b.rp
#define rph b->binary_h.rp
#define R1 b->binary_b.x1
#define X1 b->binary_b.x1[0]
#define Y1 b->binary_b.x1[1]
#define Z1 b->binary_b.x1[2]
#define R2 b->binary_b.x2
#define X2 b->binary_b.x2[0]
#define Y2 b->binary_b.x2[1]
#define Z2 b->binary_b.x2[2]
#define R3 b->binary_h.x1
#define X3 b->binary_h.x1[0]
#define Y3 b->binary_h.x1[1]
#define Z3 b->binary_h.x1[2]
#define R4 b->binary_h.x2
#define X4 b->binary_h.x2[0]
#define Y4 b->binary_h.x2[1]
#define Z4 b->binary_h.x2[2]
#define V1 b->binary_b.v1
#define VX1 b->binary_b.v1[0]
#define VY1 b->binary_b.v1[1]
#define VZ1 b->binary_b.v1[2]
#define V2 b->binary_b.v2
#define VX2 b->binary_b.v2[0]
#define VY2 b->binary_b.v2[1]
#define VZ2 b->binary_b.v2[2]
#define V3 b->binary_h.v1
#define VX3 b->binary_h.v1[0]
#define VY3 b->binary_h.v1[1]
#define VZ3 b->binary_h.v1[2]
#define V4 b->binary_h.v2
#define VX4 b->binary_h.v2[0]
#define VY4 b->binary_h.v2[1]
#define VZ4 b->binary_h.v2[2]

#define mu_b b->binary_b.mu
#define mu_h b->binary_h.mu
#define theta_b b->binary_b.true_anomaly
#define theta_h b->binary_h.true_anomaly
#define psi_b b->binary_b.eccentric_anomaly
#define psi_h b->binary_h.eccentric_anomaly
#define lb b->binary_b.angular_momentum
#define lh b->binary_h.angular_momentum
#define Eb b->binary_b.E
#define Eh b->binary_h.E

#define m4 b->binary_h.mass2
#define eh b->binary_h.ecc
#define Ph b->binary_h.P_orb
#define ah b->binary_h.a

#define rr C->r
#define vv C->v
#define EPS 1e-9

typedef struct binary{
    double mass1;
    double mass2;
    double M;               //total mass
    double mu;              //reduced mass
    double P_orb;           //period of the binary
    double ecc;           //eccentricity of the binary
    double a;          //semi major axis
    double rp;			//periapse distance
    double ao_p;        //Argument of periastron angle of the binary
    double i_a;         //inclination angle of the binary
    double la_n;        //Longitude of the Ascending Node angle for the binary
    double x1[3];        //coordinates for mass m1 (x,y,z)
    double x2[3];        //coordinates for mass m2 (x,y,z)
    double v1[3];       //velocity components for mass m1 (vx,vy,vz)
    double v2[3];       //velocity components for mass m2 (vx,vy,vz)
    double true_anomaly;	//true anomaly of the binary
    double eccentric_anomaly;	//eccentric anomaly
    double mean;		//mean anomaly
    double angular_momentum;	//angular momentum
    double E;			//internal energy
    double l[3];			//angular momentum vector
}binary;


// initial parameters struct
typedef struct iParams{
	double gamma;	// gamma value
	double beta;	// beta value
	double incA;	// inclination
	double lanA;	// longitude of ascending node
	double th0;	// initial phase

} iParams;

//BEMRI structure to be used for the program

typedef struct BEMRI{
    binary binary_b;     //little binary
    binary binary_h;    //black hole binary
    double energy;			//total system energy
    double L;			//total system ang. mom.
}BEMRI;
/*
typedef struct params{
    double pm1;
    double pm2;
    double pm3;
    double peb;
    double peh;
    double pPb;
    double pPh;
    double pab;
    double pah;
    double paop;
    double pia;
    double plan;
}params;
*/
typedef struct flags{
	unsigned int cflag:1;
	unsigned int sflag:1;
	unsigned int tflag:1;
	unsigned int ndflag:1;
	unsigned int Nflag:1;
	unsigned int soflag:1;
	unsigned int PbNflag:1;
	unsigned int th0flag:1;
	unsigned int bs2flag:1;
	unsigned int fnameflag:1;
	unsigned int betaflag:1;
	unsigned int gamflag:1;
	unsigned int ipsflag:1;
	unsigned int ehflag:1;
	unsigned int Htflag:1;
	unsigned int Hrat:1;
	unsigned int Xflag:1;
	unsigned int geoflag:1;
	unsigned int angflag:1;
	unsigned int fParflag:1;
} flags;

CHAIN *C;		//declare C, a struct pointer

//Function Declarations

void load_params(struct BEMRI * b);


//BEMRI_functions.c functions
void dvdt(double * r1, double m, double * r2, double m_2, double * out);
double two_body_energy(double * r1, double * r2, double * v1, double * v2, double M1, double M2);
double angular_momentum(double * r1, double * r2, double * v1, double * v2, double M1, double M2);
double angular_momentum_vector(double * r1, double * r2, double * v1, double * v2, double M1, double M2, double * l);
int true_anomaly(binary * b);
int calc_angles(binary * b);
void CM_values(BEMRI * b);
void advance_orbit(binary * b, double l, double theta);
double advance_orbit_set_up(BEMRI * b, double chi, double r_current, double dt1);
void calc_energies(BEMRI * b, double * E1h, double * E2h);
void calc_angular_momenta(BEMRI * b, double * l1h, double * l2h);
void calc_ecc(BEMRI * b, double * e1h, double * e2h, double E1h, double E2h, double l1h, double l2h);
void calc_total_energy(BEMRI * b);
void calc_total_angular_momentum(BEMRI * b);
double sign(double x);
double fix_angle(double x);
double peters_lifetime(double ecc, double a, double M1, double M2);
double lifetime_integrand(double ecc);
double RK4_peters(double x, double dx);
void torque(BEMRI * b, double * t);

//RK4

double N_body_main(double (*r)[3], double (*v)[3], double * m, double tmax, int N, double *dt, double dt_min, int * minned);


//main.c
void view_info(BEMRI * b);
double min(double x, double y);
double max(double x, double y);
void update_binaries(BEMRI * b, double (*r)[3], double (*v)[3]);
void update_vectors(BEMRI * b, double (*r)[3], double (*v)[3]);
void BEMRI_CM_update(BEMRI * b);
void init_params(BEMRI * b, flags f);
int print_percentage(int a, int i, double t, double tmax, int N);
void recalc_params(BEMRI * b);
void load_vectors(double (*r)[3], double (*v)[3], double * m, BEMRI * b);
void load_vectors2(double **r, double **v, double * m, BEMRI * b);
void load_special_ICs(BEMRI * b, double * t0);
double sign(double x);
int which_agent(void);
void output_time(double t, FILE * fp);
void output_positions(double (*r)[3], BEMRI * b, FILE * fp);
void output_QP(double (*Qtt)[3], FILE * fp);
void Nbody_derivs1(double x, double y[], double dydx[]);
void Nbody_derivs2(double x, double y[], double dydx[]);
void pack_y_vector(double (*r)[3], double (*v)[3], double m[], double y[]);
void unpack_y_vector(double (*r)[3], double (*v)[3], double *m, double *y);
//void setup_chain(CHAIN *C, double (*r)[3], double (*v)[3],double *m,int N);
void setup_chain(CHAIN *C, double *m,int N);
void pack_y_vector_C(CHAIN *C, double y[]);
void unpack_y_vector_C(CHAIN *C, double y[]);
void leap_derivs2(double xs, double y[], double dydx[]);

//kepler_eq.c
double kepler(double E, double e);
double solve_kepler(double M, double e);
double f_kepler(double E, double M, double e);
double bisect_kepler(double a, double b, double M, double e, double f(double E, double M, double e), double tol);

//random.c
float ran2(long *idum);
float gasdev1(long *idum);
float gasdev2();

//quadrupole.c
void point_mass_QP(double * m, double r[][3], int N, double * R, double Q[][3], double Qtt[][3]);

//args.c
int args(int argc,char * argv[],iParams *iPars, flags * f, int * N, double * frac, long * seed, double * PbN,char * fname, double * ebh, int * ipsval, double * Hrat);

#endif
