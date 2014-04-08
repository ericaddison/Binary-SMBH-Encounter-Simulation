/********************************************************************
*********************************************************************

    CHAIN algorithm attempt

    Author: Eric Addison
    Initial Date: 20Feb2011
    Affiliation: USU Astrophysics

    File: leapfrog.c
    Purpose: contains leapfrog integration method for CHAIN ala Tanikawa and Mikkola, 1999

*********************************************************************
*********************************************************************/


#include "chain.h"


/*--- leap: takes a leapfrog step ---*/
double leap(CHAIN * C, double h)
{

	double dTdP[C->N-1][3], dUdQ[C->N-1][3], du1[C->N-1][3], du2[C->N-1][3], T, a1, a0, U, Pt = -(C->E);

	//set up to calculate Q_1/2
	T = calc_T(C);
	a0 = 1/(2.0*(T + Pt));
	calc_dT(C,dTdP);
	mat_scalar_mult(C->N-1,3,dTdP,h*a0,dTdP);	//multiplied dTdP by scalar a
	mat_add(C->N-1,3,C->R,dTdP,C->R);		//this has updated C->R to equal Q_1/2

	//calculate P_1
	calc_dU1(C,du1);
	calc_dU2(C,du2);
	mat_add(C->N-1,3,du1,du2,dUdQ);
	U = calc_U(C);
	mat_scalar_mult(C->N-1,3,dUdQ, -h/U ,dUdQ);	//multiplied dUdQ by scalar -h/U
	mat_add(C->N-1,3,C->W,dUdQ,C->W);		//this has updated C->W to equal P_1

	//set up to calculate Q_1
	T = calc_T(C);
	a1 = 1/(2.0*(T + Pt));
	calc_dT(C,dTdP);
	mat_scalar_mult(C->N-1,3,dTdP,h*a1,dTdP);	//multiplied dTdP by scalar a
	mat_add(C->N-1,3,C->R,dTdP,C->R);		//this has updated C->R to equal Q_1/2


	return h/2.0 * (a1 + a0);				//returns update for time value, so use t += leap(C,h);

}

/*--- calc_T: calculates the current kinetic energy T(P) ---*/
double calc_T(CHAIN * C)
{

	double T = 0, mk, mkp1;

	for(int k=0; k<C->N-1; k++){
		mk = C->mass[C->chain[k]];
		mkp1 = C->mass[C->chain[k+1]];
		T += 0.5 * (1/mk + 1/mkp1) * dot_prod(C->W[k],C->W[k],3);
		if(k > 0)
			T -= dot_prod(C->W[k-1],C->W[k],3)/mk;
	}

	return T;
}

/*--- calc_U: calculates the current potential energy U(Q) ---*/
double calc_U(CHAIN * C)
{
	double U = 0, Rij[3], mk, mkp1, mi, mj;

	for(int k=0; k<C->N-1; k++){
		mk = C->mass[C->chain[k]];
		mkp1 = C->mass[C->chain[k+1]];
		U += G * mk * mkp1 / mag(C->R[k],3);
	}

	for(int j = 2; j < C->N; j++){
		for(int i = 0; i < j-1; i++){
			get_Rij(C,i,j,Rij);
			mi = C->mass[C->chain[i]];
			mj = C->mass[C->chain[j]];
			U += G * mi*mj / mag( Rij, 3);
		}}
	return U;
}



/*--- calc_dT: calculates the derivative dTdP ---*/
void calc_dT(CHAIN * C, double dT[][3])
{

	double mi, mip1;

	for(int i=0; i<C->N-1; i++){
		for(int k=0; k<3; k++){
			mi = C->mass[C->chain[i]];
			mip1 = C->mass[C->chain[i+1]];

			if(i == 0)
				dT[i][k] = 1/mi * (C->W[i][k]) + 1/mip1 * (C->W[i][k] - C->W[i+1][k]);
			else if(k == C->N-2)
				dT[i][k] = 1/mi * (C->W[i][k] - C->W[i-1][k]) + 1/mip1 * (C->W[i][k]);
			else
				dT[i][k] = 1/mi * (C->W[i][k] - C->W[i-1][k]) + 1/mip1 * (C->W[i][k] - C->W[i+1][k]);
	}}

		return;
}

/*--- calc_du1: calculates the first part of the derivative dUdQ ---*/
void calc_dU1(CHAIN * C, double du1[][3])
{
	double mi, mip1;

	for(int i=0; i < C->N-1; i++){
		for(int k = 0; k < 3; k++)				//calculated derivative dUdQ piece one, see 2Mar11 BEMRI notes
		{
			mi = C->mass[C->chain[i]];
			mip1 = C->mass[C->chain[i+1]];
			du1[i][k] = -C->R[i][k] * G * mi * mip1 / pow( mag( C->R[i],3),3);
		}
	}

	return;

}


/*--- calc_du2: calculates the first part of the derivative dUdQ ---*/
void calc_dU2(CHAIN * C, double du2[][3])
{
	double Rij[3], mn, mm, temp[3], a;

	for(int i=0; i<C->N-1; i++){
		du2[i][0] = 0; du2[i][1] = 0; du2[i][2] = 0;		//calculated derivative dUdQ piece two, see 2Mar11 BEMRI notes

		for(int m = 2; m < C->N; m++){
			for(int n = 0; n < m-1; n++){
				if(i >= n && i <= m-1){
					get_Rij(C,n,m,Rij);				//find new Rij vector
					mn = C->mass[C->chain[n]];
					mm = C->mass[C->chain[m]];
					a = -G*mn*mm/( pow(mag(Rij,3),3) );
					for(int l = n; l < m; l++)
					{							
						vect_mult_scalar(C->R[l],a,temp,3);
						vect_add(du2[i],temp,du2[i],3);
					}
				}
		}}}
	return;

}
