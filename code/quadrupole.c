/********************************************************************
 *********************************************************************
 
 BEMRI Tidal Disruption Simulation
 
 Author: Eric Addison
 Initial Date: 3 Dec 2010
 Affiliation: USU Astrophysics
 
 File: quadrupole.c
 Purpose: code to compute quadrupole radiation.
 
 *********************************************************************
 *********************************************************************/

#include "BEMRI.h"

/*************************************************************************/
/*************************************************************************/
/**																		**/
/** Function: point_mass_QP												**/
/** Purpose: calculate the quadrupole of a collection of point masses	**/
/**	Inputs: vector of masses, matrix of positions, number of masses,	**/
/**			CM vector (from origin), CM vector (from observer)			**/
/**			output matrix Q, output matrix Qtt (transverse traceless)	**/								
/**			                            								**/
/**	Output: populate output matrix					**/
/**									**/
/**	Note: This code is based on the quadrupole radiation 		**/
/**		scheme from the Walhquist paper				**/
/**									**/
/*************************************************************************/
/*************************************************************************/

void point_mass_QP(double * m, double r[][3], int N, double * R, double Q[][3], double Qtt[][3])
//qp should be 3x3
//the matrix r should have vectors from the coordinate origin to the masses
{

	//variable declarations
		double M=0, RR[3][3], rrr[N][N][3], rr_tensor[3][3], rr_sum[3][3], Rhat[3], RRhat[3][3];
		double I[3][3], P[3][3], rT[3],rTrT[3][3], temp[3][3], Qt[3][3];
	
	//calcuate total mass M
	for(int ii = 0; ii < N; ii++)							
		M += m[ii];
	
	//initialize rr_sum and identity
	for(int ii = 0; ii < 3; ii++)							
	{
		for(int jj = 0; jj < 3; jj++)
		{
			rr_sum[ii][jj]	= 0;
			if(ii == jj)
				I[ii][jj] = 1;
			else
				I[ii][jj] = 0;
		} }
	
	//outer product of R vectors
	outer_prod(R,R,3,RR);									
	
	//build displacement vectors between masses
	for(int jj = 0; jj < N; jj++)							
	{
		for(int ii = 0; ii < jj; ii++)				
		{
			rrr[ii][jj][0] = r[ii][0] - r[jj][0];
			rrr[ii][jj][1] = r[ii][1] - r[jj][1];
			rrr[ii][jj][2] = r[ii][2] - r[jj][2];
		} }
	
	//sum over these vectors
	for(int jj = 0; jj < N; jj++)							
	{
		for(int ii = 0; ii < jj; ii++)				
		{
			outer_prod(rrr[ii][jj],rrr[ii][jj],3,rr_tensor);	//construct current outer product
			for(int kk = 0; kk < 3; kk++)
			{
				for(int ll = 0; ll < 3; ll++)
				{
					rr_sum[kk][ll] += m[ii]*m[jj]/M * rr_tensor[kk][ll];
				} }
			
		} }
	
	//build quadrupole tensor
	for(int ii = 0; ii < 3; ii++)							
	{
		for(int jj = 0; jj < 3; jj++)				
		{
			Q[ii][jj] = M*RR[ii][jj] + rr_sum[ii][jj];
			if(log10(Q[ii][jj]) < -300)
				Q[ii][jj] = 0;
		} }
	
	//build projection tensor P
	{
		Rhat[0] = R[0] / mag(R,3);								
		Rhat[1] = R[1] / mag(R,3);
		Rhat[2] = R[2] / mag(R,3);
		outer_prod(Rhat,Rhat,3,RRhat);
		for(int ii = 0; ii < 3; ii++)						
		{
			for(int jj = 0; jj < 3; jj++)
			{
				P[ii][jj] = I[ii][jj] - RRhat[ii][jj];
			}
		}
	}

	
	//build transverse tensor Qt
	{
		
		//QT
		mat_mult(3,3,3,3,P,Q,temp);
		mat_mult(3,3,3,3,temp,P,Qt);
	}
	
	
	//build transverse traceless quadrupole tensor Qtt
	for(int ii = 0; ii < 3; ii++)						
	{
		for(int jj = 0; jj < 3; jj++)
		{
			Qtt[ii][jj] = Qt[ii][jj] - 0.5*trace(3,Qt)*P[ii][jj];
		}
	}
	
	
	//display code
	/*{
	Mdisplay(3,3,Qtt);
	
	printf("\ntrace = %.3e\n",trace(3,Qtt));
	anykey();
	}*/
	
	return;
	
}

















