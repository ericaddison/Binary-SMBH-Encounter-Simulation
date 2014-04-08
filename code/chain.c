/********************************************************************
*********************************************************************

    CHAIN algorithm attempt

    Author: Eric Addison
    Initial Date: 20Feb2011
    Affiliation: USU Astrophysics

    File: chain.c
    Purpose: Primary function file for CHAIN algorithm

*********************************************************************
*********************************************************************/

#include "chain.h"



/*--- init_chain: function to initialize chain mass array ---*/
void init_chain( CHAIN * C)
{
	for(int i = 0; i < N_MAX; i++){
		if(i < C->N)
			C->M += C->mass[i];					//find total mass M
		else
			C->mass[i] = 0;
	}

	for(int i = 0; i < C->N; i++)
		C->old_chain[i] = i;

	set_p_q(C);									//calculate center of mass positions and momenta
	set_energy(C);								//calculates initial energy terms
	build_chain(C);
	//set_R_W(C);
	return;
}

/*--- check_chain: check for chain violations and rebuild if necessary ---*/
int check_chain( CHAIN *C)
{
int f=0;
	if(check_chain1(C))
		{
			f = 1;
		}
	else if(check_chain2(C))
		{
			f = 1;
		}

	if(f == 1){
		Vcopy_int(C->N,C->chain,C->old_chain);		//first copy the current chain into old_chain
		update_positions(C);
		update_momenta(C);							//find new q's and p's
		build_chain(C);								//rebuild the chain, automatically updates R and W
		return 1;									//return 1 => chain was rebuilt
	}

	return 0;										//return 0 => chain was not rebuilt
}

/*--- check_chain1: check the first condition for chain reset ---*/
int check_chain1( CHAIN *C)
{
	//this will check all the triangles formed between two chained vectors and one non-chained vector
	double temp[3], unchained_mag=0, ave, magi, magip1;
	//try checking the triangles between the first two vectors and any others

	ave = chain_avg(C);

//Mdisplay(C->N-1,3,C->R);
	for(int i = 0; i<C->N-2; i++)
	{
		get_Rij(C,i,i+2,temp);
		unchained_mag = mag(temp,3);
		magi = mag(C->R[i],3);
		magip1 = mag(C->R[i+1],3);

//		printf("\nmag1 = %.3e",mag1);
		/* if the magnitude of the non-chained vector is greater than twice the average chain vector, don't bother */
		if( unchained_mag < 2*ave)//4*min(magi,magip1) )
		{
			if( unchained_mag < magi || unchained_mag < magip1 )		//then non-chained leg is shortest
			{
//printf("\nunchained_mag = %.3e\nR_%d = %.3e\nR_%d = %.3e\n2*ave_mag = %.3e\n",unchained_mag,i,magi,i+1,magip1,2*ave);

				return 1;			//return 1, means to rechain
			}
		}
	}

	return 0;
}


/*--- check_chain2: check the second condition for chain reset ---*/
int check_chain2( CHAIN *C)
{
	//this will check all the non-chained vectors and compare them to adjacent chained vectors
	double temp[3], mag1, Ri, Rim1;

	Ri = 0;
	for(int i = 0; i<C->N-2; i++)			//go through all but the last two points in the chain
	{
//		printf("\ni = %d\n",i);
		Rim1 = Ri;
		Ri = mag(C->R[i],3);
		//now go through the rest of the points in the chain, up to the last one
		//the very next adjacent point need not be checked, that vector is already chained
		//the point after that need not be checked either, it forms a triangle to has been checked in chain_check1
		//so the first new point to check is at j = i+3
		for(int j=i+3;j<C->N;j++)
		{
//			printf("\nchecking vector r_(%d,%d)\n",C->chain[i],C->chain[j]);
			temp[0] = 0; temp[1] = 0; temp[2] = 0;
			//find the magnitude of the vector between points j and j+3 in the chain
			//the displacement vector r_ij is equal to R_i + R_(i+1) + ... + R_(j-1)
			for(int k=i;k<j;k++){
				vect_add(C->R[k],temp,temp,3);
			}
			mag1 = mag(temp,3);
//			printf("\nmagnitude of r_(%d,%d) = %.3e\n",C->chain[i],C->chain[j],mag1);

			if(i > 0)		//if the second point is not the start of the chain
			{
//				printf("\nmag of R[i-1] = %.3e\n",Rim1);
				if(mag1 > Rim1)
					continue;
			}
//			printf("\nmag of R[i] = %.3e\n",Ri);
			if(mag1 > Ri || mag1 > mag(C->R[j-1],3))		//if mag1 is bigger than any of these vectors then it's not the smallest
				continue;

			if(j < C->N-1)		//if the second point is not the end of the chain
			{
				if(mag1 > mag(C->R[j],3) )
					continue;
			}

			return 1;		//if mag1 was not larger than any consecutive chain vector, return 1 and rebuild chain

		}
	}
		//printf("mag1 = %.3e\nmagRi = %.3e\nmagRi+1 = %.3e\n",mag1,mag(C->R[i],3),mag(C->R[i+1],3));
	return 0;
}




/*--- init_R_W: build the position chain vectors R from the absolute positions q ---*/
void init_R_W( CHAIN *C)
{
	int ci, cip1;
	double temp[3];

	for(int i=0; i<C->N-1; i++){
		ci = C->chain[i];
		cip1 = C->chain[i+1];
		//build R's
		C->R[i][0] = C->q[cip1][0] - C->q[ci][0];
		C->R[i][1] = C->q[cip1][1] - C->q[ci][1];
		C->R[i][2] = C->q[cip1][2] - C->q[ci][2];

		//build W's
		if(i==0){
			vect_mult_scalar(C->p[ci],-1,temp,3);
			Vcopy(3,temp,C->W[0]);
		}
		else{
			C->W[i][0] = C->W[i-1][0] - C->p[ci][0];
			C->W[i][1] = C->W[i-1][1] - C->p[ci][1];
			C->W[i][2] = C->W[i-1][2] - C->p[ci][2];
		}

	}

	C->init = 1;

	return;
}


/*--- set_energy: find values of kinetic and potential energy ---*/
void set_energy( CHAIN * C)
{
	double temp[3];
	C->T0 = 0;
	C->U0 = 0;
	C->E = 0;

	for(int i=0; i<C->N;i++){
		C->T0 += 1/(2*C->mass[i])*dot_prod(C->p[i],C->p[i],3);
		for(int j=0; j<C->N;j++){
			if(i < j){
				vect_mult_scalar(C->q[i],-1,temp,3);			//multiply q[i] by -1
				vect_add(temp,C->q[j],temp,3);					//calculates q[j] - q[i] vector
				C->U0 += G*C->mass[i]*C->mass[j] / mag(temp,3);	//adds term to potential energy
			}
		}
	}
	//C->T0 += 1/(2*C->M) * dot_prod(C->p0,C->p0,3);			//!!! currentlty does not include center of mass energy
	C->E = C->T0 - C->U0;

	return;
}
/*--- set_p_q: find values of p and q vectors ---*/
void set_p_q( CHAIN * C)
{
	for(int i=0; i<C->N;i++)
	{
		C->r0[0] += C->mass[i]*C->r[i][0];
		C->r0[1] += C->mass[i]*C->r[i][1];
		C->r0[2] += C->mass[i]*C->r[i][2];
		C->v0[0] += C->mass[i]*C->v[i][0];
		C->v0[1] += C->mass[i]*C->v[i][1];
		C->v0[2] += C->mass[i]*C->v[i][2];
	}
	vect_mult_scalar(C->r0, 1/C->M, C->r0, 3);
	Vcopy(3,C->v0,C->p0);
	vect_mult_scalar(C->v0, 1/C->M, C->v0, 3);
	for(int i=0; i<C->N;i++)
	{
		C->q[i][0] = C->r[i][0] - C->r0[0];
		C->q[i][1] = C->r[i][1] - C->r0[1];
		C->q[i][2] = C->r[i][2] - C->r0[2];
		C->p[i][0] = C->mass[i]*(C->v[i][0] - C->v0[0]);
		C->p[i][1] = C->mass[i]*(C->v[i][1] - C->v0[1]);
		C->p[i][2] = C->mass[i]*(C->v[i][2] - C->v0[2]);
	}
	return;
}


/*--- set_N: find the value of N, based on masses that have been added ---*/
void set_N( CHAIN * C)
{
	printf("\nset_N\n");
	int i = 0;
	while(C->mass[i] != 0){
		i++;
		C->N++;
	}

	printf("\nN=%d\n",C->N);
	return;
}

/*--- get_Rij: find the value of Rij, the displacement vector between two points ---*/
int get_Rij( CHAIN * C,int i, int j, double *Rij)
{
	int a = 1;

	if(i > (C->N)-1 || j > (C->N)-1 || i < 0 || j < 0){
		printf("\nERROR: Chain indices out of bounds\n");
		return -1;
	}

	if(i > j) {a = j; j = i; i = a; a = -1;}	//in case i > j, which means we'll find Rji and point it backwards with a = -1

	Rij[0] = 0; Rij[1] = 0; Rij[2] = 0;			//clear result vector

	for(int kk = i; kk < j; kk++)				//sum over chain variables from kk = i to j-1
	{
		Rij[0] += C->R[kk][0];
		Rij[1] += C->R[kk][1];
		Rij[2] += C->R[kk][2];
	}

	Rij[0] *= a;
	Rij[1] *= a;
	Rij[2] *= a;

	return 0;

}

/*--- update_R_W: update to new values of R and W ---*/
void update_R_W(CHAIN * C)
{

	int k0,k1,ci;
	double Rnew[C->N-1][3] , Rold[C->N-1][3] , temp[3], B;

	Mcopy(C->N-1, 3, C->R, Rold);

	//build new R's from old R's
	//method from mikkola and Aarseth CHAIN paper, 1993
	k1 = 0;


	//find initial value of k0, store it in k1
	while(C->old_chain[k1] != C->chain[0])
		k1++;

	for(int mu=0; mu<C->N-1; mu++){
		Rnew[mu][0] = 0; Rnew[mu][1] = 0; Rnew[mu][2] = 0;

		//find values k0 and k1 such that old_chain[k0] = chain[mu] and old_chain[k1] = chain[mu+1]
		k0 = k1; k1 = 0;
		//since we already have k1 from the previous iter, it is good for k0 of the current iter.
		while(C->old_chain[k1] != C->chain[mu+1])
			k1++;

		//now we should have k0 and k1 for the R_mu^new
		for(int nu = 0; nu < C->N-1; nu++){
			if(k1 > nu && k0 <= nu)
				B = 1;
			else if(k1 <= nu && k0 > nu)
				B = -1;
			else
				B = 0;

			vect_mult_scalar(Rold[nu],B,temp,3);
			vect_add(temp,Rnew[mu],Rnew[mu],3);

		}}

	//now copy Rnew into C->R

	Mcopy(C->N-1,3,Rnew,C->R);

	//rebuild W's
	for(int i=0; i<C->N-1; i++){
		ci = C->chain[i];
		if(i==0){
			C->W[i][0] = -C->p[ci][0];
			C->W[i][1] = -C->p[ci][1];
			C->W[i][2] = -C->p[ci][2];
		}
		else{
			C->W[i][0] = C->W[i-1][0] - C->p[ci][0];
			C->W[i][1] = C->W[i-1][1] - C->p[ci][1];
			C->W[i][2] = C->W[i-1][2] - C->p[ci][2];
		}
}

}

/*--- update_momenta: recover and update center of mass momenta p from chain variables R ---*/
void update_momenta(CHAIN * C)
{

	double temp[3];
	int a;

	a = C->chain[0];
	vect_mult_scalar(C->W[0],-1,C->p[a],3);
	a = C->chain[C->N-1];
	vect_mult_scalar(C->W[C->N-2],1,C->p[a],3);

	for(int i = 1; i < C->N; i++)				//recover p's
	{
		a = C->chain[i];
		vect_mult_scalar(C->W[i],-1,temp,3);
		vect_add(C->W[i-1],temp,C->p[a],3);
	}

	for(int i=0;i<C->N;i++)
	{
		for(int j=0;j<3;j++)
			C->v[i][j] = C->p[i][j]/C->mass[i];
	}
	
	return;

}

/*--- update_positions: recover and update center of mass coordinates q from chain variables R ---*/
void update_positions(CHAIN * C)
{

	double q[C->N][3], q0[3], qtemp[3];
	int a;

	q0[0] = q0[1] = q0[2] = 0;
	q[0][0] = q[0][1] = q[0][2] = 0;		//coordinate recovery scheme from Mikkola and Aarseth, page 446

	//Mdisplay(2,3,C->R);
	for(int i = 1; i < C->N; i++)				//calculate individual q's and center of mass q0
	{
		a = C->chain[i];
		vect_add(q[i-1],C->R[i-1],q[i],3);		//calculate new q's
		//Vdisplay(C->R[i-1],3);
		vect_mult_scalar(q[i],C->mass[a]/C->M,qtemp,3);	//calculate q center of mass
		vect_add(q0,qtemp,q0,3);
	}

	//Mdisplay(3,3,q);
	
	vect_mult_scalar(q0,-1,q0,3);

	for(int i = 0; i < C->N; i++)				//copy new positions into C->q
	{
		a = C->chain[i];
		/*NOTE: This stores the position of m[i] in q[i], not the position of m[chain[i]]! */
		vect_add(q0,q[i],q[i],3); 
		Vcopy(3, q[i], C->q[a]);
		//Vdisplay(C->r[i],3);
	}

	for(int i = 0; i < C->N; i++)
	{
		for(int j=0;j<3;j++)
			C->r[i][j] = C->q[i][j] + C->r0[j];
	}
	
	//Mdisplay(3,3,C->r);
	
	return;

}

/*--- output_positions: write current positions to file, center of mass coordinates ---*/
void output_positions_C(CHAIN * C, FILE * fp, double t)
{
	int a;
	double qout[C->N][3];

	for(int i = 0; i < C->N; i++)				//output positions
	{
//		a = C->chain[i];
		a=i;
		Vcopy(3, C->q[i], qout[a]);
	}

	fprintf(fp,"%.10e,",t);
	for(int i = 0; i < C->N; i++)				//output positions
		//fprintf(fp,"%.10e, %.10e, %.10e, ",q[i][0] - q0[0] + C->r0[0] ,q[i][1] - q0[1]+ C->r0[1],q[i][2] - q0[2]+ C->r0[2]);
		fprintf(fp,"%.10e, %.10e, %.10e, ",qout[i][0],qout[i][1],qout[i][2]);
	fprintf(fp,"\n");

	return;

}




/*--- build_chain: build the interparticle chain from scratch ---*/
void build_chain( CHAIN * C)
{
	double DM[C->N][C->N], DMT[C->N][C->N], start_min, end_min, x,y;
	int min1, min2, smin,xs,ys,emin,xe,ye, chained = 0, old;			//count to record how many ptcls have been chained

	if(C->N == 0)		//check if N has been set yet
		set_N(C);	if(C->N == 0)		//check if N has been set yet
			set_N(C);

	distance_matrix(C,DM);				//build distance matrix

//Mdisplay(C->N,3,C->q);

//Mdisplay(C->N,C->N,DM);
	mat_min(C->N,C->N,DM,&min1,&min2);
	/* first two entries in the chain */
	/* here the distance matrix DM is built off of actual labels, not chain labels, so min1 and min2 correspond
	* to actual mass numbers, not chain positions	*/
	C->chain[0] = min1;
	C->chain[1] = min2;
	chained = 2;
//printf("\nmin1 = %d\nmin2 = %d\n",min1,min2);
	set_infs(C->N,DM,min1,min2);			//set appropriate entries to inf
//Mdisplay(C->N,C->N,DM);
	for(int i = 2; i < C->N; i++)	//chain the rest
	{
		transpose(C->N,C->N,DM,DMT);							//compute the traspose of DM
		x = vec_min(C->N,DM[C->chain[0]],&xs);				//x is the minimum value in the row of chain[0]
		y = vec_min(C->N,DMT[C->chain[0]],&ys);				//y is the minimum value in the col of chain[0]
		if(x < y){ start_min = x; smin = xs;} 				//find ptcl closest to start of chain
		else { start_min = y; smin = ys;}

		x = vec_min(C->N,DM[C->chain[chained-1]],&xe);
		y = vec_min(C->N,DMT[C->chain[chained-1]],&ye);
		if(x < y){end_min = x; emin = xe;} 				//find ptcl closest to start of chain
		else{end_min = y; emin = ye;}

//printf("\n\ni=%d",i);
//Mdisplay(C->N,C->N,DMT);
//printf("\nhere, emin = %d\nx = %.3e\ny = %.3e\n",emin,x,y);

		if(start_min < end_min){			//then append the new ptcl at the beginning of the chain
			append_start(C->chain,smin,chained);
			old = C->chain[0];
			//set_infs(C->N,DM,smin,smin);
			set_more_infs(C->N,DM,old);
		}
		else{						//append at the end of the chain
			old = C->chain[chained-1];
			C->chain[chained] = emin;
			//set_infs(C->N,DM,emin,emin);
			set_more_infs(C->N,DM,old);
		}
		set_infs(C->N,DM,C->chain[0],C->chain[chained]);
		chained++;

	}

//Mdisplay(C->N,C->N,DM);
//Vdisplay_int(C->chain,chained);

	if(C->init)
		update_R_W(C);
	else
		init_R_W(C);

//Mdisplay(C->N-1,3,C->R);


	return;

}


/*--- set_infs: set appropriate distances to inf to avoid loops in the chain ---*/
void set_infs(int N,double (*DM)[N],int min1, int min2)
{

	DM[min1][min2] = atof("Inf");
	DM[min2][min1] = atof("Inf");
	return;
}


/*--- set_more_infs: set appropriate distances to inf to avoid loops in the chain ---*/
void set_more_infs(int N,double (*DM)[N],int old)
{
	for(int i = 0; i < N; i++){
	 DM[i][old] = atof("Inf");
	 DM[old][i] = atof("Inf");
	 }
	return;
}

/*--- append_start: add a ptcl to the beginning of the chain ---*/
void append_start(int *chain,int smin,int chained)
{
	for(int i = chained-1;i>=0;i--)		//shift chain to the right one
		chain[i+1] = chain[i];

	chain[0] = smin;					//add new ptcl

	return;
}

/*--- distance_matrix: function to calculate all the interparticle distances ---*/
void distance_matrix( CHAIN * C, double (*DM)[C->N])
{
	for(int i=0;i<C->N;i++){
		for(int j=0;j<C->N;j++){			//only calculates distances for j<i
			if(j != i) DM[i][j] = distance(C->q[i],C->q[j]);
			else DM[i][j] = atof("Inf");
//		Vdisplay(C->q[i],3);
//		Vdisplay(C->q[j],3);
//		printf("\ndistance = %.3e",distance(C->q[i],C->q[j]));
		}
	}
	return;
}

/*--- distance: function to calculate the distance between two points ---*/
double distance(double x[3], double y[3])
{
	return sqrt( pow((x[0]-y[0]),2) + pow((x[1]-y[1]),2) + pow((x[2]-y[2]),2));
}


/*--- print_masses: function to print out mass values ---*/
void print_masses( CHAIN * C)
{
	printf("\nprint_masses\n");
	if(C->N == 0)		//check if N has been set yet
		set_N(C);

	for(int i = 0; i < C->N; i++)
			printf("\nmass %d = %.3e",i,C->mass[i]);

	printf("\n");
	return;
}

/*--- print_coords: function to print out coordinates ---*/
void print_q( CHAIN * C)
{
	printf("\nprint_coords\n");
	if(C->N == 0)		//check if N has been set yet
		set_N(C);

	for(int i = 0; i < C->N; i++){
		printf("\nq[%d][0] = %.3e",i,C->q[i][0]);
		printf("\nq[%d][1] = %.3e",i,C->q[i][1]);
		printf("\nq[%d][2] = %.3e\n",i,C->q[i][2]);
	}

	printf("\n");
	return;
}

/*--- print_p: function to print out CM momenta ---*/
void print_p( CHAIN * C)
{
	printf("\nprint_coords\n");
	if(C->N == 0)		//check if N has been set yet
		set_N(C);

	for(int i = 0; i < C->N; i++){
		printf("\np[%d][0] = %.3e",i,C->p[i][0]);
		printf("\np[%d][1] = %.3e",i,C->p[i][1]);
		printf("\np[%d][2] = %.3e\n",i,C->p[i][2]);
	}

	printf("\n");
	return;
}

/*--- chain_avg: find the average chain distance ---*/
double chain_avg(CHAIN * C)
{
	double ave=0;

	for(int i = 0; i < C->N-1; i++)
	{
		ave += mag(C->R[i],3);
	//	printf("\nmagR[%d] = %.3e",i,mag(C->R[i],3));
	}


	return ave / (C->N-1);
}

