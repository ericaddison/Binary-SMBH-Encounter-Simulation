/********************************************************************
*********************************************************************

    BEMRI Tidal Disruption Simulation

    Author: Eric Addison
    Initial Date: 8 Feb 2010
    Affiliation: USU Astrophysics

    File: BEMRI_functions.c
    Purpose: Auxiliary functions needed for BEMRI simulation

**********************************************************************
*********************************************************************/

#include "BEMRI.h"

/*************************************************************************/
/*************************************************************************/
/**																		**/
/** Function: two_body_energy											**/
/** Purpose: calculate the instantaneous energy of a two body orbit		**/
/**	Inputs: struct of binary parameters             					**/
/**			                            								**/
/**	Output: double, value of energy                                     **/
/**																		**/
/*************************************************************************/
/*************************************************************************/

double two_body_energy(double * r1, double * r2, double * v1, double * v2, double M1, double M2)
{

 double x, y, z, vx, vy, vz, r, mu;
//reduced mass
mu = M1*M2/(M1+M2);

//Go into mass 1 rest frame and calculate E
x = r2[0] - r1[0];
y = r2[1] - r1[1];
z = r2[2] - r1[2];
vx = v2[0] - v1[0];
vy = v2[1] - v1[1];
vz = v2[2] - v1[2];
r = sqrt(x*x + y*y + z*z);

//total energy
return 0.5 * mu * (vx*vx + vy*vy + vz*vz) - G*M1*M2/r;

}

/*************************************************************************/
/*************************************************************************/
/**																		**/
/** Function: angular_momentum											**/
/** Purpose: calculate the instantaneous angular momentum               **/
/**          of a two body orbit		                                **/
/**	Inputs: struct of binary parameters, vector of theta values 		**/
/**			                            								**/
/**	Output: double, value of ang. mom.                                  **/
/**																		**/
/*************************************************************************/
/*************************************************************************/

double angular_momentum(double * r1, double * r2, double * v1, double * v2, double M1, double M2)
{

    double x, y, z, vx, vy, vz, r, mu, lx, ly, lz;
//reduced mass
mu = M1*M2/(M1+M2);

//Go into mass 1 rest frame and calculate l
x = r2[0] - r1[0];
y = r2[1] - r1[1];
z = r2[2] - r1[2];
vx = v2[0] - v1[0];
vy = v2[1] - v1[1];
vz = v2[2] - v1[2];
r = sqrt(x*x + y*y + z*z);

    lx = y*mu*vz - z*mu*vy;
    ly = z*mu*vx - x*mu*vz;
    lz = x*mu*vy - y*mu*vx;

    return sqrt(lx*lx + ly*ly + lz*lz);

}

/*************************************************************************/
/*************************************************************************/
/**																		**/
/** Function: angular_momentum											**/
/** Purpose: calculate the instantaneous angular momentum               **/
/**          of a two body orbit		                                **/
/**	Inputs: struct of binary parameters, vector of theta values 		**/
/**			                            								**/
/**	Output: double, value of ang. mom.                                  **/
/**																		**/
/*************************************************************************/
/*************************************************************************/

double angular_momentum_vector(double * r1, double * r2, double * v1, double * v2, double M1, double M2, double * l)
{

    double x, y, z, vx, vy, vz, r, mu, lx, ly, lz;
//reduced mass
mu = M1*M2/(M1+M2);

//Go into mass 1 rest frame and calculate l
x = r2[0] - r1[0];
y = r2[1] - r1[1];
z = r2[2] - r1[2];
vx = v2[0] - v1[0];
vy = v2[1] - v1[1];
vz = v2[2] - v1[2];
r = sqrt(x*x + y*y + z*z);

    lx = y*mu*vz - z*mu*vy;
    ly = z*mu*vx - x*mu*vz;
    lz = x*mu*vy - y*mu*vx;

    l[0] = lx;
    l[1] = ly;
    l[2] = lz;

    return sqrt(lx*lx + ly*ly + lz*lz);

}


/*************************************************************************/
/*************************************************************************/
/**																		**/
/** Function: calc_angles												**/
/** Purpose: calculate the instantaneous orbital angle parameters       **/
/**	Inputs: struct of binary parameters, output vector			 		**/
/**			                            								**/
/**	Output: vector of angles		                                    **/
/**																		**/
/*************************************************************************/
/*************************************************************************/

int calc_angles(binary * b)
{

	double x,y,z,vx,vy,vz,r,v,k,h[3],R[3],V[3],e_vec[3],vxh1[3],vxh2[3],R_temp[3],theta,edotr,ip;
	double i,z_hat[3] = {0,0,1}, ea, lan_calc, aop_calc, n[3];
	double beta;
	int ret=0;

	//First translate to mass 1 rest frame
	x = b->x2[0] - b->x1[0];
	y = b->x2[1] - b->x1[1];
	z = b->x2[2] - b->x1[2];
	vx = b->v2[0] - b->v1[0];
	vy = b->v2[1] - b->v1[1];
	vz = b->v2[2] - b->v1[2];
	r = sqrt(x*x + y*y + z*z);
	v = sqrt(vx*vx + vy*vy + vz*vz);
	k = G*b->mass1 + G*b->mass2;							//here k is the standard gravitational parameter

	//put together vectors of r and v
	R[0] = x;
	R[1] = y;
	R[2] = z;
	V[0] = vx;
	V[1] = vy;
	V[2] = vz;

	//calculate specific angular momentum vector h
	cross_prod(R,V,h);

	//calculate the eccentricity vector
	cross_prod(V,h,vxh1);
	vect_mult_scalar(vxh1, 1/k, vxh2, 3);
	vect_mult_scalar(R,-1/r,R_temp,3);

	vect_add(vxh2,R_temp,e_vec,3);

	//find true anomaly from r
	edotr = dot_prod(e_vec,R,3);

	theta = acos(edotr / ( mag(e_vec,3) * mag(R,3) ) );

	if(dot_prod(R,V,3) < 0)
	{
		theta = 2*PI - theta;
		ret = 1;
	}

	b->true_anomaly = theta;

	i = acos(h[2]/mag(h,3));

	//calculate n, unit normal vector
	vect_mult_scalar(h,1/mag(h,3),n,3);

	lan_calc = atan2(n[0],-n[1]);

	lan_calc = fix_angle(lan_calc);		//keep angles in range 0 to 2*pi


	if(i==0)				//zero inclination makes lan meaningless. Keep it at zero
		lan_calc = 0;

	beta = atan2(-cos(i) * sin(lan_calc) * x + cos(i) * cos(lan_calc) * y + sin(i) * z , cos(lan_calc) * x + sin(lan_calc) * y);
	beta = fix_angle(beta);		//keep angles in range 0 to 2*pi


	aop_calc = PI + beta - theta;
	aop_calc = fix_angle(aop_calc);		//keep angles in range 0 to 2*pi


	ea = acos((mag(e_vec,3) + cos(theta)) / (1.0 + mag(e_vec,3)*cos(theta)));

	b->i_a = i;
	b->la_n = lan_calc;
	b->ao_p = aop_calc;

	
	return ret;
}


/*************************************************************************/
/*************************************************************************/
/**																		**/
/** Function: true_anomaly												**/
/** Purpose: calculate the true anomaly of an orbit				        **/
/**	Inputs: struct of binary parameters							 		**/
/**			                            								**/
/**	Output: none, overwrite value in struct                             **/
/**																		**/
/*************************************************************************/
/*************************************************************************/

int true_anomaly(binary * b)
{
	double x,y,z,vx,vy,vz,r,v,k,h[3],R[3],V[3],e_vec[3],vxh1[3],vxh2[3],R_temp[3],theta,edotr;
	int rv=1;
	
	//First translate to mass 1 rest frame
	x = b->x2[0] - b->x1[0];
	y = b->x2[1] - b->x1[1];
	z = b->x2[2] - b->x1[2];
	vx = b->v2[0] - b->v1[0];
	vy = b->v2[1] - b->v1[1];
	vz = b->v2[2] - b->v1[2];
	r = sqrt(x*x + y*y + z*z);
	v = sqrt(vx*vx + vy*vy + vz*vz);
	k = G*b->mass1 + G*b->mass2;							//here k is the standard gravitational parameter G(m1 + m2)

	//put together vectors of r and v
	R[0] = x;
	R[1] = y;
	R[2] = z;
	V[0] = vx;
	V[1] = vy;
	V[2] = vz;

	//calculate vector h
	cross_prod(R,V,h);

	//calculate the eccentricity vector
	cross_prod(V,h,vxh1);
	vect_mult_scalar(vxh1, 1/k, vxh2, 3);
	vect_mult_scalar(R,-1/r,R_temp,3);

	vect_add(vxh2,R_temp,e_vec,3);

	//find true anomaly from r
	edotr = dot_prod(e_vec,R,3);

	theta = acos(edotr / ( mag(e_vec,3) * mag(R,3) ) );

	if(dot_prod(R,V,3) < 0)
	{
		rv = 0;
		//printf("\nold theta = %.3f\nnew theta = %.3f\n",theta,2*PI-theta);
		theta = 2*PI - theta;
	}

	b->true_anomaly = theta;

	return rv;

}


/*************************************************************************/
/*************************************************************************/
/**																		**/
/** Function: CM_values													**/
/** Purpose: calculate the position and velocity of the BEMRI CM        **/
/**	Inputs: two struct of binary parameters						 		**/
/**			                            								**/
/**	Output: none, overwrite value in struct                             **/
/**																		**/
/*************************************************************************/
/*************************************************************************/

void CM_values(BEMRI * b)
{


// first find CM coordinates
	for(int i=0;i<3;i++)
	{
		R4[i] = (R2[i]*m2 + R1[i]*m1) / (m2 + m1);
	}

// now find CM velocity
	for(int i=0;i<3;i++)
	{
		V4[i] = (V2[i]*m2 + V1[i]*m1) / (m2 + m1);
	}


	return;
}


/*************************************************************************/
/*************************************************************************/
/**																		**/
/** Function: advance_orbit												**/
/** Purpose: advance a binary orbit to a new value of true anomaly      **/
/**	Inputs: struct of binary parameters,  angular momentum, 			**/
/** 			new true anomaly value	 								**/
/**			                            								**/
/**	Output: none, overwrite values in struct                            **/
/**																		**/
/*************************************************************************/
/*************************************************************************/

void advance_orbit(binary * b, double l, double theta)
{

	double cm[3], cmv[3], x1[3], x2[3], r_new, e,a, v1[3],v2[3], px,py;
	e = b->ecc;
	a = b->a;

	
	// first find CM coordinates and velocities
		for(int i=0;i<3;i++)
		{
			cm[i] = (b->x2[i]*b->mass2 + b->x1[i]*b->mass1) / (b->mass2 + b->mass1);
			cmv[i] = (b->v2[i]*b->mass2 + b->v1[i]*b->mass1) / (b->mass2 + b->mass1);
		}

	//find new value for separation
		r_new = b->rp*(1+e) / (1 + e*cos(theta));

	//these coordinates assume the orbit is in the xy-plane with the CM at the origin
		x1[0] = (b->mass2/b->M) * r_new * cos(theta);
		x1[1] = (b->mass2/b->M) * r_new * sin(theta);
		x1[2] = 0;

		x2[0] = -(b->mass1/b->M) * r_new * cos(theta);
		x2[1] = -(b->mass1/b->M) * r_new * sin(theta);
		x2[2] = 0;

	//restore orientation using current orbital angles
		vec_rot(&x1,2,(b->ao_p));
		vec_rot(&x2,2,(b->ao_p));		//rotate by negative perihelion angle about z-axis

		vec_rot(&x1,0,(b->i_a));
		vec_rot(&x2,0,(b->i_a));		//rotate by negative inclination angle about x-axis

		vec_rot(&x1,2,(b->la_n));
		vec_rot(&x2,2,(b->la_n));		//rotate by negative lan angle about z-axis

	//restore proper position by moving out of CM frame and store results
		for(int i=0;i<3;i++)
		{
			b->x1[i] = x1[i] + cm[i];
			b->x2[i] = x2[i] + cm[i];
		}

	//now the velocities need to be adjusted

	//in the plane, the x and y momenta are given by (see my BEMRI notes, 19Aug10)

		px = -(l/r_new)*sin(theta)/(1+e*cos(theta));					//x momentum
		py = (l/r_new)*(e+cos(theta))/(1+e*cos(theta));				//y momentum


	//so the velocities are then:
		v1[0] = px/b->mass1;
		v1[1] = py/b->mass1;
		v1[2] = 0;

		v2[0] = -px/b->mass2;
		v2[1] = -py/b->mass2;
		v2[2] = 0;

	//now restore proper orientation
		vec_rot(&v1,2,(b->ao_p));
		vec_rot(&v2,2,(b->ao_p));		//rotate by negative perihelion angle about z-axis

		vec_rot(&v1,0,(b->i_a));
		vec_rot(&v2,0,(b->i_a));		//rotate by negative inclination angle about x-axis

		vec_rot(&v1,2,(b->la_n));
		vec_rot(&v2,2,(b->la_n));		//rotate by negative lan angle about z-axis

		for(int i=0;i<3;i++)
		{
			b->v1[i] = v1[i] + cmv[i];
			b->v2[i] = v2[i] + cmv[i];
		}

		return;
}

/*************************************************************************/
/*************************************************************************/
/**																		**/
/** Function: advance_orbit_set_up										**/
/** Purpose: set up for and call the advance_orbit function             **/
/**	Inputs: threshold value alpha, BEMRI structure b			 		**/
/**			                            								**/
/**	Output: double, the time elapsed during the skip			        **/
/**																		**/
/*************************************************************************/
/*************************************************************************/

double advance_orbit_set_up(BEMRI * b, double chi, double r_current, double dt1)
{
	double t_current, t_new, psi_h_target, theta_h_target, psi_b_advanced, Mb, t_elapsed, new_true;

		//find current hole eccentric anomaly
		psi_h = acos( (eh + cos(theta_h)) / (1+eh*cos(theta_h)) );

		//find current time in the orbit (since last periapse)
		t_current = Ph/(2*PI) * (psi_h - eh*sin(psi_h));

		//now find target (new) values
		psi_h_target = 2*PI - acos( (1 - chi/ah) / eh);					//hole eccentric anomaly value to skip to
		theta_h_target = fix_angle( 2*atan2(sqrt(1+eh)*sin(psi_h_target/2.0) , sqrt(1-eh)*cos(psi_h_target/2)) );	//hole true anomaly value to skip to
		t_new = Ph/(2*PI) * (psi_h_target - eh*sin(psi_h_target));		//new time value from Kepler's equation

		//then the elapsed time is given by:
		t_elapsed = t_new - t_current;

		//printf("\nt_elapsed = %.5e",t_elapsed);

//_+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++


		//force t_elapsed to equal an integer multiple of dt1

		t_elapsed = floor(t_elapsed/dt1) * dt1;		//now t_elapsed is good
		t_new = t_elapsed + t_current;				//now we're recalculating the hole orbit anomalies
		//printf("\nMh = %.15e",2*PI/Ph*t_new);
		psi_h_target = solve_kepler(2*PI/Ph*t_new,eh);	//here is the new target eccentric anomaly

		theta_h_target = fix_angle( 2*atan2(sqrt(1+eh)*sin(psi_h_target/2.0) , sqrt(1-eh)*cos(psi_h_target/2.0)) );	//new hole true anomaly value to skip to
		//printf("\npsi_h_target = %.3e\nsolving keplers eqn for t_new = %.3f",psi_h_target, Ph/(2*PI)*kepler(psi_h_target,eh));
		//printf("\nt_elapsed / dt1 = %.20f\n",t_elapsed / dt1);
//_+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++



		//now the mean anomaly of the BEMRI orbit will have advanced by:
		Mb = 2*PI/Pb*t_elapsed - 2*PI*floor(t_elapsed/Pb);	//note this is not the current mean anomaly of the BEMRI, but how much it has advanced

		//here is how much the eccentric anomaly of the BEMRI has advanced
		psi_b_advanced = solve_kepler(Mb,eb);

		//so the new value of the true anomaly of the BEMRI is given by:		(current + advanced)
		new_true = theta_b + fix_angle((2*atan2(sqrt(1+eb)*sin(psi_b_advanced/2.0),sqrt(1-eb)*cos(psi_b_advanced/2.0))) );



//_+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++


		//printf("\ncurrent hole anomaly is %.3f, target anomaly is %.3f\n",theta_h, theta_h_target);
		//printf("chi = %.3e\nah = %.3e\neh=%.3e\n",chi,ah,eh);


//_+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++


	//before advancing the BH orbit, freeze the BEMRI CM motion so the new CM motion can be added on easily
	for(int i=0; i < 3; i++)
	{
		V1[i] -= V4[i];			//remove CM motion for mass 1
		V2[i] -= V4[i];			//remove CM motion for mass 2
		R1[i] -= R4[i];
		R2[i] -= R4[i];
	}	

	    //now that we know the new true anomaly of the hole AND the BEMRI, we can advance the two orbits to the new positions
		advance_orbit(&(b->binary_b),lb,new_true);

		//now advance the BH orbit
		advance_orbit(&(b->binary_h),lh,theta_h_target);

		//and now the BEMRI CM values need to be updated:
		//BEMRI_CM_update(b);

		return t_elapsed;
}



/*************************************************************************/
/*************************************************************************/
/**																		**/
/** Function: calc_energies												**/
/** Purpose: calculate the internal energies for each orbit             **/
/**	Inputs: BEMRI structure, pointers for other energies		 		**/
/**			                            								**/
/**	Output: void						                                **/
/**																		**/
/*************************************************************************/
/*************************************************************************/
void calc_energies(BEMRI * b, double * E1h, double * E2h)
{
    Eb = two_body_energy(R1,R2,V1,V2,m1,m2);        	//calculate energy between BEMRI
    Eh = two_body_energy(R3,R4,V3,V4,m3,(m2+m1));       //CM-BH energy
    *E1h = two_body_energy(R3,R1,V3,V1,m3,m1);        	//m1-BH energy
    *E2h = two_body_energy(R3,R2,V3,V2,m3,m2);        	//m2-BH energy
    return;
}

/*************************************************************************/
/*************************************************************************/
/**																		**/
/** Function: calc_angular_momenta										**/
/** Purpose: calculate the angular momentum for each orbit             	**/
/**	Inputs: BEMRI structure, pointers for other momenta			 		**/
/**			                            								**/
/**	Output: void						                                **/
/**																		**/
/*************************************************************************/
/*************************************************************************/
void calc_angular_momenta(BEMRI * b, double * l1h, double * l2h)
{
	lb = angular_momentum(R1,R2,V1,V2,m1,m2);			//angular momentum of BEMRI
	lh = angular_momentum(R3,R4,V3,V4,m3,(m2+m1));		//CM-BH angular momentum
	*l1h = angular_momentum(R3,R1,V3,V1,m3,m1);			//m1-BH angular momentum
	*l2h = angular_momentum(R3,R2,V3,V2,m3,m2);			//m2-BH angular momentum
    return;
}


/*************************************************************************/
/*************************************************************************/
/**																		**/
/** Function: calc_ecc													**/
/** Purpose: calculate the eccentricities for each orbit             	**/
/**	Inputs: BEMRI structure, pointers for other eccentricities	 		**/
/**			                            								**/
/**	Output: void						                                **/
/**																		**/
/*************************************************************************/
/*************************************************************************/
void calc_ecc(BEMRI * b, double * e1h, double * e2h, double E1h, double E2h, double l1h, double l2h)
{
    eb   = L_E_to_e(Eb,lb,m1,m2); //sqrt(fabs(1+2*Eb*lb*lb/(mu_b*pow(G*m1*m2,2))));					//calculate BEMRI eccentricity
    eh   = L_E_to_e(Eh,lh,m3,m1+m2); //sqrt(fabs(1+2*Eh*lh*lh/(mu_h*pow(G*m3*(m2+m1),2))));				//CM-BH eccentricity
    *e1h = L_E_to_e(E1h,l1h,m1,m3); //sqrt(fabs(1+2*E1h*l1h*l1h/((m3*m1)/(m3+m1)*pow(G*m3*m1,2))));	//m1-BH eccentricity
    *e2h = L_E_to_e(E2h,l2h,m2,m3); //sqrt(fabs(1+2*E2h*l2h*l2h/((m3*m2)/(m3+m2)*pow(G*m3*m2,2))));	//m2-BH eccentricity
    return;
}


/*************************************************************************/
/*************************************************************************/
/**																		**/
/** Function: calc_total_energy											**/
/** Purpose: calculate the total energy of the system	            	**/
/**	Inputs: BEMRI structure, pointer for TE value				 		**/
/**			                            								**/
/**	Output: void						                                **/
/**																		**/
/*************************************************************************/
/*************************************************************************/
void calc_total_energy(BEMRI * b)
{
	double K=0,U=0;

	for(int i=0;i<3;i++)
		K += 0.5*( m1*V1[i]*V1[i] + m2*V2[i]*V2[i] + m3*V3[i]*V3[i] );   		//kinetic energy

	U -= G*( m1*m2/sqrt( pow((X1-X2),2) + pow((Y1-Y2),2) + pow((Z1-Z2),2)) );	//potential energy
	U -= G*( m1*m3/sqrt( pow((X1-X3),2) + pow((Y1-Y3),2) + pow((Z1-Z3),2)) );	//potential energy
	U -= G*( m3*m2/sqrt( pow((X3-X2),2) + pow((Y3-Y2),2) + pow((Z3-Z2),2)) );	//potential energy

	(*b).energy = K + U;

	//Mdisplay(3,3,V1);
	//printf("\nK = %.3e\nU = %.3e\n",K,U);
	//anykey();
	
	return;
}


/*************************************************************************/
/*************************************************************************/
/**																		**/
/** Function: calc_total_angular_momentum								**/
/** Purpose: calculate the total ang. mom. of the system            	**/
/**	Inputs: BEMRI structure, pointer for TE value				 		**/
/**			                            								**/
/**	Output: void						                                **/
/**																		**/
/*************************************************************************/
/*************************************************************************/
void calc_total_angular_momentum(BEMRI * b)
{

	double p1[3],p2[3],p3[3],l1[3],l2[3],l3[3],temp[3];

	(*b).L = 0;

		vect_mult_scalar(V1,m1,p1,3);
		vect_mult_scalar(V2,m2,p2,3);
		vect_mult_scalar(V3,m3,p3,3);

		cross_prod(R1,p1,l1);
		cross_prod(R2,p2,l2);
		cross_prod(R3,p3,l3);

		vect_add(l1,l2,temp,3);
		vect_add(l3,temp,temp,3);
		(*b).L  = mag(temp,3);

	return;
}

/*************************************************************************/
/*************************************************************************/
/**																		**/
/** Function: sign														**/
/** Purpose: determine the sign of a number                             **/
/**	Inputs: double x											 		**/
/**			                            								**/
/**	Output: 1 or -1						                                **/
/**																		**/
/*************************************************************************/
/*************************************************************************/

double sign(double x)
{
	if(x >= 0)
		return 1;
	else
		return -1;
}

/*************************************************************************/
/*************************************************************************/
/**																		**/
/** Function: fix_angle													**/
/** Purpose: keeps an angle between 0 and 2*PI                          **/
/**	Inputs: double x											 		**/
/**			                            								**/
/**	Output: new angle between 0 and 2*PI                                **/
/**																		**/
/*************************************************************************/
/*************************************************************************/

double fix_angle(double x)
{
	double theta;

	theta = x;

	while(theta < 0)
		theta += 2*PI;

	while(theta >= 2*PI)
		theta -= 2*PI;

	return theta;
}


/*************************************************************************/
/*************************************************************************/
/**																		**/
/** Function: peters_lifetime											**/
/** Purpose: compute the Peter's lifetime of a binary                   **/
/**	Inputs: binary structure									 		**/
/**			                            								**/
/**	Output: double lifetime in seconds	                                **/
/**																		**/
/*************************************************************************/
/*************************************************************************/

double peters_lifetime(double ecc, double a, double M1, double M2)
{

	double c0, beta, e, de, lifetime=0;
	int N = 1e5;

	beta = (64.0/5.0)*(G*G*G*M1 * M2 * (M1 + M2)) / pow(c,5.0);

	if(ecc < 1e-5)
		return pow(a,4)/(4.0*beta);

	c0 = a*(1.0-ecc*ecc)/pow(ecc,12.0/19.0)/pow(1.0+(121.0/304.0)*ecc*ecc , 870.0/2299.0);

	de = ecc/N;

//printf("\nc0 = %.3e\nbeta = %.3e\n\n",c0,beta);
//printf("circular lifetime a0^4/(4*beta) = %.5e\n",pow(a,4)/(4.0*beta));
	for(int ii = 0; ii <= N; ii++)
	{
		e = ii*de;
		lifetime += RK4_peters(e,de);
	}

	lifetime *= (12.0/19.0)*pow(c0,4.0)/beta;

//printf("\nLifetime = %.10e\n\n",lifetime);

//anykey();

return lifetime;

}

/*************************************************************************/
/*************************************************************************/
/**																		**/
/** Function: lifetime_integrand										**/
/** Purpose: integrand in Peter's lifetime integral	                    **/
/**	Inputs: eccentricity e										 		**/
/**			                            								**/
/**	Output: double function result		                                **/
/**																		**/
/*************************************************************************/
/*************************************************************************/

double lifetime_integrand(double ecc)
{
	return pow(ecc,29.0/19.0)*pow( 1 + (121.0/304.0)*ecc*ecc, 1181.0/2299.0) / pow(1 - ecc*ecc , 3.0/2.0);
}


/*************************************************************************/
/*************************************************************************/
/**																		**/
/** Funtion: RK4_peters                                                 **/
/** Purpose: Runge Kutta ODE stepper for lifetime integral              **/
/**	Inputs: time step dt                                                **/
/**	Output: double integral value					                    **/
/**																		**/
/*************************************************************************/
/*************************************************************************/

double RK4_peters(double x, double dx)
{
    double yk1,yk2,yk3,yk4;

    yk1 = lifetime_integrand(x);			//RK step 1

    yk2 = lifetime_integrand(x+dx/2);       //RK step 2

    yk3 = lifetime_integrand(x+dx/2);       //RK step 3

    yk4 = lifetime_integrand(x+dx);       	//RK step 4

    //return update
    return (yk1 + 2*yk2 + 2*yk3 + yk4)*dx/6.0;


}

/*************************************************************************/
/*************************************************************************/
/**																		**/
/** Funtion: torque		                                                **/
/** Purpose: calculate the torque on the BEMRI from the SMBH            **/
/**	Inputs: BEMRI structure b, torque vector                            **/
/**	Output: fill torque vector						                    **/
/**																		**/
/*************************************************************************/
/*************************************************************************/
void torque(BEMRI * b, double * t)
{
	double F1[3] = {(X3-X1),(Y3-Y1),(Z3-Z1)} , F2[3] = {(X3-X2),(Y3-Y2),(Z3-Z2)};
	double r1h = sqrt(pow((X3-X1),2) + pow((Y3-Y1),2) + pow((Z3-Z1),2));
	double r2h = sqrt(pow((X3-X2),2) + pow((Y3-Y2),2) + pow((Z3-Z2),2));
	double r1cm[3] = {(X4-X1),(Y4-Y1),(Z4-Z1)} , r2cm[3] = {(X4-X2),(Y4-Y2),(Z4-Z2)};
	double t1[3], t2[3];

	for(int ii = 0; ii < 3; ii++)
	{
		F1[ii] *= G*m1*m3/pow(r1h,3);		//compute forces on each component from SMBH
		F2[ii] *= G*m2*m3/pow(r2h,3);
	}

	cross_prod(r1cm, F1, t1);				//compute the torque from R cross F
	cross_prod(r2cm, F2, t2);

	vect_add(t1, t2, t, 3);					//add the individual torques and store in t
/*
	printf("\ntorque is:");
	Vdisplay(t,3);
	Vdisplay(F1,3);
	Vdisplay(r1cm,3);
	Vdisplay(F2,3);
	Vdisplay(r2cm,3);
	anykey();
*/
}

