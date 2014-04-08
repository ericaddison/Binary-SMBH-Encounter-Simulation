/********************************************************************
*********************************************************************

    BEMRI Tidal Disruption Simulation

    Author: Eric Addison
    Initial Date: 8 Feb 2010
    Affiliation: USU Astrophysics

    File: main.c
    Purpose: Main program execution file. Contains main() and calls
             all supporting functions.

*********************************************************************
*********************************************************************/


/********************************************************************
*********************************************************************

    STATUS:

        AS OF 10Dec2010, code has been streamlined, command line 
		arguments have been added.

*********************************************************************
*********************************************************************/
#include <stdio.h>
#include "BEMRI.h"
double max(double x, double y);
int fcount = 0;

int main(int argc, char *argv[])
{


/* Clear Screen */
	system("clear");
	printf("\nBEMRI simulator!\n");
	
/* Variable Declarations */	
    BEMRI *b;
    flags f;
    iParams iPars;
    double dt1, dt, tmax, e[3], t, t0, r_min, t_RK4 = 0.0,tau, t_node, sim_time, E0, PbN, dt_test;
    double r_bemri, r_current, r_node, r_1, r_2, theta_N, ecc_N, theta0h, r_tid, r0;
    double E1h,E2h,l1h,l2h,e1h,e2h,dAdt, dt_temp, dt_min, fmax[2], alpha = 1e-5, chi, node_phase, E_last, a1h, a2h;
    double Ph_init,eh_init, Q[3][3], Qtt[3][3], T_pm[4];
	double R[3], m[3], timer=1e300;
    double test_rad, test_angle, del_ia, del_aop, frac = 1.0, max_hard_count, eh0;
    int xed, a = 0, N=1, single, minned, status, entered_node, passed_node, times_around, max_orbits, num_params;
    int m1NegY=0, hard_count, agent, amax, astep, nvar, zeroed, angle_counter = 0, inner_runs, ipsval = 1, past;
    double t_broken, Pb0, Esys0,Lsys0, eb0, Hrat = 3, ri, r_cm, percentage_as_double, gam_now, eb_prev, ab_prev, k1, k2;
    double e1h_prev, e2h_prev, timer0;
    long seed, iseed;
	double y[22],yscal[22],dydx[22], dscale, mconv, tconv;	//these are the vectors used for BS integration, 22 = 7*N+1
	double lan0, ia0, aop0, ab0;
	char fname[80];
	


/* for BS ODE integrator setup */
	for(int i=0;i<22;i++)
		yscal[i] = 1;
	
/* Vector from observer to system center of mass */	
	R[0] = 0; R[1] = 0; R[2] = 2.45027686e20;	//distance to Sgr A*

/* timing variables */	
    time_t start;
    time_t stop;
    struct timeval t1;

/* command line argument check */
	//possible flags:
	//	-c				run to completion
	//	-s [value]		use value for seed
	//	-t [value]		run to tmax = value*Ph
	//	-nd				do not stop sim when BEMRI is disrupted
	//	-N [value]		run simulation N times, do not output position or quadrupole data
	//	-so				suppress any screen output
    //  -PbN [value]	set Pb = value*tnode
    //	-th0 [value]	set theta0 = value
    //	-RK4			use RK4 integrator instead of BS2
    //	-fname [string]	use string as file name instead of BEMRI.dat
    //	-fpars [string] use string as filename to find input parameters
    //	-beta [value]	use beta = value
	//	-gam [value]	use gamma = value
    //	-ips [value]	run [value] runs per set of parameters, sampling initial phases
    //	-eh [value]		set BH orbit eccentricity
	//	-Htest			Test Heggie values
	//	-Hrat [value]	Use [value] for rp/a ratio in Heggie test
	//	-geo			Use geometricized units -- G = c = 1
	//	-ang			Randomize binary orientation angles
	
    b = (BEMRI *) malloc( sizeof(BEMRI) );
	if(args(argc,argv,&iPars,&f,&N,&frac,&seed,&PbN,fname,&eh,&ipsval,&Hrat))
		return 1;
	

/* initialize BEMRI parameters */
    init_params(b,f);

/* file pointer declarations */	
    FILE *results_fp, *Q_fp, *pos_fp;
    if(f.fnameflag == 0)
    	results_fp = fopen("BEMRI.dat","w");
    else
    	results_fp = fopen(fname,"w");

	if(f.cflag != 1 && f.Nflag != 1)
	{
		Q_fp = fopen("QP.dat","w");
		pos_fp = fopen("pos.dat","w");
	}

/* Chain declaration */
	C = (CHAIN *) malloc( sizeof(CHAIN) );
	
/* set seed if not given on command line */
	if(f.sflag == 0)
	{gettimeofday(&t1,NULL); seed = t1.tv_usec;}
	iseed = seed;		//save initial seed to output
	
	
	
/* master loop, changes parameters */
/* this is the total number of unique (gamma,inc,lan) values that will be run */
for(int ii=0;ii<N;ii++)
    {
		
		/* randomly assign angles if -ang flag is used */
		if(f.angflag)							// angle variation
		{
			iPars.lanA = 2*PI*ran2(&seed);
			iPars.incA = acos(2*ran2(&seed)-1);
		}
		/* if not using random angles OR a parameter file, then assign angles of zero */
		else if(!f.fParflag)
			iPars.lanA = iPars.incA = 0.0;
		
		/* assign random gamma value unless certain flags are used */ 
		if(!f.fParflag && !f.betaflag && !f.gamflag){
			iPars.gamma = (0.35 + 4.65*ran2(&seed));  // randomly assigning beta from being uniform over gamma
			iPars.beta = 1.0/iPars.gamma;
		}

/* starting the loop over all theta0 initial binary phase values */
  for(angle_counter = 0; angle_counter<ipsval; angle_counter++)
    {

/* set G and c, which will be reset if -geo is used */		
	c = 299792458e0;
	G = 6.673e-11;
		
/* reset timer and m1NegY */
	timer = 1e300;
	m1NegY = 0;

/* parameter setup */
    	init_params(b,f);	// initialize BEMRI parameters
	eb = 0.0;
	aop = 0.0;			// aop is redundant for circular orbits, always set to 0	
	lan = iPars.lanA;	// set dynamic value lan
	ia = iPars.incA;	// set dynamic value ia
	if(!f.th0flag && !f.fParflag && f.ipsflag)	// assign appropriate value of th0 if using -ips flag
		iPars.th0 = 2*PI * (float)angle_counter/ipsval;		//only set theta0 if th0flag has NOT been set

	
/* Setup for testing results from Heggie and Rasio paper */			
	if(f.Htflag == 1)
		{
			G = c = 1;
			ia = 0;
			lan = 0;
			aop = 0;
			m1 = m2 = 1;
			m3 = 1;
			ab = 1;
			eb = 0.0;
			eh = 1.0;
			rph = Hrat*ab;
			r0 = 100*rph;
			theta0h = -acos(2*rph/r0 - 1);
			Pb = 2*PI*sqrt(ab*ab*ab/(G*(m1+m2)));
		}
	
/* setup for elliptical BEMRI */
	else if(eh < 1)
		theta0h = PI;
		

/* Setup for parabolic orbit evolution */	
	else if(eh >= 1)
	{
		if(f.geoflag)								// if using geometricized units
		{
			mconv = G/(c*c);						// conversion factor for masses
			tconv = c;
			G = c = 1.0;							// set constants to 1
			//ab /= dscale;
			m1 *= mconv;							// convert mass values to meters
			m2 *= mconv;
			m3 *= mconv;
			Pb *= tconv;							// recalculate the binary period with new ab and masses
		}

		r_tid = ab*pow(m3/(m1+m2),1.0/3.0);			// calculate tidal radius r_tid
		rph = r_tid/iPars.beta;						// calculate pericenter distance for BEMRI orbit
		r0 = 200.0*rph;								// calculate r0 = initial separation
		theta0h = -acos(2.0*rph/r0 - 1.0);			// calculate corresponding true anomaly
		if(isnan(theta0h) || fabs(theta0h) < 2.0)	// see if the anomaly was in acceptable range
			theta0h = -2.0;							// if not, just set to -2.0 radians
		//theta0h = -3.0;
	}


/* recalculate BEMRI parameters using new values */		
	recalc_params(b);
	
/* initial conditions and node passage time */
    lh = sqrt(rph*(1+eh)*G*(m3*m3*m4*m4)/(m3+m4));		//initial angular momentum of SMBH orbit
    advance_orbit(&(b->binary_h),lh,theta0h);			//begin cm-SMBH orbit at theta0h
   
/* timestep tau and max simulation time tmax */ 
	if(f.Htflag == 1)
	{
		tau = Pb/pow(2,7);
		tmax = 1e100;
	}
    else if(eh<1)
    {
		tau = Pb/pow(2,7);
		tmax = 1*Pb;
    }
    else
    {
    	tau = Pb/pow(2,7);
		tmax = 1e100;		  
    }


    recalc_params(b);
// update little binary positions and velocities for orbit around hole
    lb = sqrt(rpb*(1+eb)*G*(m1*m1*m2*m2)/(m1+m2));		//initial angular momentum of BEMRI
    advance_orbit(&(b->binary_b),lb,iPars.th0);						//begin BEMRI at theta0
	BEMRI_CM_update(b);
    load_vectors(rr,vv,m,b);		//load initial values into working vectors
    dt = tau;
    dt1 = dt;
			   
/* Chain setup and testing */
	setup_chain(C,m,3);
	nvar = 6*(C->N-1);
		
/* BEMRI lifetime check */			
    if ( peters_lifetime(eb,ab,m1,m2) < 100*Ph )			//BEMRI lifetime too short, abort current run
    {status = 3;times_around = -1;}

		
/* other setup */
    eh0 = eh;
    t0 = 0;
    zeroed = 0;
    t_broken = 0.0;
	xed = 0;
	minned = 0;
	t_RK4 = 0.0;
	status = 0;
	entered_node = 0;
	passed_node = 0;
	times_around = 0;
	max_orbits = 1;
    a = -1;
	if(eh < 1)
		amax = (int)(abs((tmax-t0)/dt1));
	else
		amax = 1e6;
	astep = 1;//(int)ceil((amax/200000.0));		//output management
    start = time(NULL);
	max_hard_count = 20;
	hard_count = 0;
	Pb0 = Pb;
	ri = sqrt( pow(X4-X3,2) + pow(Y4-Y3,2) + pow(Z4-Z3,2) );	// initial distance from binary cm to SMBH
	k1 = G*m3*m1;
	k2 = G*m3*m2;
		
		
/* initial energy setup */			
	calc_energies(b,&E1h,&E2h);
    E0 = Eb;
	eb0 = eb;
	ab0 = ab;
    E_last = E0;			//initialize E_last
    calc_total_energy(b);					// compute total system energy
    calc_total_angular_momentum(b);			// compute total system ang. mom.
   	Esys0 = b->energy;
	Lsys0 = b->L;
		
/* print out some diagnostic info */
	if(!f.Nflag)
	{
		printf("\nbeta = %.3f",iPars.beta);
		printf("\ngamma = %.3f",iPars.gamma);
		printf("\nt_max = %.3e",tmax);
		printf("\ntheta0h = %.3e",theta0h);
		printf("\ntheta0 = %.10e",iPars.th0);
		printf("\nab = %.4e",ab);
		printf("\nPb = %.4e",Pb);
		printf("\nr_cm0 = %.4e",ri);
		printf("\nG = %.3e\nc = %.3e",G,c);
		printf("\nia = %.3f\tlan = %.3f\n",ia*180/PI,lan*180/PI);
		
		if(angle_counter==0) anykey();
	}

		
/* main simulation loop */
	while( f.cflag*times_around < 100 && times_around >= 0 && t_RK4 <= tmax )	//something should happen before this, but if not...
		{
		/* this while loop will run under one of two conditions:
			if the -c flag is used, then completion = 1 and tmax = HUGE, so it will run until
			times_around < 100. If the -t or no flag is used, then completion = 0 and the sim
			will run until t_RK4 = tmax. */
			
			gam_now = ah*(1-eh)/(ab*pow(m3/(m1+m2),1.0/3.0));	
			if( !f.Nflag )
			{
				printf("t_RK4 = %.3e\r",t_RK4);
				//printf("Y1 = %.3e\r",Y1);
				fflush(stdout);
			}

	/* upkeep */	
		a += 1;
       	E_last = Eb;		//store old BEMRI energy
		eb_prev = eb;
		e1h_prev = e1h;
		e2h_prev = e2h;
		ab_prev = ab;
			
	/* actual integration call */
		if(!f.bs2flag)
		{
			pack_y_vector_C(C,y);
			leap_derivs2(t_RK4,y,dydx);
			for(int i=0;i<nvar;i++)
				yscal[i]=FMAX(fabs(y[i])+fabs(dydx[i]*dt)+TINY,1);
			//bs_const_step(y,nvar,&t_RK4,dt1,&dt,1e-12,1e-6,yscal,1,leap_derivs2);
			bsstep(y,nvar,&t_RK4,&dt1,1e-12,1e-9,yscal,1,leap_derivs2);
			unpack_y_vector_C(C,y);
			if(a%10 == 0)
				check_chain(C);
			update_momenta(C);		
			update_positions(C);
		}
		else{
			//t_RK4 += N_body_main(rr,vv,m,dt1,3,&dt,1e-8,&minned);	//evolves orbit from time t to t+dt1 with initial step size dt
			pack_y_vector(rr,vv,m,y);
			RK4A_const_step(y,22,&t_RK4,dt1,&dt,1e-9,1e-6,Nbody_derivs1);
			unpack_y_vector(rr,vv,m,y);
		}		
			
        update_binaries(b,rr,vv);								//copies values from v and r into the structures
        CM_values(b);										//calculate values for the center of mass
        if(dt > dt1)										//keeps integration step size at maximum of dt1
			dt = dt1;
		timer -= dt1;										// update timer value	
			
	/* calculate desired values */	
		calc_angles(&(b->binary_b));			// calculate current orbital angles
        past = true_anomaly(&(b->binary_h));	// calculate current true anomaly of BH orbit
        calc_energies(b,&E1h,&E2h);				// compute current binding energies
        calc_angular_momenta(b,&l1h,&l2h);		// compute current pairwise angular momenta
        calc_ecc(b,&e1h,&e2h,E1h,E2h,l1h,l2h);	// compute current pairwise eccentricities
        calc_total_energy(b);					// compute total system energy
        calc_total_angular_momentum(b);			// compute total system ang. mom.
        r_bemri = sqrt( pow(X1-X2,2) + pow(Y1-Y2,2) + pow(Z1-Z2,2));	//BEMRI separation
        r_1 = sqrt( pow(X1-X3,2) + pow(Y1-Y3,2) + pow(Z1-Z3,2) );		// m1-SMBH separation
        r_2 = sqrt( pow(X2-X3,2) + pow(Y2-Y3,2) + pow(Z2-Z3,2) );		// m2-SMBH separation
		r_current = min(r_1,r_2);										//current distance from SMBH to closest BEMRI component
		r_cm = sqrt( pow(X4-X3,2) + pow(Y4-Y3,2) + pow(Z4-Z3,2) );
    	ab = E_to_a(Eb,m1,m2); 					//update semi-major axis of the BEMRI
    	ah = E_to_a(Eh,m3,m1+m2);  				//update semi-major axis of the BH orbit
    	Pb = 2*PI*sqrt(pow(ab,3)/(G*(m1+m2)));	//current orbital periods
    	Ph = 2*PI*sqrt(pow(ah,3)/(G*(m1+m2+m3)));
		b->binary_b.rp = ab*(1-eb);				// compute new binary periapse
		b->binary_h.rp = ah*(1-eh);				// compute new BEMRI periapse
		point_mass_QP(m,rr,3,R,Q,Qtt);			// compute GW output
		if(eh<1) r_node = ah*(1-eh*eh);			//SMBH-BEMRI separation when BEMRI is at the node

			
	/* stopping conditions */
	if(eh0 < 1 || f.ipsflag==0)				// for elliptical BEMRIs or single runs
	{
		if( passed_node && theta_h >= PI )	//if BEMRI has passed the node AND passed apoapse
		{
			passed_node = 0;				//reset passed_node
			times_around += 1;				//increment times_around
           	if( (Eb/E_last) > 1 )			//compare current Eb to last orbit's initial Eb
          		hard_count++;				//the BEMRI has hardened, increment hard_count
           	else
           		hard_count = 0;				//the BEMRI has softened, reset hard_count
           	if( hard_count == max_hard_count)
           	{
           		status = 2;
           		break;
           	}
          	E_last = Eb;
		}
		if(r_current < r_node && entered_node == 0)	//updates values for overall while loop condition, don't want to update while BEMRI is in the node.
		{
			entered_node = 1;	//then the node has been entered
    		theta_N = theta_b;		//save the phase when BEMRI entered the node
    		ecc_N = eb;
		}
		if(r_current > r_node && entered_node == 1)	//if BEMRI has just left the node
		{
			passed_node = 1;						//set passed_node
			entered_node = 0;						//reset entered_node
		}
	}
			

	/* if m1 has gone into negative y territory, start the clock */
	if(Y1<=0 && !m1NegY){
		m1NegY = 1;
		timer0 = t_RK4;
		timer = t_RK4;
		printf("\nTimer started!\n");
	} 

	
	

	/* if timer has gone off, then break */
	if(timer<=0) 
        {
//	printf("TIMER! eb=%.3e | e1h = %.3e | e2h = %.2e\n",eb,e1h,e2h);
         	if(fabs(eb_prev-eb)/eb <1e-9 && k1/r_1 > Eb && k2/r_2 > Eb)
                {
                	status = 0;
                        printf("No change condition met\n");
                        break;
                }
		else if(fabs(e1h_prev-e1h)/e1h < 1e-9 || fabs(e2h_prev-e2h)/e2h < 1e-9)
		{
                	status = 0;
                        printf("No change condition met\n");
                        break;
                }	
        }		
			
	if(f.Htflag == 1 && r_cm > ri )	// if 
        {
           	status = 0;
           	break;
        }			
        if(r_bemri < 2*R_NS && G < 1)
        {
           	status = -1;
           	break;
        }
		
		//Energy conservation check
        if( fabs(Esys0 - b->energy) / fabs(Esys0) > 1e-6 )
        {
        	status = -3;
			printf("\nenergy violation\nEsys0 = %.10e\nEsys(t) = %.10e\n",Esys0,b->energy);
        	break;
        }
        
		// angular momentum conservation check 
		if( fabs(Lsys0 - b->L) / fabs(Lsys0) > 1e-6 )
        {
        	status = -4;
			printf("\nAng. mom violation\nLsys0 = %.10e\nLsys(t) = %.10e\n",Lsys0,b->L);
        	break;
        }


/* print progress */	
		if(f.cflag == 1 && f.soflag==0)
		{
			printf("Current Theta: %.3f\tOrbits Completed: %d\r",theta_h,times_around);
			fflush(stdout);
		}
		else if(f.soflag==0 && eh0 < 1)
		{
			if(print_percentage(a,ii,t_RK4,tmax,ipsval))
			{
				xed = 1;
				status = -2;
				fflush(stdout);
				break;
			}
		}


	/* Full Data Output */
		if(a%astep == 0 && f.cflag != 1 && f.Nflag != 1)
		{	
			output_time(t_RK4+dt1,pos_fp);
			output_positions(rr,b,pos_fp);
			fprintf(pos_fp,"%.10e,%.10e,%.10e,%.10e,%.10e,",Eb,Eh,b->energy,b->L,eb);
			fprintf(pos_fp,"\n");
			output_QP(Qtt, Q_fp);
		}


	} // END OF CURRENT SINGLE RUN

/* initial clean up */					
	stop = time(NULL);
	sim_time = difftime(stop,start)/60.0;

	a1h = E_to_a(E1h,m1,m3);
	a2h = E_to_a(E2h,m2,m3);

/* print percentage if eh0 >= 1 */		
if(f.soflag==0 && eh0 >= 1)
{
	print_percentage(0,ii,angle_counter,ipsval,N);
	printf("\nde = %.10e\n",eb - eb0);
}


/* assign proper status if came out with a zero */
if(status==0)
{
	if(Eb>0)	// binary disrupted
		status=0;
	else if(Eh < 0)	// binary survived, bound to SMBH
		status=1;
	else
		status=2;
}


/* final state output */
	if(Eb < 0)
	{ e1h = e2h = -1; }
	fprintf(results_fp,"%d,%ld,",status,iseed);
	fprintf(results_fp,"%.4e,%.4e",t_RK4,timer0);
	fprintf(results_fp,"%.10e,%.10e,%.10e,%.10e,%.10e,",eb0,eb,eh,e1h,e2h);
	fprintf(results_fp,"%.10e,%.10e,%.10e,%.10e,%.10e,",ab0,ab,ah,a1h,a2h);
	fprintf(results_fp,"%.10e,%.10e,%.10e,%.10e,%.10e,",Eh,E1h,E2h,E0,Esys0);
	fprintf(results_fp,"%.10e,%.10e,%.10e,%.10e,",lh,l1h,l2h,Lsys0);
   	fprintf(results_fp,"%.10e,%.10e,%.10e,%.10e",iPars.incA,iPars.lanA,iPars.th0,iPars.gamma);
    fprintf(results_fp,"\n");
		
/* status codes: */
	//-4: ang. mom. conservation violation
    //-3: energy conservation violation
	//-2: run manually canceled
	//-1: simulation halted due to BEMRI proximity
    //0 : good simulation, BEMRI disrupted
    //1 : good simulation, BEMRI survived but bound to SMBH
    //2 : good simulation, BEMRI survived and remained unbound
    //3 : lifetime for given parameters was too short

}	// END OF CURRENT ANGLE_COUNTER LOOP 

/* 100% print */
	if(xed == 0 && f.cflag == 0 && f.soflag==0)
	{
		fflush(stdout);
		printf("Running %d of %d... (100%%)\n\r",ii+1,N);
		printf("\n\r");
	}

}	// END OF OVERALL N LOOP
	
/* file clean up */
	fclose(results_fp);
	if(f.cflag != 1 && f.Nflag != 1)
	{
		fclose(pos_fp);
		fclose(Q_fp);
	}
	printf("\nfcount = %d\n",fcount);
	free(b);
 	free(C);	
	return 0;
}





/****************************************************************************/
/****************************************************************************/
/****************************************************************************/


void view_info(BEMRI * b)
{
	printf("\nCurrent Data for the BEMRI\n");
	printf("-----------------------------");

	printf("\nm1 \t= \t%.2e\nm2 \t= \t%.2e\nm3 \t= \t%.2e\ne_b \t= \t%.2e\ne_h \t= \t%.2e\nP_b \t= \t%.2e\nP_h \t= \t%.2e\na_b \t= \t%.2e\na_h \t= \t%.2e\nia \t= \t%.2e\naop \t= \t%.2e\nlan \t= \t%.2e\n",\
			m1,m2,m3,eb,eh,Pb,Ph,ab,ah,ia*180/PI,aop*180/PI,lan*180/PI);

	anykey();
	return;

}

/*
double min(double x, double y)
{
return x>y ? y : x;
}
*/
double max(double x, double y)
{
return x<y ? y : x;
}


void update_binaries(BEMRI * b, double (*r)[3], double (*v)[3])
{
    for(int i=0; i<3; i++)
    {
        R1[i] = r[0][i];
        R2[i] = r[1][i];
        R3[i] = r[2][i];
        V1[i] = v[0][i];
        V2[i] = v[1][i];
        V3[i] = v[2][i];
    }

}

void update_vectors(BEMRI * b, double (*r)[3], double (*v)[3])
{
    for(int i=0; i<3; i++)
    {
    	r[0][i] = R1[i];
    	r[1][i] = R2[i];
    	r[2][i] = R3[i];
    	v[0][i] = V1[i];
    	v[1][i] = V2[i];
    	v[2][i] = V3[i];
    }

}

void BEMRI_CM_update(BEMRI * b)
{

        X1 += X4;
        Y1 += Y4;
        Z1 += Z4;
        X2 += X4;
        Y2 += Y4;
        Z2 += Z4;
        VX1 += VX4;
        VY1 += VY4;
        VZ1 += VZ4;
        VX2 += VX4;
        VY2 += VY4;
        VZ2 += VZ4;

}

void init_params(BEMRI * b, flags f)
{

 m1 = 10*MSUN;
 m2 = 10*MSUN;
 m3 = 1e6*MSUN;
 eb = 0.0;
 if(!f.ehflag)	eh = 1.0;
 //ia =  0 * PI / 180.0;
 //lan = 0 * PI / 180.0;
 aop = 0;

 if(f.ehflag && eh < 1)
 {
	 Pb = 1e3;
	 Ph = 1e7;
	 ab = pow((G*(m1+m2)*(Pb*Pb)/(4*PI*PI)),0.3333333333);
	 ah = pow((G*(m1+m2+m3)*(Ph*Ph)/(4*PI*PI)),0.3333333333);
	 rpb = ab*(1-eb);
	 rph = ah*(1-eh);
 }

 if(!f.ehflag || eh >= 1.0)
 {
	 ab = 10*RSUN;
	 rpb = ab*(1-eb);
	 Ph = 2*PI*sqrt(ah*ah*ah/(G*(m1+m2+m3)));
	 Pb = 2*PI*sqrt(ab*ab*ab/(G*(m1+m2)));
 }



 b->binary_b.M = m1+m2;
 mu_b = m1*m2/(m1+m2);
 mu_h = (m3*(m2+m1))/(m3+m2+m1);
 m4 = m1+m2;
 b->binary_h.M = m3+m4;

 for(int ii = 0; ii < 3; ii++)
 {
	 R1[ii] = 0;
	 R2[ii] = 0;
	 R3[ii] = 0;
	 R4[ii] = 0;
	 V1[ii] = 0;
	 V2[ii] = 0;
	 V3[ii] = 0;
	 V4[ii] = 0;
 }

}

void recalc_params(BEMRI * b)
{

 ab = pow((G*(m1+m2)*(Pb*Pb)/(4*PI*PI)),0.3333333333);
 ah = pow((G*(m1+m2+m3)*(Ph*Ph)/(4*PI*PI)),0.3333333333);
 rpb = ab*(1-eb);
 if(eh < 1.0) rph = ah*(1-eh);
 b->binary_b.M = m1+m2;
 b->binary_b.mu = m1*m2/(m1+m2);
 mu_h = (m3*(m2+m1))/(m3+m2+m1);
 b->binary_h.mass2 = m1+m2;
 b->binary_h.M = m1+m2+m3;
}

void load_vectors(double (*r)[3], double (*v)[3], double * m, BEMRI * b)
{

    r[0][0] = X1;
    r[0][1] = Y1;
    r[0][2] = Z1;
    r[1][0] = X2;
    r[1][1] = Y2;
    r[1][2] = Z2;
    r[2][0] = X3;
    r[2][1] = Y3;
    r[2][2] = Z3;

    v[0][0] = VX1;
    v[0][1] = VY1;
    v[0][2] = VZ1;
    v[1][0] = VX2;
    v[1][1] = VY2;
    v[1][2] = VZ2;
    v[2][0] = VX3;
    v[2][1] = VY3;
    v[2][2] = VZ3;

    m[0] = m1;
    m[1] = m2;
    m[2] = m3;

}


void load_vectors2(double **r, double **v, double * m, BEMRI * b)
{
	
    r[0][0] = X1;
    r[0][1] = Y1;
    r[0][2] = Z1;
    r[1][0] = X2;
    r[1][1] = Y2;
    r[1][2] = Z2;
    r[2][0] = X3;
    r[2][1] = Y3;
    r[2][2] = Z3;
	
    v[0][0] = VX1;
    v[0][1] = VY1;
    v[0][2] = VZ1;
    v[1][0] = VX2;
    v[1][1] = VY2;
    v[1][2] = VZ2;
    v[2][0] = VX3;
    v[2][1] = VY3;
    v[2][2] = VZ3;
	
    m[0] = m1;
    m[1] = m2;
    m[2] = m3;
	
}

void load_special_ICs(BEMRI * b, double * t0)
{
    *t0 = 5.05901906e+03;

    X1 = -1.641560385613485145568847656250e+10;
    Y1 = 3.444808922980666351318359375000e+10;
    Z1 = 0;
    VX1 = -6.512512838399886339902877807617e+07;
    VY1 = 3.000371201613594219088554382324e+07;
    VZ1 = 0;

    X2 = -1.641223525744371414184570312500e+10;
    Y2 = 3.445997206149987792968750000000e+10;
    Z2 = 0;
    VX2 = -6.634873829353617876768112182617e+07;
    VY2 = 2.387156840609682723879814147949e+07;
    VZ2 = 0;

    X3 = 3.282783911358051409479230642319e+04;
    Y3 = -6.890806129135329683776944875717e+04;
    Z3 = 0;
    VX3 = 1.314738666775188846713717794046e+02;
    VY3 = -5.387528042223820534672995563596e+01;
    VZ3 = 0;
}

int print_percentage(int a, int i, double t, double tmax, int N)
{
        if((a%100)==0)
        printf("Running %d of %d... (%d%%)\r",i+1,N,(int)((100*t)/tmax));
/*
        if(kbhit())
            if(getch() == 'x')
            {
                printf("Run %d Terminated                                         \n\r",i+1);
                fflush(stdout);
                return 1;
            }
*/
        fflush(stdout);
        return 0;
}


/*	Determines which agent is running the current simulation 
	agent codes:
	0: ALTAIR
	1: ALBIREO
	2: EURISKO
	3: James-iMac
	4: greta	*/
/*
int which_agent(void)
{
	
	FILE *agent_pipe = popen("echo $HOSTNAME","r");		//run echo $HOSTNAME and assign file pointer agent_pipe to point to the output stream
	char buffer[80]; 
	fgets(buffer, sizeof(buffer), agent_pipe);			//copy output of echo $HOSTNAE to buffer
	pclose(agent_pipe);									//close the pipe

	/* Hostname checking */
/*
	if(strncmp(buffer,"ALTAIR",6) == 0)
		return 0;
	else if(strncmp(buffer,"ALBIREO",7) == 0)
		return 1;
	else if(strncmp(buffer,"EURISKO",7) == 0)
		return 2;	
	else if(strncmp(buffer,"James",5) == 0)
		return 3;
	else if(strncmp(buffer,"greta",5) == 0)
		return 4;
	else
		return -1;
}
*/
void output_time(double t, FILE * fp)
{
/* output time */
	fprintf(fp,"%.8e\t",t);	
	
	return;
}

void output_positions(double (*r)[3], BEMRI * b, FILE * fp)
{

/* output positions */
	fprintf(fp,"%.10e,%.10e,%.10e,%.10e,%.10e,%.10e,%.10e,%.10e,%.10e,%.10e,%.10e,%.10e,",r[0][0],r[0][1],r[0][2],r[1][0],r[1][1],r[1][2],r[2][0],r[2][1],r[2][2],b->binary_h.x2[0],b->binary_h.x2[1],b->binary_h.x2[2]);

	return;
}
	
void output_QP(double (*Qtt)[3], FILE * fp)
{
	
/* output quadrupole data */
	for(int ii=0;ii<3;ii++)
	{
		for(int jj=0;jj<3;jj++)
		{
			fprintf(fp,"%.10e,",Qtt[ii][jj]);
		} }
	fprintf(fp,"\n");
	return;
}

void pack_y_vector(double (*r)[3], double (*v)[3], double *m, double *y)
{
	
	int i,j;
	
	y[0] = 3;
	
	for(i=0;i<3;i++)
	{
		for(j=0;j<3;j++)
		{
			y[1+3*i+j] = r[i][j];
			y[1+3*i+j+9] = v[i][j];
		}
		y[19+i] = m[i];
	}
	
	
}


void unpack_y_vector(double (*r)[3], double (*v)[3], double *m, double *y)
{
	
	int i,j;
	
	for(i=0;i<3;i++)
	{
		for(j=0;j<3;j++)
		{
			r[i][j] = y[1+3*i+j];
			v[i][j] = y[1+3*i+j+9];
		}
	}
	
}



void Nbody_derivs1(double x, double y[], double dydx[])
//for use with standard integrators
{
	//I want the input vector to look like
	//[N x1 x2 ... v1 v2 ... m1 m2 ...]
	
	
	static int first=1, N;
	int i,j,k;
	double a;
	
	if(first==1)
		N = (int)y[0];	
	
	//build mass vector
	static double m[3];
	double d[N][N];
	if(first==1)
	{
		Vcopy(N,&y[6*N+1],m);
		first=0;
	}	
	
	//zero out output vector
	for(i=0;i<7*N+1;i++)
		dydx[i]=0;	
	
	//for N-body, position derivatives are the velocities
	for(i=1;i<=3*N;i++)
		dydx[i] = y[i+3*N];
	
	
	
	//build distance matrix
	for(i=0; i<N; i++)
	{
		for(j=0; j<i; j++)
			d[i][j] = sqrt( pow((y[3*i+1]-y[3*j+1]),2) + pow((y[3*i+2]-y[3*j+2]),2) + pow((y[3*i+3]-y[3*j+3]),2) );
	}
	
	//for simple N-body, this will calculate pairwise forces
	for(i=0; i<N; i++)
	{
		for(j=0;j<i;j++)
		{
			for(k=0;k<3;k++)
			{
				a = -G*m[j]*(y[3*i+k+1] - y[3*j+k+1]) / pow(d[i][j],3);
				dydx[3*i+k+3*N+1] += a;
				dydx[3*j+k+3*N+1] -= m[i]/m[j]*a;
			}
		}
	}
	
	
	
}

void Nbody_derivs2(double x, double y[], double dydx[])
//for use with leap frog integration
{
	
	//I want the input vector to look like
	//[N x1 x2 ... v1 v2 ... m1 m2 ...]
	
	
	static int first=1, N;
	int i,j,k,f_or_g = dydx[0];
	
	if(first==1)
		N = (int)y[0];
	
	//zero out output vector
	for(i=0;i<7*N+1;i++)
		dydx[i]=0;	
	
	
	/* Calculate f(p) */
	if(f_or_g==0)	
	{
		//for N-body, this will just return velocities
		for(i=1;i<=3*N;i++)
			dydx[i] = y[i+3*N];
	}
	
	
	/* Calculate g(y) */
	else
	{
		//build mass vector
		static double m[3];
		double d[N][N], a;
		if(first==1)
		{
			Vcopy(N,&y[6*N+1],m);
			first=0;
		}
		
		//build distance matrix
		for(i=0; i<N; i++)
		{
			for(j=0; j<i; j++)
				d[i][j] = sqrt( pow((y[3*i+1]-y[3*j+1]),2) + pow((y[3*i+2]-y[3*j+2]),2) + pow((y[3*i+3]-y[3*j+3]),2) );
		}
		
		//for simple N-body, this will calculate pairwise forces
		for(i=0; i<N; i++)
		{
			for(j=0;j<i;j++)
			{
				for(k=0;k<3;k++)
				{
					if(i!=j)
					{
						a = -G*m[j]*(y[3*i+k+1] - y[3*j+k+1]) / pow(d[i][j],3);
						dydx[3*i+k+3*N+1] += a;
						dydx[3*j+k+3*N+1] -= m[i]/m[j]*a;
					}
				}
			}
		}
	}
	

	
}


void setup_chain(CHAIN *C, double *m,int N)
{
	
	C->N = N;
	for(int i=0;i<N;i++)
	{
		C->mass[i] = m[i];
		for(int j=0;j<3;j++)
			C->R[i][j] = C->W[i][j] = 0;
	}
	C->init = 0;
	init_chain(C);
}

void pack_y_vector_C(CHAIN *C, double y[])
{

	int i,j,N=C->N-1;
	
	for(i=0;i<N;i++)
	{
		for(j=0;j<3;j++)
		{
			y[3*i+j]   = C->R[i][j];
			y[3*i+j+3*N] = C->W[i][j];
		}
	}
}

void unpack_y_vector_C(CHAIN *C, double y[])
{
	int i,j,N=C->N-1;
	
	for(i=0;i<N;i++)
	{
		for(j=0;j<3;j++)
		{
			C->R[i][j] = y[3*i+j];
			C->W[i][j] = y[3*i+j+3*N];
		}
	}
}	

/*--- Calculate vectorized leap-frog derivatives ---*/
/*--- This version does not use transformed time ---*/
void leap_derivs2(double xs, double y[], double dydx[])
{
	
	//vector y looks like [N x1 x2 ... xN v1 v2 ... vN]
	
	extern CHAIN *C;		//pick up external chain structure
	int N = C->N-1;
	
	int i,j,f_or_g = dydx[0];
	
	
	/* zero out output vector */
	for(i=0;i<6*N;i++)
		dydx[i]=0;	
	
	/* Calculate dQ part */
	if(f_or_g==0)	
	{	
		unpack_y_vector_C(C,y);
		double dT[N][3];
		calc_dT(C,dT);
		for(i=0;i<N;i++)
		{
			for(j=0;j<3;j++)
				dydx[3*i+j] = dT[i][j];
		}
	}
	
	/* Calculate dP part */
	else
	{
		unpack_y_vector_C(C,y);
		double du1[N][3], du2[N][3],U;
		calc_dU1(C,du1);
		calc_dU2(C,du2);		
		U = calc_U(C);
		for(i=0;i<N;i++)
		{
			for(j=0;j<3;j++)
				dydx[3*i+j+3*N] = (du1[i][j] + du2[i][j]);
		}
	}
	
}


/****************************************************************************/
/****************************************************************************/
/****************************************************************************/
