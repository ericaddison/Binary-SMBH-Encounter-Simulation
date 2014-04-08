
            /*
//TEST TEST TEST

			//lets reverse the orbit of the BEMRI, see if it has the same effect as ia = PI
			for(int ii = 0; ii < 3; ii++)
			{
				//V1[ii] *= -1;
				//V2[ii] *= -1;
				//V3[ii] *= -1;
				//V4[ii] *= -1;
			}

			Vdisplay(R1,3);
			Vdisplay(R2,3);
			Vdisplay(V1,3);
			Vdisplay(V2,3);
			Vdisplay(V3,3);


//TEST TEST TEST
*/

    //    	printf("\nnode time Pb ratio is: %.3e\n",t_node/Pb);

    //    	printf("\nab = %.3e m = %.3e au\n",ab,m_to_au(ab));
    //    	printf("\nah = %.3e m = %.3e au\n",ah,m_to_au(ah));
	//		anykey();





  // TEST TEST TEST TEST
    	//torque(b,tork);
    	//angular_momentum_vector(R1,R2,V1,V2,m1,m2,b->binary_b.l);
   // TEST TEST TEST TEST





















 //******************************************************************************************************************
 //******************************************************************************************************************
 //******************************************************************************************************************
 //******************************************************************************************************************
 //******************************************************************************************************************


//Code for running a single simulation


 //******************************************************************************************************************
 //												SINGLE SIMULATION
 //******************************************************************************************************************
    if(single == 1)
    {

            xed = 0;

            recalc_params(b);

            //initial conditions
            lb = sqrt(ab*(1-eb*eb)*G*(m1*m1*m2*m2)/(m1+m2));		//initial angular momentum of BEMRI
            advance_orbit(&(b->binary_b),lb,0);						//begin BEMRI at periapsis

            lh = sqrt(ah*(1-eh*eh)*G*(m3*m3*m4*m4)/(m3+m4));		//initial angular momentum of SMBH orbit
            advance_orbit(&(b->binary_h),lh,PI);					//begin hole orbit at apoapsis


            //update little binary positions and velocities for orbit around hole
            BEMRI_CM_update(b);
            t0 = 0;
            load_vectors(r,v,m,b);				//initialize vectors v and r

            dt = min(Pb,Ph)/pow(2.0,7.0);		//set initial time step dt
            dt1 = dt;                                        				//dt changes based on RK4 adaptive, dt1 stays constant and dictates time steps
            tmax = b->binary_h.P_orb*10;                  				//run for three full length orbits of SMBH
            printf("\nah = %.3e\nPh = %.3e\n",ah,Ph);
            anykey();


peters_lifetime(&(b->binary_b));

        	t_node = 2*pow(ah*ah*ah/G/(m3+m4),0.5)*( 2*atan2(sqrt(1-eh),sqrt(1+eh)) - eh*sqrt(1 - eh*eh));
        	tau = t_node/pow(2,10);

        	dt = tau;
        	dt1 = dt;

            a = -1;															//counter for staggered output

            start = time(NULL);
            alpha = 0.001;				//alpha is the parameter which ultimately dictates when to begin full integration

        t_RK4=t0;						//test to make sure time is consistent
        int amax = (int)(abs((tmax-t0)/dt1)), astep = (int)(amax/100000.0);

        printf("\ndt1 = %.3e\ntmax = %.3e\ntau = %.3e\nt0 = %.3e\n",dt1,tmax,tau,t0);
        printf("\nPh = %.3e sec",Ph);
 //       printf("\namax = %d, astep = %d",amax,astep);
		anykey();


//******************************************************************************************************************
//												BEGIN MAIN LOOP
//******************************************************************************************************************


            for(t=t0;sign(dt1)*t<tmax;t += dt1)          //stepping through orbit at intervals of dt1
            {
                a+=1;

              // printf("t = %.3f\n",t);			//display current time

                    dt_min = 1e-8;				//set value for minimum allowable value of dt

                    t_RK4 += N_body_main(r,v,m,dt1,3,&dt,dt_min,&minned);       //evolves orbit from time t to t+dt1 with initial step size dt

                    update_binaries(b,r,v);										//copies values from v and r into the structures

                    CM_values(b);												//calculate values for the center of mass

                    calc_angles(&(b->binary_b));								//calculate current orbital angles

                    r_bemri = sqrt( pow(X1-X2,2) + pow(Y1-Y2,2) + pow(Z1-Z2,2));		//BEMRI separation

    				r_current = sqrt( pow(X4-X3,2) + pow(Y4-Y3,2) + pow(Z4-Z3,2) );		//current distance from SMBH to BEMRI CM


                if(dt > dt1)            //keeps integration step size at maximum of dt1
                    dt = dt1;

                calc_energies(b,&E1h,&E2h);
                calc_angular_momenta(b,&l1h,&l2h);
                calc_ecc(b,&e1h,&e2h,E1h,E2h,l1h,l2h);
                calc_total_energy(b);



//******************************************************************************************************
//											SKIPPING CODE
//******************************************************************************************************

			if(eb < 1)		//if the BEMRI has been broken, this is all moot
			{

			//some initial updates to orbit information
                true_anomaly(&(b->binary_b));		//calculate the true anomaly of the BEMRI
                true_anomaly(&(b->binary_h));		//calculate the true anomaly of the SMBH orbit

				ab = r_bemri*(1+eb*cos(theta_b)) / (1-eb*eb);		//update semi-major axis of the BEMRI
				ah = r_current*(1+eh*cos(theta_h)) / (1-eh*eh);		//update semi-major axis of the BH orbit

				//current orbital periods
				Pb = 2*PI*sqrt(pow(ab,3)/(G*(m1+m2)));
				Ph = 2*PI*sqrt(pow(ah,3)/(G*(m1+m2+m3)));

				//threshold value chi
				chi = ab*(1+eb)*sqrt(m3/(alpha*mu_b));

				//if(a%astep == 0)
					//printf("\nr_current/chi = %.3f\n",r_current/chi);

				//if r_current is less than chi, we want to continue integrating fully
				//if r_current is greater than chi, we want to fast forward to the next point where r_current is equal to chi
/*
				if(r_current > chi)
				{

					printf("\nSkipping!\n");
					printf("\nr_current = %.3e\nchi = %.3e\nah = %.3e",r_current,chi,ah);
					t += advance_orbit_set_up(b,chi,r_current,dt1) - dt1;			//performs all advancing steps
					printf("\nt is now %.3e",t);
					update_vectors(b,r,v);
				}
*/
			}


//******************************************************************************************************************
//												DATA OUTPUT
//******************************************************************************************************************

                if(a%astep == 0)
                {
                    printf("%% done = %.3f%%\n",(t+dt1)/tmax*100);
                    //printf("r_h = %.3e\n\n",r_current);
                    //printf("t = %.3f\n",t);
                    //anykey();
                        //output time
                            fprintf(fp,"%.8e\t",t+dt1);
                        //output positions
                            //fprintf(fp,"%.30e\t%.30e\t%.30e\t%.30e\t%.30e\t%.30e\t%.30e\t%.30e\t%.30e\t%.30e\t%.30e\t%.30e\t",r[0][0],r[0][1],r[0][2],r[1][0],r[1][1],r[1][2],r[2][0],r[2][1],r[2][2],b->binary_h.x2[0],b->binary_h.x2[1],b->binary_h.x2[2]);
                        //output orbital energy and eccentricities and angular momentums
                            fprintf(fp,"%.8e\t%.8e\t%.8e\t%.8e\t%.8e\t%.8e\t%.8e\t%.8e\t%.8e\t%.8e\t%.8e\t%.8e\t%.8e\t%.8e\n",Eb,E1h,E2h,Eh,eb,e1h,e2h,eh,lb,l1h,l2h,lh,ab,ah);

                        //output momentums and total energy
                        //    fprintf(fp,"%.8e\t%.8e\t%.8e\t%.8e\t%.8e\t%.8e\t%.8e\t%.8e\n",m1*v[0][0]+m2*v[1][0]+m3*v[2][0],m1*v[0][1]+m2*v[1][1]+m3*v[2][1],m1*v[0][2]+m2*v[1][2]+m3*v[2][2],TE,r_bemri,m1*VX1,m2*VX2,m3*VX3);

                        //output velocities
                            //fprintf(fp,"%.30e\t%.30e\t%.30e\t%.30e\t%.30e\t%.30e\t%.30e\t%.30e\t%.30e\n",v[0][0],v[0][1],v[0][2],v[1][0],v[1][1],v[1][2],v[2][0],v[2][1],v[2][2]);

                        //fprintf(fp,"%.8e\t%.8e\t%.8e\t%.8e\t%.8e\t%.8e\t%.8e\t%.8e\t%.8e\t%.8e\n",m1*VX1,m1*VY1,m1*VZ1,m2*VX2,m2*VY2,m2*VZ2,m3*VX3,m3*VY3,m3*VZ3,TE);
                        //printf("\n\npx1 = %.3e\npx2 = %.3e\npx3 = %.3e\nTotal px = %.3e\n\n",m1*VX1, m2*VX2, m3*VX3, m1*VX1+m2*VX2+m3*VX3);
                        //anykey();
                      //fprintf(fp,"%.8e,%.8e,%.8e,%.8e,%.8e,%.8e,%.8e,%.8e,%.8e\n",Pb,Ph,eb,eh,ia,lan,aop,theta_b,theta_h);
                }
            }		//close main loop bracket

    }		//close if(single == 1) bracket
//******************************************************************************************************************
//												CLEAN UP
//******************************************************************************************************************

		stop = time(NULL);
    //fprintf(fp,"%.3e,%.8e,%.8e,%.8e,%.4f,%.4f,%.8e,%.8e,%.4f,%.4f,%.4f,%.8e,%.8e,%.8e,%.8e,%.8e,%.8e",t,m1,m2,m3,eb,eh,Pb,Ph,ia,aop,lan,E1,E2,E3,l1,l2,l3);

		printf("t = %.3f\n",t);
        //printf("t_RK4 = %.3f\n",t_RK4);
		printf("\nSim time = %.2f minutes\n",difftime(stop,start)/60);

    fclose(fp);

 return 0;
}


//******************************************************************************************************************
//******************************************************************************************************************
//******************************************************************************************************************
//******************************************************************************************************************














/*
 //+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
 //				TEST CODE FOR OUTLIER EXPLORATION
 //+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
 
 
 //set new parameters within a 0.5 radius of original
 
 test_rad = 0.05*ran2(&seed);
 test_angle = 2*PI*ran2(&seed);
 
 del_ia = test_rad*cos(test_angle);
 del_aop = test_rad*sin(test_angle);
 
 //point 1
 ia = 1.4022086174 + del_ia;
 aop = 1.700875351 + del_aop;
 
 //point 2
 //		ia = 1.5893396678 + del_ia;
 //		aop = 1.217821112 + del_aop;
 
 param1 = ia;
 param2 = aop;
 recalc_params(b);
 
 //+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
 //+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
 */

