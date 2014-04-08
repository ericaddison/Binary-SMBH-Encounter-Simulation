/********************************************************************
*********************************************************************

    Newtonian N-body simulator

    Author: Eric Addison
    Initial Date: 20 Feb 2010
    Affiliation: USU Astrophysics

    File: newtonian_n_body.c
    Purpose: Simulates Newtonian gravitational interactions between N bodies

*********************************************************************
*********************************************************************/


//#include <stdio.h>
//#include <stdlib.h>
//#include <math.h>
//#include "kb.h"
#include "BEMRI.h"

void derivatives(double (*r)[3], double (*v)[3], double * m, double (*dx_out)[3], double (*dv_out)[3], int N );
void RK41(double (*r)[3], double (*v)[3], double * m, int N, double (*dr)[3], double (*dv)[3], double dt);
double N_body_main(double (*r)[3], double (*v)[3], double * m, double tmax, int N, double *dt, double dt_min, int * minned);
void midpoint(double (*r)[3], double (*v)[3], double * m, int N, double dt);
double RK4_adaptive1(double (*r)[3], double (*v)[3], double * m, int N, double *dt, double tol, double dt_min);


/*************************************************************************/
/*************************************************************************/
/**																		**/
/** Funtion: N_body_main                                                **/
/** Purpose: main function for N-body sim, runs the sim until time tmax **/
/**	Inputs: vector of positions                                         **/
/**         vector of velocities                                        **/
/**         vector of masses                                            **/
/**         number of objects (N)                                       **/
/**         simulation time tmax                                        **/
/**	Output: double which will be the amount of time actually passed		**/
/** 			during the integration                                  **/
/**																		**/
/*************************************************************************/
/*************************************************************************/

double N_body_main(double (*r)[3], double (*v)[3], double * m, double tmax, int N, double *dt, double dt_min, int * minned)
{

    double tol = 1e-9, r_temp, f[2], r_min, t_passed = 0, temp, last_dt = 0;
    double dr1[N][3],dv1[N][3];
    r_min = sqrt( pow((r[0][0] - r[1][0]),2) + pow((r[0][1] - r[1][1]),2) + pow((r[0][2] - r[1][2]),2) );

    for(double t = 0; fabs(t)<=fabs(tmax); t+=last_dt)
    {
       last_dt = RK4_adaptive1(r,v,m,N,dt,tol,dt_min);
       t_passed += last_dt;

       if( fabs(t_passed)+fabs(*dt) >= fabs(tmax))
       {
    	   temp = *dt;
    	   *dt = tmax - t_passed;
    	   t_passed += *dt;
    	   RK41(r,v,m,N,dr1,dv1,*dt);
    	   *dt = temp;
    	   break;
       }

       r_temp = sqrt( pow((r[0][0] - r[1][0]),2) + pow((r[0][1] - r[1][1]),2) + pow((r[0][2] - r[1][2]),2) );
       if(r_temp < r_min)
            r_min = r_temp;
       if(r_min < 2*R_NS)
       {
       	printf("\n\nBEMRI components too close. Halting simulation\n\r");
       	fflush(stdout);
       	*minned = 1;
       	break;
       }

    }
	
    return t_passed;
}


/*************************************************************************/
/*************************************************************************/
/**																		**/
/** Funtion: RK4_adaptive                                               **/
/** Purpose: Runge Kutta ODE solver with adaptive step size control     **/
/**	Inputs: vector of positions                                         **/
/**         vector of velocities                                        **/
/**         vector of masses                                            **/
/**         number of objects (N)                                       **/
/**         time step dt                                                **/
/**         Tolerance value (fractional accuracy)                       **/
/**	Output: Double, returns the value of dt actually used in the step	**/
/**			Also fills in values of output vectors                      **/
/**																		**/
/*************************************************************************/
/*************************************************************************/

double RK4_adaptive1(double (*r)[3], double (*v)[3], double * m, int N, double *dt, double tol, double dt_min)
{

    double r_big[N][3], r_small[N][3], v_big[N][3], v_small[N][3];      //will hold results of separate steps
    double dr1[N][3], dr2[N][3], dr3[N][3], dv1[N][3], dv2[N][3], dv3[N][3];
    double delta_c, delta_i, dt_est1, dt_est = *dt, r_temp, v_temp, deriv, dt_used;
    int not_satisfied = 1, a = 0;
    //printf("RKA dt = %.3e\n",*dt);
	extern int fcount;
	fcount += 1;
	
    while(not_satisfied)
    {
        a += 1;

        not_satisfied=0;            //start by assuming tolerance will be met
        dt_est = 100 * (*dt);

        for(int i=0;i<N;i++)        //copy r and v into r_big , v_small, etc
        {
            for(int j=0;j<3;j++)
            {
                r_big[i][j] = r[i][j]; r_small[i][j] = r[i][j];
                v_big[i][j] = v[i][j]; v_small[i][j] = v[i][j];
            }
        }

        RK41(r_big,v_big,m,N,dr1,dv1,*dt);   //big step

        RK41(r_small,v_small,m,N,dr2,dv2,*dt/2.0);
        RK41(r_small,v_small,m,N,dr3,dv3,*dt/2.0);     //two half steps
		
        for(int i=0;i<N;i++)        //calculate errors and find smallest new estimated step size
        {
            for(int j=0;j<3;j++)
            {
                delta_c = fabs(r_big[i][j] - r_small[i][j]);      //delta_c is the estimated truncation error
                r_temp = (r_big[i][j] + r_small[i][j]) / 2.0;       //need a value of r(i)
                deriv = (dr1[i][j] + (dr2[i][j] + dr3[i][j])/2.0) / 2.0;      //also need value of derivatives
                delta_i = tol*fabs( (r_temp) + (*dt) * (deriv));
                dt_est1 = *dt * fabs(pow((delta_i/delta_c),0.25));
                if(fabs(dt_est1) < fabs(dt_est))
                    dt_est = dt_est1;
                if(delta_c > delta_i)
                    not_satisfied = 1;

                //now do it for the v's
                delta_c = fabs(v_big[i][j] - v_small[i][j]);
                v_temp = (v_big[i][j] + v_small[i][j]) / 2.0;       //need a value of v(i)
                deriv = (dv1[i][j] + (dv2[i][j] + dv3[i][j])/2.0) / 2.0;      //also need value of derivatives
                delta_i = tol*fabs( (v_temp) + (*dt) * (deriv));
                dt_est1 = *dt * fabs(pow((delta_i/delta_c),0.25));
                if(fabs(dt_est1) < fabs(dt_est))
                    dt_est = dt_est1;
                if(delta_c > delta_i)
                    not_satisfied = 1;
            }
        }
        	dt_used = *dt;

		
        if(S1*fabs(dt_est) > S2*fabs(*dt))        //determine new dt
            *dt *= S2;
        else if(S1*fabs(dt_est) > fabs(*dt)/S2)
            *dt /= S2;
        else
            *dt = S1 * dt_est;

        if(fabs(*dt) < dt_min)
        {
            *dt = sign(*dt)*dt_min;
            not_satisfied = 0;
        }

    }

    for(int i=0;i<N;i++)        //copy in new values
        {
            for(int j=0;j<3;j++)
            {
                r[i][j] = r_small[i][j];
                v[i][j] = v_small[i][j];
            }
        }

    return dt_used;

}



/*************************************************************************/
/*************************************************************************/
/**																		**/
/** Funtion: RK4                                                        **/
/** Purpose: Runge Kutta ODE stepper                                    **/
/**	Inputs: vector of positions                                         **/
/**         vector of velocities                                        **/
/**         vector of masses                                            **/
/**         number of objects (N)                                       **/
/**         time step dt                                                **/
/**	Output: void. Fills in values of output vectors                     **/
/**																		**/
/*************************************************************************/
/*************************************************************************/

void RK41(double (*r)[3], double (*v)[3], double * m, int N, double (*dr)[3], double (*dv)[3], double dt)
{
    double xk1[N][3], xk2[N][3], xk3[N][3], xk4[N][3];
    double vk1[N][3], vk2[N][3], vk3[N][3], vk4[N][3];
    double rtemp[N][3], vtemp[N][3];


    derivatives(r,v,m,xk1,vk1,N);               //RK step 1

    //set up for next step
    for(int i=0;i<N;i++)
    {
        for(int j=0;j<3;j++)
        {
            rtemp[i][j] = r[i][j] + xk1[i][j]*dt/2.0;
            vtemp[i][j] = v[i][j] + vk1[i][j]*dt/2.0;
        }
    }

    derivatives(rtemp,vtemp,m,xk2,vk2,N);       //RK step 2

    for(int i=0;i<N;i++)
    {
        for(int j=0;j<3;j++)
        {
            rtemp[i][j] = r[i][j] + xk2[i][j]*dt/2.0;
            vtemp[i][j] = v[i][j] + vk2[i][j]*dt/2.0;
        }
    }

    derivatives(rtemp,vtemp,m,xk3,vk3,N);       //RK step 3

    for(int i=0;i<N;i++)
    {
        for(int j=0;j<3;j++)
        {
            rtemp[i][j] = r[i][j] + xk3[i][j]*dt;
            vtemp[i][j] = v[i][j] + vk3[i][j]*dt;
        }
    }

    derivatives(rtemp,vtemp,m,xk4,vk4,N);       //RK step 4


    for(int i=0;i<N;i++)                        //update positions and velocities
    {
        for(int j=0;j<3;j++)
        {
            dr[i][j] = (xk1[i][j] + 2.0*xk2[i][j] + 2.0*xk3[i][j] + xk4[i][j])*dt/6.0;
            dv[i][j] = (vk1[i][j] + 2.0*vk2[i][j] + 2.0*vk3[i][j] + vk4[i][j])*dt/6.0;

            r[i][j] += dr[i][j];
            v[i][j] += dv[i][j];
        }
    }


    return;

}


/*************************************************************************/
/*************************************************************************/
/**																		**/
/** Funtion: midpoint                                                   **/
/** Purpose: Midpoint method ODE stepper                                **/
/**	Inputs: vector of positions                                         **/
/**         vector of velocities                                        **/
/**         vector of masses                                            **/
/**         number of objects (N)                                       **/
/**         time step dt                                                **/
/**	Output: void. Fills in values of output vectors                     **/
/**																		**/
/*************************************************************************/
/*************************************************************************/
void midpoint(double (*r)[3], double (*v)[3], double * m, int N, double dt)
{

    double dx[N][3];
    double dv[N][3];
    double rtemp[N][3], vtemp[N][3];

    derivatives(r,v,m,dx,dv,N);

    for(int i=0;i<N;i++)
    {
        for(int j=0;j<3;j++)
        {
            rtemp[i][j] = r[i][j] + dx[i][j]*dt/2.0;
            vtemp[i][j] = v[i][j] + dv[i][j]*dt/2.0;
        }
    }


    derivatives(rtemp,vtemp,m,dx,dv,N);

    for(int i=0;i<N;i++)
    {
        for(int j=0;j<3;j++)
        {
            r[i][j] += dx[i][j]*dt;
            v[i][j] += dv[i][j]*dt;
        }
    }

}


/*************************************************************************/
/*************************************************************************/
/**																		**/
/** Funtion: derivatives                                                **/
/** Purpose: calculate the derivatives for N-body sim                   **/
/**	Inputs: vector of positions                                         **/
/**         vector of velocities                                        **/
/**         vector of masses                                            **/
/**         output vector for positions                                 **/
/**         output vector for velocities                                **/
/**         number of objects (N)                                       **/
/**	Output: void. Fills in values of output vectors                     **/
/**																		**/
/*************************************************************************/
/*************************************************************************/

void derivatives(double (*r)[3], double (*v)[3], double * m, double (*dx_out)[3], double (*dv_out)[3], int N )
{

    double d[N][N], a;

    //compute distances from each object
    //initialize output vectors at the same time

    for(int i=0; i<N; i++)
    {
        for(int j=0; j<i; j++)
        {
            if(i!=j)
                d[i][j] = sqrt( pow((r[i][0]-r[j][0]),2) + pow((r[i][1]-r[j][1]),2) + pow((r[i][2]-r[j][2]),2) );
        }

        for(int j=0;j<3;j++)
        {
            dx_out[i][j] = 0;
            dv_out[i][j] = 0;
        }
    }


    //compute acceleration from each object
    for(int i=0; i<N; i++)
    {
        for(int k=0; k<3; k++)
        {
            for(int j=0; j<i; j++)
            {
                if(i!=j)
				{
					a = -G * m[j] * (r[i][k] - r[j][k]) / pow(d[i][j],3);    //kth component of acceleration on body i from body j
                    dv_out[i][k] += a;
					dv_out[j][k] -= m[i]/m[j]*a;
				}
            }
            dx_out[i][k] = v[i][k];     //derivative of position is velocity
        }
    }

    return;

}
