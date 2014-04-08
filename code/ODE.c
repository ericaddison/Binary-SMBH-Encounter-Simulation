/********************************************************************
*********************************************************************
 
	File: ODE.c
	Purpose: Main code file for ODE solvers
 
    Author: Eric Addison
    Initial Date: 12 May 2011
    Affiliation: USU Astrophysics


*********************************************************************/

#include "ODE.h"

/********************************************************************
*********************************************************************
	function:	mid
	purpose:	perform one midpoint integration step
	inputs:		double y[]	-- current value of variables
				int nvar	-- number of variables
				double xs	-- current value of independent variable
				double h	-- stepsize
				void *f		-- derivative function
	outputs:	yout		-- output vector
				h			-- output of function
********************************************************************
*********************************************************************/
double mid(double y[], int nvar, double xs, double h, 
		 double yout[],void (*f)(double, double[], double[]))
{

	double F[nvar],ytemp[nvar],dydx[nvar];
	
	/* get initial derivatives */
	f(xs,y,dydx);
	
	/* set up and perform midpoint step */
	vect_mult_scalar(dydx,h/2,ytemp,nvar);
	vect_add(y,ytemp,ytemp,nvar);
	
	f(xs+h/2,ytemp,F);
	
	/* add to previous value */
	vect_mult_scalar(F,h,F,nvar);
	vect_add(y,F,yout,nvar);
	
	return h;
}


/********************************************************************
 *********************************************************************
	function:	mmid
	purpose:	perform one modified midpoint integration step
	inputs:		double y[]	-- current value of variables
							int nvar	-- number of variables
							double xs	-- current value of independent variable
							double htot	-- total stepsize
							int nstep	-- number of substeps to take
							void *f		-- derivative function
				outputs:	yout		-- output vector
							htot		-- output of function
 ********************************************************************
 *********************************************************************/
double mmid(double y[], int nvar, double xs, double htot, int nstep, 
		   double yout[],void (*f)(double, double[], double[]))
{
	
	double h = htot/nstep, xx = xs+h, h2 = 2.0*h; 
	double dydx[nvar], ym[nvar], yn[nvar], temp[nvar];
	/* get initial derivatives */
	f(xs,y,dydx);	

	for(int i=0;i<nvar;i++)
	{
		ym[i] = y[i];
		yn[i] = y[i] + h*dydx[i];
	}
	
	f(xx,yn,yout);
	
	for(int n=2;n<=nstep;n++)
	{
		vect_mult_scalar(yout,h2,temp,nvar);
		vect_add(ym,temp,temp,nvar);
		Vcopy(nvar,yn,ym);
		Vcopy(nvar,temp,yn);
		xx += h;
		f(xx,yn,yout);
	}
	
	for(int i=0;i<nvar;i++)
		yout[i] = 0.5*(ym[i] + yn[i] + h*yout[i]);
	
	return htot;
	
}


/********************************************************************
 *********************************************************************
	function:	RK4
	purpose:	perform one Fourth Order Runge-Kutta integration step
	inputs:		double y[]	-- current value of variables
				double dy[]	-- vector to hold final change to y
				int nvar	-- number of variables
				double xs	-- current value of independent variable
				double h	-- stepsize
				void *f		-- derivative function
	outputs:	yout		-- output vector
				h			-- output of function
 ********************************************************************
 *********************************************************************/
double RK4(double y[], double dy[], int nvar, double xs, double h, 
		 double yout[],void (*f)(double, double[], double[]))
{

	double k1[nvar],k2[nvar],k3[nvar],k4[nvar];
	double x1,x2,x3,temp[nvar];
	
	x1 = xs; x2 = x1+h/2.0; x3 = x1+h;
	
	/* generate pieces k1,k2,k3,k4 */
	f(x1,y,k1);
	vect_mult_scalar(k1,h/2.0,temp,nvar);
	vect_add(y,temp,temp,nvar);	

	f(x2,temp,k2);
	vect_mult_scalar(k2,h/2.0,temp,nvar);
	vect_add(y,temp,temp,nvar);		
	
	f(x2,temp,k3);
	vect_mult_scalar(k3,h,temp,nvar);
	vect_add(y,temp,temp,nvar);	
	
	f(x3,temp,k4);
	
	/* add to previous value */
	for(int i=0;i<nvar;i++)
	{
		dy[i] = (k1[i] + 2.0*k2[i] + 2.0*k3[i] + k4[i])*h/6.0;
		yout[i] = y[i] + dy[i];
	}

	return h;
}



/********************************************************************
 *********************************************************************
	function:	RK4_adaptive
	purpose:	perform one Fourth Order Runge-Kutta integration step with adaptive step sizing
	inputs:		double y[]	-- current value of variables
				int nvar	-- number of variables
				double xs	-- current value of independent variable
				double *hh	-- last used stepsize
				double tol	-- error tolerance
				double hmin -- minimum step size
				void *f		-- derivative function
	outputs:	yout		-- output vector
				h_used		-- output of function, stepsize actually used
 ********************************************************************
 *********************************************************************/
double RK4_adaptive(double y[], int nvar, double xs, double *hh, double yout[], 
					double tol, double hmin,void (*f)(double, double[], double[])) 
{

	int not_satisfied=1, a=0;
	double ybig[nvar], ysmall[nvar], dybig[nvar],dysmall[nvar];
	double dy1[nvar], dy2[nvar], dy3[nvar], deriv, h_used;
	double delta_c, delta_i, ytemp, h_est, h_est1, h = *hh;

	extern int fcount;
	fcount += 1;
	
	
	while(not_satisfied)
	{
		a += 1;

		/* start by assuming tolerance will be met */
		not_satisfied = 0;	
		h_est = 100 * h;
		
		/* take one big step and two small steps across the interval h */
		RK4(y,dy1,nvar,xs,h,ybig,f);
		RK4(y,dy2,nvar,xs,h/2.0,ysmall,f);
		RK4(ysmall,dy3,nvar,xs+h/2.0,h/2.0,ysmall,f);
		
		/* error analysis */
		for(int i=0;i<nvar;i++)
		{
			delta_c = fabs(ybig[i]-ysmall[i]);		//estimated truncation error
			ytemp = (ybig[i] + ysmall[i])/2.0;		//average of two estimates
			deriv = (dy1[i] + (dy2[i] + dy3[i])/2.0) / 2.0;
			delta_i = tol*fabs( (ytemp) + (h)*(deriv) );
			h_est1 = h * fabs(pow((delta_i/delta_c),0.25));
		
			if(fabs(h_est1) < fabs(h_est))
				h_est = h_est1;
			if(delta_c > delta_i)
				not_satisfied = 1;
		}
							  
		h_used = h;		

		if(S1*fabs(h_est) > S2*fabs(h))
			h *= S2;
		else if(S1*fabs(h_est) > fabs(h)/S2)
			h /= S2;
		else
			h = S1 * h_est;
		
		if(fabs(h) < hmin)
		{
			h = SIGN(h)*hmin;
			not_satisfied = 0;
		}
							  
	}					

	Vcopy(nvar,ysmall,yout);
	*hh = h;
	return h_used;
							  
}


/********************************************************************
 *********************************************************************
 function:	RK4A_const_step
 purpose:	perform RK4 Adaptive integration steps up to a 
 specified value of dependent variable
 inputs:	double	y[]		-- current value of variables
			int		nvar	-- number of variables
			double *xs		-- pointer to independent variable
			double htot		-- total stepsize
			double *htry	-- stepsize to start with
			double eps		-- error tolerance
			double hmin		-- minimum step size
			void *f			-- derivative function
 ********************************************************************
 *********************************************************************/	
void RK4A_const_step(double y[], int nvar, double *xs, double htot, double *htry,double eps, 
				   double hmin, void (*f)(double, double[], double[]))
{
	int i=0;
	double xend = htot, hused, dy[nvar], t = 0;
	while(t < htot)
	{
		i++;
		hused = RK4_adaptive(y,nvar,t,htry,y,eps,hmin,f);
		t += hused;
		
		if( t+*htry > xend || i > 1000)	//if next step will land outside the range
		{
			t += RK4(y,dy,nvar,t,xend-t,y,f);
			//*htry = htot;
			break;
		}
	}
	*xs += t;
}



/********************************************************************
 *********************************************************************
 function:	bsstep
 purpose:	perform one Burlish-Stoer integration step 
 inputs:	double y[]		-- current value of variables
			int nvar		-- number of variables
			double *xs		-- pointer to independent variable
			double *htry	-- stepsize to start with
			double eps		-- error tolerance
			double hmin		-- minimum step size
			double yscal[]	-- scaling values for variables
			int method		-- integration method to use (0 => mmid, 1 => mleap)
			void *f			-- derivative function
 outputs:	h_used		-- output of function, stepsize actually used
 Source:	This is straight out of numerical recipes
 ********************************************************************
 *********************************************************************/		  
double bsstep(double y[], int nvar, double *xs, double *htry,double eps, 
			  double hmin, double yscal[], int method, void (*f)(double, double[], double[]))
{
	xBS = malloc(IMAXX*sizeof(double));
	dBS = malloc((nvar+1)*sizeof(double *));
	for(int i=0;i<=nvar;i++)	dBS[i] = malloc(IMAXX*sizeof(double));
	
	int i,iq,k,kk,km; 
	static int first=1,kmax,kopt; 
	static double epsold = -1.0,xnew; 
	double eps1,errmax,fact,h,red,scale,work,wrkmin,xest,power,xx=*xs, aa; 
	double err[KMAXX+1],yerr[nvar+1],ysav[nvar+1],yseq[nvar+1],hdid,hnext; 
	static double alf[KMAXX+1][KMAXX+1], a[IMAXX+1]; 
	static int nseq[IMAXX+1]={0,2,4,6,8,10,12,14,16,18}; 
	int reduct,exitflag=0,targetflag=0;

	
	extern int fcount;
	fcount += 1;
	
/* scale variables for error calculation */
	scale_vars(*xs,*htry,nvar,y,yscal,f);
	
	/* new tolerance is input, so reinitialize */
	if(eps != epsold)
	{
		
		hnext = -1.0e-29;
		
		xnew = -1.0e-29;
		eps1 = SAFE1*eps;		
		
		/* compute work coefficients */
		a[1] = nseq[1] + 1;	
		for(k=1;k<=KMAXX;k++)
			a[k+1] = a[k] + nseq[k+1];
		
		/* build alpha coefficients */
		for(iq = 2;iq<=KMAXX;iq++)
			for(k=1;k<iq;k++){
				power = ( a[k+1] - a[iq+1] ) / ( (a[iq+1] - a[1] + 1.0)*(2*k+1) );
				alf[k][iq] = pow(eps1,power);
			}
		epsold = eps;
		
		/* determine optimal row number for convergence */
		for(kopt = 2;kopt<(KMAXX-1);kopt++)
			if( a[kopt+1] > a[kopt]*alf[kopt-1][kopt] )
				break;
		
		kmax = kopt;
	}
	

	h = *htry;
	
	
	
	/* save starting values */
	Vcopy(nvar,y,ysav);
	
	/* either new stepsize or new integration */
	/* re-establish order window			  */
	if(xx != xnew || h != hnext){
		first = 1;
		kopt = kmax;
	}
	reduct = 0;	
	/* main loop */
	while(1)
	{

		for(k=1;k<=kmax;k++)
		{
			xnew=(xx)+h;
			if (xnew == (xx))
			{
				printf("step size underflow in bsstep\nh = %.3e\n",h);
				return -1;
			}
			if(method == 0)
				mmid(ysav,nvar,xx,h,nseq[k],yseq,f);
			else
				mleap(ysav,nvar,xx,h,nseq[k],yseq,f);
			xest=SQR(h/nseq[k]);		
			rzextr(k,xest,yseq,y,yerr,nvar);
			if (k != 1) 
			{
				errmax=TINY;
				/* scale errors */
				for (i=1;i<=nvar;i++) 
				{
					aa = fabs(yerr[i]/yscal[i-1]);
					errmax=FMAX(errmax,aa);
				}
				errmax /= eps;
				km=k-1;
				err[km]=pow(errmax/SAFE1,1.0/(2*km+1));
			}

			if (k != 1 && (k >= kopt-1 || first))
			{
				/* converged */
				if (errmax < 1.0) 
				{
					exitflag=1;
					break;
				}
				/* check for stepsize reduction */
				if (k == kmax || k == kopt+1) 
				{
					red=SAFE2/err[km];
					break;
				}
				else if (k == kopt && alf[kopt-1][kopt] < err[km]) 
				{				
					red=1.0/err[km];
					break;
				}
				else if (kopt == kmax && alf[km][kmax-1] < err[km]) 
				{				
					red=alf[km][kmax-1]*SAFE2/err[km];
					break;
				}
				else if (alf[km][kopt] < err[km]) 
				{					
					red=alf[km][kopt-1]/err[km];
					break;
				}
			}
		}
		if (exitflag) 
		{
			break;		
		}
		red=FMIN(red,REDMIN);
		red=FMAX(red,REDMAX);
		
		if(h*red < hmin)
		{
			printf("\nBS-Integration: Minimum Stepsize Reached\n");
			h = hmin;
			break;
		}
		
		h *= red;
		reduct=1;
		
	}
	
	
	xx=xnew; 
	hdid=h; 
	first=0; 
	wrkmin=1.0e35;
	
	for (kk=1;kk<=km;kk++) {
		fact=FMAX(err[kk],SCALMX);
		work=fact*a[kk+1]; 
		if (work < wrkmin) 
		{
			scale=fact; 
			wrkmin=work; 
			kopt=kk+1;
		}
	}
	
	hnext=h/scale;
	
	if (kopt >= k && kopt != kmax && !reduct) 
	{
		/* Check for possible order increase, but not if stepsize was just reduced. */
		fact=FMAX(scale/alf[kopt-1][kopt],SCALMX); 
		if (a[kopt+1]*fact <= wrkmin) 
		{
			hnext=h/fact; 
			kopt++;
		}
	}

	*xs += hdid;
	*htry = hnext;
	
	
	free(xBS);
	for(int i=0;i<=nvar;i++)	free(dBS[i]);
	free(dBS);
	
	
	return hdid;
	
}


/********************************************************************
 *********************************************************************
 function:	bs_const_step
 purpose:	perform Burlish-Stoer integration steps up to a 
			specified value of dependent variable
 inputs:		double y[]		-- current value of variables
				int	nvar		-- number of variables
				double *xs		-- pointer to independent variable
				double htot		-- total stepsize
				double *htry	-- stepsize to start with
				double eps		-- error tolerance
				double hmin		-- minimum step size
				double yscal[]	-- scaling values for variables
				void *f			-- derivative function
 ********************************************************************
 *********************************************************************/	
void bs_const_step(double y[], int nvar, double *xs, double htot, double *htry,double eps, 
				   double hmin, double yscal[], int method, void (*f)(double, double[], double[]))
{
	double x=0;
	
	while(x < htot)
	{		
		if( x+*htry > htot )	//if next step will land outside the range
		{
			if(method == 0)
				x += mmid(y,nvar,x,htot-x,8,y,f);
			else
				x += mleap(y,nvar,x,htot-x,8,y,f);
			//*htry = htot;
			break;
		}
		bsstep(y,nvar,&x,htry,eps,hmin,yscal,method,f);
	}

	*xs += htot;

}


/********************************************************************
 *********************************************************************
 function:	pzextr
 purpose:	Used by bsstep to extrapolate a function to x=0 by 
			polynomial extrapolation
 inputs:		int iest		-- current step in the extrapolation process
				double xest		-- dependent variable value
				double yest[]	-- current value for dependent variables
				double yz[]		-- output vector for extrapolated value
				double dy[]		-- output vector for estimated errors
				int nvar		-- number of variables
 ********************************************************************
 *********************************************************************/
void pzextr(int iest, double xest, double yest[], double yz[], double dy[], int nvar)
{
	extern double **dBS,*xBS;
	int k1,j; 
	double q,f2,f1,delta,c[nvar+1];
	
	/* Save current independent variable. */
	xBS[iest]=xest;	 
	for (j=0;j<nvar;j++) 
		dy[j]=yz[j]=yest[j]; 
	
	/* Store first estimate in first column. */
	if (iest == 1) {	
		for (j=1;j<=nvar;j++) 
			dBS[j][1]=yest[j-1];
		
	} 
	else {
		for (j=1;j<=nvar;j++)
			c[j]=yest[j-1]; 
		for (k1=1;k1<iest;k1++) {
			delta=1.0/(xBS[iest-k1]-xest); 
			f1=xest*delta; 
			f2=xBS[iest-k1]*delta;
			for (j=1;j<=nvar;j++) {
				q=dBS[j][k1]; 
				dBS[j][k1]=dy[j]; 
				delta=c[j]-q; 
				dy[j]=f1*delta; 
				c[j]=f2*delta; 
				yz[j-1] += dy[j];
			}
		}
		for (j=1;j<=nvar;j++) 
			dBS[j][iest]=dy[j];
	}
}

/********************************************************************
 *********************************************************************
 function:	rzextr
 purpose:	Used by bsstep to extrapolate a function to x=0 by 
			rational function extrapolation
 inputs:		int iest		-- current step in the extrapolation process
				double xest		-- dependent variable value
				double yest[]	-- current value for dependent variables
				double yz[]		-- output vector for extrapolated value
				double dy[]		-- output vector for estimated errors
				int nvar		-- number of variables
 ********************************************************************
 *********************************************************************/
void rzextr(int iest, double xest, double yest[], double yz[], double dy[], int nvar)

{
	extern double **dBS,*xBS;
	int k,j=1; 
	double yy,v,ddy,c,b1,b,fx[iest+1];
	
	xBS[iest]=xest;
	
	/* Store first estimate in first column. */
	if (iest == 1)
	{
		for (j=1;j<=nvar;j++) 
		{ 
			yz[j-1]=yest[j-1];
			dBS[j][1]=yest[j-1];
		}
	} 
	else 
	{
		dy[j]=yest[j-1];
		for (k=1;k<iest;k++) 
			fx[k+1]=xBS[iest-k]/xest;
		for (j=1;j<=nvar;j++) 
		{ 
			v=dBS[j][1];
			dBS[j][1]=yy=c=yest[j-1]; 
			for (k=2;k<=iest;k++) 
			{
				b1=fx[k]*v; 
				b=b1-c; 
				if (b) 
				{
					b=(c-v)/b; 
					ddy=c*b; 
					c=b1*b;
				} 
				else 
					ddy=v;
				/* Evaluate next diagonal in tableau. */
				if (k != iest) 
					v=dBS[j][k]; 
				dBS[j][k]=ddy; 
				yy += ddy;
			} 
			dy[j]=ddy; yz[j-1]=yy;
		}
	}
}



/********************************************************************
*********************************************************************
 function:	mleap
 purpose:	perform one modified leapfrog integration step
 inputs:		double y[]	-- current value of variables
				double dy[]	-- vector to hold final change to y
				int nvar	-- number of variables
				double xs	-- current value of independent variable
				double h	-- stepsize
				void *f		-- derivative function
 outputs:		yout		-- output vector
				h			-- output of function
 note:		this method can be used when a dynamical system has a
			Hamiltonian that can be separated, i.e.
			y_dot = f(p), p_dot = g(y)
			where y is position and p is momentum
********************************************************************
*********************************************************************/	
double mleap(double y[], int nvar, double xs, double htot, double nstep, 
			double yout[],void (*F)(double, double[], double[]))
{

	
	//the derivatives function F is slightly different than for 
	//previous methods. Now it takes an additional int argument
	//which specifies whether to evaluate f(p) or g(y), which
	//will be hidden in dydx. I choose argument 0 means do f(p), 
	//otherwise do g(y).

	
	double dydx[nvar], ysav[nvar], h = htot/nstep;	
	Vcopy(nvar,y,ysav);

/* initial caluclulation */
	//calculate y_0.5
	dydx[0]=0;				//the first value in dydx tells the derivative function whether to calculate f(p) or g(y)
	F(xs,y,dydx);
	vect_mult_scalar(dydx,h/2.0,dydx,nvar);
	vect_add(y,dydx,y,nvar);	

	//calculate p_1
	dydx[0]=1;
	F(xs,y,dydx);
	vect_mult_scalar(dydx,h,dydx,nvar);
	vect_add(y,dydx,y,nvar);	

	
/* loop over the rest */
	for(int i=1;i<nstep;i++)
	{
	//calculate y_(m-0.5)
		dydx[0]=0;
		F(xs,y,dydx);
		vect_mult_scalar(dydx,h,dydx,nvar);
		vect_add(y,dydx,y,nvar);
	
	//calculate p_m
		dydx[0]=1;
		F(xs,y,dydx);
		vect_mult_scalar(dydx,h,dydx,nvar);
		vect_add(y,dydx,y,nvar);
	}

	
/* final y calculation */
	//calculate y_n
	dydx[0]=0;
	F(xs,y,dydx);
	vect_mult_scalar(dydx,h/2.0,dydx,nvar);
	vect_add(y,dydx,y,nvar);	
	
	Vcopy(nvar,y,yout);
	Vcopy(nvar,ysav,y);
	
	return htot;
}



/********************************************************************
 *********************************************************************
 function:	scale_vars
 purpose:	computes scaling factors for variable error calculation
 inputs:	double xs	-- current value of independent variable
 double h	-- stepsize
 int nvar	-- number of variables
 double y[]	-- current value of variables
 double yscal[] -- array to hold scale values
 void *f		-- derivative function
 outputs:	fills up yscal
 ********************************************************************
 *********************************************************************/	
void scale_vars(double x, double h, int nvar, double y[], double yscal[], 
				void (*F)(double, double[], double[]))

{
	
	double dydx[nvar], a;
	
	dydx[0] = 2;
	F(x,y,dydx);
	for(int i=0;i<nvar;i++)
	{
		a = fabs(y[i])+fabs(dydx[i]*h);
		yscal[i]=FMAX(a,1);
	}

	
	
}







