/*******************************************************************
********************************************************************

    Kepler Equation and Solver

    Author: Eric Addison
    Initial Date: 14 July 2010
    Affiliation: USU Astrophysics

    File: kepler_eq.c
    Purpose: Includes function to evaluate and solve Kepler's equation
			 via bisection method.

********************************************************************
********************************************************************/

#include "BEMRI.h"

/*
#include <stdio.h>
#include <math.h>
double kepler(double E, double e);
double solve_kepler(double M, double e);
double f_kepler(double E, double M, double e);
double bisect_kepler(double a, double b, double M, double e, double f(double E, double M, double e), double tol);
#define PI 3.1415926
*/

/*************************************************************************/
/*************************************************************************/
/**																		**/
/** Funtion: kepler	                                                    **/
/** Purpose: evaluate Kepler's Equation									**/
/**	Inputs: eccentric anomaly E	                                        **/
/**			eccentricity												**/
/**         				                                            **/
/**	Output: Mean anomaly as a double                                    **/
/**																		**/
/*************************************************************************/
/*************************************************************************/
double kepler(double E, double e)
{
	return E - e*sin(E);
}

/*************************************************************************/
/*************************************************************************/
/**																		**/
/** Funtion: solve_kepler                                               **/
/** Purpose: solves Kepler's Equation for eccentric anomaly 			**/
/**	Inputs: mean anomaly M		                                        **/
/**			eccentricity												**/
/**         				                                            **/
/**	Output: eccentric anomaly as a double                               **/
/**																		**/
/*************************************************************************/
/*************************************************************************/
double solve_kepler(double M, double e)
{

	return bisect_kepler(0,2*PI,M,e,f_kepler,1e-9);

}


/*************************************************************************/
/*************************************************************************/
/**																		**/
/** Funtion: f_kepler	                                                **/
/** Purpose: function to pass to bisect for solving Kepler's eqn 		**/
/**	Inputs: mean anomaly M		                                        **/
/**			eccentricity												**/
/**         				                                            **/
/**	Output: eccentric anomaly as a double                               **/
/**																		**/
/*************************************************************************/
/*************************************************************************/
double f_kepler(double E, double M, double e)
{
	return E - e*sin(E) - M;
}


/*************************************************************************/
/*************************************************************************/
/**																		**/
/** Funtion: bisect_kepler												**/
/** Purpose: solve kepler's equation via bisection method    			**/
/**	Inputs: bracket values a and b, Mean anomaly M, eccentricity e		**/
/**			function pointer f, tolerance								**/
/**	Output: eccentric anomaly E, double									**/
/**																		**/
/*************************************************************************/
/*************************************************************************/

double bisect_kepler(double a, double b, double M, double e, double f(double E, double M, double e), double tol)
{

	double c1;		//c will hold the midpoint value

	int n;

	//first a little error checking, f(a)*f(b) < 0 is required
	if(((*f)(b,M,e))*((*f)(a,M,e)) >= 0)
	{
		printf("\nERROR: Invalid interval endpoints: f(a)*f(b) > 0\n");
		return 0;
	}



	//find the number of iterations needed

	n = ceil( (log(b - a) - log(tol)) / log(2) );

	for(int j=0;j<n;j++)
	{
		c1 = a + 0.5*(b-a);		//compute midpoint

		if( ((*f)(c1,M,e))*((*f)(a,M,e)) < 0 )	//then root lies between a and c
			b = c1;					//set new upper bound

		else if( ((*f)(c1,M,e))*((*f)(b,M,e)) < 0 )	//then root lies between b and c
			a = c1;						//set new lower bound
	}

	return c1;

}

/*
//test function
int main()
{
	double e = 0.995, M = 0.0001, E;
	printf("\n\nEccentricity e = 0.5\nMean Anomaly M = 1\n\n");

	E = solve_kepler(M,e);

	printf("Eccentric Anomaly is then E = %.3e\n",E);

	printf("Kepler's equation evaluated with that value: M = %.3e\n\n",kepler(E,e));

	return 0;
}
*/

