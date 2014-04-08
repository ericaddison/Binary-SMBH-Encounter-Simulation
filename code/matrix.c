/*
	Eric Addison
	Utah State University
	Summer 2009

	This file contains functions for vector and matrix operations
*/



#include "matrix.h"
//void Mdisplay(int rows, int cols, double X[][cols]);

/*************************************************************************/
/*************************************************************************/
/**																		**/
/** Funtion: vect_mult_element											**/
/** Purpose: performs element by element multiplication of two vectors	**/
/**	Inputs: three arrays, two to multiply and one to hold the product	**/
/**			also need the size of the array passed in					**/
/**	Output: fill the product array										**/
/**																		**/
/*************************************************************************/
/*************************************************************************/

void vect_mult_element(double * x, double * y, double * prod, int size)
{
	for(int i = 0; i < size; i++)
		prod[i] = x[i]*y[i];
}

	//can't check the length of an array that has been passed into a function,
	//so the user will have to make sure dimensions match... No error checking


/*************************************************************************/
/*************************************************************************/
/**																		**/
/** Funtion: vect_add													**/
/** Purpose: performs simple vector addition							**/
/**	Inputs: three arrays, two to add and one to hold the sum			**/
/**			and the size												**/
/**	Output: fill the sum array											**/
/**																		**/
/*************************************************************************/
/*************************************************************************/

void vect_add(double * x, double * y, double * sum, int size)
{
	for(int i = 0; i < size; i++)
		sum[i] = x[i]+y[i];
}

/*************************************************************************/
/*************************************************************************/
/**																		**/
/** Funtion: vect_add2													**/
/** Purpose: performs simple vector addition, but does not overwrite    **/
/**			 values already present in sum array						**/
/**	Inputs: three arrays, two to add and one to hold the sum			**/
/**			and the size												**/
/**	Output: fill the sum array											**/
/**																		**/
/*************************************************************************/
/*************************************************************************/

void vect_add2(double * x, double * y, double * sum, int size)
{
	for(int i = 0; i < size; i++)
		sum[i] += x[i]+y[i];
}


	//can't check the length of an array that has been passed into a function,
	//so the user will have to make sure dimensions match... No error checking



/*************************************************************************/
/*************************************************************************/
/**																		**/
/** Funtion: vect_mult_scalar											**/
/** Purpose: performs scalar multiplication of a vector					**/
/**	Inputs: two arrays, the scalar, and the size						**/
/**																		**/
/**	Output: fill the result array										**/
/**																		**/
/*************************************************************************/
/*************************************************************************/

void vect_mult_scalar(double * x, double a, double * result, int size)
{
	for(int i=0;i<size;i++)
		result[i] = a*x[i];
}



/*************************************************************************/
/*************************************************************************/
/**																		**/
/** Funtion: outer_prod													**/
/** Purpose: performs an outer product of two vectors					**/
/**	Inputs: two arrays to multiply, one 2D array to store the resulting **/
/**			matrix, and the size										**/
/**																		**/
/**	Output: fill the product matrix										**/
/**																		**/
/*************************************************************************/
/*************************************************************************/

void outer_prod(double * x, double * y, int size, double prod[][size])
{
	for(int i=0;i<size;i++)
	{
		for(int j=0;j<size;j++)
			prod[i][j] = x[i]*y[j];
	}
}


/*************************************************************************/
/*************************************************************************/
/**																		**/
/** Funtion: dot_prod													**/
/** Purpose: performs a dot product of two vectors						**/
/**	Inputs: two arrays to multiply, and the size						**/
/**																		**/
/**	Output: double value of the product									**/
/**																		**/
/*************************************************************************/
/*************************************************************************/

double dot_prod(double * x, double * y, int size)
{

	double prod = 0;

	for(int i=0;i<size;i++)
		prod += x[i]*y[i];

	return prod;

}

/*************************************************************************/
/*************************************************************************/
/**																		**/
/** Funtion: cross_prod													**/
/** Purpose: vector cross product										**/
/**	Inputs: 3 component vectors A and B, and result vector X			**/
/**																		**/
/**	Output: result vector X will be filled								**/
/**																		**/
/*************************************************************************/
/*************************************************************************/

void cross_prod(double * A, double * B, double * X)
{

	X[0] = A[1]*B[2] - A[2]*B[1];
	X[1] = A[2]*B[0] - A[0]*B[2];
	X[2] = A[0]*B[1] - A[1]*B[0];

}

/*************************************************************************/
/*************************************************************************/
/**																		**/
/** Funtion: mag														**/
/** Purpose: find the magnitude of a vector								**/
/**	Inputs: vector x													**/
/**																		**/
/**	Output: double magnitude											**/
/**																		**/
/*************************************************************************/
/*************************************************************************/

double mag(double * x, int size)
{
	double sum=0;
	for(int i = 0; i < size; i++)
		sum += x[i]*x[i];

	return sqrt(sum);

}


/*************************************************************************/
/*************************************************************************/
/**																		**/
/** Funtion: mat_mult													**/
/** Purpose: performs matrix multiplication								**/
/**	Inputs: two 2D arrays to multiply, along with the dimensions of each**/
/**			also need to input pointer to product array					**/
/**	Output: fill the product matrix										**/
/**																		**/
/*************************************************************************/
/*************************************************************************/

void mat_mult(int rowsA, int colsA, int rowsB, int colsB, double A[][colsA], double B[][colsB], double prod[][colsB])
{

	//some immediate error checking, make sure inner dimensions agree
		if(colsA != rowsB)
		{
			printf("\nERROR: Inner matrix dimensions must agree\n");
			return;
		}


	double Bt[colsB][rowsB];



	//transpose B for easier dot products
	transpose(rowsB,colsB,B,Bt);


	for(int i=0;i<rowsA;i++)
	{
		for(int j=0;j<colsB;j++)
		{
			//dot product the ith row of A with the jth column of B
			prod[i][j] = dot_prod(A[i],Bt[j],colsA);
		}
	}

	return;
}



/*************************************************************************/
/*************************************************************************/
/**																		**/
/** Funtion: mat_add													**/
/** Purpose: add two matrices											**/
/**	Inputs: two 2D arrays to add, along with the dimensions				**/
/**			also need to input pointer to sum array						**/
/**	Output: fill the product matrix										**/
/**																		**/
/*************************************************************************/
/*************************************************************************/

void mat_add(int rows, int cols, double A[][cols], double B[][cols], double sum[][cols])
{
	for(int i=0;i<rows;i++)
	{
		for(int j=0;j<cols;j++)
		{
			sum[i][j] = A[i][j] + B[i][j];
		}
	}

	return;
}

/*************************************************************************/
/*************************************************************************/
/**																		**/
/** Funtion: mat_scalar_mult											**/
/** Purpose: performs scalar*matrix multiplication						**/
/**	Inputs: one 2d array, one double,									**/
/**			also need to input pointer to product array					**/
/**	Output: fill the product matrix										**/
/**																		**/
/*************************************************************************/
/*************************************************************************/

void mat_scalar_mult(int rows, int cols, double A[][cols], double a, double prod[][cols])
{

	for(int i=0;i<rows;i++)
	{
		for(int j=0;j<cols;j++)
		{
			prod[i][j] = a*A[i][j];
		}
	}

	return;
}



/*************************************************************************/
/*************************************************************************/
/**																		**/
/** Funtion: mat_vec_mult												**/
/** Purpose: performs matrix*vector multiplication						**/
/**	Inputs: one 2d array, one vector array,								**/
/**			also need to input pointer to product array					**/
/**	Output: fill the product matrix										**/
/**																		**/
/*************************************************************************/
/*************************************************************************/

void mat_vec_mult(int rows, int cols, double A[][cols], double * v, double * prod)
{
	
	for(int i=0;i<rows;i++)
	{
		//dot product the ith row of A with the jth column of B
		prod[i] = dot_prod(A[i],v,cols);
	}
	
	return;
}

/*************************************************************************/
/*************************************************************************/
/**																		**/
/** Funtion: vec_mat_mult												**/
/** Purpose: performs vector*matrix multiplication						**/
/**	Inputs: one 2d array, one vector array,								**/
/**			also need to input pointer to product array					**/
/**	Output: fill the product matrix										**/
/**																		**/
/*************************************************************************/
/*************************************************************************/

void vec_mat_mult(int rows, int cols, double A[][cols], double * v, double * prod)
{
	
	double At[cols][rows];
	//transpose B for easier dot products
	transpose(rows,cols,A,At);
	
	for(int i=0;i<cols;i++)
	{
		//dot product the ith row of A with the jth column of B
		prod[i] = dot_prod(At[i],v,cols);
	}
	
	return;
}


/*************************************************************************/
/*************************************************************************/
/**																		**/
/** Funtion: transpose													**/
/** Purpose: performs matrix transpose									**/
/**	Inputs: matrix to transpose, result matrix, sizes					**/
/**	Output: fill the result matrix										**/
/**																		**/
/*************************************************************************/
/*************************************************************************/

void transpose(int rows, int cols, double A[][cols], double At[][rows])
{
	for(int i=0;i<rows;i++)
		for(int j=0;j<cols;j++)
			At[j][i] = A[i][j];

	return;
}

/*************************************************************************/
/*************************************************************************/
/**																		**/
/** Funtion: trace														**/
/** Purpose: find the trace of a square matrix							**/
/**	Inputs: nxn matrix, size n											**/
/**	Output: double value of the trace									**/
/**																		**/
/*************************************************************************/
/*************************************************************************/

double trace(int n, double A[][n])
{
	double t=0;
	
	for(int i=0;i<n;i++)
		t += A[i][i];
	
	return t;
}


/*************************************************************************/
/*************************************************************************/
/**																		**/
/** Funtion: MInv														**/
/** Purpose: Finds the inverse of a matrix using LU decomposition		**/
/**	Inputs: matrix A, size, and another matrix to store the result		**/
/**	Output: fill the result matrix										**/
/**																		**/
/**	Translated from Numerical Recipes in Fortran 77, 2nd Ed, pg 40		**/
/**																		**/
/*************************************************************************/
/*************************************************************************/

void MInv(int size, double A[][size], double Y[][size])
{

	int d;
	double LU[size][size], Yt[size][size], indx[size];

	//initialize Y as identity
	for(int i=0;i<size;i++)
	{
		for(int j=0;j<size;j++)
			Yt[i][j] = 0;

		Yt[i][i] = 1;
	}

	d = LUdcmp(size,A,LU,indx);		//decompose matrix


	//Y is identity, so columns = rows

	for(int j=0;j<size;j++)			//find inverse row by row, then transpose
		LUbksub(size,LU,Yt[j],indx,Yt[j]);


	transpose(size,size,Yt,Y);


	return;
}


/*************************************************************************/
/*************************************************************************/
/**																		**/
/** Funtion: LUdcmp														**/
/** Purpose: performs an LU decomposition of a square matrix A			**/
/**				using Crout's method									**/
/**	Inputs: matrix A, size, array used for recording row permutations,	**/
/**			  and another matrix to store the result					**/
/**	Output: fill the result matrix, output +-1 depending on # of pivots	**/
/**																		**/
/**	Translated from Numerical Recipes in Fortran 77, 2nd Ed, pg 38		**/
/**																		**/
/*************************************************************************/
/*************************************************************************/

int LUdcmp(int size, double A[][size], double LU[][size], double * indx)
{
	int d=1,imax;			//d will be returned, indicates even or odd number or row interchanges

	double max,scale[size],sum,dum;
	double TINY = 1e-20;

	//copy A into LU since the algorithm is set up to change the existing matrix
	for(int i=0;i<size;i++)
		for(int j=0;j<size;j++)
			LU[i][j] = A[i][j];


	//find and save the largest element in each row, used for scaling
	for(int i=0;i<size;i++)
	{
		max=0;

		for(int j=0;j<size;j++)
			if( fabs(LU[i][j]) > max)
				max = fabs(LU[i][j]);

		if( max == 0 )
			{ printf("\nERROR: Singular Matrix\n"); return 0; }

		scale[i] = 1.0/max;		//save the scaling factor for each row
	}


	//loop over columns for Crout's method

	for(int j=0;j<size;j++)
	{

		for(int i=0;i<j;i++)				//computes equation (2.3.12) in NR, except for i=j
		{
			sum=LU[i][j];
			for(int k=0;k<i;k++)			//calculates betas from NR, except diagonal beta
				sum -= LU[i][k]*LU[k][j];
			LU[i][j] = sum;
		}

		max = 0.0;							//initialize search for largest pivot element

		for(int i=j;i<size;i++)				// the i=j part of (2.3.12) and i=j+1,... of (2.3.13)
		{
			sum = LU[i][j];
			for(int k=0;k<j;k++)
				sum -= LU[i][k]*LU[k][j];
			LU[i][j] = sum;

			dum = scale[i]*fabs(sum);		//figure of merit (???) for the pivot

			if( dum >= max )				//is it better than the best so far?
				{ imax = i; max = dum; }

		}

	//	printf("\n\n\ninterim LU at j = %d\n",j);
	//	printf("\nimax = %d\n",imax);
	//	Mdisplay(size,size,LU);

		if( j != imax )						//do we need to interchange rows?
		{
			for(int k=0;k<size;k++)
			{
				dum = LU[imax][k];
				LU[imax][k] = LU[j][k];
				LU[j][k] = dum;
			}
			d = -d;
			scale[imax] = scale[j];
		}

		indx[j] = imax;

		if( LU[j][j] == 0.0 )
			LU[j][j] = TINY;

		if( j != size-1 )					//Now divide by the pivot element
		{
			dum = 1/LU[j][j];				//this is the 1/beta(j,j) in NR eqn (2.3.13)
			for(int i=j+1;i<size;i++)
				LU[i][j] *= dum;
		}

	//	printf("\n\nlatest LU at j = %d\n",j);
	//	Mdisplay(size,size,LU);

	}


	return d;

}



/*************************************************************************/
/*************************************************************************/
/**																		**/
/** Funtion: LUbksub													**/
/** Purpose: performs back substitution of an LU decomposed matrix		**/
/**				in order to solve the matrix eqn Ax = b					**/
/**	Inputs: matrix LU, decomposition of A, size, indx array from LUdcmp	**/
/**			  	and an array to store the result						**/
/**				the b input vector is the RHS of the matrix eqn			**/
/**	Output: fill the result array										**/
/**																		**/
/**	Translated from Numerical Recipes in Fortran 77, 2nd Ed, pg 38		**/
/**																		**/
/*************************************************************************/
/*************************************************************************/

void LUbksub(int size, double LU[][size], double * b, double * indx, double * res)
{
	int ii=-1,ll;
	double sum;

	//copy b into res, algorithm destroys the vector as it solves
	for(int i=0;i<size;i++)
		res[i] = b[i];


	//forward substitution part
	for(int i=0;i<size;i++)
	{
		ll = indx[i];					//ll = row that row i switched with
		sum = res[ll];					//sum starts as that value of b = b(ll)
		res[ll] = res[i];				//switch so b(ll) = b(i), i.e. scramble b vector to match scarambled LU

		if(ii != -1)
			for(int j=ii;j<i;j++)		//first time though does nothing, ii=0
				sum -= LU[i][j]*res[j];

		else if(sum != 0.0)
			ii = 0;						//this was the problem! ii=0 isntead of ii=1!

		res[i] = sum;					//element i determined
	}


	//backward substitution part
	for(int i=size-1;i>=0;i--)
	{
		sum = res[i];
		for(int j=i+1;j<size;j++)
			sum -= LU[i][j]*res[j];

		res[i] = sum/LU[i][i];
	}

	return;
}





/*************************************************************************/
/*************************************************************************/
/**																		**/
/** Funtion: Vdisplay													**/
/** Purpose: prints out a vector										**/
/**	Inputs: the array to be printed and it's size						**/
/**																		**/
/**	Output: screen output												**/
/**																		**/
/*************************************************************************/
/*************************************************************************/

void Vdisplay(double * x, int size)
{
	for(int k = 0;k<size;k++)
	{
		if(k==0)
			printf("\n[ %.1e ",x[k]);
		else if(k==(size-1))
			printf("%.1e ]\n",x[k]);
		else
			printf("%.1e ",x[k]);
	}
}


/*************************************************************************/
/*************************************************************************/
/**																		**/
/** Funtion: Vdisplay_int												**/
/** Purpose: prints out a vector of integers							**/
/**	Inputs: the array to be printed and it's size						**/
/**																		**/
/**	Output: screen output												**/
/**																		**/
/*************************************************************************/
/*************************************************************************/

void Vdisplay_int(int * x, int size)
{
	for(int k = 0;k<size;k++)
	{
		if(k==0)
			printf("\n[ %d ",x[k]);
		else if(k==(size-1))
			printf("%d ]\n",x[k]);
		else
			printf("%d ",x[k]);
	}
}


/*************************************************************************/
/*************************************************************************/
/**																		**/
/** Funtion: Mdisplay													**/
/** Purpose: prints out a matrix										**/
/**	Inputs: the matrix to be printed and it's size						**/
/**																		**/
/**	Output: screen output												**/
/**																		**/
/*************************************************************************/
/*************************************************************************/

void Mdisplay(int rows, int cols, double x[][cols])
{
	for(int k = 0;k<rows;k++)
	{
		for(int j = 0;j<cols;j++)
		{
			if(j==0)
				printf("\n[ %.2e ",x[k][j]);
			else if (j==(cols-1))
				printf(" %.2e ]",x[k][j]);
			else
				printf(" %.2e ",x[k][j]);
		}
	}
	printf("\n\n");
}


/*************************************************************************/
/*************************************************************************/
/**																		**/
/** Funtion: Mcopy  													**/
/** Purpose: copies one matrix into another								**/
/**	Inputs: the input matrix and the target     						**/
/**																		**/
/**	Output: none        												**/
/**																		**/
/*************************************************************************/
/*************************************************************************/

void Mcopy(int rows, int cols, double A[][cols], double B[][cols])
{
	for(int i=0;i<rows;i++){
		for(int j=0;j<cols;j++)
			B[i][j] = A[i][j];
	}
	return;
}

/*************************************************************************/
/*************************************************************************/
/**																		**/
/** Funtion: Vcopy  													**/
/** Purpose: copies one vector into another								**/
/**	Inputs: the input vector and the target     						**/
/**																		**/
/**	Output: none        												**/
/**																		**/
/*************************************************************************/
/*************************************************************************/

void Vcopy(int size, double A[size], double B[size])
{
	for(int i=0;i<size;i++)
			B[i] = A[i];

	return;
}

/*************************************************************************/
/*************************************************************************/
/**																		**/
/** Funtion: Vcopy_int 													**/
/** Purpose: copies one vector of ints into another						**/
/**	Inputs: the input vector and the target     						**/
/**																		**/
/**	Output: none        												**/
/**																		**/
/*************************************************************************/
/*************************************************************************/

void Vcopy_int(int size, int A[size], int B[size])
{
	for(int i=0;i<size;i++)
			B[i] = A[i];

	return;
}



/*************************************************************************/
/*************************************************************************/
/**																		**/
/** Funtion: vec_rot  													**/
/** Purpose: rotate a vector											**/
/**	Inputs: vector to rotate, number indicating axis to rotate about	**/
/**			rotation angle												**/
/**																		**/
/**	Output: none        												**/
/**																		**/
/*************************************************************************/
/*************************************************************************/

void vec_rot(double x[][3], int axis, double theta)
{

	double prod[1][3];

	if(axis == 0)		//then rotate about the x-axis
	{
		double M[3][3] = { {1,0,0} , {0,cos(theta),sin(theta)} , {0,-sin(theta),cos(theta)} };
	    mat_mult(1,3,3,3,x,M,prod);
	    Mcopy(1,3,prod,x);
	    return;
	}

	if(axis == 1)		//then rotate about the y-axis
	{
		double M[3][3] = { {cos(theta),0,-sin(theta)}, {0,1,0} , {sin(theta),0,cos(theta)} };
	    mat_mult(1,3,3,3,x,M,prod);
	    Mcopy(1,3,prod,x);
	    return;
	}


	if(axis == 2)		//then rotate about the z-axis
	{
		double M[3][3] = { {cos(theta),sin(theta),0} , {-sin(theta),cos(theta),0} , {0,0,1} };
	    mat_mult(1,3,3,3,x,M,prod);
	    Mcopy(1,3,prod,x);
	    return;
	}

}


/*************************************************************************/
/*************************************************************************/
/**																		**/
/** Funtion: mat_min  													**/
/** Purpose: find the minimum value in a matrix							**/
/**	Inputs: rows, cols, matrix, ints to hold min coords					**/
/**																		**/
/**	Output: min value      												**/
/**																		**/
/*************************************************************************/
/*************************************************************************/

double mat_min(int rows, int cols, double A[rows][cols], int *min1, int*min2)
{
	double min = atof("Inf");

	for(int i=0;i<rows;i++){
		for(int j=0;j<cols;j++){

			if( A[i][j] < min){
				min = A[i][j];
				*min1 = i;
				*min2 = j;
			}
	}
	}
	return min;

}
/*************************************************************************/
/*************************************************************************/
/**																		**/
/** Funtion: vec_min  													**/
/** Purpose: find the minimum value in a vector							**/
/**	Inputs: length, vector, int to hold min position					**/
/**																		**/
/**	Output: min value      												**/
/**																		**/
/*************************************************************************/
/*************************************************************************/

double vec_min(int s, double *A, int *min1)
{
	double min = atof("Inf");

	for(int i=0;i<s;i++){
			if( A[i] < min){
				min = A[i];
				*min1 = i;
			}
	}
	return min;

}

double min(double x, double y)
{
	if(x > y)
		return y;

	return x;

}
/*

//test function
void main()
{

	double x[3] = {1,1,0};

	Vdisplay(x,3);

	vec_rot(x,2,PI/4);

	Vdisplay(x,3);

}

*/
