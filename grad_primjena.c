#include <stdio.h>
#include <stdlib.h>
#include "f2c.h"
#include "fblaswr.h"
#include "clapack.h"
#include "math.h"


int cg( doublereal *A, doublereal *x, doublereal *b, integer n, doublereal eps )
{
	int i, j;
	int brojac = 0;
	doublereal kriterij = 1.0;
	integer m = 1;
 
	//kopirali smo vektor b u vektore rez i d
	doublereal *rez = malloc( n*sizeof(doublereal) );
	doublereal *d = malloc( n*sizeof(doublereal) );  
	char uplo;
	dlacpy_( &uplo, &m, &n, b, &m, rez, &m );
	dlacpy_( &uplo, &m, &n, b, &m, d, &m );
	
	//racunamo b - A*x
	char trans = 'n';
	doublereal alpha = -1.0;
	doublereal beta = 1.0;
	integer inc = 1;
	dgemv_( &trans, &n, &n, &alpha, A, &n, x, &inc, &beta, rez, &inc );//rez = b - A*x
	dgemv_( &trans, &n, &n, &alpha, A, &n, x, &inc, &beta, d, &inc );//d = b - A*x

	//norma od b
	doublereal norma_b = dnrm2_( &n, b, &inc );
	
	while( kriterij > eps )
	{	
			
		brojac++;
		
		doublereal pom[ n ];
		doublereal suma1 = 0.0, suma2 = 0.0, suma3 = 0.0;
		
		for( i = 0; i < n; i++ )
			suma1 += rez[ i ]*rez[ i ];		
		
		//pom = A*d
		alpha = 1.0;
		beta = 0.0;
		dgemv_( &trans, &n, &n, &alpha, A, &n, d, &inc, &beta, pom, &inc );
		
		for( i = 0; i < n; i++ )
			suma2 += pom[ i ]*d[ i ];
		
		doublereal a = suma1/suma2; 
		
		//x = a*d + x
		daxpy_( &n, &a, d, &inc, x, &inc );
		
		//rez = -a*pom + rez
		a = -a;
		daxpy_( &n, &a, pom, &inc, rez, &inc );
		
		for( i = 0; i < n; i++ )
			suma3 += rez[ i ]*rez[ i ];
		
		doublereal c = suma3/suma1;
		
		for( i = 0; i < n; i++ )
			d[ i ] = rez[ i ] +  c*d[ i ];	
					
		kriterij = (sqrt(suma3))/norma_b;

		//printf( "Relativna norma reziduala u %d. koraku: %.16f\n", brojac, kriterij );	

			
	}
	return brojac;

}
