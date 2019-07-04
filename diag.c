#include <stdio.h>
#include <stdlib.h>
#include "f2c.h"
#include "fblaswr.h"
#include "clapack.h"
#include "time.h"
#include "math.h"

int pcg( doublereal *A, doublereal *x, doublereal *b, integer n, doublereal eps )
{
		
		int i, j, brojac = 0;
		integer m = 1, info;
		doublereal kriterij = 1.0;
		char uplo, trans, diag = 'n';
		
		doublereal *M = malloc( n*n*sizeof(doublereal) );
		doublereal *rez = malloc( n*sizeof(doublereal) );
		doublereal *p = malloc( n*sizeof(doublereal) );
		doublereal *d = malloc( n*sizeof(doublereal) );
		integer ipiv[ n ];
		
		//M je dijagonalna matrica, M^-1 = D^2, D = sqrt(1/diag(A))
        for( i = 0; i < n*n; i++ )
        	M[ i ] = 0.0;
  
        for( i = 0; i < n; i++ )
        	M[ i + i*n ] = A[ i + i*n ] ;

		//kopiramo b u rez 
		dlacpy_( &uplo, &m, &n, b, &m, rez, &m );
		
		//rez = b - A*x
		trans = 'n';
		doublereal alpha = -1.0;
		doublereal beta = 1.0;
		integer inc = 1;
		dgemv_( &trans, &n, &n, &alpha, A, &n, x, &inc, &beta, rez, &inc );
		
		//kopiramo rez u p
		dlacpy_( &uplo, &m, &n, rez, &m, p, &m );
		
		//M*p = rez;
		uplo = 'u';
		dtrsv_( &uplo, &trans, &diag, &n, M, &n, p, &inc );

		//kopiramo p u d
		dlacpy_( &uplo, &m, &n, p, &m, d, &m );

		//norma od b
		doublereal norma_b = dnrm2_( &n, b, &inc );
		
		while( kriterij > eps )
		{
				brojac++;
				
				doublereal pom[ n ];
				doublereal suma1 = 0.0, suma2 = 0.0, suma3 = 0.0, norma_r;
				
				suma1 = ddot_( &n, p, &inc, rez, &inc );
				
				//pom = A*d
				alpha = 1.0;
				beta = 0.0;
				dgemv_( &trans, &n, &n, &alpha, A, &n, d, &inc, &beta, pom, &inc );
		
				suma2 = ddot_( &n, pom, &inc, d, &inc );
		
				doublereal a = suma1/suma2; 
				
				//x = a*d + x
				daxpy_( &n, &a, d, &inc, x, &inc );
		
				//rez = -a*pom + rez
				a = -a;
				daxpy_( &n, &a, pom, &inc, rez, &inc );
				
				//kopiramo rez u p
				dlacpy_( &uplo, &m, &n, rez, &m, p, &m );
		
				//M*p = rez;
				uplo = 'u';
				dtrsv_( &uplo, &trans, &diag, &n, M, &n, p, &inc );

				suma3 = ddot_( &n, p, &inc, rez, &inc ); 
					
				doublereal b = suma3/suma1;
				
				for( i = 0; i <n; i++ )
					d[ i ] = p[ i ] + b*d[ i ];	
				
				//norma reziduala
				norma_r = dnrm2_( &n, rez, &inc );
				
				kriterij = norma_r/norma_b;
				printf( "Relativna norma reziduala u %d. koraku je %.16f\n", brojac, kriterij );

		}

		free( M );
		free( rez );
		free( p );
		free( d );
		return brojac;	
}
						
		
int main()
{
		int i;
		integer n = 100;
		doublereal eps = 1e-8;
		
		FILE *f;
		doublereal *A = malloc( n*n*sizeof(doublereal) );
		
		f = fopen( "stieltjes_matr.txt", "r" );
		for( i = 0; i < n*n; i++ )
			fscanf( f, "%lf", A + i );
		fclose(f);
		
		doublereal *b = malloc( n*sizeof(doublereal) );
		doublereal *x = malloc( n*sizeof(doublereal) );
		doublereal *x0 = malloc( n*sizeof(doublereal) );
		
		for( i = 0; i < n; i++ )
		{
			x[ i ] = 1.0;
			x0[ i ] = 0.0;
		}
		
		//mnozimo A i x da dobijemo egzaktni b
		char trans = 'n';
		doublereal alpha = 1.0;
		doublereal beta = 0.0;
		integer inc = 1;
		dgemv_( &trans, &n, &n, &alpha, A, &n, x, &inc, &beta, b, &inc );
		
		int broj_iteracija = pcg( A, x0, b, n, eps );
		printf( "Potrebno je %d iteracija.\n", broj_iteracija );

		free( A );
		free( b );
		free( x );
		free( x0 );
		
		return 0;
		
}
