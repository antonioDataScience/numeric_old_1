#include <stdio.h>
#include <stdlib.h>
#include "f2c.h"
#include "fblaswr.h"
#include "clapack.h"
#include "time.h"
#include "math.h"

int sor_rjesavac( doublereal *A, doublereal *x, doublereal *b, doublereal eps, integer m, integer n )
{
			doublereal y[ n ];
			doublereal omega = 1.25;
			doublereal kriterij = 1.0;
			integer inc = 1;
			doublereal norma, norma_b;
			int brojac = 0;
			int i, j;
	
			//racunamo normu od b
			norma_b = dnrm2_( &n, b, &inc ); 
			printf("Norma vektora b iznosi %f.\n", norma_b );
	
			while( kriterij > eps )
			{
				for( i = 0; i < n; i++ )
				{
					x[ i ] = ( 1 - omega ) * x [ i ];
					doublereal pom = b[ i ];
					for( j = 0; j < i; j++ )
						pom = pom - A[ i + j*n ]*x[ j ];
					for( j = i + 1; j < n; j++ )
						pom = pom - A[ i + j*n ]*x[ j ];
					x[ i ] = x[ i ] + pom * omega / A[ i + i*n ];
				} 
		
				brojac++;
		
				//vektor b kopiramo u vektor y
				char uplo;
				dlacpy_( &uplo, &m, &n, b, &m, y, &m );
		
				//racunamo y - A*x i to spremimo u y
				char trans = 'n';
				doublereal alpha = -1.0;
				doublereal beta = 1.0;
				dgemv_( &trans, &n, &n, &alpha, A, &n, x, &inc, &beta, y, &inc );
		
				//racunamo normu od reziduala
				norma = dnrm2_( &n, y, &inc );
		
				kriterij = norma/norma_b;
			}
	
			for( i = 0; i < n; i++ )
			{
				printf("%.16f", x[ i ] );
				printf("\n");
			}
	
			return brojac;
}	
main(integer argc, char *argv[])
{
		int i, k;
		integer m = 1, n = 100;
		doublereal eps = 1e-5;
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
		/*integer idist = 3;
		integer iseed[] = { 3, 78, 1016, 4095 };
		dlarnv_( &idist, iseed, &n, x0 );*/
		
		//mnozimo A i x da dobijemo egzaktni b
		char trans = 'n';
		doublereal alpha = 1.0;
		doublereal beta = 0.0;
		integer inc = 1;
		dgemv_( &trans, &n, &n, &alpha, A, &n, x, &inc, &beta, b, &inc );
		
			
		/*for( i = 0; i < n; i++ )
			{
				printf("%f", b[ i ] );
				printf("\n");
			}*/
			
		k = sor_rjesavac( A, x0, b, eps, m, n );
		printf("Potrebno je %d iteracija.\n", k );
		
		free( A );
		free( x );
		free( b );
		free( x0 );
}
		
		

