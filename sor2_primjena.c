#include <stdio.h>
#include <stdlib.h>
#include "f2c.h"
#include "fblaswr.h"
#include "clapack.h"
#include "time.h"
#include "math.h"

int sor_rjesavac( doublereal *A, doublereal *x, doublereal *b, doublereal omega, doublereal eps, integer m, integer n )
{

			doublereal y[ n ];
			doublereal kriterij = 1.0;
			integer inc = 1;
			doublereal norma, norma_b;
			int brojac = 0;
			int i, j;
	
			//racunamo normu od b
			norma_b = dnrm2_( &n, b, &inc ); 
			//printf("Norma vektora b iznosi %f.\n", norma_b );
	
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
				char uplo = 'n';
				dlacpy_( &uplo, &m, &n, b, &m, y, &m );
	
				//racunamo y - A*x i to spremimo u y
				char trans = 'N';
				doublereal alpha = -1.0;
				doublereal beta = 1.0;
				dgemv_( &trans, &n, &n, &alpha, A, &n, x, &inc, &beta, y, &inc );

				//racunamo normu od reziduala
				norma = dnrm2_( &n, y, &inc );
		
				kriterij = norma/norma_b;
			}
			
			/*printf( "Rjesenje:\n" );
			for( i = 0; i < n; i++ )
				printf(" %f ", x[ i ] );
			printf("\n\n");*/
		
	
			return brojac;
}	

