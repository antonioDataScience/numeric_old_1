#include <stdio.h>
#include <stdlib.h>
#include "f2c.h"
#include "fblaswr.h"
#include "clapack.h"
#include "math.h"

void ispis_matrice( doublereal *A, integer n )
{
	int i, j;
	printf( "\n\n" );
	for( i = 0; i < n; i++ )
	{
		for( j = 0; j < n; j++ )
			printf( " %f ", A[ i + j*n ] );
		printf( "\n" );
	}
	printf( "\n\n" );
}

int main()
{
	integer n = 4;
	integer N = n*n;
	int i, j;

	doublereal M[ ] = { 2, 0, 0, 0, 0, 5, 0, 0, 0, 0, 3, 0, 0, 0, 0, 6 };
	doublereal K[ ] = { 24, -9, -5, 0, -9, 22, -8, -5, -5, -8, 25, -7, 0, -5, -7, 18 };

	//racunamo D = M^(-1/2)
	doublereal D[ N ];
	for( i = 0; i < n; i++ )
		for( j = 0; j < n; j++ )
			if( i == j )
				D[ i +j*n ] = 1/sqrt( M[ i + j*n ] );
			else
				D[ i + j*n ] = M[ i + j*n ];

	//A=D*K*D
	doublereal A[ N ];
	doublereal pom[ N ];
	char trans = 'n';
	doublereal alpha = 1.0;
	doublereal beta = 0.0;
	dgemm_( &trans, &trans, &n, &n, &n, &alpha, D, &n, K, &n, &beta, pom, &n );
	dgemm_( &trans, &trans, &n, &n, &n, &alpha, pom, &n, D, &n, &beta, A, &n );

	//ispis_matrice( A, n );

	//trazimo svojstvene parove od A; stupci matrice A su svojstveni vektori
	char jobz = 'v';
	char uplo = 'l';
	doublereal polje_sv[ n ];
	integer ldw = 3*n - 1;;
	doublereal work[ ldw ];
	integer info;
	dsyev_( &jobz, &uplo, &n, A, &n, polje_sv, work, &ldw, &info );

	//ispis svojstvenih vektora
	printf("Svojstveni vektori su stupci ove matrice:\n" );
	ispis_matrice( A, n );

	//ispis svojstvenih vrijednosti
	printf("Odgovarajuce svojstvene vrijednosti:\n\n" );
	for( i = 0; i < n; i++ )
		printf( " %f ", polje_sv[ i ] );

	//matrica X = D*A
	doublereal X[ N ];
	dgemm_( &trans, &trans, &n, &n, &n, &alpha, D, &n, A, &n, &beta, X, &n );
	ispis_matrice( X, n );

	return 0;
}


