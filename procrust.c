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
		integer m = 4, n = 2;
		integer N = n*n;
		int i;
		doublereal A[ ] = { 1.2, 2.9, 5.2, 6.8, 2.1, 4.3, 6.1, 8.1 };
		doublereal B[ ] = { 1, 3, 5, 7, 2, 4, 6, 8 };
		
		//C = A^T*B
		doublereal C[ N ];
		char transa = 't';
		char transb = 'n';
		doublereal alpha = 1.0;
		doublereal beta = 0.0;
		dgemm_( &transa, &transb, &n, &n, &m, &alpha, A, &m, B, &m, &beta, C, &n );
		
		//U^T*C*V = S, S je matrica singularnih vrijednosti
		char jobu = 'a';
		char jobvt = 'a';
		doublereal sing[ n ];
		doublereal U[ N ];
		doublereal VT[ N ];
		integer ldw = 5*n;
		doublereal work[ ldw ];
		integer info;
		dgesvd_( &jobu, &jobvt, &n, &n, C, &n, sing, U, &n, VT, &n, work, &ldw, &info );
		
		//ispis singularnih vrijednosti
		printf( "Singularne vrijednosti su:" );
		for( i = 0; i < n; i++ )
			printf( " %f ", sing[ i ] );
		
		//Q = U*V^T
		doublereal Q[ N ];
		char trans = 'n';
		dgemm_( &trans, &trans, &n, &n, &n, &alpha, U, &n, VT, &n, &beta, Q, &n );
		ispis_matrice( Q, n );
		
		//Frobeniusova norma matrice AQ-B
		beta = -1.0;
		dgemm_( &trans, &trans, &m, &n, &n, &alpha, A, &m, Q, &n, &beta, B, &m );
		char norm = 'f';
		doublereal norma = dlange_( &norm, &m, &n, B, &m, work );
		printf( "Frobeniusova norma matrice AQ-B: %f\n", norma );
		return 0;
}
