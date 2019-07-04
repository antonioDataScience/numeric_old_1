#include <stdio.h>
#include <stdlib.h>
#include "f2c.h"
#include "fblaswr.h"
#include "clapack.h"
#include "math.h"

void ispis_matrice( doublereal *A, integer m, integer n )
{
	int i, j;
	printf( "\n\n" );
	for( i = 0; i < m; i++ )
	{
		for( j = 0; j < n; j++ )
			printf( " %f ", A[ i + j*m ] );
		printf( "\n" );
	}
	printf( "\n\n" );
}

int main()
{
		integer m = 4, n = 2;
		integer N = n*n;
		int i;
		
		doublereal A[ ] = { 1, 0, 0, 0, 1, 1, 1, 1 };
		doublereal B[ ] = { 1, -1, 1, -1, 0, 1, 0, 1 };
	
		//QR faktorizacija od A i B
		doublereal *TAUA = malloc( n*sizeof(doublereal) );
		doublereal *TAUB = malloc( n*sizeof(doublereal) );
		doublereal *work = malloc( N*sizeof(doublereal) );
		integer info;
		dgeqrf_( &m, &n, A, &m, TAUA, work, &N, &info );
		dgeqrf_( &m, &n, B, &m, TAUB, work, &N, &info );
		
		//u A spremimo matricu Q iz QR faktorizacije od A
		dorgqr_( &m, &n, &n, A, &m, TAUA, work, &N, &info );
		
		//ispis_matrice( A, m, n );
	
		//u B spremimo matricu Q iz QR faktorizacije od B
		dorgqr_( &m, &n, &n, B, &m, TAUB, work, &N, &info );
		
		
		//ispis_matrice( B, m, n );
		
		//C = A^T*B
		doublereal C[ N ];
		char transa = 't';
		char transb = 'n';
		doublereal alpha = 1.0;
		doublereal beta = 0.0;
		dgemm_( &transa, &transb, &n, &n, &m, &alpha, A, &m, B, &m, &beta, C, &n );
		
		//ispis_matrice( C, n, n );
		
		//U^T*C*V = S, S je matrica singularnih vrijednosti
		char jobu = 'a';
		char jobvt = 'a';
		doublereal sing[ n ];
		doublereal U[ N ];
		doublereal VT[ N ];
		integer ldw = 5*n;
		doublereal work1[ ldw ];
		dgesvd_( &jobu, &jobvt, &n, &n, C, &n, sing, U, &n, VT, &n, work1, &ldw, &info );
		
		//ispis kuteva
		printf( "Kosinusi kuteva iznose: " );
		for( i = 0; i < n; i++ )
			printf( " %f ", sing[ i ] );
		printf("\n\n");
			
		//X=A*U, Y=B*V
		doublereal X[ m*n ], Y[ m*n ];
		char trans = 'n';
		dgemm_( &trans, &trans, &m, &n, &n, &alpha, A, &m, U, &n, &beta, X, &m );
		dgemm_( &trans, &transa, &m, &n, &n, &alpha, B, &m, VT, &n, &beta, Y, &m );
		
		printf("Vektori x:");
		ispis_matrice( X, m, n );
		printf("Vektori y:");
		ispis_matrice( Y, m, n );
		

		return 0;
}
		
