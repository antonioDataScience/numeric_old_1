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


double norm_out( doublereal *A, integer n )
{
	int i,j;
	double suma = 0.0;
	
	for( i = 0; i < n; i++ )
		for( j = 0; j < n; j++ )
			if( i != j )
				suma += pow( A[ i + j*n ], 2 );
				
	suma = sqrt( suma );
	
	return suma;
}


int sgn( doublereal x )
{
	if( x >= 0 )
		return 1;
	else
		return -1;
}	


void jacobi_sd( doublereal *A, integer n, doublereal tol )
{
	int p, q, k;
	doublereal tau, t, c, s;
	doublereal app, apq, aqq, pom;

	
	doublereal kriterij = norm_out( A, n );
	
	while( kriterij > tol )
	{
		for( p = 0; p < n - 1; p++ )
		{
			for( q = p + 1; q < n; q++ )
			{
				if( A[ p + q*n ] != 0 )
				{
					tau = ( A[ q + q*n ] - A[ p + p*n ] )/( 2*A[ p + q*n ] );
					t = sgn( tau )/( fabs( tau ) + sqrt( 1 + pow( tau, 2 ) ) );
					c = 1/( sqrt( 1 + pow( t, 2 ) ) );
					s = t*c;
				}
				else
				{
					c = 1;
					s = 0;
				}
				
				app = A[ p + p*n ];
				apq = A[ p + q*n ];
				aqq = A[ q + q*n ];
				
				app = app - t*apq;
				aqq = aqq + t*apq;
				
				for( k = 0; k < n; k++ )
				{
					pom = A[ k + p*n ];
					A[ k + p*n ] = c*pom - s*A[ k + q*n ];
					A[ k + q*n ] = s*pom + c*A[ k + q*n ];
					
					A[ p + k*n ] = A[ k + p*n ];
					A[ q + k*n ] = A[ k + q*n ];
				}
				
				A[ p + q*n ] = 0;
				A[ q + p*n ] = 0;
				A[ p + p*n ] = app;
				A[ q + q*n ] = aqq;
			}
		}
		kriterij = norm_out( A, n );
	
	}
}

int main()
{
		integer n = 4;
		integer N = n*n;
		doublereal eps = 4e-16;
		
		
		//generiranje slucajne matrice V	
		doublereal *V = malloc( N*sizeof( doublereal ) );
		integer idist = 1;
		integer iseed[] = { 67, 928, 1555, 13 };
		dlarnv_( &idist, iseed, &N, V );
			
		//racunanje QR faktorizacije od V
		doublereal *TAU = malloc( n*sizeof(doublereal) );
		doublereal *work = malloc( N*sizeof(doublereal) );
		integer info;
		dgeqrf_( &n, &n, V, &n, TAU, work, &N, &info );
		
		//u V spremimo matricu Q iz QR faktorizacije od V
		dorgqr_( &n, &n, &n, V, &n, TAU, work, &N, &info );
		
		//dijagonalna matrica
		doublereal *D = malloc( N*sizeof( doublereal ) );
		int i;
		for( i = 0; i < N; i++ )
			D[ i ] = 0.0;
		D[ 0 ] = -10;
		D[ 1 + 1*n ] = -5;
		D[ 2 + 2*n ] = 0.1;
		D[ 3 + 3*n ] = 0.2;
		
		//A = VDV^T
		doublereal *POM = malloc( N*sizeof( doublereal ) );
		doublereal *A = malloc( N*sizeof( doublereal ) );
		char trans = 'n', transt = 't';
		doublereal alpha = 1.0;
		doublereal beta = 0.0;
		dgemm_( &trans, &trans, &n, &n, &n, &alpha, V, &n, D, &n, &beta, POM, &n );
		dgemm_( &trans, &transt, &n, &n, &n, &alpha, POM, &n, V, &n, &beta, A, &n );
				
	
		//Frobeniusova norma matrice A
		char norm = 'f';
		doublereal norma = dlange_( &norm, &n, &n, A, &n, work );
		
		//tolerancija
		doublereal tol = eps*norma;
		
		//pozivamo Jacobijevu metodu
		jacobi_sd( A, n, tol );
		
		ispis_matrice( A, n );
		
		//ispis svojstvenih vrijednosti
		/*for( i = 0; i < n; i++ )
			printf( "%d. %f\n", i + 1, A[ i + i*n ] );
		printf( "\n" );*/
		
		free( V );
		free( TAU );
		free( work );
		free( D );
		free( POM );
		free( A );
		
		return 0;
}
		
