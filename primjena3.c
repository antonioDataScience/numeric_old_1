#include <stdio.h>
#include <stdlib.h>
#include "f2c.h"
#include "fblaswr.h"
#include "clapack.h"
#include "math.h"
#include "sor2_primjena.c"
#include "grad_primjena.c"

void ispis_matrice( doublereal *A, integer n )
{
	int i, j;
	for( i = 0; i < n; i++ )
	{
		for( j = 0; j < n; j++ )
			printf( " %f ", A[ i + j*n ] );
		printf( "\n" );
	}
	printf( "\n\n\n" );	
}

int main()
{

	integer n = 81;
	integer N = n*n;
	integer m = 1;
	int i, j, k;
	
	doublereal omega;
	doublereal eps = 1e-8;
	
	//generiramo matricu sustava
	doublereal *A = malloc( N*sizeof( doublereal ) );
	
	char uplo;
	doublereal alpha = 0.0;
	doublereal beta = 400;
	dlaset_( &uplo, &n, &n, &alpha, &beta, A, &n );
	
	//elementi ispod dijagonale
	for( i = 0; i < n - 1; i++ )
		if( ( i + 1 ) % 9 != 0  )
			A[ i + 1 + i*n ] = -100;
	
	//elementi iznad dijagonale	
	for( i = 1; i < n; i++ )
		if( i % 9 != 0 )
			A[ i - 1 + i*n ] = -100;		
	
	//jedinicne matrice
	for( i = 0; i < n - 9; i++ )
		 A[ i + 9 + i*n ] = -100;
	for( i = 9; i < n; i++ )
		 A[ i - 9 + i*n ] = -100;
		 
	//ispis_matrice( A, n );	
	
	//desna strana sustava
	doublereal *b = malloc( n*sizeof( doublereal ) );
	for( i = 0; i < n; i++ )
		b[ i ] = 0.0;
	b[ 40 ] = 10000;
	
	//pocetna iteracija
	doublereal *u0 = malloc( n*sizeof( doublereal ) );
	doublereal *u1 = malloc( n*sizeof( doublereal ) );
	doublereal *u2 = malloc( n*sizeof( doublereal ) );
	for( i = 0; i < n; i++ )
	{
		u0[ i ] = 0.0;
		u1[ i ] = 0.0;
		u2[ i ] = 0.0;
	}
	
	
	//Gauss-Seidel
	omega = 1;
	int k0 = sor_rjesavac( A, u0, b, omega, eps, m, n );
	printf("Za Gauss-Seidelovu metodu potrebno je %d iteracija.\n\n", k0 );
	
	ispis_matrice( u0, 9 );


	//SOR 
	omega = 1.53;
	int k1 = sor_rjesavac( A, u1, b, omega, eps, m, n );
	printf("Za SOR metodu potrebno je %d iteracija.\n\n", k1 );	
	
	ispis_matrice( u1, 9 );

		
	//konjugirani gradijenti
	int k2 = cg( A, u2, b, n, eps );
	printf("Za metodu konjugiranih gradijenata je potrebno %d iteracija.\n\n", k2 );
	
	ispis_matrice( u2, 9 );
	
	
	free( A );
	free( b );
	free( u0 );
	free( u1 );
	free( u2 );
	
	return 0;
}
