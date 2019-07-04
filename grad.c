#include <stdio.h>
#include <stdlib.h>
#include "f2c.h"
#include "fblaswr.h"
#include "clapack.h"
#include "time.h"
#include "math.h"

void ispis_matrice( doublereal *A, integer n )
{
	int i, j;
	for( i = 0; i < n; i++ )
	{
		for( j = 0; j < n; j++ )
			printf( " %f ", A[ i + j*n ] );
		printf( "\n" );
	}
}

void ispis_vektora( doublereal *x, integer n )
{
	int i;
	for( i = 0; i < n; i++ )
	{
		printf( "%f", x[ i ] );
		printf( "\n" );
	}
}

int gc( doublereal *A, doublereal *x, doublereal *b, integer n, doublereal eps )
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

		// rez = -a*pom + rez
		a = -a;
		daxpy_( &n, &a, pom, &inc, rez, &inc );

		for( i = 0; i < n; i++ )
			suma3 += rez[ i ]*rez[ i ];

		doublereal c = suma3/suma1;

		for( i = 0; i < n; i++ )
			d[ i ] = rez[ i ] +  c*d[ i ];

		kriterij = (sqrt(suma3))/norma_b;

		printf( "Relativna norma reziduala u %d. koraku: %.16f\n", brojac, kriterij );


	}
	return brojac;

}



main( integer argc, char *argv[] )
{
	int i;
	integer n = 100;
	integer N = n*n;
	doublereal eps = 1e-8;
	printf("%.16f\n", eps );
	doublereal kapa = 1e4;
	doublereal *A = malloc( n*n*sizeof( doublereal ) );
	doublereal *S = malloc( n*n*sizeof( doublereal ) );
	doublereal *Q = malloc( n*n*sizeof( doublereal ) );

	//generiranje slucajne matrice A
	integer idist = 2;
	integer iseed[] = { 67, 928, 1555, 13 };
	dlarnv_( &idist, iseed, &N, A );

	//ispis_matrice( A, n );

	doublereal *b = malloc( n*sizeof(doublereal) );
	doublereal *x = malloc( n*sizeof(doublereal) );
	doublereal *x0 = malloc( n*sizeof(doublereal) );

	//pocetna iteracija je nul-vektor
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

	//ispis_vektora( b, n );

	//dijagonalna matrica svojstvenih vrijednosti
	for( i = 0; i < N; i++ )
		S[ i ] = 0.0;
	for( i = 0; i < n; i++ )
		S[ i + i*n ] = (i + 1)*(i + 1);

	//ispis_matrice( S, n );

	//racunanje QR faktorizacije od A
	doublereal *TAU = malloc( n*sizeof(doublereal) );
	doublereal *work = malloc( N*sizeof(doublereal) );
	char info;
	dgeqrf_( &n, &n, A, &n, TAU, work, &N, &info );

	//racunamo S = Q*S*Q^T
	char side1 = 'l';
	char trans1 = 'n';
	dormqr_( &side1, &trans1, &n, &n, &n, A, &n, TAU, S, &n, work, &N, &info );

	char side2 = 'r';
	char trans2 = 't';
	dormqr_( &side2, &trans2, &n, &n, &n, A, &n, TAU, S, &n, work, &N, &info );

	//ispis_matrice( S, n );

	int broj_iteracija = gc( S, x0, b, n, eps );
	printf( "Potrebno je %d iteracija.\n", broj_iteracija );

}



