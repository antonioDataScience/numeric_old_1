#include <stdio.h>
#include <stdlib.h>
#include "f2c.h"
#include "fblaswr.h"
#include "clapack.h"
#include "time.h"
#include "math.h"


doublereal sor_norma( doublereal *A, integer n, doublereal omega )
{
	int i, j, N;
	N = n*n;
	
	char uplo = 'l';
	char uplo_1 = 'u';
	doublereal *L = malloc( N*sizeof(doublereal) );
	doublereal *R = malloc( N*sizeof(doublereal) );
	
	for( i = 0; i < N; i++ )
	{
		L[ i ] = 0.0;
		R[ i ] = 0.0;
	}
		
	dlacpy_( &uplo, &n, &n, A, &n, L, &n );//L je donji trokut
	dlacpy_( &uplo_1, &n, &n, A, &n, R, &n );//R je gornji trokut
	
	for( i = 0; i < n; i++ )
	{
		L[ i + i*n ] = 0.0; //postavljanje dijagonale na 0
		R[ i + i*n ] = 0.0; //postavljanje dijagonale na 0
	}

	doublereal *D = malloc( N*sizeof(doublereal) );
	for( i = 0; i < N; i++ )
		D[ i ] = 0.0;
	for( i = 0; i < n; i++ )
		D[ i + i*n ] = A[ i + i*n ];	
	
	//racunamo D + omega*L; rezultat se spremi u matricu L!!!

	for( i = 0; i < n; i++ )
	{
		for( j = 0; j < i; j++ )
			L[ i + j*n ] = omega*L[ i + j*n ];
		L[ i + i*n ] = D[ i + i*n ];
	}	

	//racunamo (1-omega)*D - omega*R; rezultat se spremi u matricu R!!!
	
	for( i = 0; i < n; i++ )
	{
		for( j = i + 1; j < n; j++ )
			R[ i + j*n ] = - omega*R[ i + j*n ];
		R[ i + i*n ] = ( 1 - omega )*D[ i + i*n ];
	}	

	//racunamo L^(-1)*R, tj. rjesavamo sustav LT = R
	
	char side = 'l';
	char uplo_2 = 'l';
	char transa = 'n';
	char diag = 'n';
	doublereal alpha = 1.0;
	
	dtrsm_( &side, &uplo_2, &transa, &diag, &n, &n, &alpha, L, &n, R, &n );
	
	//racunamo normu dobivene matrice T koja je zapravo spremljena u R!!!
	
	char norm = 'i';
	doublereal *work = malloc( n*sizeof(doublereal) );
	
	doublereal norma = dlange_( &norm, &n, &n, R, &n, work );
	
	printf("%f\n",norma );

	return norma;

}

integer sor_konvergencija( doublereal *A, integer n, doublereal omega )
{
	if ( sor_norma( A, n, omega ) < 1 )
		return 1;
	else
		return 0;
}


int kriterij( doublereal *A,  doublereal *x, doublereal *b, doublereal omega, doublereal eps, integer n )
{
	int i, j;
	integer m = 1;
	char uplo;
	doublereal x0[n];
	dlacpy_( &uplo, &m, &n, x, &m, x0, &m);//kopiramo x u x0
	
	//racunamo prvu iteraciju
	
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

	//oduzmemo x0 od x i to spremimo u x
	
	doublereal alpha = -1.0;
	integer incx = 1, incy  = 1;
	daxpy_(&n, &alpha, x0, &incx, x, &incy);
	
	//racunamo normu beskonacno razlike

	char norm = 'm';	
	doublereal *work = malloc( n*sizeof(doublereal) );
	doublereal max = dlange_(&norm, &m, &n, x, &m, work);
		
    doublereal norma = sor_norma( A, n, omega );
	
	doublereal rezultat = log( (eps*(1 - norma))/max ) / log( norma );

	int k = ceil(rezultat);
	
	return k;
}

void sor_rjesavac( doublereal *A, doublereal *x, doublereal *b, doublereal omega, integer n, integer k )
{	
	int brojac = 0;
	int i, j;
	while( brojac <= k )
	{ 
			brojac++;
			for( i = 0; i < n; i++ )
			{
				x[ i ] = ( 1 - omega ) * x [ i ];
				doublereal pom = b[ i ];
				for( j = 0; j < i; j++ )
					pom -= A[ i + j*n ]*x[ j ];
				for( j = i + 1; j < n; j++ )
					pom -= A[ i + j*n ]*x[ j ];
				x[ i ] += pom * omega/A[ i + i*n ];
			}
	}
	
	printf("Rjesenje dobiveno SOR metodom:\n");
	for( i = 0; i < n; i++ )
		printf( " %f ", x[ i ] );
		printf("\n");
}	


main(integer argc, char *argv[])
{
	doublereal A[] = { 101, -4, 8, 12,   -4, 20, -7, 3,   8, -7, 78, 32,   12, 3, 32, 113 };
	doublereal b[] = { 117, 12, 111, 160 };

	integer n = 4;
	
	/*doublereal *x0 = malloc( n*sizeof(doublereal) );
	
	integer idist = 3;
	integer iseed[] = { 3, 78, 1016, 4095 };
	dlarnv_( &idist, iseed, &n, x0 );*/
	
	doublereal x0[] = { 0, 0, 0, 0 };

	doublereal omega = 1.05;
	doublereal eps = 1e-5;
	
	if ( sor_konvergencija( A, n, omega ) )
	{
		int k = kriterij( A, x0, b, omega, eps, n );
		printf("Za SOR metodu potrebno je %d iteracija.\n", k );
		
		sor_rjesavac( A, x0, b, omega, n, k );
		printf("Kraj!\n");
	}
	else
		printf("Ne konvergira!!!\n");
	
}
	
	
