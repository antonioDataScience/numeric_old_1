#include <stdio.h>
#include <stdlib.h>
#include "f2c.h"
#include "fblaswr.h"
#include "clapack.h"
#include "math.h"
#include "sor2_primjena.c"
#include "grad_primjena.c"

int main()
{
	integer m = 1;
	integer n = 99;
	integer N = n*n;
	doublereal h = 0.01;
	doublereal eps = 1e-8;
	doublereal omega;
	int i;
	
	doublereal *A = malloc( N*sizeof( doublereal ) );
	
	for( i = 0; i < N; i++ )
		A[ i ] = 0.0;
	for( i = 0; i < n; i++ )
		A[ i + i*n ] = 1.9999;
		
	for( i = 0; i < n-1; i++ )
		A[ i + 1 + i*n ] = -1;
		
	for( i = 1; i < n; i++ )
		A[ i - 1 + i*n ] = -1;

	
	/*//kopiramo matricu A u matricu B
	doublereal *B = malloc( N*sizeof( doublereal ) );
	char uplo;
	dlacpy_( &uplo, &n, &n, A, &n, B, &n );

	//pomocu faktorizacije Choleskog provjerimo da je matrica A pozitivno definitna
	uplo = 'u';
	integer info;
	dpotrf_( &uplo, &n, B, &n, &info );
	
	if( info == 0 )
		printf( "Matrica A je pozitivno definitna.\n\n" );
	else
		printf( "Matrica A nije pozitivno definitna.\n\n" );*/
	
	
	//desna strana sustava
	doublereal *b = malloc( n*sizeof( doublereal ) );
	for( i = 0; i < n - 1; i++ )
		b[ i ] = 2*h*h*sin( (i + 1)*h );
	b[ 98 ] = 2*h*h*sin( 99*h ) + cos( 1 );
	
	//pocetna iteracija
	doublereal *x0 = malloc( n*sizeof( doublereal ) );
	doublereal *x1 = malloc( n*sizeof( doublereal ) );
	doublereal *x2 = malloc( n*sizeof( doublereal ) );
	for( i = 0; i < n; i++ )
	{
		x0[ i ] = 0.0;
		x1[ i ] = 0.0;
		x2[ i ] = 0.0;
	}
	
	//egzaktno rjesenje
	doublereal *x = malloc( n*sizeof( doublereal ) );
	for( i = 0; i < n; i++ )
		x[ i ] = ( i + 1 )*h*cos( (i + 1)*h );	


	//Gauss-Seidel
	omega = 1;
	int k0 = sor_rjesavac( A, x0, b, omega, eps, m, n );
	printf("Za Gauss-Seidelovu metodu potrebno je %d iteracija.\n\n", k0 );
	
	for( i = 0; i < n; i++ )
		printf( "%f\n", x0[ i ] - x[ i ] );
	printf("\n");


	//SOR 
	omega = 1.95;
	int k1 = sor_rjesavac( A, x1, b, omega, eps, m, n );
	printf("Za SOR metodu potrebno je %d iteracija.\n\n", k1 );	
	
	for( i = 0; i < n; i++ )
		printf( "%f\n", x1[ i ] - x[ i ] );
	printf("\n");
	
		
	//konjugirani gradijenti
	int k2 = cg( A, x2, b, n, eps );
	printf("Za metodu konjugiranih gradijenata je potrebno %d iteracija.\n\n", k2 );
	
	for( i = 0; i < n; i++ )
		printf( "%f\n", x2[ i ] - x[ i ] );
	printf("\n");

	
	free( A );
	free( b );
	free( x0 );
	free( x1 );
	free( x2 );
	free( x );
	
    return 0;
}
		
	
