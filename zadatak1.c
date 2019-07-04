#include <stdio.h>
#include <stdlib.h>
#include "f2c.h"
#include "fblaswr.h"
#include "clapack.h"
#include "time.h"

void mojmv( int n, double *A, double *x )
{
	int i, j;
	double *rez = malloc( n*sizeof(double) );
	for( i = 0; i < n; i++ )
	{
		rez[ i ] = 0.0;
		for( j = 0; j < n; j++ )
			rez[ i ] += A[ j + i*n ]*x[ j ];
 	}	
}
			


main(int argc, char *argv[])
{

	int n;
	printf("Upisite dimenziju matrice: ");
	scanf( "%d", &n );

	//generiranje slučajnog vektora

	double *x = malloc( n*sizeof(double) );
	int idist = 3;
        int iseed[] = { 3, 78, 469, 3999 };
	dlarnv_( &idist, iseed, &n, x );

	//generiranje slučajne matrice 

	double *A = malloc( n*n*sizeof(double) );
	int N = n*n;
	dlarnv_( &idist, iseed, &N, A );

	//pozvat ćemo dgemv_ i izmjeriti vrijeme

	char trans = 'n';
	double alpha = 1.0;
	double beta = 0.0;
	int incx = 1;
	int incy = 1;
	double *y = malloc( n*sizeof(double) );

	double t1_prije = clock();
	dgemv_( &trans, &n, &n, &alpha, A, &n, x, &incx, &beta, y, &incy );
	double t1_poslije = clock();

	double t1 = t1_poslije - t1_prije;
	double s1 = ((double)t1)/CLOCKS_PER_SEC;
	printf("Vrijeme od dgemv_: %f sekundi.", s1 );

	//ispis njihovog množenja

	int i;
	printf("Rezultat koji je dobio dgemv_: \n");
	for( i = 0; i < n; i++ )
		printf(" %f ", y[ i ] );

	//pozvat ćemo mojmv i izmjeriti vrijeme

	double t2_prije = clock();
	mojmv( n, A, x );
	double t2_poslije = clock();

	double t2 = t2_poslije - t2_prije;
	double s2 = ((double)t2)/CLOCKS_PER_SEC;
	printf("Moje vrijeme je: %f sekundi.", s2 );
}


	



	
	
	
	
