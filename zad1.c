#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include "f2c.h"
#include "fblaswr.h"
#include "clapack.h"

	
void mojmv_( doublereal *A, doublereal *x, int n )
{
		doublereal *rez = malloc( n*sizeof( doublereal ) );
		int i, j;
		
		for( i = 0; i < n; i++ )
		{
			rez[ i ] = 0.0;
			for( j = 0; j <n; j++ )
				rez[ i ] += A[ i + j*n ] * x[ j ];
		}
	
		/*for( i = 0; i < n; i++ )
			printf(" %f ", rez[ i ] );
		printf("\n");*/	
}

main(integer argc, char *argv[])
{	
	integer n, N, i ,j;
	scanf("%d", &n );
	N = n*n;
	
	doublereal *x = malloc( n*sizeof( doublereal ));
	doublereal *A = malloc( N*sizeof( doublereal ));
	
	//generiranje slucajnog vektora((double)t1)/CLOCKS_PER_SEC 
	integer idist = 3;
	integer iseed[] = { 3, 78, 1016, 4095 };
	dlarnv_( &idist, iseed, &n, x );

	//generiranje matrice
	dlarnv_( &idist, iseed, &N, A );
    /*for( i = 0; i < n; i++ )
	{
		for( j = 0; j < n; j++ )
			printf(" %f ", A[ i + j*n ] );
		printf("\n");
	}*/
	
	//pozivanje dgemv_ i mjerenje vremena
	char trans = 'n';
	doublereal a = 1.0;
    doublereal b = 0.0;
	integer incx = 1;
	integer incy = 1;
	doublereal *y = malloc( n*sizeof( doublereal ));

   	int t1_prije = clock();
	printf("%d\n", t1_prije);
	
	dgemv_( &trans, &n, &n, &a, A, &n, x, &incx, &b, y, &incy );
	/*for( i = 0; i < n; i++ )
		printf(" %f ", y[ i ] );
	printf("\n");*/	
			 
	int t1_poslije = clock();
	printf("%d\n", t1_poslije);
	
	int t1 = t1_poslije - t1_prije;
	printf( "Vrijeme od dgemv_ u sekundama : %.10f \n", ((double)t1)/CLOCKS_PER_SEC );

	int t2_prije = clock();
	printf("%d\n", t2_prije);
	
	mojmv_( A, x, n );
	
	int t2_poslije = clock();
	printf("%d\n", t2_poslije);
	
	int t2 = t2_poslije - t2_prije;
	printf( "Vrijeme od mojmv_ u sekundama : %.10f \n", ((double)t2)/CLOCKS_PER_SEC );
					
}
