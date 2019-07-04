#include <stdio.h>
#include <stdlib.h>
#include "f2c.h"
#include "fblaswr.h"
#include "clapack.h"
#include "time.h"
#include "GMRES.c"
#include "sor2_primjena.c"


doublereal A[] = { 11,-20,0,0,0,-2,  -5,41,-3,0,-3,0,  0,-15,7,-1,0,0,  0,0,-4,2,-10,0,  0,-6,0,-1,28,-15,  -1,0,0,0,-15,47 };
integer n = 6;


void matvec( doublereal *alpha, doublereal *x, doublereal *beta, doublereal *y )
{
		char trans = 'n';
		integer inc = 1;
		dgemv_( &trans, &n, &n, alpha, A, &n, x, &inc, beta, y, &inc );
}	

void psolve( doublereal *x, doublereal *b )
{
		integer inc = 1;
		dcopy_( &n, b, &inc, x, &inc );
}



int main()
{

	doublereal b[] = { 500, 0, 0, 0, 0, 0 };
	doublereal x0[] = { 0, 0, 0, 0, 0, 0 };
	doublereal x1[] = { 0, 0, 0, 0, 0, 0 };
	doublereal x2[] = { 0, 0, 0, 0, 0, 0 };
	
	integer m = 1;
	doublereal eps = 1e-8;
	doublereal omega;
	
	//Gauss-Seidel
	omega = 1;
	int k0 = sor_rjesavac( A, x0, b, omega, eps, m, n );
	printf("Za Gauss-Seidelovu metodu potrebno je %d iteracija.\n\n", k0 );


	//SOR 
	omega = 1.35;
	int k1 = sor_rjesavac( A, x1, b, omega, eps, m, n );
	printf("Za SOR metodu potrebno je %d iteracija.\n\n", k1 );
	
	
	//GMRES
	integer restrt = n, ldh = n + 1, iter = n;
	doublereal *work = malloc( n*( n + 4 )*sizeof( doublereal ) );
	doublereal *H = malloc( ( n + 1 )*( n + 2 )*sizeof( doublereal ) );
	integer info;
	gmres_( &n, b, x2, &restrt, work, &n, H, &ldh, &iter, &eps, matvec, psolve, &info );
	
	printf("Rjesenje dobiveno GMRES metodom:\n");
	int i;
	for( i = 0; i < n; i++ )
		printf( " %f ", x2[ i ] );
	printf("\n\n");
	
	printf("Za GMRES metodu je potrebno %d iteracija.\n\n", iter );
	
	return 0;
}
	
	
