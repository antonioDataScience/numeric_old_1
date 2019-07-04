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
    integer n = 7;
    integer N = n*n;
    int i, j;

    doublereal W[] = {0,2,3,4,0,0,0,2,0,0,7,1,0,0,3,0,0,3,0,2,1,4,7,3,0,0,0,0,0,1,0,0,0,7,3,0,0,2,0,7,0,5,0,0,1,0,3,5};

    //sumiramo retke od W i spremimo ih u vektor w
    doublereal w[ n ];
    for( i = 0; i < n; i++ )
    {
        w[ i ] = 0.0;
        for( j = 0; j < n; j++ )
            w[ i ] += W[ i + j*n ];
    }
	
    //D je dijagonalna matrica kojoj su na dijagonali tezine vrhova
    doublereal D[ N ];
    for( i = 0; i < N; i++ )
        D[ i ] = 0.0;
    for( i = 0; i < n; i++ )
        D[ i + i*n ] = w[ i ];

    //L = D - W
    doublereal L[ N ];
    for( i = 0; i < n; i++ )
        for( j = 0; j < n; j++ )
            L[ i + j*n ] = D[ i + j*n ] - W[ i + j*n ];

    //S = D^(-1/2)
    doublereal S[ N ];
    for( i = 0; i < N; i++ )
        S[ i ] = 0.0;
    for( i = 0; i < n; i++ )
        S[ i + i*n ] = 1/sqrt( D[ i + i*n ] );
   
    //LN = S1*L*S1
    doublereal LN[ N ];
	doublereal pom[ N ];
	char trans = 'n';
	doublereal alpha = 1.0;
	doublereal beta = 0.0;
	dgemm_( &trans, &trans, &n, &n, &n, &alpha, S, &n, L, &n, &beta, pom, &n );
	dgemm_( &trans, &trans, &n, &n, &n, &alpha, pom, &n, S, &n, &beta, LN, &n );

    //racunamo drugi svojstveni par od L i LN
    char jobz = 'v';
    char range = 'i';
    char uplo = 'u';
    doublereal vl, vu;
    integer il = 2, iu = 2;
    char cmach = 's';
    doublereal abstol = 2*dlamch_( &cmach );
    integer m = iu - il + 1;
    doublereal w1[ n ], w2[ n ];
    doublereal z1[ n*m ], z2[ n*m ];
    integer ldw = 8*n;
    doublereal work[ ldw ];
    doublereal iwork[ 5*n ];
    integer ifail;
    integer info;

    dsyevx_( &jobz, &range, &uplo, &n, L, &n, &vl, &vu, &il, &iu, &abstol, &m, w1, z1, &n, work, &ldw, iwork, &ifail, &info );
    dsyevx_( &jobz, &range, &uplo, &n, LN, &n, &vl, &vu, &il, &iu, &abstol, &m, w2, z2, &n, work, &ldw, iwork, &ifail, &info );


    //z1 je drugi svojstveni vektor od L
    for( i = 0; i < n; i++ )
        printf( "%f\n", z1[ i ] );
    printf("\n");
    printf("Druga najmanja svojstvena vrijednost od L je %f.\n", w1[ 0 ] );

    //z2 je drugi svojstveni vektor od LN
    for( i = 0; i < n; i++ )
        printf( "%f\n", z2[ i ] );
    printf("\n");
    printf("Druga najmanja svojstvena vrijednost od LN je %f.\n", w2[ 0 ] );

    return 0;
}
