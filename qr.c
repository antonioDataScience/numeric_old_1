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
	integer n = 21;
	integer N = n*n;
	integer m = 1;
	
	int i, j;
	doublereal x[ n ];
	doublereal y[ n ];
	doublereal A[ N ];
	doublereal b[ n ];
	
	
	for( i = 0; i < n; i++ )
	{
		x[ i ] = -5 + 0.5*i;
		y[ i ] = -5 + 0.5*i;
	}
	
	for( i = 0; i < n; i++ )
		for( j = 0; j < n; j++ )
			A[ i + j*n ] = ( x[ i ]*x[ i ]*y [ j ] - x[ i ]*x[ i ] - y[ j ]*y[ j ] + 175 )/250;
		
	//ispis_matrice( A, n );
	
	b[ 0 ] = 10;
	for( i = 1; i < n; i++ )
		b [ i ] = 1;
	
	//QR faktorizacija s pivotiranjem od A
	integer jpvt[ n ];
	for( i = 0; i < n; i++ )
		jpvt[ i ] = 0.0;
	doublereal TAU[ n ];
	integer lwork = 3*n + 1;
	doublereal work[ lwork ];
	integer info;
	dgeqp3_( &n, &n, A, &n, jpvt, TAU, work, &lwork, &info );
	

	//racunanje ranga
	integer r;
	for( i = 0; i < n; i++ )
		if( A[ i + i*n ] <= 21e-16 )
		{
			r = i;
			break;
		}	
	printf("Rang je: %d\n", r );	
		
	//uzmemo prvih r redaka od A i to je matrica R
	doublereal *R = calloc( r*n, sizeof( doublereal ) );
	for( i = 0; i < r; i++ )
		for( j = i; j < n; j++ )
		    R[ i + j*r ] = A[ i + j*n ];
	
	//generiramo matricu Q
	dorgqr_( &n, &n, &n, A, &n, TAU, work, &lwork, &info );
	doublereal Q[ N ];
	char uplo;
    dlacpy_( &uplo, &n, &n, A, &n, Q, &n );

	//racunamo Q^T*b
	char trans = 't';
	doublereal alpha = 1.0;
	doublereal beta = 0.0;
	doublereal z[ n ];
	dgemv_( &trans, &n, &n, &alpha, Q, &n, b, &m, &beta, z, &m );
	
	//uzmemo samo prva dva retka
	doublereal w[ r ];
	for( i = 0; i < r; i++ )
		w[ i ] = z[ i ];
	
	//LQ faktorizacija od R
	doublereal TAU1[ r ];
	integer lwork1 = r*n;
	doublereal work1[ lwork1 ];
	dgelqf_( &r, &n, R, &r, TAU1, work1, &lwork1, &info );

	//L*x = w, tj. x = L^-1*w, L je donji trokut od R, rjesenje je spremljeno u w
	char side = 'l';
	uplo = 'l';
	trans = 'n';
	char diag = 'n';
	dtrsm_( &side, &uplo, &trans, &diag, &r, &m, &alpha, R, &r, w, &r );
	
	//u matricu R spremimo Q iz LQ faktorizacije od R, to je zapravo matrica Z iz zadatka
	dorglq_( &r, &n, &r, R, &r, TAU1, work1, &lwork1, &info );
		    
	//mnozimo R^T*w 
	doublereal rj[ n ];
	trans = 't';
	dgemv_( &trans, &r, &n, &alpha, R, &r, w, &m, &beta, rj, &m );
	
	//permutacija rjesenja
	doublereal konacno[ n ];	
	for( i = 0; i < n; i++ )
		for( j = 0; j < n; j++ )
			if( jpvt[ j ] == i + 1 )
				konacno[ i ] = rj[ j ];
				
	for( i = 0; i < n; i++ )
		printf("%f\n", konacno[ i ] );		
	
	return 0;
}		
			
	
