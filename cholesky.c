#include <stdio.h>
#include <stdlib.h>
#include "f2c.h"
#include "fblaswr.h"
#include "clapack.h"
#include "time.h"
#include "math.h"


void ic (integer n, doublereal *a)
{

    int i,j,k;
    for (i=0; i<n; i++){
        for (k=0; k<i; k++)
            a[i+i*n]-=a[k+i*n]*a[k+i*n];
        a[i+i*n]=sqrt(a[i+i*n]);
        for (j=i+1; j<n; j++){
            if (a[i+j*n]!=0){
                for (k=0; k<i; k++)
                    a[i+j*n]-=a[k+i*n]*a[k+j*n];
                a[i+j*n]/=a[i+i*n];
            }
            
        }
    }

}

int pcg( doublereal *A, doublereal *x, doublereal *b, integer n, doublereal eps )
{
		
		int i, j, brojac = 0;
		integer m = 1;
		doublereal kriterij = 1.0;
		char uplo, trans, diag;
		
		doublereal *R = malloc( n*n*sizeof(doublereal) );
		doublereal *rez = malloc( n*sizeof(doublereal) );
		doublereal *p = malloc( n*sizeof(doublereal) );
		doublereal *s = malloc( n*sizeof(doublereal) );
		doublereal *d = malloc( n*sizeof(doublereal) );
		
		//kopiramo A u R
		dlacpy_( &uplo, &n, &n, A, &n, R, &n );
				
		//R je nekompletni faktor Choleskog
		ic( n, R );
		
		//kopiramo b u rez 
		dlacpy_( &uplo, &m, &n, b, &m, rez, &m );
		
		//rez = b - A*x
		trans = 'n';
		doublereal alpha = -1.0;
		doublereal beta = 1.0;
		integer inc = 1;
		dgemv_( &trans, &n, &n, &alpha, A, &n, x, &inc, &beta, rez, &inc );
		
		//kopiramo rez u s, s je pomocni vektor
		dlacpy_( &uplo, &m, &n,rez, &m, s, &m );

		//R^T*s = rez
		uplo = 'u';
		trans = 't';
		diag = 'n';
		dtrsv_( &uplo, &trans, &diag, &n, R, &n, s, &inc );
		
		//R*p = s
		trans = 'n';
		dtrsv_( &uplo, &trans, &diag, &n, R, &n, s, &inc );
		
		//kopiramo s u p, p je rjesenje sustava M*p = rez
		dlacpy_( &uplo, &m, &n, s, &m, p, &m );		

		//kopiramo p u d
		dlacpy_( &uplo, &m, &n, p, &m, d, &m );	
		
		//norma od b
		doublereal norma_b = dnrm2_( &n, b, &inc );
		
		while( kriterij > eps )
		{
				brojac++;
				
				doublereal pom[ n ];
				doublereal suma1 = 0.0, suma2 = 0.0, suma3 = 0.0, norma_r;
				
				suma1 = ddot_( &n, p, &inc, rez, &inc );
				
				//pom = A*d
				alpha = 1.0;
				beta = 0.0;
				dgemv_( &trans, &n, &n, &alpha, A, &n, d, &inc, &beta, pom, &inc );
		
				suma2 = ddot_( &n, pom, &inc, d, &inc );
		
				doublereal a = suma1/suma2; 
				
				//x = a*d + x
				daxpy_( &n, &a, d, &inc, x, &inc );
		
				//rez = -a*pom + rez
				a = -a;
				daxpy_( &n, &a, pom, &inc, rez, &inc );
				
				//kopiramo rez u s, s je pomocni vektor
				dlacpy_( &uplo, &m, &n,rez, &m, s, &m );

				//R^T*s = rez
				uplo = 'u';
				trans = 't';
				diag = 'n';
				dtrsv_( &uplo, &trans, &diag, &n, R, &n, s, &inc );
		
				//R*p = s
				trans = 'n';
				dtrsv_( &uplo, &trans, &diag, &n, R, &n, s, &inc );
		
				//kopiramo s u p, p je rjesenje sustava M*p = rez
				dlacpy_( &uplo, &m, &n, s, &m, p, &m );		

				suma3 = ddot_( &n, p, &inc, rez, &inc ); 
					
				doublereal b = suma3/suma1;
				
				for( i = 0; i <n; i++ )
					d[ i ] = p[ i ] + b*d[ i ];	
				
				//norma reziduala
				norma_r = dnrm2_( &n, rez, &inc );
				
				kriterij = norma_r/norma_b;
				printf( "Relativna norma reziduala u %d. koraku je %.16f\n", brojac, kriterij );

		}
		free( R );
		free( rez );
		free( p );
		free( s );
		free( d );
		return brojac;	
}
						
		
int main()
{
		int i;
		integer n = 100;
		doublereal eps = 1e-8;
		
		FILE *f;
		doublereal *A = malloc( n*n*sizeof(doublereal) );
		
		f = fopen( "stieltjes_matr.txt", "r" );
		for( i = 0; i < n*n; i++ )
			fscanf( f, "%lf", A + i );
		fclose(f);
		
		doublereal *b = malloc( n*sizeof(doublereal) );
		doublereal *x = malloc( n*sizeof(doublereal) );
		doublereal *x0 = malloc( n*sizeof(doublereal) );
		
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
		
		int broj_iteracija = pcg( A, x0, b, n, eps );
		printf( "Potrebno je %d iteracija.\n", broj_iteracija );

		free( A );
		free( b );
		free( x );
		free( x0 );
		
		return 0;
		
}
