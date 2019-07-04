#include <stdio.h>
#include <stdlib.h>
#include "f2c.h"
#include "fblaswr.h"
#include "clapack.h"
#include "time.h"
#include "math.h"


integer n = 100;
doublereal *A;

void ispis_matrice( doublereal *A, integer n )
{
	int i, j;
	for( i = 0; i < n; i++ )
		for( j = 0; j < n; j++ )
			printf( " %f ", A[ i + j*n ] );
		printf( "\n" );
	printf( "\n\n\n" );	
}

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

/* Table of constant values */

static integer c__1 = 1;
static doublereal c_b7 = -1.;
static doublereal c_b8 = 1.;
static doublereal c_b20 = 0.;

/*  -- Iterative template routine --
*     Univ. of Tennessee and Oak Ridge National Laboratory
*     October 1, 1993
*     Details of this algorithm are described in "Templates for the
*     Solution of Linear Systems: Building Blocks for Iterative
*     Methods", Barrett, Berry, Chan, Demmel, Donato, Dongarra,
*     Eijkhout, Pozo, Romine, and van der Vorst, SIAM Publications,
*     1993. (ftp netlib2.cs.utk.edu; cd linalg; get templates.ps).
*
*  Purpose
*  =======
*
*  GMRES solves the linear system Ax = b using the
*  Generalized Minimal Residual iterative method with preconditioning.
*
*  Convergence test: ( norm( b - A*x ) / norm( b ) ) < TOL.
*  For other measures, see the above reference.
*
*  Arguments
*  =========
*
*  N       (input) INTEGER.
*          On entry, the dimension of the matrix.
*          Unchanged on exit.
*
*  B       (input) DOUBLE PRECISION array, dimension N.
*          On entry, right hand side vector B.
*          Unchanged on exit.
*
*  X       (input/output) DOUBLE PRECISION array, dimension N.
*          On input, the initial guess; on exit, the iterated solution.
*
*  RESTRT  (input) INTEGER
*          Restart parameter, <= N. This parameter controls the amount
*          of memory required for matrix H (see WORK and H).
*
*  WORK    (workspace) DOUBLE PRECISION array, dimension (LDW,RESTRT+4).
*
*  LDW     (input) INTEGER
*          The leading dimension of the array WORK. LDW >= max(1,N).
*
*  H       (workspace) DOUBLE PRECISION array, dimension (LDH,RESTRT+2).
*          This workspace is used for constructing and storing the
*          upper Hessenberg matrix. The two extra columns are used to
*          store the Givens rotation matrices.
*
*  LDH    (input) INTEGER
*          The leading dimension of the array H. LDH >= max(1,RESTRT+1).
*
*  ITER    (input/output) INTEGER
*          On input, the maximum iterations to be performed.
*          On output, actual number of iterations performed.
*
*  RESID   (input/output) DOUBLE PRECISION
*          On input, the allowable convergence measure for
*          norm( b - A*x ) / norm( b ).
*          On output, the final value of this measure.
*
*  MATVEC  (external subroutine)
*          The user must provide a subroutine to perform the
*          matrix-vector product
*
*               y := alpha*A*x + beta*y,
*
*          where alpha and beta are scalars, x and y are vectors,
*          and A is a matrix. Vector x must remain unchanged.
*          The solution is over-written on vector y.
*
*          The call is:
*
*             CALL MATVEC( ALPHA, X, BETA, Y )
*
*          The matrix is passed into the routine in a common block.
*
*  PSOLVE  (external subroutine)
*          The user must provide a subroutine to perform the
*          preconditioner solve routine for the linear system
*
*               M*x = b,
*
*          where x and b are vectors, and M a matrix. Vector b must
*          remain unchanged.
*          The solution is over-written on vector x.
*
*          The call is:
*
*             CALL PSOLVE( X, B )
*
*          The preconditioner is passed into the routine in a common
*          block.
*
*  INFO    (output) INTEGER
*
*          =  0: Successful exit. Iterated approximate solution returned.
*
*          >  0: Convergence to tolerance not achieved. This will be
*                set to the number of iterations performed.
*
*          <  0: Illegal input parameter.
*
*                   -1: matrix dimension N < 0
*                   -2: LDW < N
*                   -3: Maximum number of iterations ITER <= 0.
*                   -4: LDH < RESTRT
*
*  BLAS CALLS:   DAXPY, DCOPY, DDOT, DNRM2, DROT, DROTG, DSCAL
*  ============================================================
*/

int gmres_(n, b, x, restrt, work, ldw, h, ldh, iter, resid, matvec, psolve, 
           info)
   integer *n, *restrt, *ldw, *ldh, *iter, *info;
   doublereal *b, *x, *work, *h, *resid;
   int (*matvec) (), (*psolve) ();
{
    /* System generated locals */
    integer work_dim1, work_offset, h_dim1, h_offset, i__1;
    doublereal d__1;

    /* Local variables */
    extern /* Subroutine */ int drot_();
    static doublereal bnrm2;
    extern doublereal dnrm2_();
    static integer i, k, r, s, v, w, y;
    extern /* Subroutine */ int dscal_(), basis_(), dcopy_(), drotg_();
    static integer maxit;
    static doublereal rnorm, aa, bb;
    static integer cs, av, sn;
    extern /* Subroutine */ int update_();
    static doublereal tol;

    /* Parameter adjustments */
    h_dim1 = *ldh;
    h_offset = h_dim1 + 1;
    h -= h_offset;
    work_dim1 = *ldw;
    work_offset = work_dim1 + 1;
    work -= work_offset;
    --x;
    --b;

    /* Executable Statements */
    *info = 0;

/*     Test the input parameters. */

    if (*n < 0) {
	*info = -1;
    } else if (*ldw < max(1,*n)) {
	*info = -2;
    } else if (*iter <= 0) {
	*info = -3;
    } else if (*ldh < *restrt + 1) {
	*info = -4;
    }
    if (*info != 0) {
	return 0;
    }

    maxit = *iter;
    tol = *resid;

/*     Alias workspace columns. */

    r = 1;
    s = r + 1;
    w = s + 1;
    y = w;
    av = y;
    v = av + 1;

/*     Store the Givens parameters in matrix H. */

    cs = *restrt + 1;
    sn = cs + 1;

/*     Set initial residual (AV is temporary workspace here). */

    dcopy_(n, &b[1], &c__1, &work[av * work_dim1 + 1], &c__1);
    if (dnrm2_(n, &x[1], &c__1) != 0.) {

/*        AV is temporary workspace here. */

	dcopy_(n, &b[1], &c__1, &work[av * work_dim1 + 1], &c__1);
	(*matvec)(&c_b7, &x[1], &c_b8, &work[av * work_dim1 + 1]);
    }
    (*psolve)(&work[r * work_dim1 + 1], &work[av * work_dim1 + 1]);
    bnrm2 = dnrm2_(n, &b[1], &c__1);
    if (bnrm2 == 0.) {
	bnrm2 = 1.;
    }
    if (dnrm2_(n, &work[r * work_dim1 + 1], &c__1) / bnrm2 < tol) {
	goto L70;
    }

    *iter = 0;

L10:

    i = 0;

/*        Construct the first column of V. */

    dcopy_(n, &work[r * work_dim1 + 1], &c__1, &work[v * work_dim1 + 1], &
	    c__1);
    rnorm = dnrm2_(n, &work[v * work_dim1 + 1], &c__1);
    d__1 = 1. / rnorm;
    dscal_(n, &d__1, &work[v * work_dim1 + 1], &c__1);

/*        Initialize S to the elementary vector E1 scaled by RNORM. */

    work[s * work_dim1 + 1] = rnorm;
    i__1 = *n;
    for (k = 2; k <= i__1; ++k) {
	work[k + s * work_dim1] = 0.;
/* L20: */
    }

L30:

    ++i;
    ++(*iter);

    (*matvec)(&c_b8, &work[(v + i - 1) * work_dim1 + 1], &c_b20, &work[av * 
	    work_dim1 + 1]);
    (*psolve)(&work[w * work_dim1 + 1], &work[av * work_dim1 + 1]);

/*           Construct I-th column of H orthnormal to the previous */
/*           I-1 columns. */

    basis_(&i, n, &h[i * h_dim1 + 1], &work[v * work_dim1 + 1], ldw, &work[w *
	     work_dim1 + 1]);

/*           Apply Givens rotations to the I-th column of H. This */
/*           "updating" of the QR factorization effectively reduces */
/*           the Hessenberg matrix to upper triangular form during */
/*           the RESTRT iterations. */

    i__1 = i - 1;
    for (k = 1; k <= i__1; ++k) {
	drot_(&c__1, &h[k + i * h_dim1], ldh, &h[k + 1 + i * h_dim1], ldh, &h[
		k + cs * h_dim1], &h[k + sn * h_dim1]);
/* L40: */
    }

/*           Construct the I-th rotation matrix, and apply it to H so that
 */
/*           H(I+1,I) = 0. */

    aa = h[i + i * h_dim1];
    bb = h[i + 1 + i * h_dim1];
    drotg_(&aa, &bb, &h[i + cs * h_dim1], &h[i + sn * h_dim1]);
    drot_(&c__1, &h[i + i * h_dim1], ldh, &h[i + 1 + i * h_dim1], ldh, &h[i + 
	    cs * h_dim1], &h[i + sn * h_dim1]);

/*           Apply the I-th rotation matrix to [ S(I), S(I+1) ]'. This */
/*           gives an approximation of the residual norm. If less than */
/*           tolerance, update the approximation vector X and quit. */

    drot_(&c__1, &work[i + s * work_dim1], ldw, &work[i + 1 + s * work_dim1], 
	    ldw, &h[i + cs * h_dim1], &h[i + sn * h_dim1]);
    *resid = (d__1 = work[i + 1 + s * work_dim1], abs(d__1)) / bnrm2;
    printf( "%d: %f\n", *iter, *resid );
    if (*resid <= tol) {
	update_(&i, n, &x[1], &h[h_offset], ldh, &work[y * work_dim1 + 1], &
		work[s * work_dim1 + 1], &work[v * work_dim1 + 1], ldw);
	goto L70;
    }
    if (*iter == maxit) {
	goto L50;
    }
    if (i < *restrt) {
	goto L30;
    }

L50:

/*        Compute current solution vector X. */

    update_(restrt, n, &x[1], &h[h_offset], ldh, &work[y * work_dim1 + 1], &
	    work[s * work_dim1 + 1], &work[v * work_dim1 + 1], ldw);

/*        Compute residual vector R, find norm, then check for tolerance. 
*/
/*        (AV is temporary workspace here.) */

    dcopy_(n, &b[1], &c__1, &work[av * work_dim1 + 1], &c__1);
    (*matvec)(&c_b7, &x[1], &c_b8, &work[av * work_dim1 + 1]);
    (*psolve)(&work[r * work_dim1 + 1], &work[av * work_dim1 + 1]);
    work[i + 1 + s * work_dim1] = dnrm2_(n, &work[r * work_dim1 + 1], &c__1);
    *resid = work[i + 1 + s * work_dim1] / bnrm2;
    if (*resid <= tol) {
	goto L70;
    }
    if (*iter == maxit) {
	goto L60;
    }

/*        Restart. */

    goto L10;

L60:

/*     Iteration fails. */

    *info = 1;
    return 0;

L70:

/*     Iteration successful; return. */

    return 0;

/*     End of GMRES */

} /* gmres_ */


/*     =============================================================== */
/* Subroutine */ int update_(i, n, x, h, ldh, y, s, v, ldv)
integer *i, *n;
doublereal *x, *h;
integer *ldh;
doublereal *y, *s, *v;
integer *ldv;
{
    /* System generated locals */
    integer h_dim1, h_offset, v_dim1, v_offset;

    /* Local variables */
    extern /* Subroutine */ int dgemv_(), dcopy_(), dtrsv_();



/*     This routine updates the GMRES iterated solution approximation. */


/*     .. Executable Statements .. */

/*     Solve H*Y = S for upper triangualar H. */

    /* Parameter adjustments */
    v_dim1 = *ldv;
    v_offset = v_dim1 + 1;
    v -= v_offset;
    --s;
    --y;
    h_dim1 = *ldh;
    h_offset = h_dim1 + 1;
    h -= h_offset;
    --x;

    /* Function Body */
    dcopy_(i, &s[1], &c__1, &y[1], &c__1);
    dtrsv_("UPPER", "NOTRANS", "NONUNIT", i, &h[h_offset], ldh, &y[1], &c__1 );

/*     Compute current solution vector X = X + V*Y. */

    dgemv_("NOTRANS", n, i, &c_b8, &v[v_offset], ldv, &y[1], &c__1, &c_b8, &x[
	    1], &c__1 );

    return 0;

} /* update_ */


/*     ========================================================= */
/* Subroutine */ int basis_(i, n, h, v, ldv, w)
integer *i, *n;
doublereal *h, *v;
integer *ldv;
doublereal *w;
{
    /* System generated locals */
    integer v_dim1, v_offset, i__1;
    doublereal d__1;

    /* Local variables */
    extern doublereal ddot_(), dnrm2_();
    static integer k;
    extern /* Subroutine */ int dscal_(), dcopy_(), daxpy_();



/*     Construct the I-th column of the upper Hessenberg matrix H */
/*     using the Gram-Schmidt process on V and W. */


    /* Parameter adjustments */
    --w;
    v_dim1 = *ldv;
    v_offset = v_dim1 + 1;
    v -= v_offset;
    --h;

    /* Function Body */
    i__1 = *i;
    for (k = 1; k <= i__1; ++k) {
	h[k] = ddot_(n, &w[1], &c__1, &v[k * v_dim1 + 1], &c__1);
	d__1 = -h[k];
	daxpy_(n, &d__1, &v[k * v_dim1 + 1], &c__1, &w[1], &c__1);
/* L10: */
    }
    h[*i + 1] = dnrm2_(n, &w[1], &c__1);
    dcopy_(n, &w[1], &c__1, &v[(*i + 1) * v_dim1 + 1], &c__1);
    d__1 = 1. / h[*i + 1];
    dscal_(n, &d__1, &v[(*i + 1) * v_dim1 + 1], &c__1);

    return 0;

}


	

int main( )
{
		int i, j;
		integer m = 1;
		integer N = n*n; 
		
		A = malloc( N*sizeof( doublereal ) );
		
		doublereal *x0 = malloc( n*sizeof( doublereal ) );
		doublereal *g = malloc( n*sizeof( doublereal ) );
		
		for( i = 0; i < n; i++ )
			x0[ i ] = 0.0;
			
		for( i = 0; i < n; i++ )
			g[ i ] = sqrt( 200 - 2*i - 1 );
			
		//generiranje slucajne matrice V	
		doublereal *V = malloc( N*sizeof( doublereal ) );
		integer idist = 1;
		integer iseed[] = { 67, 928, 1555, 13 };
		dlarnv_( &idist, iseed, &N, V );
		
		//ispis_matrice( V, n );
			
		//racunanje QR faktorizacije od V
		doublereal *TAU = malloc( n*sizeof(doublereal) );
		doublereal *work = malloc( N*sizeof(doublereal) );
		integer info;
		dgeqrf_( &n, &n, V, &n, TAU, work, &N, &info );
		
		//ispis_matrice( V, n );
		
		//u V spremimo matricu Q iz QR faktorizacije od V
		dorgqr_( &n, &n, &n, V, &n, TAU, work, &N, &info );
		
		//ispis_matrice( V, n );
		
		//b = V*g
		doublereal *b = malloc( n*sizeof( doublereal ) );
		char trans = 'n';
		doublereal alpha = 1.0;
		doublereal beta = 0.0;
		integer inc = 1;
		dgemv_( &trans, &n, &n, &alpha, V, &n, g, &inc, &beta, b, &inc );	
		printf( "Norma vektora b: %f\n", dnrm2_( &n, b, &inc ) );
		
		
		//B definiramo t.d. joj je prvi stupac b, a ostali stupci su stupci od Q( bez zadnjeg )
		doublereal *B = malloc( N*sizeof( doublereal ) );
		for( i = 0; i < n; i++ )
			for( j = 0; j < n; j++ )
				if( i == 0 )
					B[ j + i*n ] = b[ j ];
				else
					B[ j + i*n ] = V[ j + (i - 1)*n ];
		
		//ispis_matrice( B, n );			
		
		//definiramo matricu M 
		doublereal *M = malloc( N*sizeof( doublereal ) );
		for( i = 0; i < N; i++ )
			M[ i ] = 0.0;
		M[ ( n - 1 )*n ] = 1.0;	
		for( i = 0; i < n; i++ )
			M[ i + 1 + i*n ] = 1.0;
		
		//ispis_matrice( M, n );	
		
		//C=B*M
		alpha = 1.0;
		beta = 0.0;
		doublereal *C = malloc( N*sizeof( doublereal ) );
		dgemm_( &trans, &trans, &n, &n, &n, &alpha, B, &n, M, &n, &beta, C, &n );
		
		//LU faktorizacija od B
		integer *ipiv = malloc( n*sizeof( integer ) );
		dgetrf_( &n, &n, B, &n, ipiv, &info );
		
		//racunanje B^-1 pomocu LU faktorizacije
		doublereal *work1 = malloc( n*sizeof( doublereal ) );
		dgetri_( &n, B, &n, ipiv, work1, &n, &info );
		
		//A = C*B, u B je spremljen inverz od B
		dgemm_( &trans, &trans, &n, &n, &n, &alpha, C, &n, B, &n, &beta, A, &n );
		
		
		//definiramo varijable koje su nam potrebne za poziv GMRES-a
		doublereal *work2 = malloc( n*( n + 4 )*sizeof( doublereal ) );
		doublereal *H = malloc( ( n + 1 )*( n + 2 )*sizeof( doublereal ) );
		integer restrt = n, ldh = n + 1, iter = n;
		doublereal resid = 1e-5;
		
		gmres_( &n, b, x0, &restrt, work2, &n, H, &ldh, &iter, &resid, matvec, psolve, &info );	
		
		//svojstvene vrijednosti
		char jobvl = 'n';
		char jobvr = 'n';
		doublereal *wr = malloc( n*sizeof( doublereal ) );
		doublereal *wi = malloc( n*sizeof( doublereal ) );
		dgeev_( &jobvl, &jobvr, &n, A, &n, wr, wi, NULL, &n, NULL, &n, work, &N, &info );
		
		doublereal *fi = malloc( n*sizeof( doublereal ) );
        printf( "Svojstvene vrijednosti su: " );
        for( i = 0; i < n; i++ )
        {
        	fi[ i ] = atan2( wi[ i ], wr[ i ] );
        	printf( "%f + %fi", cos( 100*fi[ i ] ), sin( 100*fi[ i ] ) );
        	printf( "\n" );
        }

		
		free( x0 );		 
		free( g );
		free( V );
		free( b );
		free( B );
		free( M );
	
		free( TAU );
		free( work );
		free( C );
		free( ipiv );
		free( work1 );
		free( work2 );
		free( H );
		free( wr );
		free( wi );
		free( fi );
		
		return 0;
}

