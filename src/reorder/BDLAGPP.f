      SUBROUTINE BDLAGPP( ISIDE, M, N, NB, A, LDA, NITRAF, ITRAF, 
     $                    DTRAF, WORK )
	IMPLICIT NONE
*
*     .. Scalar Arguments ..
      INTEGER            ISIDE, LDA, M, N, NB, NITRAF
*     ..
*     .. Array Arguments ..
      INTEGER            ITRAF( * )
      DOUBLE PRECISION   A( LDA, * ), DTRAF( * ), WORK( * )
*
*  ??? Check whether BDLAAPP and BDLAGPP could be merged.
*
*  Purpose
*  =======
*
*  BDLAAPP computes 
*                            
*          B = Q**T * A       or       B = A * Q,
*
*  where A is an M-by-N matrix and Q is an orthogonal matrix represented
*  by the parameters in the arrays ITRAF and DTRAF as described in
*  BDTGEXC.
*
*  This is an auxiliary routine called by BDTGSEN.
*
*  Arguments
*  =========
*
*  ISIDE   (input) INTEGER
*          Specifies whether Q multiplies A from the left or right as
*          follows: 
*          = 0: compute B = Q**T * A;
*          = 1: compute B = A * Q.
*
*  M       (input) INTEGER
*          The number of rows of A.
*
*  N       (input) INTEGER
*          The number of columns of A.
*
*  NB      (input) INTEGER
*          If ISIDE = 0, the Q is applied block column-wise to the rows
*          of A and NB specifies the maximal width of the block columns.
*          If ISIDE = 1, this variable is not referenced.
*
*  A       (input/output) DOUBLE PRECISION array, dimension (LDA,N)
*          On entry, the matrix A.
*          On exit, A is overwritten by B.
*
*  LDA     (input) INTEGER
*          The leading dimension of the array A. LDA >= max(1,N).
*
*  NITRAF  (input) INTEGER
*          Length of the array ITRAF. NITRAF >= 0.
*
*  ITRAF   (input) INTEGER array, length ITRAF
*          List of parameters for representing the transformation
*          matrix Q, see BDTGEXC.
*
*  DTRAF   (output) DOUBLE PRECISION array, length k, where
*          List of parameters for representing the transformation
*          matrix Q, see BDTGEXC.
*
*  WORK    (workspace) DOUBLE PRECISION array, dimension (4*N)
*
*  =====================================================================
*

*     .. Parameters ..
      DOUBLE PRECISION   ZERO, ONE
      PARAMETER          ( ZERO = 0.0D+0, ONE = 1.0D+0 )
*     ..
*     .. Local Scalars ..
      INTEGER            I, IT, J, K, NNB, PD
*     ..
*     .. External Subroutines ..
      EXTERNAL           DLARFX, DROT
*     .. Intrinsic Functions ..
      INTRINSIC          MIN
*     ..
*     .. Executable Statements ..
*
*     Quick return if possible
*
      IF( M.LE.0 .OR. N.LE.0 )
     $   RETURN
*
      IF( ISIDE.EQ.0 ) THEN
*
*        Apply Q from left.
*
         DO 20 J = 1, N, NB
            PD = 1
            NNB = MIN( NB, N - J + 1 )
C	      PRINT*, 'NITRAF = ', NITRAF
C	      PRINT*, '--------------'
            DO 10 I = 1, NITRAF
	         IT = ITRAF(I)
	         IF( IT.LE.M ) THEN
*
*                 Apply Givens rotation.
*
	            CALL DROT( NNB, A(IT,J), LDA, A(IT+1,J), LDA,
     $                       DTRAF(PD), DTRAF(PD+1) )
	            PD = PD + 2
               ELSE
*
*                 Apply 3x3 or 4x4 matrix.
*
                  IF( IT.LE.2*M ) THEN
	               K = 3
	               IT = IT - M
	            ELSE
	               K = 4
	               IT = IT - 2*M
                  END IF
C	            PRINT*, 'IT = ', IT
C	            PRINT*, 'J = ', J
                  CALL DGEMM( 'T', 'N', K, NNB, K, ONE, DTRAF(PD), K,
     $                        A( IT, J ), LDA, ZERO, WORK, K )
                  CALL DLACPY( 'Full', K, NNB, WORK, K, A( IT, J ),
     $                         LDA )
	            PD = PD + K*K
	         END IF
   10       CONTINUE
   20    CONTINUE
      ELSE
	   PD = 1
	   DO 30 I = 1, NITRAF
	      IT = ITRAF(I)
	      IF( IT.LE.N ) THEN
*
*              Apply Givens rotation.
*
	         CALL DROT( M, A(1,IT), 1, A(1,IT+1), 1, DTRAF(PD),
     $                    DTRAF(PD+1) )
	         PD = PD + 2
            ELSE
*
*              Apply 3x3 or 4x4 matrix.
*
               IF( IT.LE.2*N ) THEN
	            K = 3
	            IT = IT - N
	         ELSE
	            K = 4
	            IT = IT - 2*N
               END IF
               CALL DGEMM( 'N', 'N', M, K, K, ONE, A( 1, IT ), LDA, 
     $                     DTRAF(PD), K, ZERO, WORK, M )
               CALL DLACPY( 'Full', M, K, WORK, M, A( 1, IT ), LDA ) 
	         PD = PD + K*K
	      END IF
   30    CONTINUE
	END IF
*
      RETURN
*
*     End of DLAAPP
*
      END
