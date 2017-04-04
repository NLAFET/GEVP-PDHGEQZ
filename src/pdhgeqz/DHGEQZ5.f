***********************************************************************
*                                                                     *
*     DHGEQZ5.f:                                                      *
*         Auxillary routine in the package PDHGEQZ.                   *
*                                                                     *
*     Contributors: Bjorn Adlerborn                                   *
*                   Bo Kagstrom                                       *
*                   Daniel Kressner                                   *
*                                                                     *
*     Department of Computing Science and HPC2N, Umea University      *
*     MATHICSE ANCHP, EPF Lausanne                                    *
*                                                                     * 
***********************************************************************
      SUBROUTINE DHGEQZ5( N, A, LDA, B, LDB, IR, NCOL, WSIZE, ESHIFT,
     $   TALPHAR, TALPHAI, TBETA, Q, Z, LDQZ, ATOL )
      IMPLICIT NONE
*     ..
*     .. Scalar Arguments ..
*     ..
      INTEGER             LDA, LDB, N, NCOL, WSIZE, ESHIFT, IR
      DOUBLE PRECISION    ATOL
*     ..
*     .. Array Arguments ..
*     ..
      DOUBLE PRECISION    A( LDA, * ), B( LDB, * ), TALPHAI(*), TALPHAR(
     $   *), TBETA(*), Q( LDQZ, * ), Z( LDQZ, * )
*     
*     Purpose
*     =======
*     
*     DHGEQZ5 chases a bulge down the diagonal of H and T, while keeping
*     orthogonal Q and Z updates. If a bulge collapses due to a vigialent 
*     deflation it is recreated.
*     
*     All the inputs are assumed to be valid without checking
*     
*     Arguments 
*     ========= 
*     
*     N       (input) INTEGER
*     Number of steps the bulge should be chased.
*     
*     A       (input/output) DOUBLE PRECISION array of dimension ( LDA, * )
*     Hessenberg matrix with a bulge to be chased.
*     
*     LDA     (input) INTEGER
*     Leading dimension of the matrix A.
*     
*     B       (input/output) DOUBLE PRECISION array of dimension (LDB,*)
*     Triangular matrix with a bulge to be chased.
*     
*     LDB     (input) INTEGER
*     Leading dimension of the matrix B.
*     
*     IR      (input) INTEGER
*     Row index of where to start chasing the bulge.
*     
*     NCOL    (input) INTEGER
*     Number of steps to not chase the bulge - N is descreased with 
*     this value in loop 10. Either 3, for a full chase, or 2 when 
*     matrix border has been reached and the bulge will be pushed of
*     afterwards.
*     
*     WSIZE   (input) INTEGER
*     Size of the matrices A, B, Q, and Z.
*     
*     ESHIFT  (input) INTEGER
*     Index for where shift information is stored in TALPHAR, TALPHAI,
*     and TBETA.
*     
*     TALPHAR (input) DOUBLE PRECSION array of dimension (WSIZE)
*     TALPHAI (input) DOUBLE PRECSION array of dimension (WSIZE)
*     TBETA   (input) DOUBLE PRECSION array of dimension (WSIZE)a
*     TALPHAR( ESHIFT : ESHIFT + 1 ) contains the real parts and 
*     TALPHAI( ESHIFT : ESHIFT + 1 ) contains the imaginary
*     parts and TBETA( ESHIFT : ESHIFT + 1) contains the scale parts of 
*     the shifts that was used to create the bulge which is about
*     to be to chase.
*     
*     Q       (input/output) DOUBLE PRECISION array of dimension ( LDQZ, * )
*     Q is multiplied, from right, with all transformations
*     applied from left on H and T.          
*     
*     Z       (input/output) DOUBLE PRECISION array of dimension ( LDQZ, * )
*     Z is multiplied, from right, with all transformations
*     applied from right on H and T.   
*     
*     LDQZ    (input) INTEGER 
*     Leading dimension of the matrices Q and Z.
*     
*     ATOL    (input) DOUBLE PRECISION
*     Holds number where safe to set a value of A to zero, used 
*     for vigialent deflations.
*     
*     ==========================================================

*     ..
*     .. Parameters ..
*     ..
      DOUBLE PRECISION    ZERO, ONE
      PARAMETER           ( ZERO = 0.0D+0, ONE = 1.0D+0 )
      DOUBLE PRECISION    SMLNUM, ULP, TFNORM, T1NORM, HFNORM, HTOL,
     $   TTOL
      COMMON /PREC/       SMLNUM, ULP, TFNORM, T1NORM, HFNORM, HTOL,
     $   TTOL

*     ..
*     .. Local Scalars ..
*     ..
      INTEGER             J, JC, JR, INFO, LDQZ
      DOUBLE PRECISION    TAU, TAU1, TEMP, REFSUM, TST1,
     $   TST2, SCL, H12, H21, H11, H22

*     ..
*     .. Local Arrays ..
*     ..
      DOUBLE PRECISION    V( 4 ), V1( 4 )

*     ..
*     .. Externals 
*     ..      
      
*     ..
*     .. Intrinsic Functions ..
*     ..
      INTRINSIC          ABS, DBLE, MAX, MIN, SQRT

      DO 10 J = 2 + IR, N - NCOL + IR
*        
*        Use 3x3 Householder transforms
*        to zero (j-1)st column of A
         V( 1 ) = A( J, J - 1)
         V( 2 ) = A( J+1, J - 1 )
         V( 3 ) = A( J+2, J - 1 )
         CALL DLARFG( 3, V( 1 ), V( 2 ), 1, TAU )   

*        ==== A Bulge may collapse because of vigilant
*        .    deflation or destructive underflow.  In the
*        .    underflow case, try the two-small-subdiagonals
*        .    trick to try to reinflate the bulge.  ====
         IF( A( J+2, J-1 ).NE.ZERO .OR. A( J+2, J ).NE.
     $      ZERO .OR. A( J+2, J+1 ).EQ.ZERO ) THEN 
*           ==== Typical case: not collapsed (yet). ====
            A( J, J - 1) = V( 1 )
            A( J+1, J-1 ) = ZERO
            A( J+2, J-1 ) = ZERO
            
         ELSE
*           ==== Atypical case: collapsed.  Attempt to
*           .    reintroduce ignoring A(J,J-1) and H(J+1,J-1).
*           .    If the fill resulting from the new
*           .    reflector is too large, then abandon it.
*           .    Otherwise, use the new one. ====

     $         
            CALL QZFCOL(3, 2, TALPHAR( ESHIFT ), TALPHAI( ESHIFT ),
     $         TBETA( ESHIFT ), A( J, J ), LDA, B( J , J ), LDB, V1,
     $         INFO ) 

            CALL DLARFG( 3, V1( 1 ), V1( 2 ), 1, TAU1 )

            
            A( J+1, J-1 ) = ZERO
            A( J+2, J-1 ) = ZERO              
            REFSUM = TAU1*( A( J, J - 1 ) + V1( 2 )*
     $         A( J + 1, J - 1 ) )
            
            IF( ABS( A( J + 1, J - 1 ) - REFSUM * V1( 2 ) )+ ABS( REFSUM
     $         * V1( 3 ) ).GT.ULP * ( ABS( A( J - 1, J - 1 ) ) + ABS( A(
     $         J, J ) ) + ABS( A( J + 1, J + 1 ) ) ) ) THEN

     $         
*              ==== Starting a new bulge here would
*              .    create non-negligible fill.  Use
*              .    the old one with trepidation. ====
               A( J, J-1 ) = V(1)
            ELSE
*              ==== Stating a new bulge here would
*              .    create only negligible fill.
*              .    Replace the old reflector with
*              .    the new one. ====
               A( J, J - 1 ) = A( J, J - 1 ) - REFSUM
               TAU = TAU1
               V( 2 ) = V1( 2 )
               V( 3 ) = V1( 3 )         
            END IF
            V( 1 ) = ONE            
            
         END IF

*        Apply transformations from the left.
         DO 20 JC = J, WSIZE
*           On A.
            
            TEMP = TAU * ( A( J, JC ) + V( 2 ) * A( J + 1, JC ) + V( 3 )
     $         * A( J + 2, JC ) )
            A( J, JC ) = A( J, JC ) - TEMP
            A( J + 1, JC ) = A( J + 1, JC ) - TEMP * V( 2 )
            A( J + 2, JC ) = A( J + 2, JC ) - TEMP * V( 3 )
*           On B.
            TEMP = TAU * ( B( J, JC ) + V( 2 ) * B( J + 1, JC ) + V( 3 )
     $         * B( J + 2, JC ) )
            B( J, JC ) = B( J, JC ) - TEMP
            B( J + 1, JC ) = B( J + 1, JC ) - TEMP * V( 2 )
            B( J + 2, JC ) = B( J + 2, JC ) - TEMP * V( 3 )  
 20      CONTINUE

*        Apply left transformations on Q (apply from right)
         DO 30 JR = 1, WSIZE
            TEMP = TAU * ( Q( JR, J ) + V( 2 ) * Q ( JR, J + 1 ) + 
     $         V( 3 ) * Q( JR, J + 2 ) )
            Q( JR, J ) = Q( JR, J ) - TEMP
            Q( JR, J + 1 ) = Q( JR, J + 1 ) - TEMP * V( 2 )
            Q( JR, J + 2 ) = Q( JR, J + 2 ) - TEMP * V( 3 )
 30      CONTINUE
         

*        Zero j-th column of B.
         CALL INVHSE( 3, B( J, J ), LDB, TAU, V, INFO )

*        Apply transformations from the right on A.
         DO 40 JR = 1, MIN(J + 3, WSIZE )
            TEMP = TAU * ( A( JR, J ) + V( 2 ) * A( JR, J + 1 ) + V( 3 )
     $         * A( JR, J + 2 ) )
            A( JR, J ) = A( JR, J ) - TEMP
            A( JR, J + 1 ) = A( JR, J + 1 ) - TEMP * V( 2 )
            A( JR, J + 2 ) = A( JR, J + 2 ) - TEMP * V( 3 )
 40      CONTINUE

*        Apply transformations from the right on B.
         DO 50 JR = 1, MIN(J + 2, WSIZE)
            TEMP = TAU *( B( JR, J ) + V( 2 ) * B( JR, J + 1 ) + V( 3 )
     $         * B( JR, J + 2 ) )
            B( JR, J ) = B( JR, J ) - TEMP
            B( JR, J + 1 ) = B( JR, J + 1 ) - TEMP * V( 2 )
            B( JR, J + 2 ) = B( JR, J + 2 ) - TEMP * V( 3 )
 50      CONTINUE

*        The elements below are annihilated so set to zero.
         B( J + 1, J ) = ZERO
         B( J + 2, J ) = ZERO

*        Apply transformations from the right on Z.
         DO 60 JR = 1, WSIZE
            TEMP = TAU * ( Z( JR, J ) + V( 2 ) * Z( JR, J + 1 ) + V( 3 )
     $         * Z( JR, J + 2 ) )
            Z( JR, J ) = Z( JR, J ) - TEMP
            Z( JR, J + 1 ) = Z( JR, J + 1 ) - TEMP * V( 2 )
            Z( JR, J + 2 ) = Z( JR, J + 2 ) - TEMP * V( 3 )
 60      CONTINUE
 10   CONTINUE 
      
*     Vigilant deflation check
      DO 100 J = 2 + IR, N - NCOL + IR
         IF( A( J, J - 1 ).NE.ZERO ) THEN
            TST1 = ABS( A( J - 1, J - 1 ) ) + ABS( A( J, J ) )
            IF( ABS( A( J, J - 1 ) ).LE.MAX( SMLNUM, ULP * TST1 ) )THEN 
               H12 = MAX( ABS( A( J, J - 1 ) ), ABS( A( J - 1, J ) ) )
               H21 = MIN( ABS( A( J, J - 1 ) ), ABS( A( J - 1, J ) ) )
               H11 = MAX( ABS( A( J, J ) ),
     $            ABS( A( J - 1, J - 1 ) - A( J, J ) ) )
               H22 = MIN( ABS( A( J, J ) ),
     $            ABS( A( J - 1, J - 1 ) - A( J, J ) ) )
               SCL = H11 + H12
               TST2 = H22 * ( H11 / SCL )
*              
               IF( TST2 .EQ. ZERO .OR. H21 * ( H12 / SCL ).LE.
     $            MAX( SMLNUM, ULP * TST2 ) ) THEN
                  A( J, J - 1 ) = ZERO
c                  A( J + 1, J - 1 ) = ZERO
c                  A( J + 2, J - 1 ) = ZERO
               END IF
            END IF         
         END IF
 100  CONTINUE
      RETURN
      END
     

