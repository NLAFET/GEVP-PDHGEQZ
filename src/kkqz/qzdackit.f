      SUBROUTINE QZDACKIT( LSCHUR, LCOMPQ, LCOMPZ, LEXC, N, ILO, IHI,
     $   ESHIFT, A, LDA, B, LDB, Q, LDQ, Z, LDZ,
     $   WORK, LWORK, NB, ANORM, BNORM, NSHF )
      IMPLICIT NONE
*     
*     This is an isolated piece of code of the routine QZDACK only
*     containing one step of the QZ iteration.
*     The extra parameter NSHFT returns the number of applied shifts.
*     
*     .. Scalar Arguments ..
      LOGICAL            LSCHUR, LCOMPQ, LCOMPZ, LEXC
      INTEGER            IHI, ILO, LDA, LDB, LDQ, LDZ, LWORK, N
      INTEGER            NB, NSHF
      DOUBLE PRECISION   ANORM, BNORM, ESHIFT
*     ..
*     .. Array Arguments ..
      DOUBLE PRECISION   A( LDA, * ), 
     $   B( LDB, * ), Q( LDQ, * ), WORK( * ),
     $   Z( LDZ, * )
      DOUBLE PRECISION   D1, D2, D3, D4, D5

*     .. Parameters ..
*     $                     SAFETY = 1.0E+0 )
      DOUBLE PRECISION   HALF, ZERO, ONE, SAFETY
      PARAMETER          ( HALF = 0.5D+0, ZERO = 0.0D+0, ONE = 1.0D+0,
     $   SAFETY = 1.0D+2 )
C      INCLUDE           'defcfg.f'
C      INCLUDE           'infcfg.f'
*     ..
*     .. Local Scalars ..
      LOGICAL            ILPIVT, NDONE
      INTEGER            J, JC, JR, NNB, RT, UPDLO, UPDHI, NCOLS, 
     $   JJ, NREF, NBL
      DOUBLE PRECISION   AD11, AD11L, AD12, AD12L, AD21, AD21L, AD22, 
     $   AD22L, AD32L, ASCALE, BSCALE, C, S, S1, S2, SAFMAX, SAFMIN, 
     $   SCALE, T, TAU, TEMP, TEMP2, TEMPR, U1, U12, U12L, U2, VS, 
     $   W11, W12, W21, W22, WI, WR, WR2
*     ..
*     .. Local Arrays ..
      DOUBLE PRECISION   V( 3 )
*     ..
*     .. External Functions ..
      LOGICAL            LSAME
      DOUBLE PRECISION   DLAMCH, DLANHS, DLANTR, DLAPY2, DLAPY3
      EXTERNAL           LSAME, DLAMCH, DLANHS, DLANTR, DLAPY2, DLAPY3
      INTEGER            ILAENV
      EXTERNAL           ILAENV
*     ..
*     .. External Subroutines ..
      EXTERNAL           DLAG2, DLARFG, DLARTG, DLASET, DLASV2, DROT,
     $   XERBLA
*     ..
*     .. Intrinsic Functions ..
      INTRINSIC          ABS, DBLE, MAX, MIN, SQRT
*     ..
*     .. Executable Statements ..
*     
      RT = 3*N+1
      SAFMIN = DLAMCH( 'S' )
      SAFMAX = ONE / SAFMIN
      ASCALE = ONE / MAX( SAFMIN, ANORM )
      BSCALE = ONE / MAX( SAFMIN, BNORM )
*     
      IF ( LSCHUR ) THEN
         UPDLO = 1
         UPDHI = N
      ELSE
         UPDLO = ILO
         UPDHI = IHI
      END IF
      IF ( LEXC ) THEN
*        
*        Exceptional shift.  Chosen for no particularly good reason.
*        (Single shift only.)
*        
         IF( ( DBLE( 30*N )*SAFMIN )*ABS( A( IHI-1, IHI ) ).LT.
     $      ABS( B( IHI-1, IHI-1 ) ) ) THEN
            ESHIFT = ESHIFT + A( IHI-1, IHI ) /
     $         B( IHI-1, IHI-1 )
         ELSE
            ESHIFT = ESHIFT + ONE / ( SAFMIN*DBLE( 30*N ) )
         END IF
         S1 = ONE
         WR = ESHIFT 
      ELSE
*        
*        Shifts based on the generalized eigenvalues of the
*        bottom-right 2x2 block of A and B. The first eigenvalue
*        returned by DLAG2 is the Wilkinson shift (AEP p.512),
*        
         CALL DLAG2( A( IHI-1, IHI-1 ), LDA, B( IHI-1, IHI-1 ), LDB,
     $      SAFMIN*SAFETY, S1, S2, WR, WR2, WI )
*        
         TEMP = MAX( S1, SAFMIN*MAX( ONE, ABS( WR ), ABS( WI ) ) )
         IF( WI.NE.ZERO )
     $      GO TO 200
*        
      END IF


*     
*     Fiddle with shift to avoid overflow
*     
      TEMP = MIN( ASCALE, ONE )*( HALF*SAFMAX )
      IF( S1.GT.TEMP ) THEN
         SCALE = TEMP / S1
      ELSE
         SCALE = ONE
      END IF
*     
      TEMP = MIN( BSCALE, ONE )*( HALF*SAFMAX )
      IF( ABS( WR ).GT.TEMP )
     $   SCALE = MIN( SCALE, TEMP / ABS( WR ) )
      S1 = SCALE*S1
      WR = SCALE*WR
*     

*     
*     Do an implicit single-shift QZ sweep.
*     
*     Initial Q
*     
      TEMP = S1*A( ILO, ILO ) - WR*B( ILO, ILO )
      TEMP2 = S1*A( ILO+1, ILO )
      CALL DLARTG( TEMP, TEMP2, C, S, TEMPR )
*     
*     Sweep
*     
      J = ILO
*     
      DO 140 JC = J, UPDHI
         TEMP = C*A( J, JC ) + S*A( J+1, JC )
         A( J+1, JC ) = -S*A( J, JC ) + C*A( J+1, JC )
         A( J, JC ) = TEMP
         TEMP2 = C*B( J, JC ) + S*B( J+1, JC )
         B( J+1, JC ) = -S*B( J, JC ) + C*B( J+1, JC )
         B( J, JC ) = TEMP2
 140  CONTINUE
      IF( LCOMPQ ) THEN
         DO 150 JR = 1, N
            TEMP = C*Q( JR, J ) + S*Q( JR, J+1 )
            Q( JR, J+1 ) = -S*Q( JR, J ) + C*Q( JR, J+1 )
            Q( JR, J ) = TEMP
 150     CONTINUE
      ENDIF
*     
      TEMP = B( J+1, J+1 )
      CALL DLARTG( TEMP, B( J+1, J ), C, S, B( J+1, J+1 ) )
      B( J+1, J ) = ZERO
*     
      DO 160 JR = UPDLO, MIN( J+2, IHI )
         TEMP = C*A( JR, J+1 ) + S*A( JR, J )
         A( JR, J ) = -S*A( JR, J+1 ) + C*A( JR, J )
         A( JR, J+1 ) = TEMP
 160  CONTINUE
      DO 170 JR = UPDLO, J
         TEMP = C*B( JR, J+1 ) + S*B( JR, J )
         B( JR, J ) = -S*B( JR, J+1 ) + C*B( JR, J )
         B( JR, J+1 ) = TEMP
 170  CONTINUE
      IF( LCOMPZ ) THEN
         DO 180 JR = 1, N
            TEMP = C*Z( JR, J+1 ) + S*Z( JR, J )
            Z( JR, J ) = -S*Z( JR, J+1 ) + C*Z( JR, J )
            Z( JR, J+1 ) = TEMP
 180     CONTINUE
      END IF
*     
*     New Blocked SIngle Code Starts Here
*     
      JJ = ILO
      NDONE = .TRUE.
      NREF = 1
 190  CONTINUE
      NNB = MIN(NB, IHI-JJ+1)
      NCOLS = 2
      NBL = NNB - 3
*     
*     Move bulge and accumulate transformations Q and Z
*     
*     Last block?
      IF ((JJ + NNB -1) .EQ. IHI) THEN
         NCOLS = 1
         NDONE = .FALSE.
         NBL = UPDHI - (JJ + 2)
      ENDIF
*     
*     Transformations from the left
*     
      IF (JJ-ILO .GT. 0) THEN
         CALL MYROT('Rows', A(1,JJ+3),
     $      LDA, ILO+1, 1, 1, JJ-ILO,
     $      1, NBL, WORK)
         CALL MYROT('Rows', B(1,JJ+3),
     $      LDB, ILO+1, 1, 1, JJ-ILO,
     $      1, NBL, WORK)
      ENDIF
      CALL GIV(NNB, A(JJ, JJ), LDA, B(JJ, JJ), LDB,
     $   WORK(NREF), WORK(RT), NCOLS)
*     
*     Transformations from the left (last iteration)
*     
      IF (.NOT. NDONE) THEN
         CALL MYROT('Rows', A(JJ,JJ+NNB), LDA,
     $      2, 1, 1, NNB-(NCOLS+1), 1, UPDHI-JJ-NNB+1,
     $      WORK(NREF) )
         CALL MYROT('Rows', B(JJ,JJ+NNB), LDB,
     $      2, 2, 1, NNB-(NCOLS+1), 1, UPDHI-JJ-NNB+1,
     $      WORK(NREF))
      ENDIF
      
      IF (LCOMPQ) THEN
         CALL MYROT('X', Q(1,JJ), LDQ, 
     $      1, 2, 1, NNB-(NCOLS+1),
     $      1, N, WORK(NREF) )
      ENDIF
*     
*     Transformations from the right
*     
      
      IF (JJ-UPDLO .GT. 0) THEN
         CALL MYROT('Cols', A(UPDLO,JJ), LDA,
     $      1, 2, 1, NNB-(NCOLS+1), 1, JJ-UPDLO,
     $      WORK(RT))
         CALL MYROT('Cols', B(UPDLO,JJ), LDB,
     $      1, 2, 1, NNB-(NCOLS+1), 1, JJ-UPDLO,
     $      WORK(RT))
      ENDIF
      
      IF (LCOMPZ) THEN
         CALL MYROT('Cols', Z(1,JJ), LDZ, 1, 2, 1, 
     $      NNB-(NCOLS+1), 1, N, WORK(RT))
      ENDIF
      NREF = NREF + (NNB-(NCOLS+1)) * 2
      JJ = JJ + NB - 3
      IF (NDONE)  GOTO 190
*     
*     New Blocked Code Ends Here
*     
      NSHF = 1
      GO TO 350
*     
 200  CONTINUE

*     
*     Usual case: 3x3 or larger block, using Francis implicit
*     double-shift
*     
*     2
*     Eigenvalue equation is  w  - c w + d = 0,
*     
*     -1 2        -1
*     so compute 1st column of  (A B  )  - c A B   + d
*     using the formula in QZIT (from EISPACK)
*     
*     We assume that the block is at least 3x3
*     
      AD11 = ( ASCALE*A( IHI-1, IHI-1 ) ) /
     $   ( BSCALE*B( IHI-1, IHI-1 ) )
      AD21 = ( ASCALE*A( IHI, IHI-1 ) ) /
     $   ( BSCALE*B( IHI-1, IHI-1 ) )
      AD12 = ( ASCALE*A( IHI-1, IHI ) ) /
     $   ( BSCALE*B( IHI, IHI ) )
      AD22 = ( ASCALE*A( IHI, IHI ) ) /
     $   ( BSCALE*B( IHI, IHI ) )
      U12 = B( IHI-1, IHI ) / B( IHI, IHI )
      AD11L = ( ASCALE*A( ILO, ILO ) ) /
     $   ( BSCALE*B( ILO, ILO ) )
      AD21L = ( ASCALE*A( ILO+1, ILO ) ) /
     $   ( BSCALE*B( ILO, ILO ) )
      AD12L = ( ASCALE*A( ILO, ILO+1 ) ) /
     $   ( BSCALE*B( ILO+1, ILO+1 ) )
      AD22L = ( ASCALE*A( ILO+1, ILO+1 ) ) /
     $   ( BSCALE*B( ILO+1, ILO+1 ) )
      AD32L = ( ASCALE*A( ILO+2, ILO+1 ) ) /
     $   ( BSCALE*B( ILO+1, ILO+1 ) )
      U12L = B( ILO, ILO+1 ) / B( ILO+1, ILO+1 )
*     
      V( 1 ) = ( AD11-AD11L )*( AD22-AD11L ) - AD12*AD21 +
     $   AD21*U12*AD11L + ( AD12L-AD11L*U12L )*AD21L
      V( 2 ) = ( ( AD22L-AD11L )-AD21L*U12L-( AD11-AD11L )-
     $   ( AD22-AD11L )+AD21*U12 )*AD21L
      V( 3 ) = AD32L*AD21L
*     
      CALL DLARFG( 3, V( 1 ), V( 2 ), 1, TAU )
      V( 1 ) = ONE
*     
*     Sweep XXX
*     
      J = ILO

*     
*     All but last elements: use 3x3 Householder transforms.
*     
*     Zero (j-1)st column of A
*     
      DO 230 JC = J, UPDHI
         TEMP = TAU*( A( J, JC )+V( 2 )*A( J+1, JC )+V( 3 )*
     $      A( J+2, JC ) )
         A( J, JC ) = A( J, JC ) - TEMP
         A( J+1, JC ) = A( J+1, JC ) - TEMP*V( 2 )
         A( J+2, JC ) = A( J+2, JC ) - TEMP*V( 3 )
         TEMP2 = TAU*( B( J, JC )+V( 2 )*B( J+1, JC )+V( 3 )*
     $      B( J+2, JC ) )
         B( J, JC ) = B( J, JC ) - TEMP2
         B( J+1, JC ) = B( J+1, JC ) - TEMP2*V( 2 )
         B( J+2, JC ) = B( J+2, JC ) - TEMP2*V( 3 )
 230  CONTINUE
      IF( LCOMPQ ) THEN
         DO 240 JR = 1, N
            TEMP = TAU*( Q( JR, J )+V( 2 )*Q( JR, J+1 )+V( 3 )*
     $         Q( JR, J+2 ) )
            Q( JR, J ) = Q( JR, J ) - TEMP
            Q( JR, J+1 ) = Q( JR, J+1 ) - TEMP*V( 2 )
            Q( JR, J+2 ) = Q( JR, J+2 ) - TEMP*V( 3 )
 240     CONTINUE
      END IF
*     
*     Zero j-th column of B (see DLAGBC for details)
*     
*     Swap rows to pivot
*     
      ILPIVT = .FALSE.
      TEMP = MAX( ABS( B( J+1, J+1 ) ), ABS( B( J+1, J+2 ) ) )
      TEMP2 = MAX( ABS( B( J+2, J+1 ) ), ABS( B( J+2, J+2 ) ) )
      IF( MAX( TEMP, TEMP2 ).LT.SAFMIN ) THEN
         SCALE = ZERO
         U1 = ONE
         U2 = ZERO
         GO TO 250
      ELSE IF( TEMP.GE.TEMP2 ) THEN
         W11 = B( J+1, J+1 )
         W21 = B( J+2, J+1 )
         W12 = B( J+1, J+2 )
         W22 = B( J+2, J+2 )
         U1 = B( J+1, J )
         U2 = B( J+2, J )
      ELSE
         W21 = B( J+1, J+1 )
         W11 = B( J+2, J+1 )
         W22 = B( J+1, J+2 )
         W12 = B( J+2, J+2 )
         U2 = B( J+1, J )
         U1 = B( J+2, J )
      END IF
*     
*     Swap columns if nec.
*     
      IF( ABS( W12 ).GT.ABS( W11 ) ) THEN
         ILPIVT = .TRUE.
         TEMP = W12
         TEMP2 = W22
         W12 = W11
         W22 = W21
         W11 = TEMP
         W21 = TEMP2
      END IF
*     
*     LU-factor
*     
      TEMP = W21 / W11
      U2 = U2 - TEMP*U1
      W22 = W22 - TEMP*W12
      W21 = ZERO
*     
*     Compute SCALE
*     
      SCALE = ONE
      IF( ABS( W22 ).LT.SAFMIN ) THEN
         SCALE = ZERO
         U2 = ONE
         U1 = -W12 / W11
         GO TO 250
      END IF
      IF( ABS( W22 ).LT.ABS( U2 ) )
     $   SCALE = ABS( W22 / U2 )
      IF( ABS( W11 ).LT.ABS( U1 ) )
     $   SCALE = MIN( SCALE, ABS( W11 / U1 ) )
*     
*     Solve
*     
      U2 = ( SCALE*U2 ) / W22
      U1 = ( SCALE*U1-W12*U2 ) / W11
*     
 250  CONTINUE
      IF( ILPIVT ) THEN
         TEMP = U2
         U2 = U1
         U1 = TEMP
      END IF
*     
*     Compute Householder Vector
*     
      T = SQRT( SCALE**2+U1**2+U2**2 )
      TAU = ONE + SCALE / T
      VS = -ONE / ( SCALE+T )
      V( 1 ) = ONE
      V( 2 ) = VS*U1
      V( 3 ) = VS*U2
*     
*     Apply transformations from the right.
*     
      DO 260 JR = UPDLO, MIN( J+3, IHI )
         TEMP = TAU*( A( JR, J )+V( 2 )*A( JR, J+1 )+V( 3 )*
     $      A( JR, J+2 ) )
         A( JR, J ) = A( JR, J ) - TEMP
         A( JR, J+1 ) = A( JR, J+1 ) - TEMP*V( 2 )
         A( JR, J+2 ) = A( JR, J+2 ) - TEMP*V( 3 )
 260  CONTINUE
      DO 270 JR = UPDLO, J + 2
         TEMP = TAU*( B( JR, J )+V( 2 )*B( JR, J+1 )+V( 3 )*
     $      B( JR, J+2 ) )
         B( JR, J ) = B( JR, J ) - TEMP
         B( JR, J+1 ) = B( JR, J+1 ) - TEMP*V( 2 )
         B( JR, J+2 ) = B( JR, J+2 ) - TEMP*V( 3 )
 270  CONTINUE
      IF( LCOMPZ ) THEN
         DO 280 JR = 1, N
            TEMP = TAU*( Z( JR, J )+V( 2 )*Z( JR, J+1 )+V( 3 )*
     $         Z( JR, J+2 ) )
            Z( JR, J ) = Z( JR, J ) - TEMP
            Z( JR, J+1 ) = Z( JR, J+1 ) - TEMP*V( 2 )
            Z( JR, J+2 ) = Z( JR, J+2 ) - TEMP*V( 3 )
 280     CONTINUE
      END IF
      B( J+1, J ) = ZERO
      B( J+2, J ) = ZERO
*     
*     New Blocked Code Starts Here
*     
      JJ = ILO
      NDONE = .TRUE.
      NREF = 1
 1000 CONTINUE
      NNB = MIN(NB, IHI-JJ+1)
      NCOLS = 3
      NBL = NNB - 4
*     
*     Move bulge and accumulate transformations Q and Z
*     
*     
*     Last block? 
      IF ((JJ + NNB -1) .EQ. IHI) THEN
         NCOLS = 2
         NDONE = .FALSE.
         NBL = UPDHI - (JJ + 3) 
      ENDIF                             
*     
*     Transformations from the left
*     
      IF (JJ-ILO .GT. 0) THEN
         CALL DLAREF('Rows', A(1,JJ+4), 
     $      LDA, .FALSE., Q,
     $      LDQ, .TRUE., ILO+1, 1, 1, JJ-ILO, 
     $      1, NBL, 1, 1, WORK, 
     $      D1, D2, D3, D4, D5)
         CALL DLAREF('Rows', B(1,JJ+4), 
     $      LDB, .FALSE., Q,
     $      LDQ, .TRUE., ILO+1, 1, 1, JJ-ILO, 
     $      1, NBL, 1, 1, WORK, 
     $      D1, D2, D3, D4, D5)
      ENDIF

      CALL QZ2(NNB, A(JJ, JJ), LDA, B(JJ, JJ), LDB, 
     $   WORK(NREF), WORK(RT), NCOLS)
      IF (.NOT. NDONE) THEN
         CALL DLAREF('Rows', A(JJ,JJ+NNB), LDA, .FALSE., Q,
     $      LDQ, .TRUE., 2, 1, 1, NNB-(NCOLS+1),
     $      1, UPDHI-JJ-NNB+1,
     $      1, 1, WORK(NREF), D1, D2, D3, D4, D5)
         CALL DLAREF('Rows', B(JJ,JJ+NNB), LDB, .FALSE., Q,
     $      LDQ, .TRUE., 2, 2, 1, NNB-(NCOLS+1),
     $      1, UPDHI-JJ-NNB+1,
     $      1, 1, WORK(NREF), D1, D2, D3, D4, D5)
      ENDIF

      IF (LCOMPQ) THEN
         CALL DLAREF('Cols', Q(1,JJ), LDQ, .FALSE., Q,
     $      LDQ, .TRUE., 1, 2, 1, NNB-(NCOLS+1), 
     $      1, N, 1, 1, WORK(NREF), 
     $      D1, D2, D3, D4, D5)
      ENDIF
*     
*     Transformations from the right
*     

      IF (JJ-UPDLO .GT. 0) THEN
         CALL DLAREF('Cols', A(UPDLO,JJ), LDA, .FALSE., Q,
     $      LDQ, .TRUE., 1, 2, 1, NNB-(NCOLS+1), 1, 
     $      JJ-UPDLO,
     $      1, 1, WORK(RT), D1, D2, D3, D4, D5)
         CALL DLAREF('Cols', B(UPDLO,JJ), LDB, .FALSE., Q,
     $      LDQ, .TRUE., 1, 2, 1, NNB-(NCOLS+1), 1, JJ-
     $      UPDLO,
     $      1, 1, WORK(RT), D1, D2, D3, D4, D5)
      ENDIF

      IF (LCOMPZ) THEN
         CALL DLAREF('Cols', Z(1,JJ), LDZ, .FALSE., Z,
     $      LDQ, .TRUE., 1, 2, 1, NNB-(NCOLS+1), 
     $      1, N,
     $      1, 1, WORK(RT), D1, D2, D3, D4, D5)
      ENDIF
      NREF = NREF + (NNB-(NCOLS+1)) * 3
      JJ = JJ + NB - 4
      IF (NDONE)  GOTO 1000
*     
*     New Blocked Code Ends Here
*     
*     
*     Last elements: Use Givens rotations
*     
*     Rotations from the left
*     
      J = IHI - 1
      TEMP = A( J, J-1 )
      CALL DLARTG( TEMP, A( J+1, J-1 ), C, S, A( J, J-1 ) )
      A( J+1, J-1 ) = ZERO
*     
      DO 300 JC = J, UPDHI
         TEMP = C*A( J, JC ) + S*A( J+1, JC )
         A( J+1, JC ) = -S*A( J, JC ) + C*A( J+1, JC )
         A( J, JC ) = TEMP
         TEMP2 = C*B( J, JC ) + S*B( J+1, JC )
         B( J+1, JC ) = -S*B( J, JC ) + C*B( J+1, JC )
         B( J, JC ) = TEMP2
 300  CONTINUE
      IF( LCOMPQ ) THEN
         DO 310 JR = 1, N
            TEMP = C*Q( JR, J ) + S*Q( JR, J+1 )
            Q( JR, J+1 ) = -S*Q( JR, J ) + C*Q( JR, J+1 )
            Q( JR, J ) = TEMP
 310     CONTINUE
      END IF
*     
*     Rotations from the right.
*     
      TEMP = B( J+1, J+1 )
      CALL DLARTG( TEMP, B( J+1, J ), C, S, B( J+1, J+1 ) )
      B( J+1, J ) = ZERO
*     
      DO 320 JR = UPDLO, IHI
         TEMP = C*A( JR, J+1 ) + S*A( JR, J )
         A( JR, J ) = -S*A( JR, J+1 ) + C*A( JR, J )
         A( JR, J+1 ) = TEMP
 320  CONTINUE
      DO 330 JR = UPDLO, IHI - 1
         TEMP = C*B( JR, J+1 ) + S*B( JR, J )
         B( JR, J ) = -S*B( JR, J+1 ) + C*B( JR, J )
         B( JR, J+1 ) = TEMP
 330  CONTINUE
      IF( LCOMPZ ) THEN
         DO 340 JR = 1, N
            TEMP = C*Z( JR, J+1 ) + S*Z( JR, J )
            Z( JR, J ) = -S*Z( JR, J+1 ) + C*Z( JR, J )
            Z( JR, J+1 ) = TEMP
 340     CONTINUE
      END IF
      NSHF = 2
*     
*     End of Double-Shift code
*     
      GO TO 350
*     
*     End of iteration loop
*     
 350  CONTINUE


      END
