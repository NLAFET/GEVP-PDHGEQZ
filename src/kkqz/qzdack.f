      SUBROUTINE QZDACK( JOB, COMPQ, COMPZ, N, ILO, IHI, 
     $   A, LDA, B, LDB,
     $   ALPHAR, ALPHAI, BETA, Q, LDQ, Z, LDZ, WORK,
     $   LWORK, INFO, NB )
      implicit none
*     
*     .. Scalar Arguments ..
      CHARACTER          COMPQ, COMPZ, JOB
      INTEGER            IHI, ILO, INFO, LDA, LDB, LDQ, LDZ, LWORK, N
      INTEGER            NB
*     ..
*     .. Array Arguments ..
      DOUBLE PRECISION   A( LDA, * ), ALPHAI( * ), ALPHAR( * ),
     $   B( LDB, * ), BETA( * ), Q( LDQ, * ), WORK( * ),
     $   Z( LDZ, * )
      DOUBLE PRECISION   D1, D2, D3, D4, D5
*     ..
*     
*     Purpose
*     =======
*     
*     KDHGEQZ implements a single-/double-shift version of the QZ method for
*     finding the generalized eigenvalues
*     
*     w(j)=(ALPHAR(j) + i*ALPHAI(j))/BETAR(j)   of the equation
*     
*     det( A - w(i) B ) = 0
*     
*     In addition, the pair A,B may be reduced to generalized Schur form:
*     B is upper triangular, and A is block upper triangular, where the
*     diagonal blocks are either 1-by-1 or 2-by-2, the 2-by-2 blocks having
*     complex generalized eigenvalues (see the description of the argument
*     JOB.)
*     
*     If JOB='S', then the pair (A,B) is simultaneously reduced to Schur
*     form by applying one orthogonal tranformation (usually called Q) on
*     the left and another (usually called Z) on the right.  The 2-by-2
*     upper-triangular diagonal blocks of B corresponding to 2-by-2 blocks
*     of A will be reduced to positive diagonal matrices.  (I.e.,
*     if A(j+1,j) is non-zero, then B(j+1,j)=B(j,j+1)=0 and B(j,j) and
*     B(j+1,j+1) will be positive.)
*     
*     If JOB='E', then at each iteration, the same transformations
*     are computed, but they are only applied to those parts of A and B
*     which are needed to compute ALPHAR, ALPHAI, and BETAR.
*     
*     If JOB='S' and COMPQ and COMPZ are 'V' or 'I', then the orthogonal
*     transformations used to reduce (A,B) are accumulated into the arrays
*     Q and Z s.t.:
*     
*     Q(in) A(in) Z(in)* = Q(out) A(out) Z(out)*
*     Q(in) B(in) Z(in)* = Q(out) B(out) Z(out)*
*     
*     Ref: C.B. Moler & G.W. Stewart, "An Algorithm for Generalized Matrix
*     Eigenvalue Problems", SIAM J. Numer. Anal., 10(1973),
*     pp. 241--256.
*     
*     Arguments
*     =========
*     
*     JOB     (input) CHARACTER*1
*     = 'E': compute only ALPHAR, ALPHAI, and BETA.  A and B will
*     not necessarily be put into generalized Schur form.
*     = 'S': put A and B into generalized Schur form, as well
*     as computing ALPHAR, ALPHAI, and BETA.
*     
*     COMPQ   (input) CHARACTER*1
*     = 'N': do not modify Q.
*     = 'V': multiply the array Q on the right by the transpose of
*     the orthogonal tranformation that is applied to the
*     left side of A and B to reduce them to Schur form.
*     = 'I': like COMPQ='V', except that Q will be initialized to
*     the identity first.
*     
*     COMPZ   (input) CHARACTER*1
*     = 'N': do not modify Z.
*     = 'V': multiply the array Z on the right by the orthogonal
*     tranformation that is applied to the right side of
*     A and B to reduce them to Schur form.
*     = 'I': like COMPZ='V', except that Z will be initialized to
*     the identity first.
*     
*     N       (input) INTEGER
*     The order of the matrices A, B, Q, and Z.  N >= 0.
*     
*     ILO     (input) INTEGER
*     IHI     (input) INTEGER
*     It is assumed that A is already upper triangular in rows and
*     columns 1:ILO-1 and IHI+1:N.
*     1 <= ILO <= IHI <= N, if N > 0; ILO=1 and IHI=0, if N=0.
*     
*     A       (input/output) DOUBLE PRECISION array, dimension (LDA, N)
*     On entry, the N-by-N upper Hessenberg matrix A.  Elements
*     below the subdiagonal must be zero.
*     If JOB='S', then on exit A and B will have been
*     simultaneously reduced to generalized Schur form.
*     If JOB='E', then on exit A will have been destroyed.
*     The diagonal blocks will be correct, but the off-diagonal
*     portion will be meaningless.
*     
*     LDA     (input) INTEGER
*     The leading dimension of the array A.  LDA >= max( 1, N ).
*     
*     B       (input/output) DOUBLE PRECISION array, dimension (LDB, N)
*     On entry, the N-by-N upper triangular matrix B.  Elements
*     below the diagonal must be zero.  2-by-2 blocks in B
*     corresponding to 2-by-2 blocks in A will be reduced to
*     positive diagonal form.  (I.e., if A(j+1,j) is non-zero,
*     then B(j+1,j)=B(j,j+1)=0 and B(j,j) and B(j+1,j+1) will be
*     positive.)
*     If JOB='S', then on exit A and B will have been
*     simultaneously reduced to Schur form.
*     If JOB='E', then on exit B will have been destroyed.
*     Elements corresponding to diagonal blocks of A will be
*     correct, but the off-diagonal portion will be meaningless.
*     
*     LDB     (input) INTEGER
*     The leading dimension of the array B.  LDB >= max( 1, N ).
*     
*     ALPHAR  (output) DOUBLE PRECISION array, dimension (N)
*     ALPHAR(1:N) will be set to real parts of the diagonal
*     elements of A that would result from reducing A and B to
*     Schur form and then further reducing them both to triangular
*     form using unitary transformations s.t. the diagonal of B
*     was non-negative real.  Thus, if A(j,j) is in a 1-by-1 block
*     (i.e., A(j+1,j)=A(j,j+1)=0), then ALPHAR(j)=A(j,j).
*     Note that the (real or complex) values
*     (ALPHAR(j) + i*ALPHAI(j))/BETA(j), j=1,...,N, are the
*     generalized eigenvalues of the matrix pencil A - wB.
*     
*     ALPHAI  (output) DOUBLE PRECISION array, dimension (N)
*     ALPHAI(1:N) will be set to imaginary parts of the diagonal
*     elements of A that would result from reducing A and B to
*     Schur form and then further reducing them both to triangular
*     form using unitary transformations s.t. the diagonal of B
*     was non-negative real.  Thus, if A(j,j) is in a 1-by-1 block
*     (i.e., A(j+1,j)=A(j,j+1)=0), then ALPHAR(j)=0.
*     Note that the (real or complex) values
*     (ALPHAR(j) + i*ALPHAI(j))/BETA(j), j=1,...,N, are the
*     generalized eigenvalues of the matrix pencil A - wB.
*     
*     BETA    (output) DOUBLE PRECISION array, dimension (N)
*     BETA(1:N) will be set to the (real) diagonal elements of B
*     that would result from reducing A and B to Schur form and
*     then further reducing them both to triangular form using
*     unitary transformations s.t. the diagonal of B was
*     non-negative real.  Thus, if A(j,j) is in a 1-by-1 block
*     (i.e., A(j+1,j)=A(j,j+1)=0), then BETA(j)=B(j,j).
*     Note that the (real or complex) values
*     (ALPHAR(j) + i*ALPHAI(j))/BETA(j), j=1,...,N, are the
*     generalized eigenvalues of the matrix pencil A - wB.
*     (Note that BETA(1:N) will always be non-negative, and no
*     BETAI is necessary.)
*     
*     Q       (input/output) DOUBLE PRECISION array, dimension (LDQ, N)
*     If COMPQ='N', then Q will not be referenced.
*     If COMPQ='V' or 'I', then the transpose of the orthogonal
*     transformations which are applied to A and B on the left
*     will be applied to the array Q on the right.
*     
*     LDQ     (input) INTEGER
*     The leading dimension of the array Q.  LDQ >= 1.
*     If COMPQ='V' or 'I', then LDQ >= N.
*     
*     Z       (input/output) DOUBLE PRECISION array, dimension (LDZ, N)
*     If COMPZ='N', then Z will not be referenced.
*     If COMPZ='V' or 'I', then the orthogonal transformations
*     which are applied to A and B on the right will be applied
*     to the array Z on the right.
*     
*     LDZ     (input) INTEGER
*     The leading dimension of the array Z.  LDZ >= 1.
*     If COMPZ='V' or 'I', then LDZ >= N.
*     
*     WORK    (workspace/output) DOUBLE PRECISION array, dimension (LWORK)
*     On exit, if INFO >= 0, WORK(1) returns the optimal LWORK.
*     
*     LWORK   (input) INTEGER
*     The dimension of the array WORK.  LWORK >= 4*N.
*     
*     INFO    (output) INTEGER
*     = 0: successful exit
*     < 0: if INFO = -i, the i-th argument had an illegal value
*     = 1,...,N: the QZ iteration did not converge.  (A,B) is not
*     in Schur form, but ALPHAR(i), ALPHAI(i), and
*     BETA(i), i=INFO+1,...,N should be correct.
*     = N+1,...,2*N: the shift calculation failed.  (A,B) is not
*     in Schur form, but ALPHAR(i), ALPHAI(i), and
*     BETA(i), i=INFO-N+1,...,N should be correct.
*     > 2*N:     various "impossible" errors.
*     
*     Further Details
*     ===============
*     
*     Iteration counters:
*     
*     JITER  -- counts iterations.
*     IITER  -- counts iterations run since ILAST was last
*     changed.  This is therefore reset only when a 1-by-1 or
*     2-by-2 block deflates off the bottom.
*     
*     =====================================================================
*     
*     .. Parameters ..
*     $                     SAFETY = 1.0E+0 )
      DOUBLE PRECISION   HALF, ZERO, ONE, SAFETY
      DOUBLE PRECISION   TIME1, TM2, TIME2
      PARAMETER          ( HALF = 0.5D+0, ZERO = 0.0D+0, ONE = 1.0D+0,
     $   SAFETY = 1.0D+2 )
      COMMON /ooeigdef/    NINF 
      INTEGER              NINF 

      INTEGER           PARA_AGGCRIT, PARA_AGGSHF, PARA_AGGMIN,
     $   PARA_AGGWIN, PARA_HDEF
      
      LOGICAL           PARA_AGGRQUP
      DOUBLE PRECISION  PARA_AGGREP   
      
      COMMON /PARA_DEF/ PARA_AGGCRIT,
     $   PARA_AGGSHF, PARA_AGGMIN, PARA_AGGWIN, PARA_HDEF, PARA_AGGRQUP,
     $   PARA_AGGREP 
      
      LOGICAL           PARA_INFALL, PARA_INFBLOCK
      INTEGER           PARA_INFDEF, PARA_INFDIR, PARA_INFNZ,
     $   PARA_INFWDW
      DOUBLE PRECISION  PARA_INFTOL
      COMMON /PARA_INF/ PARA_INFALL, PARA_INFBLOCK, PARA_INFDEF, 
     $   PARA_INFDIR, PARA_INFNZ, PARA_INFWDW, PARA_INFTOL
      
*     ..
*     .. Local Scalars ..
      LOGICAL            ILAZR2, ILAZRO, ILPIVT, ILQ, ILSCHR, ILZ
      LOGICAL            NDONE
      INTEGER            ICOMPQ, ICOMPZ, IFIRST, IFRSTM, IITER, ILAST,
     $   ILASTM, IN, ISCHUR, ISTART, J, JC, JCH, JITER,
     $   JR, MAXIT, DT, NNB, DT2, RT, NCOLS, JJ, NREF, NBL
      DOUBLE PRECISION   A11, A12, A1I, A1R, A21, A22, A2I, A2R, AD11,
     $   AD11L, AD12, AD12L, AD21, AD21L, AD22, AD22L,
     $   AD32L, AN, ANORM, ASCALE, ATOL, B11, B1A, B1I,
     $   B1R, B22, B2A, B2I, B2R, BN, BNORM, BSCALE,
     $   BTOL, C, C11I, C11R, C12, C21, C22I, C22R, CL,
     $   CQ, CR, CZ, ESHIFT, S, S1, S1INV, S2, SAFMAX,
     $   SAFMIN, SCALE, SL, SQI, SQR, SR, SZI, SZR, T,
     $   TAU, TEMP, TEMP2, TEMPI, TEMPR, U1, U12, U12L,
     $   U2, VS, W11, W12, W21, W22, WABS, WI, WR,
     $   WR2

      DOUBLE PRECISION    SMLNUM, ULP, TFNORM, T1NORM, HFNORM, HTOL,
     $   TTOL
      COMMON /PREC/       SMLNUM, ULP, TFNORM, T1NORM, HFNORM, HTOL,
     $   TTOL

*     ..
*     .. Local Arrays ..
      DOUBLE PRECISION   V( 3 )
*     ..
*     .. External Functions ..
      LOGICAL            LSAME
      DOUBLE PRECISION   DLAMCH, DLANHS, DLANTR, DLAPY2, DLAPY3, DSECND
      EXTERNAL           LSAME, DLAMCH, DLANHS, DLANTR, DLAPY2, DLAPY3,
     $   DSECND
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
*     NB = ILAENV( 1, '   ', ' ', N, N, -1, -1 )
*     write(*,*)'NB =',NB
*     
*     Decode JOB, COMPQ, COMPZ
*     
      RT = 3*N+1
      DT = 0
      DT2 = 0
      TIME1 = 0.d0
      TIME2 = 0.d0
      IF( LSAME( JOB, 'E' ) ) THEN
         ILSCHR = .FALSE.
         ISCHUR = 1
      ELSE IF( LSAME( JOB, 'S' ) ) THEN
         ILSCHR = .TRUE.
         ISCHUR = 2
      ELSE
         ISCHUR = 0
      END IF
*     
      IF( LSAME( COMPQ, 'N' ) ) THEN
         ILQ = .FALSE.
         ICOMPQ = 1
      ELSE IF( LSAME( COMPQ, 'V' ) ) THEN
         ILQ = .TRUE.
         ICOMPQ = 2
      ELSE IF( LSAME( COMPQ, 'I' ) ) THEN
         ILQ = .TRUE.
         ICOMPQ = 3
      ELSE
         ICOMPQ = 0
      END IF
*     
      IF( LSAME( COMPZ, 'N' ) ) THEN
         ILZ = .FALSE.
         ICOMPZ = 1
      ELSE IF( LSAME( COMPZ, 'V' ) ) THEN
         ILZ = .TRUE.
         ICOMPZ = 2
      ELSE IF( LSAME( COMPZ, 'I' ) ) THEN
         ILZ = .TRUE.
         ICOMPZ = 3
      ELSE
         ICOMPZ = 0
      END IF
*     
*     Check Argument Values
*     
      INFO = 0
      IF( ISCHUR.EQ.0 ) THEN
         INFO = -1
      ELSE IF( ICOMPQ.EQ.0 ) THEN
         INFO = -2
      ELSE IF( ICOMPZ.EQ.0 ) THEN
         INFO = -3
      ELSE IF( N.LT.0 ) THEN
         INFO = -4
      ELSE IF( ILO.LT.1 ) THEN
         INFO = -5
      ELSE IF( IHI.GT.N .OR. IHI.LT.ILO-1 ) THEN
         INFO = -6
      ELSE IF( LDA.LT.N ) THEN
         INFO = -8
      ELSE IF( LDB.LT.N ) THEN
         INFO = -10
      ELSE IF( LDQ.LT.1 .OR. ( ILQ .AND. LDQ.LT.N ) ) THEN
         INFO = -15
      ELSE IF( LDZ.LT.1 .OR. ( ILZ .AND. LDZ.LT.N ) ) THEN
         INFO = -17
      ELSE IF( LWORK.LT.MAX( 1, 4*N ) ) THEN
         INFO = -19
      ELSE IF (NB .LT. 5) THEN 
         INFO = -20
      END IF
      IF( INFO.NE.0 ) THEN
         CALL XERBLA( 'DHGEQZ', -INFO )
         RETURN
      END IF
*     
*     Quick return if possible
*     
      IF( N.LE.0 ) THEN
         WORK( 1 ) = DBLE( 1 )
         RETURN
      END IF
*     
*     Initialize Q and Z
*     
      IF( ICOMPQ.EQ.3 )
     $   CALL DLASET( 'Full', N, N, ZERO, ONE, Q, LDQ )
      IF( ICOMPZ.EQ.3 )
     $   CALL DLASET( 'Full', N, N, ZERO, ONE, Z, LDZ )
*     
*     Machine Constants
*     
      IN = IHI + 1 - ILO
      SAFMIN = DLAMCH( 'S' )
      SAFMAX = ONE / SAFMIN
c     ULP = DLAMCH( 'E' )*DLAMCH( 'B' )
      ANORM = HFNORM
      BNORM = TFNORM
      ATOL = MAX( SAFMIN, ULP*ANORM )
      BTOL = MAX( SAFMIN, ULP*BNORM )
      ASCALE = ONE / MAX( SAFMIN, ANORM )
      BSCALE = ONE / MAX( SAFMIN, BNORM )
*     
*     Set Eigenvalues IHI+1:N
*     
C     DO 30 J = IHI + 1, N
C     IF( B( J, J ).LT.ZERO ) THEN
C     IF( ILSCHR ) THEN
C     DO 10 JR = 1, J
C     A( JR, J ) = -A( JR, J )
C     B( JR, J ) = -B( JR, J )
C     10          CONTINUE
C     ELSE
C     A( J, J ) = -A( J, J )
C     B( J, J ) = -B( J, J )
C     END IF
C     IF( ILZ ) THEN
C     DO 20 JR = 1, N
C     Z( JR, J ) = -Z( JR, J )
C     20          CONTINUE
C     END IF
C     END IF
C     ALPHAR( J ) = A( J, J )
C     ALPHAI( J ) = ZERO
C     BETA( J ) = B( J, J )
C     30 CONTINUE
*     
*     If IHI < ILO, skip QZ steps
*     
      IF( IHI.LT.ILO )
     $   GO TO 380
*     
*     MAIN QZ ITERATION LOOP
*     
*     Initialize dynamic indices
*     
*     Eigenvalues ILAST+1:N have been found.
*     Column operations modify rows IFRSTM:whatever.
*     Row operations modify columns whatever:ILASTM.
*     
*     If only eigenvalues are being computed, then
*     IFRSTM is the row of the last splitting row above row ILAST;
*     this is always at least ILO.
*     IITER counts iterations since the last eigenvalue was found,
*     to tell when to use an extraordinary shift.
*     MAXIT is the maximum number of QZ sweeps allowed.
*     
      ILAST = IHI
      IF( ILSCHR ) THEN
         IFRSTM = 1
         ILASTM = N
      ELSE
         IFRSTM = ILO
         ILASTM = IHI
      END IF
      IITER = 0
      ESHIFT = ZERO
      MAXIT = 30*( IHI-ILO+1 )
*     
      DO 360 JITER = 1, MAXIT
*        
*        Split the matrix if possible.
*        
*        Two tests:
*        1: A(j,j-1)=0  or  j=ILO
*        2: B(j,j)=0
*       
         
         IF( ILAST.EQ.ILO ) THEN
*           
*           Special case: j=ILAST
*           
            GO TO 80
         ELSE
            IF( ABS( A( ILAST, ILAST-1 ) ).LE.ATOL ) THEN
               A( ILAST, ILAST-1 ) = ZERO
               GO TO 80
            END IF
         END IF

         IF( ABS( B( ILAST, ILAST ) ).LE.BTOL ) THEN
            B( ILAST, ILAST ) = ZERO
            NINF = NINF + 1
            GO TO 70
         END IF
*        
*        General case: j<ILAST
*        
         DO 60 J = ILAST - 1, ILO, -1
*           
*           Test 1: for A(j,j-1)=0 or j=ILO
*           
            IF( J.EQ.ILO ) THEN
               ILAZRO = .TRUE.
            ELSE
               IF( ABS( A( J, J-1 ) ).LE.ATOL ) THEN
                  A( J, J-1 ) = ZERO
                  ILAZRO = .TRUE.
               ELSE
                  ILAZRO = .FALSE.
               END IF
            END IF
*           
*           Test 2: for B(j,j)=0
*           
*           
            IF( ABS( B( J, J ) ).LT.BTOL ) THEN
               B( J, J ) = ZERO
               NINF = NINF + 1
*              
*              Test 1a: Check for 2 consecutive small subdiagonals in A
*              Modification: This test is skipped.
*              
               ILAZR2 = .FALSE.
*              
*              If both tests pass (1 & 2), i.e., the leading diagonal
*              element of B in the block is zero, split a 1x1 block off
*              at the top. (I.e., at the J-th row/column) The leading
*              diagonal element of the remainder can also be zero, so
*              this may have to be done repeatedly.
*              
               IF( ILAZRO .OR. ILAZR2 ) THEN
                  DO 40 JCH = J, ILAST - 1
                     TEMP = A( JCH, JCH )
                     CALL DLARTG( TEMP, A( JCH+1, JCH ), C, S,
     $                  A( JCH, JCH ) )
                     A( JCH+1, JCH ) = ZERO
                     CALL DROT( ILASTM-JCH, A( JCH, JCH+1 ), LDA,
     $                  A( JCH+1, JCH+1 ), LDA, C, S )
                     CALL DROT( ILASTM-JCH, B( JCH, JCH+1 ), LDB,
     $                  B( JCH+1, JCH+1 ), LDB, C, S )
                     IF( ILQ )
     $                  CALL DROT( N, Q( 1, JCH ), 1, Q( 1, JCH+1 ), 1,
     $                  C, S )
                     IF( ILAZR2 )
     $                  A( JCH, JCH-1 ) = A( JCH, JCH-1 )*C
                     ILAZR2 = .FALSE.
                     IF( ABS( B( JCH+1, JCH+1 ) ).GE.BTOL ) THEN
                        IF( JCH+1.GE.ILAST ) THEN
                           GO TO 80
                        ELSE
                           IFIRST = JCH + 1
                           GO TO 110
                        END IF
                     END IF
                     B( JCH+1, JCH+1 ) = ZERO
 40               CONTINUE
                  GO TO 70
               ELSE
*                 
*                 Only test 2 passed -- chase the zero to B(ILAST,ILAST)
*                 Then process as in the case B(ILAST,ILAST)=0
*                 
                  DO 50 JCH = J, ILAST - 1
                     TEMP = B( JCH, JCH+1 )
                     CALL DLARTG( TEMP, B( JCH+1, JCH+1 ), C, S,
     $                  B( JCH, JCH+1 ) )
                     B( JCH+1, JCH+1 ) = ZERO
                     IF( JCH.LT.ILASTM-1 )
     $                  CALL DROT( ILASTM-JCH-1, B( JCH, JCH+2 ), LDB,
     $                  B( JCH+1, JCH+2 ), LDB, C, S )
                     CALL DROT( ILASTM-JCH+2, A( JCH, JCH-1 ), LDA,
     $                  A( JCH+1, JCH-1 ), LDA, C, S )
                     IF( ILQ )
     $                  CALL DROT( N, Q( 1, JCH ), 1, Q( 1, JCH+1 ), 1,
     $                  C, S )
                     TEMP = A( JCH+1, JCH )
                     CALL DLARTG( TEMP, A( JCH+1, JCH-1 ), C, S,
     $                  A( JCH+1, JCH ) )
                     A( JCH+1, JCH-1 ) = ZERO
                     CALL DROT( JCH+1-IFRSTM, A( IFRSTM, JCH ), 1,
     $                  A( IFRSTM, JCH-1 ), 1, C, S )
                     CALL DROT( JCH-IFRSTM, B( IFRSTM, JCH ), 1,
     $                  B( IFRSTM, JCH-1 ), 1, C, S )
                     IF( ILZ )
     $                  CALL DROT( N, Z( 1, JCH ), 1, Z( 1, JCH-1 ), 1,
     $                  C, S )
 50               CONTINUE
                  GO TO 70
               END IF
            ELSE IF( ILAZRO ) THEN
*              
*              Only test 1 passed -- work on J:ILAST
*              
               IFIRST = J
               GO TO 110
            END IF
*           
*           Neither test passed -- try next J
*           
 60      CONTINUE
*        
*        (Drop-through is "impossible")
*        
         INFO = N + 1
         GO TO 420
*        
*        B(ILAST,ILAST)=0 -- clear A(ILAST,ILAST-1) to split off a
*        1x1 block.
*        
 70      CONTINUE
         TEMP = A( ILAST, ILAST )
         CALL DLARTG( TEMP, A( ILAST, ILAST-1 ), C, S,
     $      A( ILAST, ILAST ) )
         A( ILAST, ILAST-1 ) = ZERO
         CALL DROT( ILAST-IFRSTM, A( IFRSTM, ILAST ), 1,
     $      A( IFRSTM, ILAST-1 ), 1, C, S )
         CALL DROT( ILAST-IFRSTM, B( IFRSTM, ILAST ), 1,
     $      B( IFRSTM, ILAST-1 ), 1, C, S )
         IF( ILZ )
     $      CALL DROT( N, Z( 1, ILAST ), 1, Z( 1, ILAST-1 ), 1, C, S )
*        
*        A(ILAST,ILAST-1)=0 -- Standardize B, set ALPHAR, ALPHAI,
*        and BETA
*        
 80      CONTINUE
         IF( B( ILAST, ILAST ).LT.ZERO ) THEN
            IF( ILSCHR ) THEN
               DO 90 J = IFRSTM, ILAST
                  A( J, ILAST ) = -A( J, ILAST )
                  B( J, ILAST ) = -B( J, ILAST )
 90            CONTINUE
            ELSE
               A( ILAST, ILAST ) = -A( ILAST, ILAST )
               B( ILAST, ILAST ) = -B( ILAST, ILAST )
            END IF
            IF( ILZ ) THEN
               DO 100 J = 1, N
                  Z( J, ILAST ) = -Z( J, ILAST )
 100           CONTINUE
            END IF
         END IF
         ALPHAR( ILAST ) = A( ILAST, ILAST )
         ALPHAI( ILAST ) = ZERO
         BETA( ILAST ) = B( ILAST, ILAST )
*        
*        Go to next block -- exit if finished.
*        
         ILAST = ILAST - 1
         IF( ILAST.LT.ILO )
     $      GO TO 380
*        
*        Reset counters
*        
         IITER = 0
         ESHIFT = ZERO
         IF( .NOT.ILSCHR ) THEN
            ILASTM = ILAST
            IF( IFRSTM.GT.ILAST )
     $         IFRSTM = ILO
         END IF
         GO TO 350
*        
*        QZ step
*        
*        This iteration only involves rows/columns IFIRST:ILAST. We
*        assume IFIRST < ILAST, and that the diagonal of B is non-zero.
*        
 110     CONTINUE
         IITER = IITER + 1
         DT = DT + 1
         IF( .NOT.ILSCHR ) THEN
            IFRSTM = IFIRST
         END IF
*        
*        Compute single shifts.
*        
*        At this point, IFIRST < ILAST, and the diagonal elements of
*        B(IFIRST:ILAST,IFIRST,ILAST) are larger than BTOL (in
*        magnitude)
*        
         IF( ( IITER / 10 )*10.EQ.IITER ) THEN
*           
*           Exceptional shift.  Chosen for no particularly good reason.
*           (Single shift only.)
*           
            IF( ( DBLE( MAXIT )*SAFMIN )*ABS( A( ILAST-1, ILAST ) ).LT.
     $         ABS( B( ILAST-1, ILAST-1 ) ) ) THEN
               ESHIFT = ESHIFT + A( ILAST-1, ILAST ) /
     $            B( ILAST-1, ILAST-1 )
            ELSE
               ESHIFT = ESHIFT + ONE / ( SAFMIN*DBLE( MAXIT ) )
            END IF
            S1 = ONE
            WR = ESHIFT
*           
         ELSE
*           
*           Shifts based on the generalized eigenvalues of the
*           bottom-right 2x2 block of A and B. The first eigenvalue
*           returned by DLAG2 is the Wilkinson shift (AEP p.512),
*           
            CALL DLAG2( A( ILAST-1, ILAST-1 ), LDA,
     $         B( ILAST-1, ILAST-1 ), LDB, SAFMIN*SAFETY, S1,
     $         S2, WR, WR2, WI )
*           
            TEMP = MAX( S1, SAFMIN*MAX( ONE, ABS( WR ), ABS( WI ) ) )
            IF( WI.NE.ZERO )
     $         GO TO 200
         END IF
*        
*        Fiddle with shift to avoid overflow
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
     $      SCALE = MIN( SCALE, TEMP / ABS( WR ) )
         S1 = SCALE*S1
         WR = SCALE*WR
*        
*        Now check for two consecutive small subdiagonals.
*        
         DO 120 J = ILAST - 1, IFIRST + 1, -1
            ISTART = J
            TEMP = ABS( S1*A( J, J-1 ) )
            TEMP2 = ABS( S1*A( J, J )-WR*B( J, J ) )
            TEMPR = MAX( TEMP, TEMP2 )
            IF( TEMPR.LT.ONE .AND. TEMPR.NE.ZERO ) THEN
               TEMP = TEMP / TEMPR
               TEMP2 = TEMP2 / TEMPR
            END IF
            IF( ABS( ( ASCALE*A( J+1, J ) )*TEMP ).LE.( ASCALE*ATOL )*
     $         TEMP2 )GO TO 130
 120     CONTINUE
*        
         ISTART = IFIRST
 130     CONTINUE
*        
*        Do an implicit single-shift QZ sweep.
*        
*        Initial Q
*        
         TEMP = S1*A( ISTART, ISTART ) - WR*B( ISTART, ISTART )
         TEMP2 = S1*A( ISTART+1, ISTART )
         CALL DLARTG( TEMP, TEMP2, C, S, TEMPR )
*        
*        Sweep
*        
         TM2 = DSECND()
         J = ISTART
*        
         DO 140 JC = J, ILASTM
            TEMP = C*A( J, JC ) + S*A( J+1, JC )
            A( J+1, JC ) = -S*A( J, JC ) + C*A( J+1, JC )
            A( J, JC ) = TEMP
            TEMP2 = C*B( J, JC ) + S*B( J+1, JC )
            B( J+1, JC ) = -S*B( J, JC ) + C*B( J+1, JC )
            B( J, JC ) = TEMP2
 140     CONTINUE
         IF( ILQ ) THEN
            DO 150 JR = 1, N
               TEMP = C*Q( JR, J ) + S*Q( JR, J+1 )
               Q( JR, J+1 ) = -S*Q( JR, J ) + C*Q( JR, J+1 )
               Q( JR, J ) = TEMP
 150        CONTINUE
         ENDIF
*        
         TEMP = B( J+1, J+1 )
         CALL DLARTG( TEMP, B( J+1, J ), C, S, B( J+1, J+1 ) )
         B( J+1, J ) = ZERO
*        
         DO 160 JR = IFRSTM, MIN( J+2, ILAST )
            TEMP = C*A( JR, J+1 ) + S*A( JR, J )
            A( JR, J ) = -S*A( JR, J+1 ) + C*A( JR, J )
            A( JR, J+1 ) = TEMP
 160     CONTINUE
         DO 170 JR = IFRSTM, J
            TEMP = C*B( JR, J+1 ) + S*B( JR, J )
            B( JR, J ) = -S*B( JR, J+1 ) + C*B( JR, J )
            B( JR, J+1 ) = TEMP
 170     CONTINUE
         IF( ILZ ) THEN
            DO 180 JR = 1, N
               TEMP = C*Z( JR, J+1 ) + S*Z( JR, J )
               Z( JR, J ) = -S*Z( JR, J+1 ) + C*Z( JR, J )
               Z( JR, J+1 ) = TEMP
 180        CONTINUE
         END IF
*        
*        New Blocked SIngle Code Starts Here
*        
         JJ = IFIRST
         NDONE = .TRUE.
         NREF = 1
 190     CONTINUE
         NNB = MIN(NB, ILAST-JJ+1)
         NCOLS = 2
         NBL = NNB - 3
*        
*        Move bulge and accumulate transformations Q and Z
*        
*        Last block?
         IF ((JJ + NNB -1) .EQ. ILAST) THEN
            NCOLS = 1
            NDONE = .FALSE.
            NBL = ILASTM - (JJ + 2)
         ENDIF
*        
*        Transformations from the left
*        
         IF (JJ-IFIRST .GT. 0) THEN
            CALL MYROT('Rows', A(1,JJ+3),
     $         LDA, IFIRST+1, 1, 1, JJ-IFIRST,
     $         1, NBL, WORK)
            CALL MYROT('Rows', B(1,JJ+3),
     $         LDB, IFIRST+1, 1, 1, JJ-IFIRST,
     $         1, NBL, WORK)
         ENDIF
         CALL GIV(NNB, A(JJ, JJ), LDA, B(JJ, JJ), LDB,
     $      WORK(NREF), WORK(RT), NCOLS)
*        
*        Transformations from the left (last iteration)
*        
         IF (.NOT. NDONE) THEN
            CALL MYROT('Rows', A(JJ,JJ+NNB), LDA,
     $         2, 1, 1, NNB-(NCOLS+1), 1, ILASTM-JJ-NNB+1,
     $         WORK(NREF) )
            CALL MYROT('Rows', B(JJ,JJ+NNB), LDB,
     $         2, 2, 1, NNB-(NCOLS+1), 1, ILASTM-JJ-NNB+1,
     $         WORK(NREF))
         ENDIF
         
         IF (ILQ) THEN
            CALL MYROT('X', Q(1,JJ), LDQ, 
     $         1, 2, 1, NNB-(NCOLS+1),
     $         1, N, WORK(NREF) )
         ENDIF
*        
*        Transformations from the right
*        
         
         IF (JJ-IFRSTM .GT. 0) THEN
            CALL MYROT('Cols', A(IFRSTM,JJ), LDA,
     $         1, 2, 1, NNB-(NCOLS+1), 1, JJ-IFRSTM,
     $         WORK(RT))
            CALL MYROT('Cols', B(IFRSTM,JJ), LDB,
     $         1, 2, 1, NNB-(NCOLS+1), 1, JJ-IFRSTM,
     $         WORK(RT))
         ENDIF
         
         IF (ILZ) THEN
            CALL MYROT('Cols', Z(1,JJ), LDZ, 1, 2, 1, 
     $         NNB-(NCOLS+1), 1, N, WORK(RT))
         ENDIF
         NREF = NREF + (NNB-(NCOLS+1)) * 2
         JJ = JJ + NB - 3
         IF (NDONE)  GOTO 190
*        
*        New Blocked Code Ends Here
*        
         TIME2 = TIME2 +  DSECND() - TM2
*        
         GO TO 350
*        
*        Use Francis double-shift
*        
*        Note: the Francis double-shift should work with real shifts,
*        but only if the block is at least 3x3.
*        This code may break if this point is reached with
*        a 2x2 block with real eigenvalues.
*        
 200     CONTINUE
         IF( IFIRST+1.EQ.ILAST ) THEN
*           
*           Special case -- 2x2 block with complex eigenvectors
*           
*           Step 1: Standardize, that is, rotate so that
*           
*           ( B11  0  )
*           B = (         )  with B11 non-negative.
*           (  0  B22 )
*           
            CALL DLASV2( B( ILAST-1, ILAST-1 ), B( ILAST-1, ILAST ),
     $         B( ILAST, ILAST ), B22, B11, SR, CR, SL, CL )
*           
            IF( B11.LT.ZERO ) THEN
               CR = -CR
               SR = -SR
               B11 = -B11
               B22 = -B22
            END IF
*           
            CALL DROT( ILASTM+1-IFIRST, A( ILAST-1, ILAST-1 ), LDA,
     $         A( ILAST, ILAST-1 ), LDA, CL, SL )
            CALL DROT( ILAST+1-IFRSTM, A( IFRSTM, ILAST-1 ), 1,
     $         A( IFRSTM, ILAST ), 1, CR, SR )
*           
            IF( ILAST.LT.ILASTM )
     $         CALL DROT( ILASTM-ILAST, B( ILAST-1, ILAST+1 ), LDB,
     $         B( ILAST, ILAST+1 ), LDA, CL, SL )
            IF( IFRSTM.LT.ILAST-1 )
     $         CALL DROT( IFIRST-IFRSTM, B( IFRSTM, ILAST-1 ), 1,
     $         B( IFRSTM, ILAST ), 1, CR, SR )
*           
            IF( ILQ )
     $         CALL DROT( N, Q( 1, ILAST-1 ), 1, Q( 1, ILAST ), 1, CL,
     $         SL )
            IF( ILZ )
     $         CALL DROT( N, Z( 1, ILAST-1 ), 1, Z( 1, ILAST ), 1, CR,
     $         SR )
*           
            B( ILAST-1, ILAST-1 ) = B11
            B( ILAST-1, ILAST ) = ZERO
            B( ILAST, ILAST-1 ) = ZERO
            B( ILAST, ILAST ) = B22
*           
*           If B22 is negative, negate column ILAST
*           
            IF( B22.LT.ZERO ) THEN
               DO 210 J = IFRSTM, ILAST
                  A( J, ILAST ) = -A( J, ILAST )
                  B( J, ILAST ) = -B( J, ILAST )
 210           CONTINUE
*              
               IF( ILZ ) THEN
                  DO 220 J = 1, N
                     Z( J, ILAST ) = -Z( J, ILAST )
 220              CONTINUE
               END IF
            END IF
*           
*           Step 2: Compute ALPHAR, ALPHAI, and BETA (see refs.)
*           
*           Recompute shift
*           
            CALL DLAG2( A( ILAST-1, ILAST-1 ), LDA,
     $         B( ILAST-1, ILAST-1 ), LDB, SAFMIN*SAFETY, S1,
     $         TEMP, WR, TEMP2, WI )
*           
*           If standardization has perturbed the shift onto real line,
*           do another (real single-shift) QR step.
*           
            IF( WI.EQ.ZERO )
     $         GO TO 350
            S1INV = ONE / S1
*           
*           Do EISPACK (QZVAL) computation of alpha and beta
*           
            A11 = A( ILAST-1, ILAST-1 )
            A21 = A( ILAST, ILAST-1 )
            A12 = A( ILAST-1, ILAST )
            A22 = A( ILAST, ILAST )
*           
*           Compute complex Givens rotation on right
*           (Assume some element of C = (sA - wB) > unfl )
*           __
*           (sA - wB) ( CZ   -SZ )
*           ( SZ    CZ )
*           
            C11R = S1*A11 - WR*B11
            C11I = -WI*B11
            C12 = S1*A12
            C21 = S1*A21
            C22R = S1*A22 - WR*B22
            C22I = -WI*B22
*           
            IF( ABS( C11R )+ABS( C11I )+ABS( C12 ).GT.ABS( C21 )+
     $         ABS( C22R )+ABS( C22I ) ) THEN
               T = DLAPY3( C12, C11R, C11I )
               CZ = C12 / T
               SZR = -C11R / T
               SZI = -C11I / T
            ELSE
               CZ = DLAPY2( C22R, C22I )
               IF( CZ.LE.SAFMIN ) THEN
                  CZ = ZERO
                  SZR = ONE
                  SZI = ZERO
               ELSE
                  TEMPR = C22R / CZ
                  TEMPI = C22I / CZ
                  T = DLAPY2( CZ, C21 )
                  CZ = CZ / T
                  SZR = -C21*TEMPR / T
                  SZI = C21*TEMPI / T
               END IF
            END IF
*           
*           Compute Givens rotation on left
*           
*           (  CQ   SQ )
*           (  __      )  A or B
*           ( -SQ   CQ )
*           
            AN = ABS( A11 ) + ABS( A12 ) + ABS( A21 ) + ABS( A22 )
            BN = ABS( B11 ) + ABS( B22 )
            WABS = ABS( WR ) + ABS( WI )
            IF( S1*AN.GT.WABS*BN ) THEN
               CQ = CZ*B11
               SQR = SZR*B22
               SQI = -SZI*B22
            ELSE
               A1R = CZ*A11 + SZR*A12
               A1I = SZI*A12
               A2R = CZ*A21 + SZR*A22
               A2I = SZI*A22
               CQ = DLAPY2( A1R, A1I )
               IF( CQ.LE.SAFMIN ) THEN
                  CQ = ZERO
                  SQR = ONE
                  SQI = ZERO
               ELSE
                  TEMPR = A1R / CQ
                  TEMPI = A1I / CQ
                  SQR = TEMPR*A2R + TEMPI*A2I
                  SQI = TEMPI*A2R - TEMPR*A2I
               END IF
            END IF
            T = DLAPY3( CQ, SQR, SQI )
            CQ = CQ / T
            SQR = SQR / T
            SQI = SQI / T
*           
*           Compute diagonal elements of QBZ
*           
            TEMPR = SQR*SZR - SQI*SZI
            TEMPI = SQR*SZI + SQI*SZR
            B1R = CQ*CZ*B11 + TEMPR*B22
            B1I = TEMPI*B22
            B1A = DLAPY2( B1R, B1I )
            B2R = CQ*CZ*B22 + TEMPR*B11
            B2I = -TEMPI*B11
            B2A = DLAPY2( B2R, B2I )
*           
*           Normalize so beta > 0, and Im( alpha1 ) > 0
*           
            BETA( ILAST-1 ) = B1A
            BETA( ILAST ) = B2A
            ALPHAR( ILAST-1 ) = ( WR*B1A )*S1INV
            ALPHAI( ILAST-1 ) = ( WI*B1A )*S1INV
            ALPHAR( ILAST ) = ( WR*B2A )*S1INV
            ALPHAI( ILAST ) = -( WI*B2A )*S1INV
*           
*           Step 3: Go to next block -- exit if finished.
*           
            ILAST = IFIRST - 1
            IF( ILAST.LT.ILO )
     $         GO TO 380
*           
*           Reset counters
*           
            IITER = 0
            ESHIFT = ZERO
            IF( .NOT.ILSCHR ) THEN
               ILASTM = ILAST
               IF( IFRSTM.GT.ILAST )
     $            IFRSTM = ILO
            END IF
            GO TO 350
         ELSE
            TM2 = DSECND()
            DT2 = DT2 + 1
*           
*           Usual case: 3x3 or larger block, using Francis implicit
*           double-shift
*           
*           2
*           Eigenvalue equation is  w  - c w + d = 0,
*           
*           -1 2        -1
*           so compute 1st column of  (A B  )  - c A B   + d
*           using the formula in QZIT (from EISPACK)
*           
*           We assume that the block is at least 3x3
*           
            AD11 = ( ASCALE*A( ILAST-1, ILAST-1 ) ) /
     $         ( BSCALE*B( ILAST-1, ILAST-1 ) )
            AD21 = ( ASCALE*A( ILAST, ILAST-1 ) ) /
     $         ( BSCALE*B( ILAST-1, ILAST-1 ) )
            AD12 = ( ASCALE*A( ILAST-1, ILAST ) ) /
     $         ( BSCALE*B( ILAST, ILAST ) )
            AD22 = ( ASCALE*A( ILAST, ILAST ) ) /
     $         ( BSCALE*B( ILAST, ILAST ) )
            U12 = B( ILAST-1, ILAST ) / B( ILAST, ILAST )
            AD11L = ( ASCALE*A( IFIRST, IFIRST ) ) /
     $         ( BSCALE*B( IFIRST, IFIRST ) )
            AD21L = ( ASCALE*A( IFIRST+1, IFIRST ) ) /
     $         ( BSCALE*B( IFIRST, IFIRST ) )
            AD12L = ( ASCALE*A( IFIRST, IFIRST+1 ) ) /
     $         ( BSCALE*B( IFIRST+1, IFIRST+1 ) )
            AD22L = ( ASCALE*A( IFIRST+1, IFIRST+1 ) ) /
     $         ( BSCALE*B( IFIRST+1, IFIRST+1 ) )
            AD32L = ( ASCALE*A( IFIRST+2, IFIRST+1 ) ) /
     $         ( BSCALE*B( IFIRST+1, IFIRST+1 ) )
            U12L = B( IFIRST, IFIRST+1 ) / B( IFIRST+1, IFIRST+1 )
*           
            V( 1 ) = ( AD11-AD11L )*( AD22-AD11L ) - AD12*AD21 +
     $         AD21*U12*AD11L + ( AD12L-AD11L*U12L )*AD21L
            V( 2 ) = ( ( AD22L-AD11L )-AD21L*U12L-( AD11-AD11L )-
     $         ( AD22-AD11L )+AD21*U12 )*AD21L
            V( 3 ) = AD32L*AD21L
*           
            ISTART = IFIRST
*           
            CALL DLARFG( 3, V( 1 ), V( 2 ), 1, TAU )
            V( 1 ) = ONE
*           
*           Sweep XXX
*           
            J = ISTART

*           
*           All but last elements: use 3x3 Householder transforms.
*           
*           Zero (j-1)st column of A
*           
            DO 230 JC = J, ILASTM
               TEMP = TAU*( A( J, JC )+V( 2 )*A( J+1, JC )+V( 3 )*
     $            A( J+2, JC ) )
               A( J, JC ) = A( J, JC ) - TEMP
               A( J+1, JC ) = A( J+1, JC ) - TEMP*V( 2 )
               A( J+2, JC ) = A( J+2, JC ) - TEMP*V( 3 )
               TEMP2 = TAU*( B( J, JC )+V( 2 )*B( J+1, JC )+V( 3 )*
     $            B( J+2, JC ) )
               B( J, JC ) = B( J, JC ) - TEMP2
               B( J+1, JC ) = B( J+1, JC ) - TEMP2*V( 2 )
               B( J+2, JC ) = B( J+2, JC ) - TEMP2*V( 3 )
 230        CONTINUE
            IF( ILQ ) THEN
               DO 240 JR = 1, N
                  TEMP = TAU*( Q( JR, J )+V( 2 )*Q( JR, J+1 )+V( 3 )*
     $               Q( JR, J+2 ) )
                  Q( JR, J ) = Q( JR, J ) - TEMP
                  Q( JR, J+1 ) = Q( JR, J+1 ) - TEMP*V( 2 )
                  Q( JR, J+2 ) = Q( JR, J+2 ) - TEMP*V( 3 )
 240           CONTINUE
            END IF
*           
*           Zero j-th column of B (see DLAGBC for details)
*           
*           Swap rows to pivot
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
*           Swap columns if nec.
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
*           LU-factor
*           
            TEMP = W21 / W11
            U2 = U2 - TEMP*U1
            W22 = W22 - TEMP*W12
            W21 = ZERO
*           
*           Compute SCALE
*           
            SCALE = ONE
            IF( ABS( W22 ).LT.SAFMIN ) THEN
               SCALE = ZERO
               U2 = ONE
               U1 = -W12 / W11
               GO TO 250
            END IF
            IF( ABS( W22 ).LT.ABS( U2 ) )
     $         SCALE = ABS( W22 / U2 )
            IF( ABS( W11 ).LT.ABS( U1 ) )
     $         SCALE = MIN( SCALE, ABS( W11 / U1 ) )
*           
*           Solve
*           
            U2 = ( SCALE*U2 ) / W22
            U1 = ( SCALE*U1-W12*U2 ) / W11
*           
 250        CONTINUE
            IF( ILPIVT ) THEN
               TEMP = U2
               U2 = U1
               U1 = TEMP
            END IF
*           
*           Compute Householder Vector
*           
            T = SQRT( SCALE**2+U1**2+U2**2 )
            TAU = ONE + SCALE / T
            VS = -ONE / ( SCALE+T )
            V( 1 ) = ONE
            V( 2 ) = VS*U1
            V( 3 ) = VS*U2
*           
*           Apply transformations from the right.
*           
            DO 260 JR = IFRSTM, MIN( J+3, ILAST )
               TEMP = TAU*( A( JR, J )+V( 2 )*A( JR, J+1 )+V( 3 )*
     $            A( JR, J+2 ) )
               A( JR, J ) = A( JR, J ) - TEMP
               A( JR, J+1 ) = A( JR, J+1 ) - TEMP*V( 2 )
               A( JR, J+2 ) = A( JR, J+2 ) - TEMP*V( 3 )
 260        CONTINUE
            DO 270 JR = IFRSTM, J + 2
               TEMP = TAU*( B( JR, J )+V( 2 )*B( JR, J+1 )+V( 3 )*
     $            B( JR, J+2 ) )
               B( JR, J ) = B( JR, J ) - TEMP
               B( JR, J+1 ) = B( JR, J+1 ) - TEMP*V( 2 )
               B( JR, J+2 ) = B( JR, J+2 ) - TEMP*V( 3 )
 270        CONTINUE
            IF( ILZ ) THEN
               DO 280 JR = 1, N
                  TEMP = TAU*( Z( JR, J )+V( 2 )*Z( JR, J+1 )+V( 3 )*
     $               Z( JR, J+2 ) )
                  Z( JR, J ) = Z( JR, J ) - TEMP
                  Z( JR, J+1 ) = Z( JR, J+1 ) - TEMP*V( 2 )
                  Z( JR, J+2 ) = Z( JR, J+2 ) - TEMP*V( 3 )
 280           CONTINUE
            END IF
            B( J+1, J ) = ZERO
            B( J+2, J ) = ZERO
*           
*           New Blocked Code Starts Here
*           
            JJ = IFIRST
            NDONE = .TRUE.
            NREF = 1
 1000       CONTINUE
            NNB = MIN(NB, ILAST-JJ+1)
            NCOLS = 3
            NBL = NNB - 4
*           
*           Move bulge and accumulate transformations Q and Z
*           
*           
*           Last block? 
            IF ((JJ + NNB -1) .EQ. ILAST) THEN
               NCOLS = 2
               NDONE = .FALSE.
               NBL = ILASTM - (JJ + 3) 
            ENDIF                             
*           
*           Transformations from the left
*           
            IF (JJ-IFIRST .GT. 0) THEN
               CALL DLAREF('Rows', A(1,JJ+4), 
     $            LDA, .FALSE., Q,
     $            LDQ, .TRUE., IFIRST+1, 1, 1, JJ-IFIRST, 
     $            1, NBL, 1, 1, WORK, 
     $            D1, D2, D3, D4, D5)
               CALL DLAREF('Rows', B(1,JJ+4), 
     $            LDB, .FALSE., Q,
     $            LDQ, .TRUE., IFIRST+1, 1, 1, JJ-IFIRST, 
     $            1, NBL, 1, 1, WORK, 
     $            D1, D2, D3, D4, D5)
            ENDIF
            CALL QZ2(NNB, A(JJ, JJ), LDA, B(JJ, JJ), LDB, 
     $         WORK(NREF), WORK(RT), NCOLS)
            IF (.NOT. NDONE) THEN
               CALL DLAREF('Rows', A(JJ,JJ+NNB), LDA, .FALSE., Q,
     $            LDQ, .TRUE., 2, 1, 1, NNB-(NCOLS+1),
     $            1, ILASTM-JJ-NNB+1,
     $            1, 1, WORK(NREF), D1, D2, D3, D4, D5)
               CALL DLAREF('Rows', B(JJ,JJ+NNB), LDB, .FALSE., Q,
     $            LDQ, .TRUE., 2, 2, 1, NNB-(NCOLS+1),
     $            1, ILASTM-JJ-NNB+1,
     $            1, 1, WORK(NREF), D1, D2, D3, D4, D5)
            ENDIF

            IF (ILQ) THEN
               CALL DLAREF('Cols', Q(1,JJ), LDQ, .FALSE., Q,
     $            LDQ, .TRUE., 1, 2, 1, NNB-(NCOLS+1), 
     $            1, N, 1, 1, WORK(NREF), 
     $            D1, D2, D3, D4, D5)
            ENDIF
*           
*           Transformations from the right
*           

            IF (JJ-IFRSTM .GT. 0) THEN
               CALL DLAREF('Cols', A(IFRSTM,JJ), LDA, .FALSE., Q,
     $            LDQ, .TRUE., 1, 2, 1, NNB-(NCOLS+1), 1, 
     $            JJ-IFRSTM,
     $            1, 1, WORK(RT), D1, D2, D3, D4, D5)
               CALL DLAREF('Cols', B(IFRSTM,JJ), LDB, .FALSE., Q,
     $            LDQ, .TRUE., 1, 2, 1, NNB-(NCOLS+1), 1, JJ-
     $            IFRSTM,
     $            1, 1, WORK(RT), D1, D2, D3, D4, D5)
            ENDIF

            IF (ILZ) THEN
               CALL DLAREF('Cols', Z(1,JJ), LDZ, .FALSE., Z,
     $            LDQ, .TRUE., 1, 2, 1, NNB-(NCOLS+1), 
     $            1, N,
     $            1, 1, WORK(RT), D1, D2, D3, D4, D5)
            ENDIF
            NREF = NREF + (NNB-(NCOLS+1)) * 3
            JJ = JJ + NB - 4
            IF (NDONE)  GOTO 1000
*           
*           New Blocked Code Ends Here
*           
*           
*           Last elements: Use Givens rotations
*           
*           Rotations from the left
*           
            J = ILAST - 1
            TEMP = A( J, J-1 )
            CALL DLARTG( TEMP, A( J+1, J-1 ), C, S, A( J, J-1 ) )
            A( J+1, J-1 ) = ZERO
*           
            DO 300 JC = J, ILASTM
               TEMP = C*A( J, JC ) + S*A( J+1, JC )
               A( J+1, JC ) = -S*A( J, JC ) + C*A( J+1, JC )
               A( J, JC ) = TEMP
               TEMP2 = C*B( J, JC ) + S*B( J+1, JC )
               B( J+1, JC ) = -S*B( J, JC ) + C*B( J+1, JC )
               B( J, JC ) = TEMP2
 300        CONTINUE
            IF( ILQ ) THEN
               DO 310 JR = 1, N
                  TEMP = C*Q( JR, J ) + S*Q( JR, J+1 )
                  Q( JR, J+1 ) = -S*Q( JR, J ) + C*Q( JR, J+1 )
                  Q( JR, J ) = TEMP
 310           CONTINUE
            END IF
*           
*           Rotations from the right.
*           
            TEMP = B( J+1, J+1 )
            CALL DLARTG( TEMP, B( J+1, J ), C, S, B( J+1, J+1 ) )
            B( J+1, J ) = ZERO
*           
            DO 320 JR = IFRSTM, ILAST
               TEMP = C*A( JR, J+1 ) + S*A( JR, J )
               A( JR, J ) = -S*A( JR, J+1 ) + C*A( JR, J )
               A( JR, J+1 ) = TEMP
 320        CONTINUE
            DO 330 JR = IFRSTM, ILAST - 1
               TEMP = C*B( JR, J+1 ) + S*B( JR, J )
               B( JR, J ) = -S*B( JR, J+1 ) + C*B( JR, J )
               B( JR, J+1 ) = TEMP
 330        CONTINUE
            IF( ILZ ) THEN
               DO 340 JR = 1, N
                  TEMP = C*Z( JR, J+1 ) + S*Z( JR, J )
                  Z( JR, J ) = -S*Z( JR, J+1 ) + C*Z( JR, J )
                  Z( JR, J+1 ) = TEMP
 340           CONTINUE
            END IF
            TIME1 = TIME1 + DSECND() - TM2 
*           
*           End of Double-Shift code
*           
         END IF
*        
         GO TO 350
*        
*        End of iteration loop
*        
 350     CONTINUE
 360  CONTINUE
*     
*     Drop-through = non-convergence
*     

      INFO = ILAST
      GO TO 420
*     
*     Successful completion of all QZ steps
*     
 380  CONTINUE
*     
*     Set Eigenvalues 1:ILO-1
*     
C     DO 410 J = 1, ILO - 1
C     IF( B( J, J ).LT.ZERO ) THEN
C     IF( ILSCHR ) THEN
C     DO 390 JR = 1, J
C     A( JR, J ) = -A( JR, J )
C     B( JR, J ) = -B( JR, J )
C     390          CONTINUE
C     ELSE
C     A( J, J ) = -A( J, J )
C     B( J, J ) = -B( J, J )
C     END IF
C     IF( ILZ ) THEN
C     DO 400 JR = 1, N
C     Z( JR, J ) = -Z( JR, J )
C     400          CONTINUE
C     END IF
C     END IF
C     ALPHAR( J ) = A( J, J )
C     ALPHAI( J ) = ZERO
C     BETA( J ) = B( J, J )
C     410 CONTINUE
*     
*     Normal Termination
*     
      INFO = 0
      
*     
*     Exit (other than argument error) -- return optimal workspace size
*     
 420  CONTINUE
      WORK( 1 ) = DBLE( N )
      RETURN
*     
*     End of DHGEQZ
*     
      END



      SUBROUTINE DLAREF( TYPE, A, LDA, WANTZ, Z, LDZ, BLOCK, IROW1,
     $                    ICOL1, ISTART, ISTOP, ITMP1, ITMP2, LILOZ,
     $                    LIHIZ, VECS, V2, V3, T1, T2, T3 )
      IMPLICIT NONE
*

*  -- ScaLAPACK routine (version 1.4 ALPHA) --
*     University of Tennessee, Knoxville, Oak Ridge National Laboratory,
*     and University of California, Berkeley.
*     November 17, 1996
*
*     .. Scalar Arguments ..
      LOGICAL            BLOCK, WANTZ
      CHARACTER          TYPE
      INTEGER            ICOL1, IROW1, ISTART, ISTOP, ITMP1, ITMP2, LDA,
     $                   LDZ, LIHIZ, LILOZ
      DOUBLE PRECISION   T1, T2, T3, V2, V3
*     ..
*     .. Array Arguments ..
      DOUBLE PRECISION   A( LDA, * ), VECS( * ), Z( LDZ, * )
*     ..
*
*  Purpose
*  =======
*
*  DLAREF applies one or several Householder reflectors of size 3
*     to one or two matrices (if column is specified) on either their
*     rows or columns.
*
*  Arguments
*  =========
*
*  TYPE    (global input) CHARACTER*1
*          If 'R': Apply reflectors to the rows of the matrix
*              (apply from left)
*          Otherwise: Apply reflectors to the columns of the matrix
*          Unchanged on exit.
*
*  A       (global input/output) DOUBLE PRECISION array, (LDA,*)
*          On entry, the matrix to receive the reflections.
*          The updated matrix on exit.
*
*  LDA     (local input) INTEGER
*          On entry, the leading dimension of A.  Unchanged on exit.
*
*  WANTZ   (global input) LOGICAL
*          If .TRUE., then apply any column reflections to Z as well.
*          If .FALSE., then do no additional work on Z.
*
*  Z       (global input/output) DOUBLE PRECISION array, (LDZ,*)
*          On entry, the second matrix to receive column reflections.
*          This is changed only if WANTZ is set.
*
*  LDZ     (local input) INTEGER
*          On entry, the leading dimension of Z.  Unchanged on exit.
*
*  BLOCK   (global input) LOGICAL
*          If .TRUE., then apply several reflectors at once and read
*             their data from the VECS array.
*          If .FALSE., apply the single reflector given by V2, V3,
*             T1, T2, and T3.
*
*  IROW1   (local input/output) INTEGER
*          On entry, the local row element of A.
*          Undefined on output.
*
*
*  ICOL1   (local input/output) INTEGER
*          On entry, the local column element of A.
*          Undefined on output.
*
*  ISTART  (global input) INTEGER
*          Specifies the "number" of the first reflector.  This is
*              used as an index into VECS if BLOCK is set.
*              ISTART is ignored if BLOCK is .FALSE..
*
*  ISTOP   (global input) INTEGER
*          Specifies the "number" of the last reflector.  This is
*              used as an index into VECS if BLOCK is set.
*              ISTOP is ignored if BLOCK is .FALSE..
*
*  ITMP1   (local input) INTEGER
*          Starting range into A.  For rows, this is the local
*              first column.  For columns, this is the local first row.
*
*  ITMP2   (local input) INTEGER
*          Ending range into A.  For rows, this is the local last
*              column.  For columns, this is the local last row.
*
*  LILOZ
*  LIHIZ   (local input) INTEGER
*          These serve the same purpose as ITMP1,ITMP2 but for Z
*              when WANTZ is set.
*
*  VECS    (global input) DOUBLE PRECISION array of size 3*N (matrix
*                                                             size)
*          This holds the size 3 reflectors one after another and this
*              is only accessed when BLOCK is .TRUE.
*
*  V2
*  V3
*  T1
*  T2
*  T3      (global input/output) DOUBLE PRECISION
*          This holds information on a single size 3 Householder
*              reflector and is read when BLOCK is .FALSE., and
*              overwritten when BLOCK is .TRUE.
*
*  Implemented by:  G. Henry, November 17, 1996
*
*  =====================================================================
*
*     .. Local Scalars ..
      INTEGER            J, K, IR1, IC1
      DOUBLE PRECISION   SUM
*     ..
*     .. External Functions ..
      LOGICAL            LSAME
      EXTERNAL           LSAME
*     ..
*     .. Intrinsic Functions ..
      INTRINSIC          MOD
*     ..
*     .. Executable Statements ..
*
      IF (( ISTART .GT. ISTOP ) .OR. (ITMP1 .GT. ITMP2 )) THEN
         RETURN
      END IF

      IR1 = IROW1 
      IC1 = ICOL1
      IF( LSAME( TYPE, 'R' ) ) THEN
         IF( BLOCK ) THEN
               DO 30 J = ITMP1, ITMP2
            IR1 = IROW1
            DO 40 K = ISTART, ISTOP
               V2 = VECS( ( K-1 )*3+1 )
               V3 = VECS( ( K-1 )*3+2 )
               T1 = VECS( ( K-1 )*3+3 )
               T2 = T1*V2
               T3 = T1*V3
                  SUM = A( IR1, J ) + V2*A( IR1+1, J ) +
     $                  V3*A( IR1+2, J )
                  A( IR1, J ) = A( IR1, J ) - SUM*T1
                  A( IR1+1, J ) = A( IR1+1, J ) - SUM*T2
                  A( IR1+2, J ) = A( IR1+2, J ) - SUM*T3
               IR1 = IR1 + 1
   40       CONTINUE
   30          CONTINUE
         ELSE
            DO 50 J = ITMP1, ITMP2
               SUM = A( IR1, J ) + V2*A( IR1+1, J ) +
     $               V3*A( IROW1+2, J )
               A( IR1, J ) = A( IROW1, J ) - SUM*T1
               A( IR1+1, J ) = A( IR1+1, J ) - SUM*T2
               A( IR1+2, J ) = A( IR1+2, J ) - SUM*T3
   50       CONTINUE
         END IF
      ELSE
*
*        Do column transforms
*
         IF( BLOCK ) THEN
            DO 110 K = ISTART, ISTOP
               V2 = VECS( ( K-1 )*3+1 )
               V3 = VECS( ( K-1 )*3+2 )
               T1 = VECS( ( K-1 )*3+3 )
               T2 = T1*V2
               T3 = T1*V3
               DO 90 J = ITMP1, ITMP2
                  SUM = A( J, IC1 ) + V2*A( J, IC1+1 ) +
     $                  V3*A( J, IC1+2 )
                  A( J, IC1 ) = A( J, IC1 ) - SUM*T1
                  A( J, IC1+1 ) = A( J, IC1+1 ) - SUM*T2
                  A( J, IC1+2 ) = A( J, IC1+2 ) - SUM*T3
   90          CONTINUE
               IF( WANTZ ) THEN
                  DO 100 J = LILOZ, LIHIZ
                     SUM = Z( J, IC1 ) + V2*Z( J, IC1+1 ) +
     $                     V3*Z( J, IC1+2 )
                     Z( J, IC1 ) = Z( J, IC1 ) - SUM*T1
                     Z( J, IC1+1 ) = Z( J, IC1+1 ) - SUM*T2
                     Z( J, IC1+2 ) = Z( J, IC1+2 ) - SUM*T3
  100             CONTINUE
               END IF
               IC1 = IC1 + 1
  110       CONTINUE
         ELSE
            DO 120 J = ITMP1, ITMP2
               SUM = A( J, IC1 ) + V2*A( J, IC1+1 ) +
     $               V3*A( J, IC1+2 )
               A( J, IC1 ) = A( J, IC1 ) - SUM*T1
               A( J, IC1+1 ) = A( J, IC1+1 ) - SUM*T2
               A( J, IC1+2 ) = A( J, IC1+2 ) - SUM*T3
  120       CONTINUE
         END IF
      END IF
      RETURN
*
*     End of DLAREF
*
      END
*

      SUBROUTINE MYROT( TYPE, A, LDA, IROW1, ICOL1, ISTART,
     $                   ISTOP, ITMP1, ITMP2, CS)
      implicit none
*     .. Scalar Arguments ..
      CHARACTER          TYPE
      INTEGER            ICOL1, IROW1, ISTART, ISTOP, ITMP1, ITMP2, LDA

*     ..
*     .. Array Arguments ..
      DOUBLE PRECISION   A( LDA, * ), CS( * )
*     ..
*
*  Purpose
*  =======
*
*  MYROT applies one or several Givens rotations to a matrix
*     on either their rows or columns.
*
*  Arguments
*  =========
*
*  TYPE    (global input) CHARACTER*1
*          If 'R': Apply reflectors to the rows of the matrix
*              (apply from left)
*          Otherwise: Apply reflectors to the columns of the matrix
*          Unchanged on exit.
*
*  A       (global input/output) DOUBLE PRECISION array, (LDA,*)
*          On entry, the matrix to receive the reflections.
*          The updated matrix on exit.
*
*  LDA     (local input) INTEGER
*          On entry, the leading dimension of A.  Unchanged on exit.
*
*  IROW1   (local input/output) INTEGER
*          On entry, the local row element of A.
*          Undefined on output.
*
*
*  ICOL1   (local input/output) INTEGER
*          On entry, the local column element of A.
*          Undefined on output.
*
*  ISTART  (global input) INTEGER
*          Specifies the "number" of the first cspair.  This is
*              used as an index into VECS.
*
*  ISTOP   (global input) INTEGER
*          Specifies the "number" of the last reflector.
*
*  ITMP1   (local input) INTEGER
*          Starting range into A.  For rows, this is the local
*              first column.  For columns, this is the local first row.
*
*  ITMP2   (local input) INTEGER
*          Ending range into A.  For rows, this is the local last
*              column.  For columns, this is the local last row.
*
*  CS    (global input) DOUBLE PRECISION array of size 2*N (matrix
*                                                             size)
*          This holds the C and S one after another and this
*
*  =====================================================================
*
*     .. Local Scalars ..
      INTEGER            J, K, JC, JR
      DOUBLE PRECISION   TEMP, C, S
*     ..
*     .. External Functions ..
      LOGICAL            LSAME
      EXTERNAL           LSAME
*     ..

*     .. Intrinsic Functions ..
      INTRINSIC          MOD
*     ..
*     .. Executable Statements ..
*
*     
      IF (( ISTART .GT. ISTOP ) .OR. (ITMP1 .GT. ITMP2 )) THEN
         RETURN
      END IF
      IF( LSAME( TYPE, 'R' ) ) THEN
            DO 10 JC = ITMP1, ITMP2
         J = IROW1
         DO 20 K = ISTART, ISTOP
            C = CS( ( K-1 )*2+1 )
            S = CS( ( K-1 )*2+2 )
               TEMP = C*A( J, JC ) + S*A( J+1, JC )
               A( J+1, JC ) = -S*A( J, JC ) + C*A( J+1, JC )
               A( J, JC ) = TEMP
            J = J + 1
20       CONTINUE
10          CONTINUE
      ENDIF
      IF( LSAME( TYPE, 'C' ) ) THEN
         J = ICOL1
         DO 40 K = ISTART, ISTOP
            C = CS( ( K-1 )*2+1 )
            S = CS( ( K-1 )*2+2 )
            DO 30 JR = ITMP1, ITMP2
               TEMP = C*A( JR, J+1 ) + S*A( JR, J )
               A( JR, J ) = -S*A( JR, J+1 ) + C*A( JR, J )
               A( JR, J+1 ) = TEMP
30          CONTINUE
            J = J + 1
40       CONTINUE
      ENDIF
      IF( LSAME( TYPE, 'X' ) ) THEN
        J = ICOL1
        DO 80 K = ISTART, ISTOP        
            C = CS( ( K-1 )*2+1 )
            S = CS( ( K-1 )*2+2 )
            DO 90 JR = ITMP1, ITMP2
               TEMP = C*A( JR, J ) + S*A( JR, J+1 )
               A( JR, J+1 ) = -S*A( JR, J ) + C*A( JR, J+1 )
               A( JR, J ) = TEMP
90          CONTINUE
            J = J + 1
80      CONTINUE
      END IF
      RETURN
*
*     End of MYROT
*
      END
      
      SUBROUTINE QZ2( N, A, LDA, B, LDB, VL, VR, NCOL)
      IMPLICIT NONE
*     
*     
*     .. Scalar Arguments ..
      INTEGER            LDA, LDB, N, NCOL
*     ..
*     .. Array Arguments ..
      DOUBLE PRECISION   A( LDA, * ), B( LDB, * ), VL( * ), VR( * )
*     ..
*     
*     .. Parameters ..
      DOUBLE PRECISION   ZERO, ONE
      PARAMETER          ( ZERO = 0.0D+0, ONE = 1.0D+0 )
*     ..
*     .. Local Scalars ..
      LOGICAL            ILPIVT
      INTEGER            J, JC, JR, JMP
      DOUBLE PRECISION   SCALE, T, VS, U2, SAFMIN, SAFMAX, TAU, TEMP, 
     $   TEMP2, U1, W11, W12, W21, W22

*     ..
*     .. Local Arrays ..
      DOUBLE PRECISION   V( 3 )
*     ..
*     .. External Functions ..
      LOGICAL            LSAME
      DOUBLE PRECISION   DLAMCH, DLANHS, DLAPY2, DLAPY3
      EXTERNAL           LSAME, DLAMCH, DLANHS, DLAPY2, DLAPY3

*     ..
*     .. Intrinsic Functions ..
      INTRINSIC          ABS, DBLE, MAX, MIN, SQRT
*     ..
*     .. Executable Statements ..
*     
*     Machine Constants
*     
      SAFMIN = DLAMCH( 'S' )
      SAFMAX = ONE / SAFMIN
      JMP = 0
      DO 290 J = 2, N - NCOL 
*        
*        Use 3x3 Householder transforms.
*        
*        Zero (j-1)st column of A
*        
         V( 1 ) = A( J, J-1 )
         V( 2 ) = A( J+1, J-1 )
         V( 3 ) = A( J+2, J-1 )
*        
         CALL DLARFG( 3, A( J, J-1 ), V( 2 ), 1, TAU )
         V( 1 ) = ONE
         A( J+1, J-1 ) = ZERO
         A( J+2, J-1 ) = ZERO
*        
*        Save Householder transforms.
*        
         VL( JMP + 1) = V( 2 )
         VL( JMP + 2) = V( 3 )
         VL( JMP + 3) = TAU 
*        
         DO 230 JC = J, N
            TEMP = TAU*( A( J, JC )+V( 2 )*A( J+1, JC )+V( 3 )*
     $         A( J+2, JC ) )
            A( J, JC ) = A( J, JC ) - TEMP
            A( J+1, JC ) = A( J+1, JC ) - TEMP*V( 2 )
            A( J+2, JC ) = A( J+2, JC ) - TEMP*V( 3 )
            TEMP2 = TAU*( B( J, JC )+V( 2 )*B( J+1, JC )+V( 3 )*
     $         B( J+2, JC ) )
            B( J, JC ) = B( J, JC ) - TEMP2
            B( J+1, JC ) = B( J+1, JC ) - TEMP2*V( 2 )
            B( J+2, JC ) = B( J+2, JC ) - TEMP2*V( 3 )
 230     CONTINUE
*        
*        Zero j-th column of B (see DLAGBC for details)
*        
*        Swap rows to pivot
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
*        Swap columns if nec.
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
*        LU-factor
*        
         TEMP = W21 / W11
         U2 = U2 - TEMP*U1
         W22 = W22 - TEMP*W12
         W21 = ZERO
*        
*        Compute SCALE
*        
         SCALE = ONE
         IF( ABS( W22 ).LT.SAFMIN ) THEN
            SCALE = ZERO
            U2 = ONE
            U1 = -W12 / W11
            GO TO 250
         END IF
         IF( ABS( W22 ).LT.ABS( U2 ) )
     $      SCALE = ABS( W22 / U2 )
         IF( ABS( W11 ).LT.ABS( U1 ) )
     $      SCALE = MIN( SCALE, ABS( W11 / U1 ) )
*        
*        Solve
*        
         U2 = ( SCALE*U2 ) / W22
         U1 = ( SCALE*U1-W12*U2 ) / W11
*        
 250     CONTINUE
         IF( ILPIVT ) THEN
            TEMP = U2
            U2 = U1
            U1 = TEMP
         END IF
*        
*        Compute Householder Vector
*        
         T = SQRT( SCALE**2+U1**2+U2**2 )
         TAU = ONE + SCALE / T
         VS = -ONE / ( SCALE+T )
         V( 1 ) = ONE
         V( 2 ) = VS*U1
         V( 3 ) = VS*U2
*        
*        Save Householder transforms.
*        
         VR( JMP + 1) = V( 2 )
         VR( JMP + 2) = V( 3 )
         VR( JMP + 3) = TAU 
*        
*        Apply transformations from the right.
*        
         DO 260 JR = 1, MIN( J+3, N )
            TEMP = TAU*( A( JR, J )+V( 2 )*A( JR, J+1 )+V( 3 )*
     $         A( JR, J+2 ) )
            A( JR, J ) = A( JR, J ) - TEMP
            A( JR, J+1 ) = A( JR, J+1 ) - TEMP*V( 2 )
            A( JR, J+2 ) = A( JR, J+2 ) - TEMP*V( 3 )
 260     CONTINUE
         DO 270 JR = 1, J + 2
            TEMP = TAU*( B( JR, J )+V( 2 )*B( JR, J+1 )+V( 3 )*
     $         B( JR, J+2 ) )
            B( JR, J ) = B( JR, J ) - TEMP
            B( JR, J+1 ) = B( JR, J+1 ) - TEMP*V( 2 )
            B( JR, J+2 ) = B( JR, J+2 ) - TEMP*V( 3 )
 270     CONTINUE
         B( J+1, J ) = ZERO
         B( J+2, J ) = ZERO
         JMP = JMP + 3
 290  CONTINUE
      RETURN
*     
*     End of QZ2
*     
      END

      SUBROUTINE GIV(N, A, LDA, B, LDB, LCS, RCS, NCOLS)
      IMPLICIT NONE
      INTEGER N, LDA, LDB, NCOLS
      DOUBLE PRECISION A(LDA, *), B(LDB, *), LCS( * ), RCS( * )
*
*     Moves a bulge using Givens Rotations
*
      INTEGER            J, JC, JR, I 
      DOUBLE PRECISION   TEMP, TEMP2, ZERO
      PARAMETER          ( ZERO = 0.0D+0 )

      I = 1
      DO 190 J = 2, N - NCOLS
         TEMP = A( J, J-1 )
         CALL DLARTG( TEMP, A( J+1, J-1 ), LCS(I), LCS(I+1), A(J,J-1))
         A( J+1, J-1 ) = ZERO
*
         DO 140 JC = J, N
            TEMP = LCS(I)*A( J, JC ) + LCS(I+1)*A( J+1, JC )
            A( J+1, JC ) = -LCS(I+1)*A( J, JC ) + LCS(I)*A( J+1, JC )
            A( J, JC ) = TEMP
            TEMP2 = LCS(I)*B( J, JC ) + LCS(I+1)*B( J+1, JC )
            B( J+1, JC ) = -LCS(I+1)*B( J, JC ) + LCS(I)*B( J+1, JC )
            B( J, JC ) = TEMP2
  140    CONTINUE
*
         TEMP = B( J+1, J+1 )
         CALL DLARTG( TEMP, B( J+1,J), RCS(I), RCS(I+1), B( J+1, J+1))
         B( J+1, J ) = ZERO
*
         DO 160 JR = 1, MIN( J+2, N )
            TEMP = RCS(I)*A( JR, J+1 ) + RCS(I+1)*A( JR, J )
            A( JR, J ) = -RCS(I+1)*A( JR, J+1 ) + RCS(I)*A( JR, J )
            A( JR, J+1 ) = TEMP
  160    CONTINUE
         DO 170 JR = 1, J
            TEMP = RCS(I)*B( JR, J+1 ) + RCS(I+1)*B( JR, J )
            B( JR, J ) = -RCS(I+1)*B( JR, J+1 ) + RCS(I)*B( JR, J )
            B( JR, J+1 ) = TEMP
  170    CONTINUE
         I = I + 2
  190 CONTINUE
      RETURN
      END
