      SUBROUTINE QZEARLY( N, S, HNORM, H, LDH, T, LDT, KBOT, SR, SI,
     $   SBETA, Q, LDQ, Z, LDZ, DWORK, LDWORK, INFO )
      IMPLICIT NONE
C     
C     PURPOSE
C     
C     To perform aggressive early deflation on the n-by-n
C     Hessenberg-triangular pencil (H,T). This routine computes
C     a possibly reduced Hessenberg-triangular pencil (HR,TR) and
C     orthogonal matrices Q, Z so that
C     
C     (1) H*Z = Q*HR,  T*Z = Q*TR,
C     
C     (2) the trailing principal subpencil
C     ( HR(KBOT:N,KBOT:N), TR(KBOT:N,KBOT:N) )
C     is in generalized Schur form, normalized as
C     returned by the LAPACK routine DHGEQZ,
C     
C     (3) the components of S*Q(1,KBOT:N,1) are small enough
C     to be neglected (that also depends on the parameters
C     in PARA, see TTQZ), and if KBOT is returned with a
C     value of zero, then S*Q(1,1) is also negligible.
C     
C     The procedure is designed for deflating eigenvalues during the QZ
C     iteration.
C     
C     Things NOT done:
C     - Implementation of different deflation strategies
C     - use of modified QZ for small-sized problem
C     - proper LDWORK checking
C     
C     ARGUMENTS
C     
C     Input/Output Parameters
C     
C     N       (input) INTEGER
C     The order of the matrices H and T.  N >= 0.
C     
C     S       (input) DOUBLE PRECISION
C     Scaling factor.
C     
C     HNORM   (input) DOUBLE PRECISION
C     Frobenius norm of H.
C     
C     H       (input/output) DOUBLE PRECISION array, dimension (LDH,N)
C     On entry, the leading N-by-N part of this array must
C     contain the upper Hessenberg matrix H.
C     On exit, the leading N-by-N part of this array contains
C     the (possibly reduced) upper Hessenberg matrix HR.
C     
C     LDH     (input) INTEGER
C     The leading dimension of the array H. LDH >= max(1,N).
C     
C     T       (input/output) DOUBLE PRECISION array, dimension (LDT,N)
C     On entry, the leading N-by-N part of this array must
C     contain the upper triangular matrix T.
C     On exit, the leading N-by-N part of this array contains
C     the upper triangular matrix TR.
C     
C     LDT     (input) INTEGER
C     The leading dimension of the array T. LDT >= max(1,N).
C     
C     KBOT    (output) INTEGER
C     On exit, KBOT specifies that the pencil
C     (HR(KBOT:N,KBOT:N),TR(KBOT:N,KBOT:N)) is in generalized
C     Schur form. If the entire vector S*Q(1,:) is negligible,
C     then KBOT is returned with a value of zero.
C     
C     SR      (output) DOUBLE PRECISION array, dimension (N)
C     SI      (output) DOUBLE PRECISION array, dimension (N)
C     SBETA   (output) DOUBLE PRECISION array, dimension (N)
C     On exit, the leading N elements of the arrays SR, SI and
C     SBETA contain the eigenvalues of the pencil (HR,TR) as
C     returned by the LAPACK routine DHGEQZ. The last N-KBOT+1
C     elements correspond to the eigenvalues of the subpencil
C     (HR(KBOT:N,KBOT:N),TR(KBOT:N,KBOT:N)).
C     
C     Q       (input/output) DOUBLE PRECISION array, dimension (LDQ,N)
C     On entry, the leading N-by-N part of this array must
C     contain an orthogonal matrix Q.
C     On exit, the leading N-by-N part of this array contains
C     the matrix Q post-multiplied by the transformations which
C     are applied to H and T on the left.
C     
C     LDQ     (input) INTEGER
C     The leading dimension of the array Q.  LDQ >= MAX(1,N).
C     
C     Z       (input/output) DOUBLE PRECISION array, dimension (LDZ,N)
C     On entry, the leading N-by-N part of this array must
C     contain an orthogonal matrix Z.
C     On exit, the leading N-by-N part of this array contains
C     the matrix Z post-multiplied by the transformations which
C     are applied to H and T on the right.
C     
C     LDZ     (input) INTEGER
C     The leading dimension of the array Z.  LDZ >= MAX(1,N).
C     
C     Workspace
C     
C     DWORK   DOUBLE PRECISION array, dimension (LDWORK)
C     On exit, if INFO = 0, DWORK(1) returns the optimal
C     value of LDWORK, ????.
C     On exit, if  INFO = -16, DWORK(1) returns the minimum
C     value of LDWORK.
C     
C     LDWORK  INTEGER
C     The length of the array DWORK.  LDWORK >= ???.
C     
C     Error Indicator
C     
C     INFO    INTEGER
C     = 0:  successful exit;
C     < 0:  if INFO = -i, the i-th argument had an illegal 
C     value;
C     > 0:  the QZ algorithm failed to compute the first INFO
C     eigenvalues of (H,T). KBOT is set to N and the
C     not converged eigenvalues are replaced by the
C     diagonal elements of H and T.
C     
C     CONTRIBUTOR
C     
C     D. Kressner, Univ. Umea, Sweden, January 2005.
C     
C     ******************************************************************
C     
C     .. Parameters ..
      DOUBLE PRECISION  ZERO, ONE, TWO
      PARAMETER         ( ZERO = 0.0D0, ONE = 1.0D0, TWO = 2.0D0 )
      INTEGER           PARA_AGGCRIT, PARA_AGGSHF, PARA_AGGMIN,
     $   PARA_AGGWIN, PARA_HDEF
      
      LOGICAL           PARA_AGGRQUP
      DOUBLE PRECISION  PARA_AGGREP   
      
      COMMON /PARA_DEF/ PARA_AGGCRIT,
     $   PARA_AGGSHF, PARA_AGGMIN, PARA_AGGWIN, PARA_HDEF, PARA_AGGRQUP,
     $   PARA_AGGREP  


      LOGICAL           PARA_ITTRI
      INTEGER           PARA_ITASH, PARA_ITBULG, PARA_ITDACK, 
     $   PARA_ITDANB, PARA_ITEXC, PARA_ITEXT, PARA_ITKK, PARA_ITMAX, 
     $   PARA_ITNS, PARA_ITOPP
      COMMON /PARA_IT/  PARA_ITASH, PARA_ITBULG, PARA_ITDACK, 
     $   PARA_ITDANB, PARA_ITEXC, PARA_ITEXT, PARA_ITKK, PARA_ITMAX, 
     $   PARA_ITNS, PARA_ITOPP, PARA_ITTRI
      
      INTEGER           MSGLVL
      COMMON /MESSAGES/ MSGLVL
C     
C     .. Scalar Arguments ..
      INTEGER           INFO, KBOT, LDH, LDQ, LDT, LDWORK, LDZ, N
      DOUBLE PRECISION  HNORM,S


C     .. Array Arguments ..
      DOUBLE PRECISION  DWORK(LDWORK), H(LDH,*), Q(LDQ,*), SBETA(*),
     $   SI(*), SR(*), T(LDT,*), Z(LDZ,*)
C     .. Local Scalars ..
      DOUBLE PRECISION  ALPHA, RTDET, TAU, TB1,
     $   TB2, TEMP, TI1, TI2, TR1, TR2
      INTEGER           I, IERR, IFST, ILST, KNT, WRKMIN
      INTEGER           OILST, OIFST
      LOGICAL           BULGE, DFLATE

      DOUBLE PRECISION    SMLNUM, ULP, TFNORM, T1NORM, HFNORM, HTOL,
     $   TTOL
      COMMON /PREC/       SMLNUM, ULP, TFNORM, T1NORM, HFNORM, HTOL, 
     $   TTOL

C     .. External Functions ..
      DOUBLE PRECISION  DLAMCH, DLANHS, DLAPY2
      EXTERNAL          DLAMCH, DLANHS, DLAPY2
C     .. External Subroutines ..
      EXTERNAL          DAXPY, DCOPY, DGGHRD, DHGEQZ, DLABAD, DLARF,
     $   DLARFG, DLARTG, DLASET, DROT, XERBLA, DTGEX2, DETGEXC
C     .. Intrinsic Functions ..
      INTRINSIC         ABS, DBLE, INT, MAX, SQRT
C     
C     .. Executable Statements ..
C     
C     Decode and check the scalar input parameters.
C     
      INFO = 0
      WRKMIN = 4*N + 16
C     
      IF ( N.LT.0 ) THEN
         INFO = -1
      ELSE IF ( LDH.LT.MAX( 1, N ) ) THEN
         INFO = -4
      ELSE IF ( LDT.LT.MAX( 1, N ) ) THEN
         INFO = -6
      ELSE IF ( LDQ.LT.MAX( 1, N ) ) THEN
         INFO = -12
      ELSE IF ( LDZ.LT.MAX( 1, N ) ) THEN
         INFO = -14
      ELSE IF ( LDWORK.LT.WRKMIN ) THEN
         INFO = -16
         DWORK(1) = DBLE( WRKMIN )
      END IF
C     
C     Return if there were illegal values.
C     
      IF ( INFO.NE.0 ) THEN
         CALL XERBLA( 'QZEARLY', -INFO )
         RETURN
      END IF
C     
C     Quick return if N = 0.
C     
      IF ( N.EQ.0 ) THEN
         DWORK(1) = ONE
         RETURN
      END IF
C     
C     Set machine constants for stopping criterion.
C     


C     
C     Convert to (spike-)triangular form.
C     
c     CALL DHGEQZ( 'S',  'I', 'I', N, 1, N,
c     $        H, LDH, T, LDT,
c     $        SR, SI, SBETA, Q, LDQ, Z, LDZ,
c     $        DWORK, LDWORK, INFO)

      IF ( N.GT.PARA_ITDACK ) THEN
         CALL QZDACK( 'Schur form', 'V', 'V', N, 1, N,  H, LDH, T, LDT,
     $      SR, SI, SBETA, Q, LDQ, Z, LDZ, DWORK, LDWORK,
     $      INFO, PARA_ITDANB )
      ELSE
         CALL QZLAP( 'Schur form', 'V', 'V', N, 1, N, H, LDH, T, LDT,
     $      SR, SI, SBETA, Q, LDQ, Z, LDZ, DWORK, LDWORK,
     $      INFO )
      END IF
      IF ( INFO.GT.1 ) THEN
C        
C        Set not converged eigenvalues and leave.
C        
         IF ( MSGLVL.GE.1 ) THEN
            WRITE(*,*) 'WARNING: QZ DID NOT CONVERGE AT A SUBPROBLEM'
         END IF
         DO 10 I = 1, INFO
            SR(I) = H(I,I)
            SI(I) = ZERO
            SBETA(I) = T(I,I)
 10      CONTINUE
         KBOT = N
         RETURN
      END IF
C     
      ILST = 1
      KBOT = N
      KNT = 0
C     
C     Check spike for deflated eigenvalues.
C     while ( knt < n )
C     
 20   CONTINUE
      IF ( KNT.LT.N ) THEN
C        
C        BULGE signals whether we are in a two-by-two block.
C        
         IF ( KBOT.EQ.1 ) THEN
            BULGE = .FALSE.
         ELSE
            BULGE = H(KBOT, KBOT-1).NE.ZERO
         END IF
C        
         DFLATE = .FALSE.
         IF ( .NOT.BULGE ) THEN
            KNT = KNT + 1
            TEMP = ABS( S ) + ABS( H(KBOT,KBOT) )
            DFLATE = ABS( S*Q(1,KBOT)).LE.MAX( ULP*TEMP, SMLNUM )
         ELSE
            KNT = KNT + 2
            IF ( ( H(KBOT-1,KBOT-1 ).EQ.ZERO )  .AND.
     $         ( H(KBOT,KBOT-1).EQ.ZERO ) ) THEN
               RTDET = ZERO
            ELSE IF ( ABS( H(KBOT-1,KBOT-1) ).GE.
     $            ABS( H(KBOT,KBOT-1) ) ) THEN
               RTDET = SQRT( ABS( H(KBOT-1,KBOT-1) ) )*
     $            SQRT( ABS( H(KBOT,KBOT) - ( H(KBOT,KBOT-1)/
     $            H(KBOT-1,KBOT-1) )*H(KBOT-1,KBOT) ) )
            ELSE
               RTDET = SQRT( ABS( H(KBOT,KBOT-1) ) )*
     $            SQRT( ABS( H(KBOT-1,KBOT) -
     $            ( H(KBOT-1,KBOT-1)/
     $            H(KBOT,KBOT-1) )*H(KBOT,KBOT) ) )
            END IF
            TEMP = ABS(S) + RTDET
            DFLATE = MAX( ABS( S*Q(1,KBOT) ), ABS( S*Q(1,KBOT-1) ) )
     $         .LE.MAX( ULP*TEMP, SMLNUM )
         END IF

         IF ( DFLATE ) THEN
C           
C           Deflation.
C           
            IF ( BULGE ) THEN
               KBOT = KBOT - 2
            ELSE
               KBOT = KBOT - 1
            END IF
         ELSE
C           
C           Move undeflatable eigenvalue up.
C           
            IFST = KBOT
            OIFST = IFST
            OILST = ILST


            CALL DTGEXC( .TRUE., .TRUE., N, H, LDH, T, LDT, Q, LDQ, Z,
     $         LDZ, IFST, ILST, DWORK, LDWORK, IERR )
C           
C           Go directly to reduction to Hessenberg-triangular form if 
C           eigenvalue swapping failed. This might not be a good idea! 
C           
            IF ( IERR.NE.0 ) THEN 
c              PRINT*, '% WARNING: SWAPPING FAILED', OIFST, OILST, 
c              $            IFST, ILST,
c              $            BULGE, N, IERR
               GO TO 42 
            END IF 
C           
            
C           
            IF ( .NOT.BULGE ) THEN
               TR1 = SR(KBOT)
               TI1 = SI(KBOT)
               TB1 = SBETA(KBOT)
               DO 30 I = KBOT - 1, ILST, -1
                  SR(I+1) = SR(I)
                  SI(I+1) = SI(I)
                  SBETA(I+1) = SBETA(I)
 30            CONTINUE
               SR(ILST) = TR1
               SI(ILST) = TI1
               SBETA(ILST) = TB1
               ILST   = ILST + 1
            ELSE
               TR1 = SR(KBOT-1)
               TI1 = SI(KBOT-1)
               TB1 = SBETA(KBOT-1)
               TR2 = SR(KBOT)
               TI2 = SI(KBOT)
               TB2 = SBETA(KBOT)
               DO 40 I = KBOT - 2, ILST, -1
                  SR(I+2) = SR(I)
                  SI(I+2) = SI(I)
                  SBETA(I+2) = SBETA(I)
 40            CONTINUE
               SR(ILST) = TR1
               SI(ILST) = TI1
               SBETA(ILST) = TB1
               SR(ILST+1) =  TR2
               SI(ILST+1) =  TI2
               SBETA(ILST+1) = TB2
               ILST = ILST + 2
            END IF
         END IF
C        
         GO TO 20
C        
      END IF
C     
C     Return to Hessenberg form.
C     
 42   CONTINUE
C     
      IF ( KBOT.GT.0 ) THEN
         CALL DCOPY( KBOT, Q, LDQ, DWORK, 1 )
         ALPHA = DWORK(1)
         CALL DLARFG( KBOT, ALPHA, DWORK(2), 1, TAU )
         DWORK(1) = ONE
         CALL DLARF( 'Left', KBOT, N, DWORK, 1, TAU, H, LDH,
     $      DWORK(N+1) )
         CALL DLARF( 'Right', N, KBOT, DWORK, 1, TAU, Q, LDQ,
     $      DWORK(N+1) )
C        
         IF ( KBOT.GT.1 .AND. .NOT.PARA_AGGRQUP ) THEN
            CALL DLARF( 'Left', KBOT, N, DWORK(1), 1, TAU, T, LDT,
     $         DWORK(N+1) )
            CALL DGERQF( KBOT, KBOT, T, LDT, DWORK(1), DWORK(KBOT+1),
     $         LDWORK-KBOT, IERR )
            CALL DORMRQ( 'Right', 'Transpose', KBOT, KBOT, KBOT, T, LDT,
     $         DWORK(1), H, LDH, DWORK(KBOT+1), LDWORK-KBOT,
     $         IERR )
            CALL DORMRQ( 'Right', 'Transpose', N, KBOT, KBOT, T, LDT,
     $         DWORK(1), Z, LDZ, DWORK(KBOT+1), LDWORK-KBOT,
     $         IERR )
            CALL DLASET( 'Lower', KBOT-1, KBOT-1, ZERO, ZERO, T(2,1),
     $         LDT )
         ELSE IF ( KBOT.GT.1 ) THEN
C           
C           The updated matrix (I-tau*u*u') T = T - tau*u*(T'*u)' is a
C           rank-1 perturbation of an upper triangular matrix. Thus,
C           an RQ-update consisting of 2*KBOT-2 Givens rotations can be
C           employed.
C           
            IF ( N.GT.KBOT ) THEN
               CALL DLARF( 'Left', KBOT, N-KBOT, DWORK(1), 1, TAU,
     $            T(1,KBOT+1), LDT, DWORK(N+1) )
            END IF
            CALL DGEMV( 'Transpose', KBOT, KBOT, -TAU, T, LDT, DWORK, 1,
     $         ZERO, DWORK(KBOT+1), 1 )
C           
            DO 17 I = 1, KBOT-1
               ALPHA = DWORK(KBOT+I+1)
               CALL DLARTG( ALPHA, DWORK(KBOT+I), DWORK(2*KBOT+I),
     $            DWORK(3*KBOT+I), DWORK(KBOT+I+1) )
               DWORK(3*KBOT+I) = -DWORK(3*KBOT+I)
               CALL DROT( I+1, T(1,I), 1, T(1,I+1), 1, DWORK(2*KBOT+I),
     $            DWORK(3*KBOT+I) )
 17         CONTINUE
            CALL DAXPY( KBOT, DWORK(2*KBOT), DWORK, 1, T(1,KBOT), 1 )
            DO 18 I = 1, KBOT-1
               CALL DROT( KBOT, H(1,I), 1, H(1,I+1), 1, DWORK(2*KBOT+I),
     $            DWORK(3*KBOT+I) )
 18         CONTINUE
            DO 19 I = 1, KBOT-1
               CALL DROT( N, Z(1,I), 1, Z(1,I+1), 1, DWORK(2*KBOT+I),
     $            DWORK(3*KBOT+I) )
 19         CONTINUE
C           
            DO 28 I = KBOT-1, 1, -1
               ALPHA = T(I+1,I+1)
               CALL DLARTG( ALPHA, T(I+1,I), DWORK(I), DWORK(KBOT+I),
     $            T(I+1,I+1) )
               DWORK(KBOT+I) = -DWORK(KBOT+I)
               T(I+1,I) = ZERO
               CALL DROT( I, T(1,I), 1, T(1,I+1), 1, DWORK(I),
     $            DWORK(KBOT+I) )
 28         CONTINUE
            DO 21 I = KBOT-1, 1, -1
               CALL DROT( KBOT, H(1,I), 1, H(1,I+1), 1, DWORK(I),
     $            DWORK(KBOT+I) )
 21         CONTINUE
            DO 22 I = KBOT-1, 1, -1
               CALL DROT( N, Z(1,I), 1, Z(1,I+1), 1, DWORK(I),
     $            DWORK(KBOT+I) )
 22         CONTINUE
         END IF
C        
C        Hessenberg-triangular reduction. 
C        
         CALL DGGHRD( 'V', 'V', N, 1, KBOT, H, LDH, T, LDT, Q,
     $      LDQ, Z, LDZ, IERR ) 
      END IF
C     
C     Last line of QZEARLY
C     
      END

