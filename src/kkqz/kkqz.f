      SUBROUTINE KKQZ( JOB, COMPQ, COMPZ, N, ILO, IHI, H, LDH, T, LDT,
     $   ALPHAR, ALPHAI, BETA, Q, LDQ, Z, LDZ, DWORK,
     $   LDWORK, IWORK, LIWORK, INFO )
      IMPLICIT NONE
C     
C     This is the final test implementation of the QZ algorithm as
C     described in [Kagstrom and Kressner''05]. It differs in some
C     aspects from what would be called a proper LAPACK-like
C     implementation:
C     
C     - too much input checking (H and T are tested for being upper
C     Hessenberg/triangular, respectively);
C     - no proper error handling (an error message is printed and
C     INFO = 1 is returned regardless of the nature of the error);
C     - too many subroutines, too many parameters;
C     - parameters must be set via the subroutine QZCONF;
C     - some global parameters are used for convenience;
C     - timing;
C     - SLICOT-like documentation;
C     - IWORK is required because of the block algorithm for deflating
C     infinite eigenvalues, this part can possibly be removed.
C     
C     ARGUMENTS
C     
C     Input/Output Parameters
C     
C     JOB     (input) CHARACTER*1
C     = 'E': compute only ALPHAR, ALPHAI, and BETA.  A and B will
C     not necessarily be put into generalized Schur form.
C     = 'S': put A and B into generalized Schur form, as well
C     as computing ALPHAR, ALPHAI, and BETA.
C     
C     COMPQ   (input) CHARACTER*1
C     = 'N': do not modify Q.
C     = 'V': multiply the array Q on the right by the transpose
C     of the orthogonal tranformation that is applied to
C     the left side of H and T to reduce them to Schur
C     form.
C     = 'I': like COMPQ='V', except that Q will be initialized
C     to the identity first.
C     
C     COMPZ   (input) CHARACTER*1
C     = 'N': do not modify Z.
C     = 'V': multiply the array Z on the right by the orthogonal
C     tranformation that is applied to the right side of
C     H and T to reduce them to Schur form.
C     = 'I': like COMPZ='V', except that Z will be initialized
C     to the identity first.
C     
C     N       (input) INTEGER
C     The order of the matrices H and T.  N >= 0.
C     
C     ILO     (input) INTEGER
C     IHI     (input) INTEGER
C     It is assumed that H and T are already upper triangular in
C     rows and columns 1:ILO-1 and IHI+1:N. ILO and IHI are
C     normally set by a previous call to DGGBAL. Otherwise ILO
C     and IHI should be set to 1 and N, respectively.
C     1 <= ILO <= IHI <= N, if N > 0; ILO=1 and IHI=0, if N=0. 
C     
C     H       (input/output) DOUBLE PRECISION array, dimension (LDH,N)
C     On entry, the leading n-by-n part of this array contains
C     the upper Hessenberg matrix H.
C     On exit, the leading n-by-n part of this array contains
C     the upper quasi-triangular matrix H from the generalized
C     Schur decomposition.
C     
C     LDH     (input) INTEGER
C     The leading dimension of the array H. LDH >= max(1,N).
C     
C     T       (input/output) DOUBLE PRECISION array, dimension (LDT,N)
C     On entry, the leading n-by-n part of this array contains
C     the upper triangular matrix T.
C     On exit, the leading n-by-n part of this array contains
C     the upper triangular matrix T from the generalized Schur
C     decomposition. 2-by-2 blocks in T corresponding to 2-by-2
C     blocks in H will be reduced to positive diagonal form. 
C     (I.e., if H(j+1,j) is non-zero, then T(j+1,j)=T(j,j+1)=0
C     and T(j,j) and T(j+1,j+1) will be positive.) 
C     
C     LDT     (input) INTEGER
C     The leading dimension of the array T. LDT >= max(1,N).
C     
C     ALPHAR  (output) DOUBLE PRECISION array, dimension (N)
C     ALPHAR(1:N) will be set to the real parts of the diagonal
C     elements of H that would result from reducing H and T to
C     Schur form and then further reducing them both to
C     triangular form using unitary transformations s.t. the
C     diagonal of T was non-negative real.  Thus, if H(j,j) is
C     in a 1-by-1 block, then ALPHAR(j)=H(j,j).
C     Note that the (real or complex) values
C     (ALPHAR(j) + i*ALPHAI(j))/BETA(j), j=1,...,N, are the
C     generalized eigenvalues of the matrix pencil H - wT.
C     
C     ALPHAI  (output) DOUBLE PRECISION array, dimension (N)
C     ALPHAI(1:N) will be set to imaginary parts of the diagonal
C     elements of H that would result from reducing H and T to
C     Schur form and then further reducing them both to
C     triangular form using unitary transformations s.t. the
C     diagonal of T was non-negative real.  
C     
C     BETA    (output) DOUBLE PRECISION array, dimension (N)
C     BETA(1:N) will be set to the (real) diagonal elements of
C     T that would result from reducing H and T to Schur form
C     and then further reducing them both to triangular form
C     using unitary transformations s.t. the diagonal of T was
C     non-negative real.
C     
C     Q       (input/output) DOUBLE PRECISION array, dimension (LDQ,N)
C     On entry, the leading N-by-N part of this array must
C     contain an orthogonal matrix Q.
C     On exit, the leading N-by-N part of this array contains
C     the updated orthogonal matrix Q.
C     
C     LDQ     (input) INTEGER
C     The leading dimension of the array Q. LDQ >= max(1,N).
C     
C     Z       (input/output) DOUBLE PRECISION array, dimension (LDZ,N)
C     On entry, the leading N-by-N part of this array must
C     contain an orthogonal matrix Z.
C     On exit, the leading N-by-N part of this array contains
C     the updated orthogonal matrix Z.
C     
C     LDZ     (input) INTEGER
C     The leading dimension of the array Z. LDZ >= max(1,N).
C     
C     Workspace
C     
C     DWORK   DOUBLE PRECISION array, dimension (LDWORK)
C     On exit, if INFO = 0,  DWORK(1)  returns the optimal
C     value of LDWORK, ????.
C     On exit, if  INFO = -??,  DWORK(1)  returns the minimum
C     value of LDWORK.
C     
C     LDWORK  INTEGER
C     The length of the array DWORK.  LDWORK >= ???.
C     
C     IWORK   DOUBLE PRECISION array, dimension (LIWORK)
C     On exit, if  INFO = -??,  IWORK(1)  returns the minimum
C     value of LIWORK.
C     
C     LDWORK  INTEGER
C     The length of the array IWORK.  LIWORK >= ???.
C     
C     Error Indicator
C     
C     INFO    INTEGER
C     = 0:  successful exit;
C     < 0:  if INFO = -i, the i-th argument had an illegal 
C     value.
C     > 0:  if INFO = i, KKQZ failed to compute all of the
C     eigenvalues in a total of 30*n iterations;
C     elements i+1:n of ALPHAR, ALPHAI and BETA contain
C     those eigenvalues which have been successfully
C     computed.
C     
C     CONTRIBUTORS
C     
C     B. Kagstrom and D. Kressner, Univ. Umea, Sweden, February 2005.
C     
C     ******************************************************************
C     
C     .. Parameters ..
      DOUBLE PRECISION  ZERO, ONE, TWO
      PARAMETER         ( ZERO = 0.0D0, ONE = 1.0D0, TWO = 2.0D0 )
      INTEGER           BULGMX
      PARAMETER         ( BULGMX = 12 )
      INTEGER           PARA_LDADD
      COMMON /PARA_MIS/ PARA_LDADD
      LOGICAL           PARA_ITTRI
      INTEGER           PARA_ITASH, PARA_ITBULG, PARA_ITDACK, 
     $   PARA_ITDANB, PARA_ITEXC, PARA_ITEXT, PARA_ITKK, PARA_ITMAX, 
     $   PARA_ITNS, PARA_ITOPP
      COMMON /PARA_IT/  PARA_ITASH, PARA_ITBULG, PARA_ITDACK, 
     $   PARA_ITDANB, PARA_ITEXC, PARA_ITEXT, PARA_ITKK, PARA_ITMAX, 
     $   PARA_ITNS, PARA_ITOPP, PARA_ITTRI
      INTEGER           PARA_AGGCRIT, PARA_AGGSHF, PARA_AGGMIN,
     $   PARA_AGGWIN, PARA_HDEF
      
      LOGICAL           PARA_AGGRQUP
      DOUBLE PRECISION  PARA_AGGREP   
      
      COMMON /PARA_DEF/ PARA_AGGCRIT,
     $   PARA_AGGSHF, PARA_AGGMIN, PARA_AGGWIN, PARA_HDEF, PARA_AGGRQUP,
     $   PARA_AGGREP  


C     .. Scalar Arguments ..
      CHARACTER         COMPQ, COMPZ, JOB
      INTEGER           IHI, ILO, INFO, LDH, LDQ, LDT, LDWORK, LIWORK,
     $   LDZ, N
C     .. Array Arguments ..
      INTEGER           IWORK(*)
      DOUBLE PRECISION  ALPHAI(*), ALPHAR(*), BETA(*), DWORK(*),
     $   H(LDH,*), Q(LDQ,*), T(LDT,*), Z(LDZ,*)
C     .. Local Scalars ..
      LOGICAL           DOAGGR,
     $   LCOMPQ, LCOMPZ, LEXC, LINITQ, LINITZ, LSCHUR
      CHARACTER         UPDQ, UPDZ
      INTEGER           ACTHI, ACTLO, AGGSHF, CNTIT, DIM1, DIM2,
     $   EXCIT, EROW, F, G, I, IERR, J, K, KBOT, L, LEN, NBUMPS, 
     $   NCHASE, NH, NSHF, OFF, OLDHI, OLDLO, OVRSHF, P, PDW, 
     $   PHESS, PLV, PRV, PTRIU, PU, PU11, PU12, PU21, PU22,
     $   PV, PV11, PV12, PV21, PV22, SROW, UPDHI, UPDLO, USIZE,
     $   VPOS, WINSIZ
      INTEGER           TOTNDINF, NDINF1, NDINF2
      DOUBLE PRECISION  ESHIFT, SCAL, SUM, T1, T2, TAU, TEMP, HCMP

      DOUBLE PRECISION    SMLNUM, ULP, TFNORM, T1NORM, HFNORM, HTOL,
     $   TTOL
      COMMON /PREC/       SMLNUM, ULP, TFNORM, T1NORM, HFNORM, HTOL, 
     $   TTOL
      INTEGER IAM, NPROCS
C     .. Plenty of Timers ..
C     TIME_ALL  = Overall computational time
C     TIME_DEFL = Time spent for aggressive deflation/shift compt.
C     TIME_INF  = Time spent for deflating infinite eigenvalues
C     TIME_KK   = Time spent for block QZ iterations
C     TIME_BL3  = Time for BLAS 3 operations
C     TIME_DACK = Time spent for Dackland/Kagstrom QZ
C     TIME MOLS = Time spent for Moler/Stewart QZ
      DOUBLE PRECISION  TIME_ALL, TIME_BL3, TIME_DACK, TIME_DEFL,
     $   TIME_INF, TIME_KK, TIME_MOLS
      COMMON /TIMERS/   TIME_ALL, TIME_DACK, TIME_DEFL,
     $   TIME_INF, TIME_KK, TIME_MOLS
      INTEGER           MSGLVL
      COMMON /MESSAGES/ MSGLVL
C     .. Local Arrays ..
      DOUBLE PRECISION  V(BULGMX+1)
C     .. External Functions ..
      LOGICAL           LSAME
      DOUBLE PRECISION  DDOT, DLAMCH, DLANHS, DSECND
      EXTERNAL          DDOT, DLAMCH, DLANHS, DSECND, LSAME
C     
C     .. Executable Statements ..
C     
C     MSGLVL = 0 - only error messages are displayed
C     MSGLVL = 1 - also timings and major events are displayed
C     MSGLVL = 2 - lots of details are displayed
C     
      MSGLVL = 2     
      MSGLVL = 0
      INFO   = 0
C     
C     Decode and check input parameters.
C     
      LSCHUR = LSAME( JOB, 'S' )
      LCOMPQ = LSAME( COMPQ, 'V' ).OR.LSAME( COMPQ, 'I' )
      LINITQ = LSAME( COMPQ, 'I' )
      LCOMPZ = LSAME( COMPZ, 'V' ).OR.LSAME( COMPZ, 'I' )
      LINITZ = LSAME( COMPZ, 'I' )
      IF ( LCOMPQ )
     $   UPDQ = 'V'
      IF ( LCOMPZ )
     $   UPDZ = 'V'
      
      IF ( INFO.NE.0 ) THEN
         CALL XERBLA( 'KKQZ', -INFO )
         RETURN
      END IF
C     
C     Quick return if possible
C     
      IF ( N.LE.0 ) THEN
         DWORK(1) = ONE
         RETURN
      END IF
      TOTNDINF = 0

C     
C     Initialize Q and Z if requested.
C     
      IF ( LINITQ )
     $   CALL DLASET( 'Full', N, N, ZERO, ONE, Q, LDQ )
      IF ( LINITZ )
     $   CALL DLASET( 'Full', N, N, ZERO, ONE, Z, LDZ )
C     
C     Machine constants for neighborwise and normwise deflation
C     criteria.
C     
      NH   = IHI + 1 - ILO
      HCMP = HTOL

C     deflation is applied in the current iteration.
C     
      TIME_ALL  = DSECND()
      TIME_BL3  = ZERO
      TIME_KK   = ZERO
      TIME_DACK = ZERO
      TIME_DEFL = ZERO
      TIME_INF  = ZERO
      TIME_MOLS = ZERO
C     
      OVRSHF = 0
      AGGSHF = 0
      EXCIT  = 0
      CNTIT  = 0
      ACTLO  = ILO
      ACTHI  = IHI
      IF ( LSCHUR ) THEN
         UPDLO = 1
         UPDHI = N
      ELSE
         UPDLO = ILO
         UPDHI = IHI
      END IF
      DOAGGR = .TRUE.
      P = N+1
C     
C     ------------------------------------------------------------------
C     Begin of main loop.
C     ------------------------------------------------------------------
C     while acthi >= actlo & cntswp <= maxswp*nh do
C     
 42   IF ( ( ACTHI.GE.ACTLO ) .AND.
     $   ( OVRSHF.LE.PARA_ITMAX*NH ) ) THEN
C        
C        Search for deflations in Hessenberg matrix.
C        
         DO 70  J = ACTHI, ACTLO+1, -1
            TEMP = ABS( H(J,J) ) + ABS( H(J-1,J-1) )
            IF ( TEMP.EQ.ZERO ) THEN
               TEMP = DLANHS( '1', ACTHI-ACTLO+1, H(ACTLO,ACTLO),
     $            LDH, DWORK )
            END IF
            HCMP = MAX( SMLNUM, ULP*TEMP )
            IF ( ABS( H(J,J-1)).LE.HCMP ) THEN
               H(J,J-1) = ZERO
               ACTLO = J
               IF ( .NOT.LSCHUR ) THEN
                  UPDLO = ACTLO
               END IF
               GO TO 80
            END IF
 70      CONTINUE
C        
 80      CONTINUE

C        
C        Search for deflation of infinite eigenvalues.
C        
         OLDLO = ACTLO
         OLDHI = ACTHI
         T1 = DSECND()
         CALL QZINF( LSCHUR, LCOMPQ, LCOMPZ, N, ACTLO, ACTHI, TFNORM, H,
     $      LDH, T, LDT, Q, LDQ, Z, LDZ, DWORK, LDWORK, IWORK,
     $      LIWORK, IERR )

c        ACTLO = MIN(ACTLO, OLDHI)
c        ACTHI = MAX(ACTHI, OLDLO)

         NDINF1 = ACTLO - OLDLO
         NDINF2 = OLDHI - ACTHI
         IF (NDINF1 + NDINF2 .GT. 0) THEN
            TOTNDINF = TOTNDINF + NDINF1 + NDINF2
            WRITE(*,*)'% KKQZ2 Deflated ', NDINF1+NDINF2, 
     $         ' oo eig.values'
         END IF
         TIME_INF = TIME_INF + ( DSECND() - T1 )
         IF ( MSGLVL.GE.1 .AND. ACTLO-OLDLO.GE.1 ) THEN
            WRITE(*,'(I4,A)') ACTLO - OLDLO,
     $         '% infinite eigenvalues deflated at top.'
         END IF
         IF ( MSGLVL.GE.1 .AND. OLDHI-ACTHI.GE.1 ) THEN
            WRITE(*,'(I4,A)') OLDHI-ACTHI,
     $         '% infinite eigenvalues deflated at bottom.'
         END IF
C        
         IF ( ACTLO.NE.OLDLO .OR. ACTHI.NE.OLDHI ) THEN
            DO I = OLDLO+1, MIN(ACTLO, OLDHI)
               BETA(I) = T(I,I)
            END DO
            DO I = MAX(ACTHI,OLDLO)+1, OLDHI
               BETA(I) = T(I,I)
            END DO

C           
C           Try again to search for standard deflations after an
C           infinite eigenvalue has been found.
C           
            IF (ACTHI.LE.ACTLO) THEN
               ACTHI = OLDLO - 1
               ACTLO = ILO
            END IF

            GO TO 42
         END IF


C        
C        Try aggressive early deflation.
C        
         IF ( DOAGGR .AND. ACTHI-ACTLO+1.GE.PARA_AGGMIN ) THEN
            T1 = DSECND()
            WINSIZ = MIN( PARA_AGGWIN, ACTHI-ACTLO+1 )
            P = ACTHI-WINSIZ+1
            IF ( P.GT.ACTLO )  THEN
               SCAL = H(P,P-1)
            ELSE
               SCAL = ZERO
            END IF
C           
C           Set up workspace.
C           
C           PU    - points to matrix containing left trafos
            PU    = 1
C           PV    - points to matrix containing right trafos
            PV    = PU +    ( WINSIZ + PARA_LDADD )*WINSIZ
C           PHESS - points to Hessenberg part of reduced window pair
            PHESS = PV +    ( WINSIZ + PARA_LDADD )*WINSIZ
C           PTRIU - points to triangular part of reduced window pair
            PTRIU = PHESS + ( WINSIZ + PARA_LDADD )*WINSIZ
C           PDW   - points to remaining workspace
            PDW   = PTRIU + ( WINSIZ + PARA_LDADD )*WINSIZ
C           
            CALL DLASET( 'All', WINSIZ, WINSIZ, ZERO, ONE, DWORK(PU),
     $         WINSIZ+PARA_LDADD )
            CALL DLASET( 'All', WINSIZ, WINSIZ, ZERO, ONE, DWORK(PV),
     $         WINSIZ+PARA_LDADD )
            CALL DLACPY( 'All', WINSIZ, WINSIZ, H(P,P), LDH,
     $         DWORK(PHESS), WINSIZ+PARA_LDADD )
            CALL DLACPY( 'All', WINSIZ, WINSIZ, T(P,P), LDT,
     $         DWORK(PTRIU), WINSIZ+PARA_LDADD )
C           
C           Do aggressive early deflation.
C           
            CALL QZEARLY( WINSIZ, SCAL, HFNORM, DWORK(PHESS),
     $         WINSIZ+PARA_LDADD, DWORK(PTRIU),
     $         WINSIZ+PARA_LDADD, KBOT, ALPHAR(P), ALPHAI(P),
     $         BETA(P), DWORK(PU), WINSIZ+PARA_LDADD,
     $         DWORK(PV), WINSIZ+PARA_LDADD, DWORK(PDW),
     $         LDWORK-PDW+1, IERR )
C           
            IF ( MSGLVL.GE.2 ) THEN
               WRITE(*,'(I4,A,A)') WINSIZ-KBOT, ' eigenvalues deflated',
     $            ' by aggressive early deflation.'
            END IF
C           
            IF ( KBOT.LT.WINSIZ ) THEN
C              
C              Update H and T.
C              
               CALL DLACPY( 'All', WINSIZ, WINSIZ, DWORK(PHESS),
     $            WINSIZ+PARA_LDADD, H(P,P), LDH )
               CALL DLACPY( 'All', WINSIZ, WINSIZ, DWORK(PTRIU),
     $            WINSIZ+PARA_LDADD, T(P,P), LDT )
               PDW = PHESS
C              
               IF ( P.GT.UPDLO ) THEN
                  CALL DGEMM( 'No transpose', 'No transpose', P-UPDLO,
     $               WINSIZ, WINSIZ, ONE, H(UPDLO,P), LDH,
     $               DWORK(PV), WINSIZ+PARA_LDADD, ZERO,
     $               DWORK(PDW), P-UPDLO+PARA_LDADD )
                  CALL DLACPY( 'All', P-UPDLO, WINSIZ, DWORK(PDW),
     $               P-UPDLO+PARA_LDADD, H(UPDLO,P), LDH )
                  H(P,P-1) = SCAL*DWORK(PU)
               END IF
C              
               IF ( UPDHI.GT.ACTHI ) THEN
                  CALL DGEMM( 'Transpose', 'No Transpose', WINSIZ,
     $               UPDHI-ACTHI, WINSIZ, ONE, DWORK(PU),
     $               WINSIZ+PARA_LDADD, H(P,ACTHI+1), LDH,
     $               ZERO, DWORK(PDW), WINSIZ+PARA_LDADD )
                  CALL DLACPY( 'All', WINSIZ, UPDHI-ACTHI, DWORK(PDW),
     $               WINSIZ+PARA_LDADD, H(P,ACTHI+1), LDH )
               END IF
C              
               IF ( P.GT.UPDLO ) THEN
                  CALL DGEMM( 'No transpose', 'No transpose', P-UPDLO,
     $               WINSIZ, WINSIZ, ONE, T(UPDLO,P), LDT,
     $               DWORK(PV), WINSIZ+PARA_LDADD, ZERO,
     $               DWORK(PDW), P-UPDLO+PARA_LDADD )
                  CALL DLACPY( 'All', P-UPDLO, WINSIZ, DWORK(PDW),
     $               P-UPDLO+PARA_LDADD, T(UPDLO,P), LDT )
               END IF
C              
               IF ( UPDHI.GT.ACTHI ) THEN
                  CALL DGEMM( 'Transpose', 'No Transpose', WINSIZ,
     $               UPDHI-ACTHI, WINSIZ, ONE, DWORK(PU),
     $               WINSIZ+PARA_LDADD, T(P,ACTHI+1), LDT,
     $               ZERO, DWORK(PDW), WINSIZ+PARA_LDADD )
                  CALL DLACPY( 'All', WINSIZ, UPDHI-ACTHI, DWORK(PDW),
     $               WINSIZ+PARA_LDADD, T(P,ACTHI+1), LDT )
               END IF
C              
C              Update Q and Z.
C              
               IF ( LCOMPQ ) THEN
                  CALL DGEMM( 'No transpose', 'No transpose', N, WINSIZ,
     $               WINSIZ, ONE, Q(1,P), LDQ, DWORK(PU),
     $               WINSIZ+PARA_LDADD, ZERO, DWORK(PDW),
     $               N+PARA_LDADD )
                  CALL DLACPY( 'All', N, WINSIZ, DWORK(PDW),
     $               N+PARA_LDADD, Q(1,P), LDQ )
               END IF
               IF ( LCOMPZ ) THEN
                  CALL DGEMM( 'No transpose', 'No transpose', N, WINSIZ,
     $               WINSIZ, ONE, Z(1,P), LDZ, DWORK(PV),
     $               WINSIZ+PARA_LDADD, ZERO, DWORK(PDW),
     $               N+PARA_LDADD )
                  CALL DLACPY( 'All', N, WINSIZ, DWORK(PDW),
     $               N+PARA_LDADD, Z(1,P), LDZ )
               END IF
            END IF
C           
            TIME_DEFL = TIME_DEFL + ( DSECND() - T1 )
            DOAGGR = .FALSE.
C           
C           If sufficiently many eigenvalues have been deflated, let us
C           try another aggressive early deflation instead of applying a
C           QZ iteration. In any case, we jump back to the beginning to
C           also be able to take care of infinite eigenvalues. Note that
C           P points to the eigenvalues that will be used as shifts
C           later on.
C           
            IF ( DBLE(WINSIZ-KBOT) / DBLE(WINSIZ) .GE. PARA_AGGREP )THEN
               DOAGGR = .TRUE.
            END IF
            GO TO 42
         END IF
C        
C        Having tried all the deflations we could think of, it is time
C        to do some QZ iterations!
C        
         EXCIT = EXCIT + 1
         IF ( EXCIT.GT.PARA_ITEXC ) THEN
            LEXC = .TRUE.
            EXCIT = 0
            IF ( MSGLVL.GT.0 )
     $         WRITE(*,*) '     !!!Exceptional shift applied!!!'
         ELSE
            LEXC = .FALSE.
         END IF
C        
         IF ( ACTHI-ACTLO+1.GE.PARA_ITKK .AND.
     $      ACTHI-ACTLO+1.GT.PARA_ITNS .AND.
     $      .NOT.LEXC ) THEN
C           
C           Do multishift QZ by chasing tightly coupled chains of tiny
C           bulges.
C           
            IF ( MSGLVL.GE.2 )
     $         WRITE(*,'(A,I4,A,I4)') 'Apply KK QZ to:      ', ACTLO,
     $         ' ...', ACTHI
            T1 = DSECND()
C           
C           There are some shifts from aggressive early deflation
C           residing between P and ACTHI, lets use them if possible.
C           Otherwise, we have to recompute the shifts.
C           
            IF ( PARA_ITASH.GT.0 .AND. P.GT.ACTLO .AND.
     $         (ACTHI-P+1).GE.PARA_ITNS ) THEN
C              
C              Reverse if desired.
C              
               IF ( PARA_ITASH.EQ.2 ) THEN
                  DO 90  I = 0, ( ACTHI - P - 1 ) / 2
                     TEMP = ALPHAR(P+I)
                     ALPHAR(P+I) = ALPHAR(ACTHI-I)
                     ALPHAR(ACTHI-I) = TEMP
                     TEMP = ALPHAI(P+I)
                     ALPHAI(P+I) = ALPHAI(ACTHI-I)
                     ALPHAI(ACTHI-I) = TEMP
                     TEMP = BETA(P+I)
                     BETA(P+I) = BETA(ACTHI-I)
                     BETA(ACTHI-I) = TEMP
 90               CONTINUE
               END IF
C              
C              Shuffle real shifts into pairs.
C              
               DO 100 I = P, ACTHI-2, 2
                  IF ( ALPHAI(I).EQ.ZERO .AND. ALPHAI(I+1).NE.ZERO )
     $               THEN
                     TEMP = ALPHAR(I)
                     ALPHAR(I) = ALPHAR(I+1)
                     ALPHAR(I+1) = ALPHAR(I+2)
                     ALPHAR(I+2) = TEMP
                     TEMP = ALPHAI(I)
                     ALPHAI(I) = ALPHAI(I+1)
                     ALPHAI(I+1) = ALPHAI(I+2)
                     ALPHAI(I+2) = TEMP
                     TEMP = BETA(I)
                     BETA(I) = BETA(I+1)
                     BETA(I+1) = BETA(I+2)
                     BETA(I+2) = TEMP
                  END IF
 100           CONTINUE
            ELSE
C              
C              Apply QZ to compute some shifts.
C              
               P = ACTHI - PARA_ITNS + 1
               PHESS = 1
               PTRIU = PHESS + ( PARA_ITNS + PARA_LDADD )*PARA_ITNS
               PDW   = PTRIU + ( PARA_ITNS + PARA_LDADD )*PARA_ITNS
               CALL DLACPY( 'All', PARA_ITNS, PARA_ITNS, H(P,P), LDH,
     $            DWORK(PHESS), PARA_ITNS+PARA_LDADD )
               CALL DLACPY( 'All', PARA_ITNS, PARA_ITNS, T(P,P), LDT,
     $            DWORK(PTRIU), PARA_ITNS+PARA_LDADD )
C              
               IF ( PARA_ITNS.GE.PARA_ITDACK ) THEN
                  CALL QZLAP( 'Eigenvalues', 'N', 'N', PARA_ITNS, 1,
     $               PARA_ITNS, DWORK(PHESS),
     $               PARA_ITNS+PARA_LDADD, DWORK(PTRIU),
     $               PARA_ITNS+PARA_LDADD, ALPHAR(P),
     $               ALPHAI(P), BETA(P), DWORK(PDW), 1,
     $               DWORK(PDW), 1, DWORK(PDW), LDWORK, IERR )
               ELSE
                  CALL QZLAP( 'Eigenvalues', 'N', 'N', PARA_ITNS, 1,
     $               PARA_ITNS, DWORK(PHESS),
     $               PARA_ITNS+PARA_LDADD, DWORK(PTRIU),
     $               PARA_ITNS+PARA_LDADD, ALPHAR(P),
     $               ALPHAI(P), BETA(P), DWORK(PDW), 1,
     $               DWORK(PDW), 1, DWORK(PDW), LDWORK, IERR )
               END IF
C              
               IF ( IERR.GT.1 ) THEN
C                 
C                 Set non-converged eigenvalues.
C                 
                  IF ( MSGLVL.GE.1 ) THEN
                     WRITE(*,*) 'WARNING: QZ DID NOT CONVERGE AT SOME ',
     $                  'SUBPROBLEM'
                  END IF                  
                  DO 110 I = 0, INFO-1
                     ALPHAR(P+I) = DWORK(PHESS + I*( PARA_ITNS
     $                  + PARA_LDADD ) + I )
                     ALPHAI(P+I) = ZERO
                     BETA(P+I)   = DWORK(PTRIU + I*( PARA_ITNS
     $                  + PARA_LDADD ) + I )
 110              CONTINUE
               END IF
            END IF
C           
C           Set up the sweep.
C           
            SROW   = ACTLO
            NBUMPS = PARA_ITNS / PARA_ITBULG
            USIZE  = ( PARA_ITBULG + 1 )*NBUMPS

            EROW   = USIZE + SROW - 1
            EROW   = MIN( EROW , ACTHI )
            USIZE = MIN( EROW - SROW + 1, USIZE )
            NBUMPS = MIN( NBUMPS, USIZE / ( PARA_ITBULG + 1 ) )

            PU     = 1
            PV     = PU + ( USIZE + PARA_LDADD )*USIZE
            PDW    = PV + ( USIZE + PARA_LDADD )*USIZE
            CALL DLASET( 'All', USIZE, USIZE, ZERO, ONE, DWORK(PU),
     $         USIZE + PARA_LDADD )
            CALL DLASET( 'All', USIZE, USIZE, ZERO, ONE, DWORK(PV),
     $         USIZE + PARA_LDADD )
            DO 142 F = 1, NBUMPS
C              
C              Compute first column of shift polynomial.
C              
               CALL QZFCOL( PARA_ITBULG+1, PARA_ITBULG,
     $            ALPHAR(P+(F-1)*PARA_ITBULG),
     $            ALPHAI(P+(F-1)*PARA_ITBULG),
     $            BETA(P+(F-1)*PARA_ITBULG), H(ACTLO,ACTLO),
     $            LDH, T(ACTLO,ACTLO), LDT, V, IERR )
C              
               CALL DLARFG( PARA_ITBULG+1, V(1), V(2), 1, TAU )
               V(1) = ONE
               CALL DLARFX( 'Left', PARA_ITBULG+1, USIZE, V, TAU,
     $            H(SROW,SROW), LDH, DWORK(PDW) )
               CALL DLARFX( 'Left', PARA_ITBULG+1, USIZE, V, TAU,
     $            T(SROW,SROW), LDT, DWORK(PDW) )
               CALL DLARFX( 'Right', USIZE, PARA_ITBULG+1, V, TAU,
     $            DWORK(PU), USIZE+PARA_LDADD, DWORK(PDW) )
C              
               CALL INVHSE( PARA_ITBULG+1, T(SROW,SROW), LDT, TAU, V,
     $            IERR )
               IF ( IERR.NE.0 ) THEN
                  IF ( MSGLVL.GE.1 ) THEN
                     WRITE(*,*) 'Opposite Householder failed.'
                  END IF
                  INFO = 1
                  RETURN
               END IF
               CALL DLARFX( 'Right', PARA_ITBULG+1, PARA_ITBULG+1, V,
     $            TAU, T(SROW,SROW), LDT, DWORK(PDW) )
               DO 116 J = 1, PARA_ITBULG
                  T(SROW+J,SROW) = ZERO
 116           CONTINUE
               CALL DLARFX( 'Right', PARA_ITBULG+2, PARA_ITBULG+1, V,
     $            TAU, H(SROW,SROW), LDH, DWORK(PDW) )
               CALL DLARFX( 'Right', USIZE, PARA_ITBULG+1, V, TAU,
     $            DWORK(PV), USIZE+PARA_LDADD, DWORK(PDW) )
C              
C              Move bump down.
C              
               DO 117 G = SROW, EROW - ( PARA_ITBULG + 1 )*F
                  CALL DLARFG( PARA_ITBULG+1, H(G+1,G), H(G+2,G), 1,
     $               TAU )
                  TEMP = H(G+1,G)
                  H(G+1,G) = ONE
                  CALL DLARFX( 'Left', PARA_ITBULG+1, EROW-G, H(G+1,G),
     $               TAU, H(G+1,G+1), LDH, DWORK(PDW) )
                  CALL DLARFX( 'Left', PARA_ITBULG+1, EROW-G, H(G+1,G),
     $               TAU, T(G+1,G+1), LDT, DWORK(PDW) )
                  CALL DLARFX( 'Right', USIZE, PARA_ITBULG+1, H(G+1,G),
     $               TAU,
     $               DWORK(PU+(G-SROW+1)*(USIZE+PARA_LDADD)),
     $               USIZE+PARA_LDADD, DWORK(PDW) )
                  H(G+1,G) = TEMP
                  DO 118 J = 2, PARA_ITBULG+1
                     H(G+J,G) = ZERO
 118              CONTINUE
                  CALL INVHSE( PARA_ITBULG+1, T(G+1,G+1), LDT, TAU, V,
     $               IERR )
                  IF ( IERR.NE.0 ) THEN
                     IF ( MSGLVL.GE.1 ) THEN
                        WRITE(*,*) 'Opposite Householder failed.'
                     END IF
                     INFO = 1
                     RETURN
                  END IF
                  
                  CALL DLARFX( 'Right', G+PARA_ITBULG-SROW+2,
     $               PARA_ITBULG+1, V, TAU, T(SROW,G+1), LDT,
     $               DWORK(PDW) )
                  DO 716 J = 2, PARA_ITBULG+1
                     T(G+J,G+1) = ZERO
 716              CONTINUE
                  CALL DLARFX( 'Right', G+PARA_ITBULG-SROW+3,
     $               PARA_ITBULG+1, V, TAU, H(SROW,G+1), LDH,
     $               DWORK(PDW) )
                  CALL DLARFX( 'Right', USIZE, PARA_ITBULG+1, V, TAU,
     $               DWORK(PV+(G-SROW+1)*(USIZE+PARA_LDADD)),
     $               USIZE+PARA_LDADD, DWORK(PDW) )
 117           CONTINUE
 142        CONTINUE
C           
            T2 = DSECND()
C           
C           Update H.
C           
            IF ( UPDHI.GT.EROW ) THEN
               CALL DGEMM( 'Transpose', 'No transpose', USIZE,
     $            UPDHI-EROW, USIZE, ONE, DWORK(PU),
     $            USIZE+PARA_LDADD, H(SROW,EROW+1), LDH,
     $            ZERO, DWORK(PDW), USIZE+PARA_LDADD )
               CALL DLACPY( 'All', USIZE, UPDHI-EROW, DWORK(PDW),
     $            USIZE+PARA_LDADD, H(SROW,EROW+1), LDH )
            END IF
C           
            
            IF ( SROW.GT.UPDLO ) THEN
               IF (UPDLO.GT.LDH .OR. SROW.GT.LDH) THEN
                  CALL BLACS_PINFO(IAM, NPROCS)
                  WRITE(*,*)'Err KKQZ2:', UPDLO, LDH, SROW, F, N, ILO, 
     $               IHI, ACTLO, ACTHI, P, IAM
                  STOP

               END IF
               CALL DGEMM( 'No transpose', 'No transpose',
     $            SROW-UPDLO, USIZE, USIZE, ONE,
     $            H(UPDLO,SROW), LDH, DWORK(PV),
     $            USIZE+PARA_LDADD, ZERO, DWORK(PDW),
     $            SROW-UPDLO+PARA_LDADD )
               CALL DLACPY( 'All', SROW-UPDLO, USIZE, DWORK(PDW),
     $            SROW-UPDLO+PARA_LDADD, H(UPDLO,SROW), LDH )
            END IF
C           
C           Update T.
C           
            IF ( UPDHI.GT.EROW ) THEN
               CALL DGEMM( 'Transpose', 'No transpose', USIZE,
     $            UPDHI-EROW, USIZE, ONE, DWORK(PU),
     $            USIZE+PARA_LDADD, T(SROW,EROW+1), LDT,
     $            ZERO, DWORK(PDW), USIZE+PARA_LDADD )
               CALL DLACPY( 'All', USIZE, UPDHI-EROW, DWORK(PDW),
     $            USIZE+PARA_LDADD, T(SROW,EROW+1), LDT )
            END IF
C           
            IF ( SROW.GT.UPDLO ) THEN
               CALL DGEMM( 'No transpose', 'No transpose',
     $            SROW-UPDLO, USIZE, USIZE, ONE, T(UPDLO,SROW),
     $            LDT, DWORK(PV), USIZE+PARA_LDADD, ZERO,
     $            DWORK(PDW), SROW-UPDLO+PARA_LDADD )
               CALL DLACPY( 'All', SROW-UPDLO, USIZE, DWORK(PDW),
     $            SROW-UPDLO+PARA_LDADD, T(UPDLO,SROW), LDT )
            END IF
C           
C           Update Q.
C           
            IF ( LCOMPQ ) THEN
               CALL DGEMM( 'No transpose', 'No transpose', N, USIZE,
     $            USIZE, ONE, Q(1,SROW), LDQ, DWORK(PU),
     $            USIZE+PARA_LDADD, ZERO, DWORK(PDW),
     $            N+PARA_LDADD )
               CALL DLACPY( 'All', N, USIZE, DWORK(PDW), N+PARA_LDADD,
     $            Q(1,SROW), LDQ )
            END IF
C           
C           Update Z.
C           
            IF ( LCOMPZ ) THEN
               CALL DGEMM( 'No transpose', 'No transpose', N, USIZE,
     $            USIZE, ONE, Z(1,SROW), LDZ, DWORK(PV),
     $            USIZE+PARA_LDADD, ZERO, DWORK(PDW),
     $            N+PARA_LDADD )
               CALL DLACPY( 'All', N, USIZE, DWORK(PDW), N+PARA_LDADD,
     $            Z(1,SROW), LDZ )
            END IF
            TIME_BL3 = TIME_BL3 + ( DSECND() - T2 )
C           
C           Chase Bulges.
C           
            USIZE  = (PARA_ITBULG+1)*NBUMPS + PARA_ITEXT
            NCHASE = MAX( 0, ( ACTHI - ACTLO -
     $         (PARA_ITBULG+1)*NBUMPS ) / PARA_ITEXT )
            PU  = 1
            PV  = PU + ( USIZE + PARA_LDADD )*USIZE
            PLV = PV + ( USIZE + PARA_LDADD )*USIZE
C           ??? check whether PRV can be removed !!!
            PRV = PLV + (PARA_ITBULG+1)*PARA_ITEXT
            PDW = PRV + (PARA_ITBULG+1)*PARA_ITEXT
C           
            DO 222  K = 1, NCHASE
               EROW = SROW + USIZE - 1
               CALL DLASET( 'All', USIZE, USIZE, ZERO, ONE,
     $            DWORK(PU), USIZE+PARA_LDADD )
               CALL DLASET( 'All', USIZE, USIZE, ZERO, ONE,
     $            DWORK(PV), USIZE+PARA_LDADD )
C              
               DO 223 L = 0, EROW-SROW-(PARA_ITBULG+1)*NBUMPS
                  DO 224 F = 1, NBUMPS
                     G = SROW + (PARA_ITBULG+1)*(NBUMPS-F) + L
                     VPOS = PLV + (F-1)*(PARA_ITBULG+1)
                     CALL DLARFG( PARA_ITBULG+1, H(G+1,G), H(G+2,G), 1,
     $                  TAU )
                     DWORK(VPOS+1) = TAU
                     DO 218 J = 2, PARA_ITBULG+1
                        DWORK(VPOS+J) = H(G+J,G)
                        H(G+J,G) = ZERO
 218                 CONTINUE
 224              CONTINUE
C                 
C                 Multiply from left to H within extended window.
C                 
                  DO 225 J = SROW+L+1, EROW
                     DO 226 F = MAX( 1, 1+(SROW+(PARA_ITBULG+1)*NBUMPS+
     $                  L-J) / (PARA_ITBULG+1) ), NBUMPS
                        G = SROW + (PARA_ITBULG+1)*(NBUMPS-F) + L
                        VPOS = PLV + (F-1)*(PARA_ITBULG+1)
                        SUM = H(G+1,J) + DDOT( PARA_ITBULG, 
     $                     DWORK(VPOS+2), 1, H(G+2,J), 1 )
                        TEMP = -SUM*DWORK(VPOS+1)
                        H(G+1,J) = H(G+1,J) + TEMP
                        CALL DAXPY( PARA_ITBULG, TEMP, DWORK(VPOS+2), 1,
     $                     H(G+2,J), 1 )
 226                 CONTINUE  
 225              CONTINUE
C                 
C                 Multiply from left to T within extended window.
C                 
                  DO 725 J = SROW+L+1, EROW
                     DO 726 F = MAX( 1, 1+(SROW+(PARA_ITBULG+1)*NBUMPS+
     $                  L-J) / (PARA_ITBULG+1) ), NBUMPS
                        G = SROW + (PARA_ITBULG+1)*(NBUMPS-F) + L
                        VPOS = PLV + (F-1)*(PARA_ITBULG+1)
                        SUM = T(G+1,J) + DDOT( PARA_ITBULG, 
     $                     DWORK(VPOS+2), 1, T(G+2,J), 1 )
                        TEMP = -SUM*DWORK(VPOS+1)
                        T(G+1,J) = T(G+1,J) + TEMP
                        CALL DAXPY( PARA_ITBULG, TEMP, DWORK(VPOS+2), 1,
     $                     T(G+2,J), 1 )
 726                 CONTINUE
 725              CONTINUE
C                 
C                 Multiply from right to U.
C                 
                  DO 727 F = 1, NBUMPS
                     G = SROW + (PARA_ITBULG+1)*(NBUMPS-F) + L
                     VPOS = PLV + (F-1)*(PARA_ITBULG+1)
                     TAU = DWORK(VPOS+1)
                     DWORK(VPOS+1) = ONE
                     CALL DLARFX( 'Right', USIZE, PARA_ITBULG+1,
     $                  DWORK(VPOS+1), TAU, DWORK(PU+
     $                  (G-SROW+1)*(USIZE+PARA_LDADD)),
     $                  USIZE+PARA_LDADD, DWORK(PDW) )
 727              CONTINUE
C                 
C                 Annihilate subdiagonal entries of T.
C                 
                  DO 724 F = 1, NBUMPS
                     G = SROW + (PARA_ITBULG+1)*(NBUMPS-F) + L
                     CALL INVHSE( PARA_ITBULG+1, T(G+1,G+1), LDT, TAU,
     $                  V, IERR )
                     IF ( IERR.NE.0 ) THEN
                        IF ( MSGLVL.GE.1 ) THEN
                           WRITE(*,*) 'Opposite Householder failed.'
                        END IF
                        INFO = 1
                        RETURN
                     END IF
                     CALL DLARFX( 'Right', G+PARA_ITBULG-SROW+2,
     $                  PARA_ITBULG+1, V, TAU, T(SROW,G+1),
     $                  LDT, DWORK(PDW) )
                     DO 816 J = 2, PARA_ITBULG+1
                        T(G+J,G+1) = ZERO
 816                 CONTINUE
                     CALL DLARFX( 'Right', G-SROW+PARA_ITBULG+3,
     $                  PARA_ITBULG+1, V, TAU, H(SROW,G+1),
     $                  LDH, DWORK(PDW) )
                     CALL DLARFX( 'Right', USIZE, PARA_ITBULG+1, V, TAU,
     $                  DWORK(PV+(G-SROW+1)*(USIZE+
     $                  PARA_LDADD)), USIZE+PARA_LDADD,
     $                  DWORK(PDW) )
 724              CONTINUE
 223           CONTINUE
C              
               T2 = DSECND()
C              
C              Multiply U and V to rest of H and T, as well as Q and Z.
C              
               PDW = PV + (USIZE+PARA_LDADD)*USIZE
               IF ( PARA_ITTRI ) THEN
C                 
C                 [  1    0    0  ]
C                 Exploit U = [  0   U11  U12 ], where U12 is lower
C                 [  0   U21  U22 ]
C                 
C                 triangular and U21 is upper triangular. U11 is
C                 dim1-by-(dim2+off) and U22 is (dim2+off)-by-dim1.
C                 The same applies to V.
C                 
                  OFF  = NBUMPS-PARA_ITBULG
                  IF ( OFF.LT.0 ) THEN
                     IF ( MSGLVL.GE.1 )
     $                  WRITE(*,*) 'ERROR: NEGATIVE OFF !'
                     INFO = 1
                     RETURN
                  END IF
                  DIM1 = NBUMPS*( PARA_ITBULG + 1 ) - PARA_ITBULG
                  DIM2 = USIZE - DIM1 - OFF - 1
                  PU11 = PU + USIZE + PARA_LDADD + 1
                  PU21 = PU11 + DIM1 + OFF*( USIZE + PARA_LDADD )
                  PU12 = PU11 + ( DIM2 + OFF )*( USIZE + PARA_LDADD )
                  PU22 = PU12 + DIM1
                  PV11 = PV + USIZE + PARA_LDADD + 1
                  PV21 = PV11 + DIM1 + OFF*( USIZE + PARA_LDADD )
                  PV12 = PV11 + ( DIM2 + OFF )*( USIZE + PARA_LDADD )
                  PV22 = PV12 + DIM1
C                 
C                 Update H..
C                 Copy bottom of H (den. by H2) to top of workspace.
C                 
                  CALL DLASET( 'All', OFF, UPDHI-EROW, ZERO, ZERO, 
     $               DWORK(PDW), USIZE+PARA_LDADD-1 )
                  IF ( DIM2.GT.0 ) THEN
                     CALL DLACPY( 'All', DIM2, UPDHI-EROW,
     $                  H(SROW+DIM1+1,EROW+1), LDH,
     $                  DWORK(PDW+OFF), USIZE+PARA_LDADD-1 )
C                    
C                    Compute U21''*H2.
C                    
                     CALL DTRMM( 'Left', 'Upper', 'Transpose',
     $                  'No transpose', DIM2, UPDHI-EROW, ONE,
     $                  DWORK(PU21), USIZE+PARA_LDADD,
     $                  DWORK(PDW+OFF), USIZE+PARA_LDADD-1 )
                  END IF
C                 
C                 Add U11''*H1.
C                 
                  CALL DGEMM( 'Transpose', 'No transpose', DIM2+OFF,
     $               UPDHI-EROW, DIM1, ONE, DWORK(PU11),
     $               USIZE+PARA_LDADD, H(SROW+1,EROW+1), LDH,
     $               ONE, DWORK(PDW), USIZE+PARA_LDADD-1 )
C                 
C                 Copy top of H (den. by H1) to bottom of workspace.
C                 
                  CALL DLACPY( 'All', DIM1, UPDHI-EROW,
     $               H(SROW+1,EROW+1), LDH,
     $               DWORK(PDW+OFF+DIM2), USIZE+PARA_LDADD-1 )
C                 
C                 Compute U12''*H1.
C                 
                  CALL DTRMM( 'Left', 'Lower', 'Transpose',
     $               'No transpose', DIM1, UPDHI-EROW, ONE,
     $               DWORK(PU12), USIZE+PARA_LDADD,
     $               DWORK(PDW+OFF+DIM2), USIZE+PARA_LDADD-1 )
C                 
C                 Add U22''*H2.
C                 
                  CALL DGEMM( 'Transpose', 'No Transpose', DIM1,
     $               UPDHI-EROW, DIM2+OFF, ONE, DWORK(PU22),
     $               USIZE+PARA_LDADD, H(SROW+DIM1+1,EROW+1),
     $               LDH, ONE, DWORK(PDW+OFF+DIM2),
     $               USIZE+PARA_LDADD-1 )
                  CALL DLACPY( 'All', USIZE-1, UPDHI-EROW, DWORK(PDW),
     $               USIZE+PARA_LDADD-1, H(SROW+1,EROW+1),
     $               LDH )
C                 
                  IF ( SROW.GT.UPDLO ) THEN
C                    
C                    Copy right of H (den. by H2) to left of workspace.
C                    
                     CALL DLASET( 'All', SROW-UPDLO, OFF, ZERO, ZERO,
     $                  DWORK(PDW), SROW+PARA_LDADD-UPDLO )
                     IF ( DIM2.GT.0 ) THEN
                        CALL DLACPY( 'All', SROW-UPDLO, DIM2,
     $                     H(UPDLO,SROW+DIM1+1), LDH,
     $                     DWORK(PDW+OFF*(SROW-UPDLO
     $                     +PARA_LDADD)), SROW-UPDLO
     $                     +PARA_LDADD )
C                       
C                       Compute H2*V21.
C                       
                        CALL DTRMM( 'Right', 'Upper', 'No transpose',
     $                     'No transpose', SROW-UPDLO, DIM2,
     $                     ONE, DWORK(PV21), USIZE+PARA_LDADD,
     $                     DWORK(PDW+OFF*(SROW-UPDLO
     $                     +PARA_LDADD)),
     $                     SROW-UPDLO+PARA_LDADD )
                     END IF
C                    
C                    Add H1*V11.
C                    
                     CALL DGEMM( 'No transpose', 'No transpose',
     $                  SROW-UPDLO, DIM2+OFF, DIM1, ONE,
     $                  H(UPDLO,SROW+1), LDH, DWORK(PV11),
     $                  USIZE+PARA_LDADD, ONE, DWORK(PDW),
     $                  SROW-UPDLO+PARA_LDADD )
C                    
C                    Copy left of H (den. by H1) to right of workspace.
C                    
                     CALL DLACPY( 'All', SROW-UPDLO, DIM1,
     $                  H(UPDLO,SROW+1), LDH,
     $                  DWORK(PDW+(DIM2+OFF)*(SROW+
     $                  PARA_LDADD-UPDLO)),
     $                  SROW+PARA_LDADD-UPDLO )
C                    
C                    Compute H1*V12.
C                    
                     CALL DTRMM( 'Right', 'Lower', 'No transpose',
     $                  'No transpose', SROW-UPDLO, DIM1,
     $                  ONE, DWORK(PV12), USIZE+PARA_LDADD,
     $                  DWORK(PDW+(DIM2+OFF)*(SROW+
     $                  PARA_LDADD-UPDLO)),
     $                  SROW+PARA_LDADD-UPDLO )
C                    
C                    Add H2*V22.
C                    
                     CALL DGEMM( 'No transpose', 'No transpose', 
     $                  SROW-UPDLO, DIM1, DIM2+OFF, ONE,
     $                  H(UPDLO,SROW+DIM1+1), LDH, DWORK(PV22),
     $                  USIZE+PARA_LDADD, ONE,
     $                  DWORK(PDW+(DIM2+OFF)*
     $                  (SROW+PARA_LDADD-UPDLO)),
     $                  SROW+PARA_LDADD-UPDLO )
                     CALL DLACPY( 'All', SROW-UPDLO, USIZE-UPDLO,
     $                  DWORK(PDW), SROW-UPDLO+PARA_LDADD,
     $                  H(UPDLO,SROW+1), LDH )
                  END IF
C                 
C                 Update T..
C                 Copy bottom of T (den. by T2) to top of workspace.
C                 
                  CALL DLASET( 'All', OFF, UPDHI-EROW, ZERO, ZERO, 
     $               DWORK(PDW), USIZE+PARA_LDADD-1 )
                  IF ( DIM2.GT.0 ) THEN
                     CALL DLACPY( 'All', DIM2, UPDHI-EROW,
     $                  T(SROW+DIM1+1,EROW+1), LDT,
     $                  DWORK(PDW+OFF), USIZE+PARA_LDADD-1 )
C                    
C                    Compute U21''*T2.
C                    
                     CALL DTRMM( 'Left', 'Upper', 'Transpose',
     $                  'No transpose', DIM2, UPDHI-EROW, ONE,
     $                  DWORK(PU21), USIZE+PARA_LDADD,
     $                  DWORK(PDW+OFF), USIZE+PARA_LDADD-1 )
                  END IF
C                 
C                 Add U11''*T1.
C                 
                  CALL DGEMM( 'Transpose', 'No transpose', DIM2+OFF,
     $               UPDHI-EROW, DIM1, ONE, DWORK(PU11),
     $               USIZE+PARA_LDADD, T(SROW+1,EROW+1), LDT,
     $               ONE, DWORK(PDW), USIZE+PARA_LDADD-1 )
C                 
C                 Copy top of T (den. by T1) to bottom of workspace.
C                 
                  CALL DLACPY( 'All', DIM1, UPDHI-EROW,
     $               T(SROW+1,EROW+1), LDT,
     $               DWORK(PDW+OFF+DIM2), USIZE+PARA_LDADD-1 )
C                 
C                 Compute U12''*T1.
C                 
                  CALL DTRMM( 'Left', 'Lower', 'Transpose',
     $               'No transpose', DIM1, UPDHI-EROW, ONE,
     $               DWORK(PU12), USIZE+PARA_LDADD,
     $               DWORK(PDW+OFF+DIM2), USIZE+PARA_LDADD-1 )
C                 
C                 Add U22''*T2.
C                 
                  CALL DGEMM( 'Transpose', 'No Transpose', DIM1,
     $               UPDHI-EROW, DIM2+OFF, ONE, DWORK(PU22),
     $               USIZE+PARA_LDADD, T(SROW+DIM1+1,EROW+1),
     $               LDT, ONE, DWORK(PDW+OFF+DIM2),
     $               USIZE+PARA_LDADD-1 )
                  CALL DLACPY( 'All', USIZE-1, UPDHI-EROW, DWORK(PDW),
     $               USIZE+PARA_LDADD-1, T(SROW+1,EROW+1),
     $               LDT )
C                 
                  IF ( SROW.GT.UPDLO ) THEN
C                    
C                    Copy right of T (den. by T2) to left of workspace.
C                    
                     CALL DLASET( 'All', SROW-UPDLO, OFF, ZERO, ZERO,
     $                  DWORK(PDW), SROW+PARA_LDADD-UPDLO )
                     IF ( DIM2.GT.0 ) THEN
                        CALL DLACPY( 'All', SROW-UPDLO, DIM2,
     $                     T(UPDLO,SROW+DIM1+1), LDT,
     $                     DWORK(PDW+OFF*(SROW-UPDLO
     $                     +PARA_LDADD)),
     $                     SROW+PARA_LDADD-UPDLO )
C                       
C                       Compute T2*V21.
C                       
                        CALL DTRMM( 'Right', 'Upper', 'No transpose',
     $                     'No transpose', SROW-UPDLO, DIM2,
     $                     ONE, DWORK(PV21), USIZE+PARA_LDADD,
     $                     DWORK(PDW+OFF*(SROW+PARA_LDADD
     $                     -UPDLO)),
     $                     SROW+PARA_LDADD-UPDLO )
                     END IF
C                    
C                    Add T1*V11.
C                    
                     CALL DGEMM( 'No transpose', 'No transpose',
     $                  SROW-UPDLO, DIM2+OFF, DIM1, ONE,
     $                  T(UPDLO,SROW+1), LDT, DWORK(PV11),
     $                  USIZE+PARA_LDADD, ONE, DWORK(PDW),
     $                  SROW+PARA_LDADD-UPDLO )
C                    
C                    Copy left of T (den. by T1) to right of workspace.
C                    
                     CALL DLACPY( 'All', SROW-UPDLO, DIM1,
     $                  T(UPDLO,SROW+1), LDT,
     $                  DWORK(PDW+(DIM2+OFF)*(SROW+PARA_LDADD
     $                  -UPDLO)), SROW+PARA_LDADD-UPDLO )
C                    
C                    Compute T1*V12.
C                    
                     CALL DTRMM( 'Right', 'Lower', 'No transpose',
     $                  'No transpose', SROW-UPDLO, DIM1, ONE,
     $                  DWORK(PV12), USIZE+PARA_LDADD,
     $                  DWORK(PDW+(DIM2+OFF)*(SROW+PARA_LDADD
     $                  -UPDLO)), SROW+PARA_LDADD-UPDLO )
C                    
C                    Add T2*V22.
C                    
                     CALL DGEMM( 'No transpose', 'No transpose',
     $                  SROW-UPDLO, DIM1, DIM2+OFF, ONE,
     $                  T(UPDLO,SROW+DIM1+1), LDT, DWORK(PV22),
     $                  USIZE+PARA_LDADD, ONE,
     $                  DWORK(PDW+(DIM2+OFF)*(SROW+PARA_LDADD
     $                  -UPDLO)), SROW+PARA_LDADD-UPDLO )
                     CALL DLACPY( 'All', SROW-UPDLO, USIZE-1,
     $                  DWORK(PDW), SROW+PARA_LDADD-UPDLO,
     $                  T(UPDLO,SROW+1), LDT )
                  END IF
C                 
                  IF ( LCOMPQ ) THEN
C                    
C                    Update Q..
C                    Copy right of Q (den. by Q2) to left of workspace.
C                    
                     CALL DLASET( 'All', N, OFF, ZERO, ZERO, DWORK(PDW),
     $                  N+PARA_LDADD )
                     CALL DLACPY( 'All', N, DIM2, Q(1,SROW+DIM1+1), LDQ,
     $                  DWORK(PDW+OFF*(N+PARA_LDADD)),
     $                  N+PARA_LDADD )
                     IF ( DIM2.GT.0 ) THEN
C                       
C                       Compute Q2*U21.
C                       
                        CALL DTRMM( 'Right', 'Upper', 'No transpose',
     $                     'No transpose', N, DIM2, ONE,
     $                     DWORK(PU21), USIZE+PARA_LDADD,
     $                     DWORK(PDW+OFF*(N+PARA_LDADD)),
     $                     N+PARA_LDADD )
                     END IF
C                    
C                    Add Q1*U11.
C                    
                     CALL DGEMM( 'No transpose', 'No transpose', N,
     $                  DIM2+OFF, DIM1, ONE, Q(1,SROW+1), LDQ,
     $                  DWORK(PU11), USIZE+PARA_LDADD, ONE,
     $                  DWORK(PDW), N+PARA_LDADD )
C                    
C                    Copy left of Q (den. by Q1) to right of workspace.
C                    
                     CALL DLACPY( 'All', N, DIM1, Q(1,SROW+1), LDQ,
     $                  DWORK(PDW+(DIM2+OFF)*(N+PARA_LDADD)),
     $                  N+PARA_LDADD )
C                    
C                    Compute Q1*U12.
C                    
                     CALL DTRMM( 'Right', 'Lower', 'No transpose',
     $                  'No transpose', N, DIM1, ONE,
     $                  DWORK(PU12), USIZE+PARA_LDADD,
     $                  DWORK(PDW+(DIM2+OFF)*(N+PARA_LDADD)),
     $                  N+PARA_LDADD )
C                    
C                    Add Q2*U22.
C                    
                     CALL DGEMM( 'No transpose', 'No transpose', 
     $                  N, DIM1, DIM2+OFF, ONE,
     $                  Q(1,SROW+DIM1+1), LDQ, DWORK(PU22),
     $                  USIZE+PARA_LDADD, ONE,
     $                  DWORK(PDW+(DIM2+OFF)*(N+PARA_LDADD)),
     $                  N+PARA_LDADD )
                     CALL DLACPY( 'All', N, USIZE-1, DWORK(PDW),
     $                  N+PARA_LDADD, Q(1,SROW+1), LDQ )
                  END IF
C                 
                  IF ( LCOMPZ ) THEN
C                    
C                    Update Z..
C                    Copy right of Z (den. by Z2) to left of workspace.
C                    
                     CALL DLASET( 'All', N, OFF, ZERO, ZERO, DWORK(PDW),
     $                  N+PARA_LDADD )
                     CALL DLACPY( 'All', N, DIM2, Z(1,SROW+DIM1+1), LDZ,
     $                  DWORK(PDW+OFF*(N+PARA_LDADD)),
     $                  N+PARA_LDADD )
C                    
                     IF ( DIM2.GT.0 ) THEN
C                       
C                       Compute Z2*V21.
C                       
                        CALL DTRMM( 'Right', 'Upper', 'No transpose',
     $                     'No transpose', N, DIM2, ONE,
     $                     DWORK(PV21), USIZE+PARA_LDADD,
     $                     DWORK(PDW+OFF*(N+PARA_LDADD)),
     $                     N+PARA_LDADD )
                     END IF
C                    
C                    Add Z1*V11.
C                    
                     CALL DGEMM( 'No transpose', 'No transpose', N,
     $                  DIM2+OFF, DIM1, ONE, Z(1,SROW+1), LDZ,
     $                  DWORK(PV11), USIZE+PARA_LDADD, ONE,
     $                  DWORK(PDW), N+PARA_LDADD )
C                    
C                    Copy left of Z (den. by Z1) to right of workspace.
C                    
                     CALL DLACPY( 'All', N, DIM1, Z(1,SROW+1), LDZ,
     $                  DWORK(PDW+(DIM2+OFF)*(N+PARA_LDADD)),
     $                  N+PARA_LDADD )
C                    
C                    Compute Z1*V12.
C                    
                     CALL DTRMM( 'Right', 'Lower', 'No transpose',
     $                  'No transpose', N, DIM1, ONE,
     $                  DWORK(PV12), USIZE+PARA_LDADD,
     $                  DWORK(PDW+(DIM2+OFF)*(N+PARA_LDADD)),
     $                  N+PARA_LDADD )
C                    
C                    Add Z2*V22.
C                    
                     CALL DGEMM( 'No transpose', 'No transpose', 
     $                  N, DIM1, DIM2+OFF, ONE,
     $                  Z(1,SROW+DIM1+1), LDZ, DWORK(PV22),
     $                  USIZE+PARA_LDADD, ONE,
     $                  DWORK(PDW+(DIM2+OFF)*(N+PARA_LDADD)),
     $                  N+PARA_LDADD )
                     CALL DLACPY( 'All', N, USIZE-1, DWORK(PDW),
     $                  N+PARA_LDADD, Z(1,SROW+1), LDZ )
                  END IF
               ELSE
                  CALL DGEMM( 'Transpose', 'No transpose', USIZE-1,
     $               UPDHI-EROW, USIZE-1, ONE,
     $               DWORK(PU+USIZE+PARA_LDADD+1),
     $               USIZE+PARA_LDADD,
     $               H(SROW+1,EROW+1), LDH, ZERO,
     $               DWORK(PDW), USIZE+PARA_LDADD-1 )
                  CALL DLACPY( 'All', USIZE-1, UPDHI-EROW, DWORK(PDW),
     $               USIZE+PARA_LDADD-1, H(SROW+1,EROW+1),
     $               LDH )
                  IF ( SROW.GT.UPDLO ) THEN
                     CALL DGEMM( 'No transpose', 'No transpose',
     $                  SROW-UPDLO, USIZE-1, USIZE-1, ONE,
     $                  H(UPDLO,SROW+1), LDH,
     $                  DWORK(PV+USIZE+PARA_LDADD+1),
     $                  USIZE+PARA_LDADD, ZERO, DWORK(PDW),
     $                  SROW+PARA_LDADD-UPDLO ) 
                     CALL DLACPY( 'All', SROW-UPDLO, USIZE-1,
     $                  DWORK(PDW), SROW-UPDLO+PARA_LDADD,
     $                  H(UPDLO,SROW+1), LDH )
                  END IF
                  CALL DGEMM( 'Transpose', 'No transpose', USIZE-1,
     $               UPDHI-EROW, USIZE-1, ONE,
     $               DWORK(PU+USIZE+PARA_LDADD+1),
     $               USIZE+PARA_LDADD,
     $               T(SROW+1,EROW+1), LDT, ZERO,
     $               DWORK(PDW), USIZE+PARA_LDADD-1 )
                  CALL DLACPY( 'All', USIZE-1, UPDHI-EROW, DWORK(PDW),
     $               USIZE+PARA_LDADD-1, T(SROW+1,EROW+1),
     $               LDT )
                  IF ( SROW.GT.UPDLO ) THEN
                     CALL DGEMM( 'No transpose', 'No transpose',
     $                  SROW-UPDLO, USIZE-1, USIZE-1, ONE,
     $                  T(UPDLO,SROW+1), LDT,
     $                  DWORK(PV+USIZE+PARA_LDADD+1),
     $                  USIZE+PARA_LDADD, ZERO, DWORK(PDW),
     $                  SROW+PARA_LDADD-1 ) 
                     CALL DLACPY( 'All', SROW-UPDLO, USIZE-1,
     $                  DWORK(PDW), SROW+PARA_LDADD-1,
     $                  T(UPDLO,SROW+1), LDT )
                  END IF
C                 
                  IF ( LCOMPQ ) THEN
                     CALL DGEMM( 'No transpose', 'No transpose', N,
     $                  USIZE-1, USIZE-1, ONE, Q(1,SROW+1),
     $                  LDQ, DWORK(PU+USIZE+PARA_LDADD+1),
     $                  USIZE+PARA_LDADD, ZERO, DWORK(PDW),
     $                  N+PARA_LDADD )
                     CALL DLACPY( 'All', N, USIZE-1, DWORK(PDW),
     $                  N+PARA_LDADD, Q(1,SROW+1), LDQ )
                  END IF
C                 
                  IF ( LCOMPZ ) THEN
                     CALL DGEMM( 'No transpose', 'No transpose', N,
     $                  USIZE-1, USIZE-1, ONE, Z(1,SROW+1),
     $                  LDZ, DWORK(PV+USIZE+PARA_LDADD+1),
     $                  USIZE+PARA_LDADD, ZERO, DWORK(PDW),
     $                  N+PARA_LDADD )
                     CALL DLACPY( 'All', N, USIZE-1, DWORK(PDW),
     $                  N+PARA_LDADD, Z(1,SROW+1), LDZ )
                  END IF
               END IF
               TIME_BL3 = TIME_BL3 + ( DSECND() - T2 )
               SROW = SROW + PARA_ITEXT
 222        CONTINUE
C           
C           Step 3: Wipe Bumps off.
C           
            USIZE = ACTHI-SROW+1
            PU = 1
            PV = PU + ( USIZE+PARA_LDADD )*USIZE
            PDW = PV + ( USIZE+PARA_LDADD )*USIZE
C           
            CALL DLASET( 'All', USIZE, USIZE, ZERO, ONE, DWORK(PU),
     $         USIZE+PARA_LDADD )
            CALL DLASET( 'All', USIZE, USIZE, ZERO, ONE, DWORK(PV),
     $         USIZE+PARA_LDADD )
            DO 230 F = 1, NBUMPS
               DO 231  G = SROW+(PARA_ITBULG+1)*(NBUMPS-F), ACTHI-2
                  LEN = MIN( PARA_ITBULG+1, ACTHI-G )
                  CALL DLARFG( LEN, H(G+1,G), H(G+2,G), 1, TAU )
                  TEMP = H(G+1,G)
                  H(G+1,G) = ONE
                  CALL DLARFX( 'Left', LEN, ACTHI-G, H(G+1,G),
     $               TAU, H(G+1,G+1), LDH, DWORK(PDW) )
                  CALL DLARFX( 'Left', LEN, ACTHI-G, H(G+1,G),
     $               TAU, T(G+1,G+1), LDT, DWORK(PDW) )
                  CALL DLARFX( 'Right', USIZE, LEN, H(G+1,G),
     $               TAU, DWORK(PU+(G-SROW+1)*
     $               (USIZE+PARA_LDADD)), USIZE+PARA_LDADD,
     $               DWORK(PDW) )
                  H(G+1,G) = TEMP
                  DO 238 J = 2, LEN
                     H(G+J,G) = ZERO
 238              CONTINUE
C                 
                  CALL INVHSE( LEN, T(G+1,G+1), LDT, TAU, V, IERR )
C                 
                  CALL DLARFX( 'Right', MIN( N, G+PARA_ITBULG+1 )
     $               -SROW+1, LEN, V, TAU, T(SROW,G+1), LDT,
     $               DWORK(PDW) )
                  DO 438 J = 2, LEN
                     T(G+J,G+1) = ZERO
 438              CONTINUE
                  CALL DLARFX( 'Right', MIN( N, G+PARA_ITBULG+2 )
     $               -SROW+1, LEN, V, TAU, H(SROW,G+1), LDH,
     $               DWORK(PDW) )
                  CALL DLARFX( 'Right', USIZE, LEN, V, TAU,
     $               DWORK(PV+(G-SROW+1)*(USIZE+PARA_LDADD)),
     $               USIZE+PARA_LDADD, DWORK(PDW) )
 231           CONTINUE
 230        CONTINUE
C           
            T2 = DSECND()
C           
C           Update rest of H.
C           
            IF ( UPDHI.GT.ACTHI ) THEN
               CALL DGEMM( 'Transpose', 'No transpose', USIZE,
     $            UPDHI-ACTHI, USIZE, ONE, DWORK(PU),
     $            USIZE+PARA_LDADD, H(SROW,ACTHI+1),
     $            LDH, ZERO, DWORK(PDW), USIZE+PARA_LDADD )
               CALL DLACPY( 'All', USIZE, UPDHI-ACTHI, DWORK(PDW),
     $            USIZE+PARA_LDADD, H(SROW,ACTHI+1), LDH )
            END IF
C           
            IF ( SROW.GT.UPDLO ) THEN
               CALL DGEMM( 'No transpose', 'No transpose', SROW-UPDLO,
     $            USIZE, USIZE, ONE, H(UPDLO,SROW), LDH,
     $            DWORK(PV), USIZE+PARA_LDADD, ZERO,
     $            DWORK(PDW), SROW-UPDLO+PARA_LDADD )
               CALL DLACPY( 'All', SROW-1, USIZE, DWORK(PDW),
     $            SROW-UPDLO+PARA_LDADD, H(UPDLO,SROW), LDH )
            END IF
C           
C           Update rest of T.
C           
            IF ( UPDHI.GT.ACTHI ) THEN
               CALL DGEMM( 'Transpose', 'No transpose', USIZE,
     $            UPDHI-ACTHI, USIZE, ONE, DWORK(PU),
     $            USIZE+PARA_LDADD, T(SROW,ACTHI+1),
     $            LDT, ZERO, DWORK(PDW), USIZE+PARA_LDADD )
               CALL DLACPY( 'All', USIZE, UPDHI-ACTHI, DWORK(PDW),
     $            USIZE+PARA_LDADD, T(SROW,ACTHI+1), LDT )
            END IF
C           
            IF ( SROW.GT.UPDLO ) THEN
               CALL DGEMM( 'No transpose', 'No transpose', SROW-UPDLO,
     $            USIZE, USIZE, ONE, T(UPDLO,SROW), LDT,
     $            DWORK(PV), USIZE+PARA_LDADD, ZERO,
     $            DWORK(PDW), SROW-UPDLO+PARA_LDADD )
               CALL DLACPY( 'All', SROW-UPDLO, USIZE, DWORK(PDW),
     $            SROW-UPDLO+PARA_LDADD, T(UPDLO,SROW), LDT )
            END IF
C           
C           Update Q and Z.
C           
            IF ( LCOMPQ) THEN
               CALL DGEMM( 'No transpose', 'No transpose', N, USIZE,
     $            USIZE, ONE, Q(1,SROW), LDQ, DWORK(PU),
     $            USIZE+PARA_LDADD, ZERO, DWORK(PDW),
     $            N+PARA_LDADD )
               CALL DLACPY( 'All', N, USIZE, DWORK(PDW), N+PARA_LDADD,
     $            Q(1,SROW), LDQ )
            END IF
C           
            IF ( LCOMPZ ) THEN
               CALL DGEMM( 'No transpose', 'No transpose', N, USIZE,
     $            USIZE, ONE, Z(1,SROW), LDZ, DWORK(PV),
     $            USIZE+PARA_LDADD, ZERO, DWORK(PDW),
     $            N+PARA_LDADD )
               CALL DLACPY( 'All', N, USIZE, DWORK(PDW), N+PARA_LDADD,
     $            Z(1,SROW), LDZ )
            END IF
C           
            TIME_BL3 = TIME_BL3 + ( DSECND() - T2 )
            OVRSHF   = OVRSHF + PARA_ITNS
            AGGSHF   = AGGSHF + PARA_ITNS
            TIME_KK  = TIME_KK + ( DSECND() - T1 )

C           
         ELSE IF ( ACTHI-ACTLO+1.GE.PARA_ITDACK ) THEN
            IF ( MSGLVL.GE.2 )
     $         WRITE(*,'(A,I4,A,I4)') 'Apply DACK QZ to: ', ACTLO,
     $         ' ... ', ACTHI
            T1 = DSECND()
            CALL QZDACKIT( LSCHUR, LCOMPQ, LCOMPZ, LEXC, N, ACTLO,
     $         ACTHI, ESHIFT, H, LDH, T, LDT, Q, LDQ, Z,
     $         LDZ, DWORK, LDWORK, PARA_ITDANB, HFNORM,
     $         TFNORM, NSHF )
            TIME_DACK = TIME_DACK + ( DSECND() - T1 )
            OVRSHF = OVRSHF + NSHF
            AGGSHF = AGGSHF + NSHF
C           
         ELSE IF ( ACTHI-ACTLO+1.GT.2 ) THEN
C           
            IF ( MSGLVL.GE.2 )
     $         WRITE(*,'(A,I4,A,I4)') 'Apply LAPACK QZ to: ', ACTLO,
     $         ' ... ', ACTHI
C           
            T1 = DSECND()
            CALL QZLAPIT( LSCHUR, LCOMPQ, LCOMPZ, LEXC, N, ACTLO, ACTHI,
     $         ESHIFT, H, LDH, T, LDT, Q, LDQ, Z, LDZ, HFNORM,
     $         TFNORM, NSHF )
            TIME_MOLS = TIME_MOLS + ( DSECND() - T1 )
            OVRSHF = OVRSHF + NSHF
            AGGSHF = AGGSHF + NSHF
C           
         ELSE IF ( ACTHI-ACTLO+1.GT.1 ) THEN
C           
            IF ( MSGLVL.GE.2 )
     $         WRITE(*,'(A,I4)') 'Standardize 2x2 block at:    ',ACTLO
C           
            CALL QZLAP( JOB, UPDQ, UPDZ, N, ACTLO, ACTHI, H, LDH, T,
     $         LDT, ALPHAR, ALPHAI, BETA, Q, LDQ, Z, LDZ,
     $         DWORK, LDWORK, IERR )
C           
            IF ( IERR.NE.0 ) THEN
               IF ( MSGLVL.GE.1 )
     $            WRITE(*,*) 'Call to QZLAP failed'
               INFO = 1
               RETURN
            END IF
C           
            ACTHI = ACTLO-1
            ACTLO = ILO
            IF ( .NOT.LSCHUR ) THEN
               UPDLO = ACTLO
               UPDHI = ACTHI
            END IF
            EXCIT = 0
            ESHIFT = ZERO
C           
         ELSE
            IF ( MSGLVL.GE.2 )
     $         WRITE(*,'(A,I4)') 'Standardize eigenvalue at:', ACTLO
            IF ( T( ACTLO,ACTLO ).LT.ZERO ) THEN
               IF ( LSCHUR ) THEN
                  CALL DSCAL( ACTLO, -ONE, H(1,ACTLO), 1 )
                  CALL DSCAL( ACTLO, -ONE, T(1,ACTLO), 1 )
               ELSE
                  H(ACTLO,ACTLO) = -H(ACTLO,ACTLO)
                  T(ACTLO,ACTLO) = -T(ACTLO,ACTLO)
               END IF
               IF ( LCOMPZ )
     $            CALL DSCAL( N, -ONE, Z(1,ACTLO), 1 )
            END IF
            ALPHAR(ACTLO) = H(ACTLO,ACTLO)
            ALPHAI(ACTLO) = ZERO
            BETA(ACTLO) = T(ACTLO,ACTLO)
C           
            ACTHI = ACTLO-1
            ACTLO = ILO
            IF ( .NOT.LSCHUR ) THEN
               UPDLO = ACTLO
               UPDHI = ACTHI
            END IF
            EXCIT = 0
            ESHIFT = ZERO
C           
         END IF
C        
C        Invalidate shifts and check whether aggressive early deflation
C        shall be applied in the next step.
C        
         P = N + 1
         IF ( AGGSHF.GE.PARA_AGGSHF ) THEN
            DOAGGR = .TRUE.
            AGGSHF = 0
         END IF
C        
         GO TO 42
      END IF
C     
C     end while
C     ------------------------------------------------------------------
C     End of main loop.
C     ------------------------------------------------------------------
C     Number of iterations exceeded the maximum.
      


C     
C     
C     



      TIME_ALL = DSECND() - TIME_ALL
      IF ( MSGLVL.GE.1 ) THEN
         WRITE(*,'(A,F8.2,A)') 'Overall compt. time      : ',
     $      TIME_ALL, ' sec.'
         WRITE(*,'(A,F8.2,A)') 'Time for BLAS 3          : ',
     $      TIME_BL3, ' sec.'
         WRITE(*,'(A,F8.2,A)') 'Time for KK QZ           : ',
     $      TIME_KK, ' sec.'
         WRITE(*,'(A,F8.2,A)') 'Time for DacK QZ         : ',
     $      TIME_DACK, ' sec.'
         WRITE(*,'(A,F8.2,A)') 'Time for LAPACK QZ       : ',
     $      TIME_MOLS, ' sec.'
         WRITE(*,'(A,F8.2,A)') 'Time for early deflation : ',
     $      TIME_DEFL, ' sec.'
         WRITE(*,'(A,F8.2,A)') 'Time for inf. deflation  : ',
     $      TIME_INF, ' sec.'
         WRITE(*, '(A,F8.2)')  'Number of shifts / N     : ',
     $      DBLE( OVRSHF ) / DBLE ( N )

      END IF
C     
C     Last line of KKQZ
C     
      END
