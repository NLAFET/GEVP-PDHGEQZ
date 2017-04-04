      SUBROUTINE QZINF( UPDHT, UPDQ, UPDZ, N, ILO, IHI, TNORM, 
     $   H, LDH, T, LDT, Q, LDQ, Z, LDZ, DWORK, LDWORK, IWORK,
     $   LIWORK, INFO )
      IMPLICIT NONE
C     
C     PURPOSE
C     
C     This routine deflates infinite eigenvalues from a matrix pair
C     (H,T) in Hessenberg-triangular form.
C     
C     ARGUMENTS
C     
C     Mode Parameters
C     
C     UPDHT   LOGICAL
C     Specifies how the matrices H and T are to be updated:
C     = .F.:  update only the active submatrix pair
C     ( H(ILO:IHI,ILO:IHI), T(ILO:IHI,ILO:IHI) );
C     = .T.:  update the complete matrices H and T.
C     
C     UPDQ    LOGICAL
C     = .F.:  the orthogonal factor Q is not updated;
C     = .T.:  Q is updated.
C     
C     UPDZ    LOGICAL
C     = .F.:  the orthogonal factor Z is not updated;
C     = .T.:  Z is updated.
C     
C     Input/Output Parameters
C     
C     N       (input) INTEGER
C     The order of the matrices H and T.  N >= 0.
C     
C     ILO     (input/output) INTEGER
C     IHI     (input/output) INTEGER
C     On entry, it is assumed that H and T are already upper
C     triangular in rows and columns 1:ILO-1 and IHI+1:N. The
C     deflation of infinite eigenvalues is restricted to the
C     active submatrix pair ( H(ILO:IHI,ILO:IHI),
C     T(ILO:IHI,ILO:IHI) ).
C     1 <= ILO <= IHI <= N, if N > 0; ILO=1 and IHI=0, if N=0. 
C     
C     TNORM   (input) DOUBLE PRECISION
C     Frobenius norm of T.
C     
C     H       (input/output) DOUBLE PRECISION array, dimension (LDH,N)
C     On entry, the leading n-by-n part of this array contains
C     an upper Hessenberg matrix H.
C     On exit, the leading n-by-n part of this array contains
C     the updated upper Hessenberg matrix H.
C     
C     LDH     (input) INTEGER
C     The leading dimension of the array H. LDH >= max(1,N).
C     
C     T       (input/output) DOUBLE PRECISION array, dimension (LDT,N)
C     On entry, the leading n-by-n part of this array contains
C     an upper triangular matrix T.
C     On exit, the leading n-by-n part of this array contains
C     the updated upper triangular matrix T.
C     
C     LDT     (input) INTEGER
C     The leading dimension of the array T. LDT >= max(1,N).
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
C     Error Indicator
C     
C     INFO    INTEGER
C     = 0:  successful exit;
C     < 0:  if INFO = -i, the i-th argument had an illegal 
C     value.
C     
C     CONTRIBUTOR
C     
C     D. Kressner, Univ. Umea, Sweden, November 2004.
C     
C     ******************************************************************
C     
C     
C     .. Parameters ..
      DOUBLE PRECISION  ZERO, ONE
      PARAMETER         ( ZERO = 0.0D0, ONE = 1.0D0 )
C     .. Configuration parameters ..
      LOGICAL           PARA_INFALL, PARA_INFBLOCK
      INTEGER           PARA_INFDEF, PARA_INFDIR, PARA_INFNZ,
     $   PARA_INFWDW
      DOUBLE PRECISION  PARA_INFTOL
      COMMON /PARA_INF/ PARA_INFALL, PARA_INFBLOCK, PARA_INFDEF, 
     $   PARA_INFDIR, PARA_INFNZ, PARA_INFWDW, PARA_INFTOL

      LOGICAL           PARA_ITTRI
      INTEGER           PARA_ITASH, PARA_ITBULG, PARA_ITDACK, 
     $   PARA_ITDANB, PARA_ITEXC, PARA_ITEXT, PARA_ITKK, PARA_ITMAX, 
     $   PARA_ITNS, PARA_ITOPP
      COMMON /PARA_IT/  PARA_ITASH, PARA_ITBULG, PARA_ITDACK, 
     $   PARA_ITDANB, PARA_ITEXC, PARA_ITEXT, PARA_ITKK, PARA_ITMAX, 
     $   PARA_ITNS, PARA_ITOPP, PARA_ITTRI

C     .. Precision paramters ..
      DOUBLE PRECISION  SMLNUM
      COMMON /PREC/     SMLNUM
C     .. Scalar Arguments ..
      LOGICAL           UPDHT, UPDQ, UPDZ
      INTEGER           IHI, ILO, INFO, LDH, LDQ, LDT, LDWORK, LDZ,
     $   LIWORK, N
      DOUBLE PRECISION  TNORM
C     .. Array Arguments ..
      INTEGER           IWORK(*)
      DOUBLE PRECISION  DWORK(*), H(LDH,*), Q(LDQ,*), T(LDT,*), Z(LDZ,*)
C     .. Local Scalars ..
      INTEGER           BPOS, DEF, I, IFRST, ILAST, J, LGTHL, LGTHR,
     $   TPOS
      DOUBLE PRECISION  C, OVFL, S, TEMP, UNFL
C     .. External Functions ..
      LOGICAL           INFDEF
      DOUBLE PRECISION  DLAMCH
      EXTERNAL          DLAMCH, INFDEF

C     .. External Subroutines ..
      EXTERNAL          DLABAD, DLARTG, DROT
C     .. Intrinsic Functions ..
      INTRINSIC         ABS, MAX, MIN
C     
C     .. Executable Statements ..
C     
      INFO = 0
C     
C     Quick return if possible.
C     
      IF ( PARA_INFDEF.EQ.0 .OR. IHI-ILO+1.EQ.1 )
     $   RETURN
C     
      UNFL = DLAMCH( 'Safe minimum' )
      OVFL = ONE / UNFL
      CALL DLABAD( UNFL, OVFL )
      SMLNUM = UNFL*( ( IHI-ILO+1 ) / DLAMCH( 'Precision' ) ) 
      IFRST = 1
      ILAST = N
C     
      IF ( .NOT.PARA_INFBLOCK ) THEN
C        
C        Unblocked deflation algorithm.
C        
 10      CONTINUE
C        Do
C        
C        Try to deflate an infinite eigenvalues at top left corner
C        (search at positions ilo, ..., j).
C        
         IF ( .NOT.PARA_INFALL ) THEN
            IF ( IHI-ILO+1.LT.PARA_ITKK ) THEN
               J = MIN( IHI, ILO+1 )
            ELSE
               J = MIN( IHI, ILO +
     $            ( PARA_ITNS/PARA_ITBULG )*( PARA_ITBULG+1 ) - 2 )
            END IF
         ELSE
            IF ( PARA_INFDIR.EQ.1 ) THEN
               J = IHI
            ELSE IF ( PARA_INFDIR.EQ.2 ) THEN
               J = ILO-1
            ELSE
               J = ILO + ( IHI-ILO ) / 2
            END IF
         END IF
C        
         TPOS = 0
         DO 30  I = ILO, J
            IF ( INFDEF( IHI-ILO+1, I-ILO+1, TNORM, T(ILO,ILO),
     $         LDT ) ) THEN
               T(I,I) = ZERO
               TPOS = I
               GO TO 40
            END IF
 30      CONTINUE
C        
         GO TO 60
C        
 40      CONTINUE

C        
C        Found an infinite eigenvalue at position TPOS. Chase this
C        eigenvalue to top left corner.
C        
         IF ( .NOT.UPDHT ) THEN
            IFRST = ILO
            ILAST = IHI   
         END IF
C        
         DO 50  I = TPOS-1, ILO, -1
            TEMP = T(I,I+1)
            CALL DLARTG( TEMP, T(I,I), C, S, T(I,I+1) )
            T(I,I) = ZERO 
            CALL DROT( I-IFRST, T(IFRST,I), 1, T(IFRST,I+1), 1,
     $         C, -S )
            CALL DROT( MIN(IHI,I+2)-IFRST+1, H(IFRST,I), 1,
     $         H(IFRST,I+1), 1, C, -S )
            IF ( UPDZ )
     $         CALL DROT( N, Z(1,I), 1, Z(1,I+1), 1, C, -S )
            IF ( IHI.GE.I+2 ) THEN
               TEMP = H(I+1,I)
               CALL DLARTG( TEMP, H(I+2,I), C, S, H(I+1,I) )
               H(I+2,I) = ZERO
               CALL DROT( ILAST-I, H(I+1,I+1), LDH, H(I+2,I+1), LDH,
     $            C, S )
               CALL DROT( ILAST-I-1, T(I+1,I+2), LDT, T(I+2,I+2),
     $            LDT, C, S )
               IF ( UPDQ )
     $            CALL DROT( N, Q(1,I+1), 1, Q(1,I+2), 1, C, S )
            END IF
 50      CONTINUE
C        
C        Deflate infinite eigenvalue at top left corner.
C        
         IF (ILO.LT.IHI) THEN
            TEMP = H(ILO,ILO)
            CALL DLARTG( TEMP, H(ILO+1,ILO), C, S, H(ILO,ILO) )
            H(ILO+1,ILO) = ZERO
            CALL DROT( ILAST-ILO, H(ILO,ILO+1), LDH, 
     $         H(ILO+1,ILO+1),
     $         LDH, C, S )
            CALL DROT( ILAST-ILO, T(ILO,ILO+1), LDT, T(ILO+1,ILO+1),
     $         LDT, C, S )
            IF ( UPDQ )
     $         CALL DROT( N, Q(1,ILO), 1, Q(1,ILO+1), 1, C, S )
C           

         ELSE
            TEMP = H(IHI,IHI)
            CALL DLARTG( TEMP, H(IHI,IHI-1), C, S, H(IHI,IHI) )
            H(IHI,IHI-1 ) = ZERO
            CALL DROT( IHI-IFRST, H(IFRST,IHI), 1, H(IFRST,IHI-1), 1,
     $         C, S )
            CALL DROT( IHI-IFRST, T(IFRST,IHI), 1, T(IFRST,IHI-1), 1,
     $         C, S )
            IF ( UPDZ )
     $         CALL DROT( N, Z(1,IHI), 1, Z(1,IHI-1), 1, C, S ) 
         END IF
         ILO = ILO + 1
C        
         GO TO 90
C        
 60      CONTINUE
C        
C        Try to deflate infinite eigenvalues at bottom right
C        corner (search at positions ihi, ihi-1,..., j).
C        
         IF ( .NOT.PARA_INFALL ) THEN
            IF ( IHI-ILO+1.LT.PARA_ITKK ) THEN
               J = MAX( IHI-1, ILO )
            ELSE
               J = MAX( IHI-PARA_ITNS+1, ILO )
            END IF
         ELSE
            IF ( PARA_INFDIR.EQ.1 ) THEN
               J = IHI+1
            ELSE IF ( PARA_INFDIR.EQ.2 ) THEN
               J = ILO
            ELSE
               J = ILO + ( IHI-ILO ) / 2 + 1
            END IF
         END IF
C        
         BPOS = 0
         DO 70  I = IHI, J, -1
            IF ( INFDEF( IHI-ILO+1, I-ILO+1, TNORM, T(ILO,ILO),
     $         LDT ) ) THEN
               T(I,I) = ZERO
               BPOS = I
               GO TO 80
            END IF
 70      CONTINUE
C        
         GO TO 90
C        
 80      CONTINUE
C        
C        Found an infinite eigenvalue at position BPOS. Chase this
C        eigenvalue to bottom right corner.
C        
         IF ( .NOT.UPDHT ) THEN
            IFRST = ILO
            ILAST = IHI   
         END IF
C        
         DO 85 I = BPOS, IHI-1
            TEMP = T(I,I+1 )
            CALL DLARTG( TEMP, T(I+1,I+1 ), C, S, T(I,I+1) )
            T(I+1,I+1 ) = ZERO
            IF ( I.LT.ILAST-1 )
     $         CALL DROT( ILAST-I-1, T(I,I+2), LDT, T(I+1,I+2), LDT,
     $         C, S )
            CALL DROT( ILAST-MAX(I-1,ILO)+1, H(I,MAX(I-1,ILO)), LDH,
     $         H(I+1,MAX(I-1,ILO)), LDH, C, S )
            IF ( UPDQ )
     $         CALL DROT( N, Q(1,I), 1, Q(1,I+1), 1, C, S )
            IF ( I.GT.ILO ) THEN
               TEMP = H(I+1,I)
               CALL DLARTG( TEMP, H(I+1,I-1), C, S, H(I+1,I) )
               H(I+1,I-1 ) = ZERO
               CALL DROT( I+1-IFRST, H(IFRST,I), 1, H(IFRST,I-1), 1,
     $            C, S )
               CALL DROT( I-IFRST, T(IFRST,I), 1, T(IFRST,I-1), 1,
     $            C, S )
               IF ( UPDZ )
     $            CALL DROT( N, Z(1,I), 1, Z(1,I-1), 1, C, S )
            END IF
 85      CONTINUE 
C        
C        Deflate infinite eigenvalues at bottom right corner.
C        
         IF (IHI.GT.ILO) THEN 
            TEMP = H(IHI,IHI)
            CALL DLARTG( TEMP, H(IHI,IHI-1), C, S, H(IHI,IHI) )
            H(IHI,IHI-1 ) = ZERO
            CALL DROT( IHI-IFRST, H(IFRST,IHI), 1, H(IFRST,IHI-1), 1,
     $         C, S )
            CALL DROT( IHI-IFRST, T(IFRST,IHI), 1, T(IFRST,IHI-1), 1,
     $         C, S )
            IF ( UPDZ )
     $         CALL DROT( N, Z(1,IHI), 1, Z(1,IHI-1), 1, C, S ) 
C           
         ELSE
            TEMP = H(ILO,ILO)
            CALL DLARTG( TEMP, H(ILO+1,ILO), C, S, H(ILO,ILO) )
            H(ILO+1,ILO) = ZERO
            CALL DROT( ILAST-ILO, H(ILO,ILO+1), LDH,
     $         H(ILO+1,ILO+1),
     $         LDH, C, S )
            CALL DROT( ILAST-ILO, T(ILO,ILO+1), LDT, T(ILO+1,ILO+1),
     $         LDT, C, S )
            IF ( UPDQ )
     $         CALL DROT( N, Q(1,ILO), 1, Q(1,ILO+1), 1, C, S )
         END IF
         IHI = IHI - 1
C        
 90      CONTINUE
C        while TPOS<>0 or BPOS<>0.
         IF ( (TPOS.NE.0 .OR. BPOS.NE.0) 
     $      .AND. (IHI .GE. ILO) .AND. (ILO.LE.IHI ))  GO TO 10
      ELSE
C        
C        Blocked deflation algorithm.
C        
 110     CONTINUE
C        Do
C        
C        Try to deflate infinite eigenvalues at top left corner.
C        
         IF ( .NOT.PARA_INFALL ) THEN
            IF ( IHI-ILO+1.LT.PARA_ITKK ) THEN
               J = MIN( IHI, ILO+1 )
            ELSE
               J = MIN( IHI, ILO +
     $            ( PARA_ITNS/PARA_ITBULG )*( PARA_ITBULG+1 ) - 2 )
            END IF
         ELSE
            IF ( PARA_INFDIR.EQ.1 ) THEN
               J = IHI
            ELSE IF ( PARA_INFDIR.EQ.2 ) THEN
               J = ILO-1
            ELSE
               J = ILO + ( IHI - ILO ) / 2
            END IF
         END IF
C        
         DEF = 0
         I = ILO
C        
 120     CONTINUE
         IF ( DEF .LT. PARA_INFNZ .AND. I .LE. J ) THEN
C           while  def < PARA_INFNZ  and  i <= j  do
            IF ( INFDEF( IHI - ILO + 1, I - ILO + 1, TNORM, 
     $         T( ILO, ILO ), LDT ) ) THEN
               T( I, I ) = ZERO
               DEF = DEF + 1
               I = I + 2
            ELSE
               I = I + 1
            END IF
            GO TO 120
C           end while
         END IF
C        
         IF ( DEF.GT.0 ) THEN
C           
            BPOS = MIN( I - 1, IHI )
            TPOS = MAX( ILO, BPOS - PARA_INFWDW )
            IF ( .NOT. UPDHT ) THEN
               IFRST = ILO
               ILAST = IHI   
            END IF
C           
 130        CONTINUE
C           Do
            CALL UPWD( ( TPOS.EQ.ILO ), ( BPOS.EQ.IHI ),
     $         BPOS-TPOS+1, H(TPOS,TPOS),
     $         LDH, T(TPOS,TPOS), LDT, DEF, LGTHL,
     $         IWORK, DWORK, LGTHR,
     $         IWORK(PARA_INFNZ*PARA_INFWDW+1),
     $         DWORK(2*PARA_INFNZ*PARA_INFWDW+1) )
            CALL ROTAPPL( .FALSE., LGTHR, TPOS-IFRST,
     $         IWORK(PARA_INFNZ*PARA_INFWDW+1),
     $         DWORK(2*PARA_INFNZ*PARA_INFWDW+1),
     $         H(IFRST,TPOS), LDH )
            IF ( BPOS.LT.ILAST )
     $         CALL ROTAPPL( .TRUE., LGTHL, ILAST-BPOS, IWORK,
     $         DWORK, H(TPOS,BPOS+1), LDH )

            CALL ROTAPPL( .FALSE., LGTHR, TPOS-IFRST, 
     $         IWORK(PARA_INFNZ*PARA_INFWDW+1),
     $         DWORK(2*PARA_INFNZ*PARA_INFWDW+1),
     $         T(IFRST,TPOS), LDT )
            IF ( BPOS.LT.ILAST )
     $         CALL ROTAPPL( .TRUE., LGTHL, ILAST-BPOS, IWORK,
     $         DWORK, T(TPOS,BPOS+1), LDT )
            IF ( UPDQ )
     $         CALL ROTAPPL( .FALSE., LGTHL, N, IWORK, DWORK,
     $         Q(1,TPOS), LDQ )
            IF ( UPDZ )
     $         CALL ROTAPPL( .FALSE., LGTHR, N,
     $         IWORK(PARA_INFNZ*PARA_INFWDW+1),
     $         DWORK(2*PARA_INFNZ*PARA_INFWDW+1),
     $         Z(1,TPOS), LDZ )
            IF ( TPOS .EQ. ILO )
     $         ILO = ILO + DEF

C           WRITE(*,*)TPOS, ILO, DEF
C           
C           while TPOS >= ILO
            IF ( TPOS .GT. ILO ) THEN
               IF ( T( TPOS - 1,TPOS - 1 ).EQ.ZERO ) THEN
                  BPOS = MIN( TPOS + DEF * 2, IHI )
               ELSE
                  BPOS = MIN( TPOS + DEF * 2 - 1, IHI )
               END IF
               TPOS = MAX( ILO, BPOS - PARA_INFWDW )
C              WRITE(*,*)TPOS, BPOS, ILO, IHI, DEF
               GO TO 130
            END IF
         END IF
C        
         IF ( DEF.NE.0 .AND. ILO .LE. IHI) GO TO 190
C        
C        Try to deflate infinite eigenvalues at bottom right corner.
C        
         IF ( .NOT.PARA_INFALL ) THEN
            IF ( IHI-ILO+1.LT.PARA_ITKK ) THEN
               J = MAX( IHI-1, ILO )
            ELSE
               J = MAX( IHI-PARA_ITNS+1, ILO )
            END IF
         ELSE
            IF ( PARA_INFDIR.EQ.1 ) THEN
               J = IHI+1
            ELSE IF ( PARA_INFDIR.EQ.2 ) THEN
               J = ILO
            ELSE
               J = ILO + ( IHI-ILO ) / 2 + 1
            END IF
         END IF
C        
         DEF = 0
         I = IHI
C        
 140     CONTINUE
         IF ( DEF.LT.PARA_INFNZ .AND. I.GE.J ) THEN
C           while  def < PARA_INFNZ  and  i >= j  do
            IF ( INFDEF( IHI-ILO+1, I-ILO+1, TNORM, T(ILO,ILO),
     $         LDT ) ) THEN
               T(I,I) = ZERO
               DEF = DEF + 1
               I = I - 2
            ELSE
               I = I - 1
            END IF
            GO TO 140
C           end while
         END IF
C        
         IF ( DEF.GT.0 ) THEN
C           
            TPOS = MAX( I+1, ILO )
            BPOS = MIN( IHI, TPOS + PARA_INFWDW )
            IF ( .NOT.UPDHT ) THEN
               IFRST = ILO
               ILAST = IHI   
            END IF
C           
 150        CONTINUE
C           Do
            CALL DOWNWD( ( TPOS.EQ.ILO ), ( BPOS.EQ.IHI ),
     $         BPOS-TPOS+1, H(TPOS,TPOS),
     $         LDH, T(TPOS,TPOS), LDT, DEF, LGTHL,
     $         IWORK, DWORK, LGTHR,
     $         IWORK(PARA_INFNZ*PARA_INFWDW+1),
     $         DWORK(2*PARA_INFNZ*PARA_INFWDW+1) )
            CALL ROTAPPL( .FALSE., LGTHR, TPOS-IFRST,
     $         IWORK(PARA_INFNZ*PARA_INFWDW+1),
     $         DWORK(2*PARA_INFNZ*PARA_INFWDW+1),
     $         H(IFRST,TPOS), LDH )
            IF ( BPOS.LT.ILAST )
     $         CALL ROTAPPL( .TRUE., LGTHL, ILAST-BPOS, IWORK,
     $         DWORK, H(TPOS,BPOS+1), LDH )
            CALL ROTAPPL( .FALSE., LGTHR, TPOS-IFRST, 
     $         IWORK(PARA_INFNZ*PARA_INFWDW+1),
     $         DWORK(2*PARA_INFNZ*PARA_INFWDW+1),
     $         T(IFRST,TPOS), LDT )
            IF ( BPOS.LT.ILAST )
     $         CALL ROTAPPL( .TRUE., LGTHL, ILAST-BPOS, IWORK,
     $         DWORK, T(TPOS,BPOS+1), LDT )
            IF ( UPDQ )
     $         CALL ROTAPPL( .FALSE., LGTHL, N, IWORK, DWORK,
     $         Q(1,TPOS), LDQ )
            IF ( UPDZ )
     $         CALL ROTAPPL( .FALSE., LGTHR, N,
     $         IWORK(PARA_INFNZ*PARA_INFWDW+1),
     $         DWORK(2*PARA_INFNZ*PARA_INFWDW+1),
     $         Z(1,TPOS), LDZ )
            IF ( BPOS .EQ. IHI )
     $         IHI = IHI - DEF
C           
C           while BPOS <= IHI
            IF ( BPOS .LT. IHI ) THEN
               IF ( T(BPOS+1,BPOS+1).EQ.ZERO ) THEN
                  TPOS = MAX( BPOS - DEF * 2, ILO )
               ELSE
                  TPOS = MAX( BPOS - DEF * 2 + 1, ILO )
               END IF
               BPOS = MIN( IHI, TPOS + PARA_INFWDW )
               GO TO 150
            END IF
         END IF
 190     CONTINUE
C        PRINT*, '(ILO,IHI) = ', '(', ILO, ',', IHI, ')'
C        
C        while DEF<>0.
         IF ( DEF.NE.0 .AND. IHI .GE. ILO )  GO TO 110
C        
      END IF
C     
      RETURN
C     *** Last line of QZINF ***
      END
C     
      SUBROUTINE UPWD( ISTOP, ISBOT, N, H, LDH, T, LDT, DEF, LGTHL,
     $   IDXL, CSL, LGTHR, IDXR, CSR )
      IMPLICIT NONE
C     
C     Chases all infinite eigenvalues to the top of the submatrix pair
C     (H,T) and stores transformations from the left in IDXL(1:LGTHL),
C     CSL(1:2*LGTHL), and transformations from the right in
C     IDXR(1:LGTHR), CSR(1:2*LGTHR).
C     
C     If ISTOP, this submatrix pair resides at the top left corner
C     of the original matrix pair.
C     If ISBOT, this submatrix pair resides at the bottom right corner
C     of the original matrix pair.
C     
C     
C     .. Parameters ..
      DOUBLE PRECISION  ZERO
      PARAMETER         ( ZERO = 0.0D0 )
C     .. Configuration parameters ..
      LOGICAL           PARA_INFALL, PARA_INFBLOCK
      INTEGER           PARA_INFDEF, PARA_INFDIR, PARA_INFNZ,
     $   PARA_INFWDW
      DOUBLE PRECISION  PARA_INFTOL
      COMMON /PARA_INF/ PARA_INFALL, PARA_INFBLOCK, PARA_INFDEF, 
     $   PARA_INFDIR, PARA_INFNZ, PARA_INFWDW, PARA_INFTOL
C     .. Scalar Arguments ..
      LOGICAL           ISBOT, ISTOP
      INTEGER           LDH, LDT, LGTHL, LGTHR, N
C     .. Array Arguments ..
      DOUBLE PRECISION  CSL(*), CSR(*), H(LDH,*), T(LDT,*)
      INTEGER           IDXL(*), IDXR(*)
C     .. Local Scalars ..
      INTEGER           DEF, I, J, MAXJ, TOP
      DOUBLE PRECISION  C, S, TEMP
C     .. External Subroutines ..
      EXTERNAL          DLARTG, DROT
C     .. External Functions ..
      LOGICAL           INFDEF
      EXTERNAL          INFDEF
C     .. Intrinsic Functions ..
      INTRINSIC         MAX
C     
C     .. Executable Statements ..
C     
      LGTHL = 0
      LGTHR = 0
      IF ( ISBOT ) THEN
         MAXJ = N
      ELSE
         MAXJ = N - 1
      END IF
C     J = search position
      J = 1
C     TOP = position to chase next eigenvalue to
      TOP = 1
C     DEF = number of deflatable eigenvalues found
      DEF = 0
C     
 10   CONTINUE
C     do      
      IF ( T(J,J).EQ.ZERO ) THEN
C        
C        Chase zero to position TOP.
C        
         DEF = DEF + 1
         DO 50  I = J-1, TOP, -1
            TEMP = T(I,I+1)
            CALL DLARTG( TEMP, T(I,I), C, S, T(I,I+1) )
            T(I,I) = ZERO
            S = -S 
            CALL DROT( I-1, T(1,I), 1, T(1,I+1), 1, C, S )
            CALL DROT( MIN(N,I+2), H(1,I), 1, H(1,I+1), 1, C, S )
            LGTHR = LGTHR + 1
            IDXR(LGTHR) = I
            CSR(LGTHR*2 - 1) = C
            CSR(LGTHR*2) = S
            IF ( N.GE.I+2 ) THEN
               TEMP = H(I+1,I)
               CALL DLARTG( TEMP, H(I+2,I), C, S, H(I+1,I) )
               H(I+2,I) = ZERO
               CALL DROT( N-I, H(I+1,I+1), LDH, H(I+2,I+1), LDH,
     $            C, S )
               CALL DROT( N-I-1, T(I+1,I+2), LDT, T(I+2,I+2),
     $            LDT, C, S )
               LGTHL = LGTHL + 1
               IDXL(LGTHL) = I + 1
               CSL(LGTHL*2 - 1) = C
               CSL(LGTHL*2) = S
            END IF
 50      CONTINUE
C        
         IF ( ISTOP ) THEN
C           
C           Deflate infinite eigenvalue at top left corner.
C           
            TEMP = H(TOP,TOP)
            CALL DLARTG( TEMP, H(TOP+1,TOP), C, S, H(TOP,TOP) )
            H(TOP+1,TOP) = ZERO
            CALL DROT( N-TOP, H(TOP,TOP+1), LDH, H(TOP+1,TOP+1), LDH,
     $         C, S )
            CALL DROT( N-TOP, T(TOP,TOP+1), LDT, T(TOP+1,TOP+1), LDT,
     $         C, S )
            LGTHL = LGTHL + 1
            IDXL(LGTHL) = TOP
            CSL(LGTHL*2 - 1) = C
            CSL(LGTHL*2) = S
            TOP = TOP + 1
         ELSE
            TOP = TOP + 2
         END IF
         J = TOP
      ELSE
         J = J + 1
      END IF
C     
      
C     while ( def < PARA_INFNZ ) & ( j <= MAXJ )
      IF ( DEF.LT.PARA_INFNZ .AND. J.LE.MAXJ )  GO TO 10
      
      RETURN
C     *** Last line of UPWD ***
      END
C     
      SUBROUTINE DOWNWD( ISTOP, ISBOT, N, H, LDH, T, LDT, DEF, LGTHL,
     $     IDXL, CSL, LGTHR, IDXR, CSR )
      IMPLICIT NONE
C     
C     Chases all infinite eigenvalues to the bottom of the submatrix
C     pair (H,T) and stores transformations from the left in IDXL(1:LGTHL),
C     CSL(1:2*LGTHL), and transformations from the right in
C     IDXR(1:LGTHR), CSR(1:2*LGTHR).
C     
C     If ISTOP, this submatrix pair resides at the top left corner
C     of the original matrix pair.
C     If ISBOT, this submatrix pair resides at the bottom right corner
C     of the original matrix pair.
C     
C     .. Parameters ..
      DOUBLE PRECISION  ZERO
      PARAMETER         ( ZERO = 0.0D0 )
C     .. Configuration parameters ..

      LOGICAL           PARA_INFALL, PARA_INFBLOCK
      INTEGER           PARA_INFDEF, PARA_INFDIR, PARA_INFNZ,
     $   PARA_INFWDW
      DOUBLE PRECISION  PARA_INFTOL
      COMMON /PARA_INF/ PARA_INFALL, PARA_INFBLOCK, PARA_INFDEF, 
     $   PARA_INFDIR, PARA_INFNZ, PARA_INFWDW, PARA_INFTOL

C     .. Scalar Arguments ..
      LOGICAL           ISBOT, ISTOP
      INTEGER           LDH, LDT, LGTHL, LGTHR, N
C     .. Array Arguments ..
      DOUBLE PRECISION  CSL(*), CSR(*), H(LDH,*), T(LDT,*)
      INTEGER           IDXL(*), IDXR(*)
C     .. Local Scalars ..
      INTEGER           BOT, DEF, I, J, MINJ
      DOUBLE PRECISION  C, S, TEMP
C     .. External Subroutines ..
      EXTERNAL          DLARTG, DROT
C     .. External Functions ..
      LOGICAL           INFDEF
      EXTERNAL          INFDEF
C     .. Intrinsic Functions ..
      INTRINSIC         MAX
C     
C     .. Executable Statements ..
C     
      LGTHL = 0
      LGTHR = 0
      IF ( ISTOP ) THEN
         MINJ = 1
      ELSE
         MINJ = 2
      END IF
C     J = search position
      J = N
C     BOT = position to chase next eigenvalue to
      BOT = N
C     DEF = number of deflatable eigenvalues found
      DEF = 0
C     
 10   CONTINUE
C     do      
      IF ( T(J,J).EQ.ZERO ) THEN
C     
C     Chase zero to position BOT.
C     
         DEF = DEF + 1
C     
         DO 20 I = J, BOT-1
            TEMP = T(I,I+1 )
            CALL DLARTG( TEMP, T(I+1,I+1 ), C, S, T(I,I+1) )
            T(I+1,I+1 ) = ZERO
            IF ( I.LT.N-1 )
     $           CALL DROT( N-I-1, T(I,I+2), LDT, T(I+1,I+2), LDT,
     $           C, S )
            CALL DROT( N-MAX(I-1,1)+1, H(I,MAX(I-1,1)), LDH,
     $           H(I+1,MAX(I-1,1)), LDH, C, S )
            LGTHL = LGTHL + 1
            IDXL(LGTHL) = I
            CSL(LGTHL*2 - 1) = C
            CSL(LGTHL*2) = S
            IF ( I.GT.1 ) THEN
               TEMP = H(I+1,I)
               CALL DLARTG( TEMP, H(I+1,I-1), C, S, H(I+1,I) )
               H(I+1,I-1 ) = ZERO
               S = -S
               CALL DROT( I, H(1,I-1), 1, H(1,I), 1, C, S )
               CALL DROT( I-1, T(1,I-1), 1, T(1,I), 1, C, S )
               LGTHR = LGTHR + 1
               IDXR(LGTHR) = I-1
               CSR(LGTHR*2 - 1) = C
               CSR(LGTHR*2) = S
            END IF
 20      CONTINUE
C     
         IF ( ISBOT ) THEN
C     
C     Deflate infinite eigenvalue at bottom right corner.
C     
            TEMP = H(BOT,BOT)
            CALL DLARTG( TEMP, H(BOT,BOT-1), C, S, H(BOT,BOT) )
            H(BOT,BOT-1) = ZERO
            S = -S
            CALL DROT( BOT-1, H(1,BOT-1), 1, H(1,BOT), 1, C, S )
            CALL DROT( BOT-1, T(1,BOT-1), 1, T(1,BOT), 1, C, S )
            LGTHR = LGTHR + 1
            IDXR(LGTHR) = BOT-1
            CSR(LGTHR*2 - 1) = C
            CSR(LGTHR*2) = S
            BOT = BOT - 1
         ELSE
            BOT = BOT - 2
         END IF
         J = BOT
      ELSE
         J = J - 1
      END IF
C     while ( def < PARA_INFNZ ) & ( j >= MINJ )
      IF ( DEF.LT.PARA_INFNZ .AND. J.GE.MINJ )  GO TO 10
      
      RETURN
C     *** Last line of DOWNWD ***
      END

C
      SUBROUTINE ROTAPPL( LEFT, LGTH, M, IDX, CS, A, LDA )
C
C     If LEFT, this routine applies the transformations contained
C     in IDX, CS to M columns of A from the left.
C     If .NOT.LEFT, this routine applies the transformations contained
C     in IDX, CS to M rows of A from the right.
C
C     .. Scalar Arguments ..
      LOGICAL           LEFT
      INTEGER           LDA, LGTH, M
C     .. Array Arguments ..
      DOUBLE PRECISION  CS(*), A(LDA,*)
      INTEGER           IDX(*)
C     .. Local Scalars ..
      INTEGER           J
C     .. External Subroutines ..
      EXTERNAL          DLARTG, DROT
C     
C     .. Executable Statements ..
C
      IF ( LEFT ) THEN
         DO 10  J = 1, LGTH
            CALL DROT( M, A(IDX(J),1), LDA, A(IDX(J)+1,1), LDA,
     $         CS(J*2-1), CS(J*2) )
 10      CONTINUE
      ELSE
         DO 20  J = 1, LGTH
            CALL DROT( M, A(1,IDX(J)), 1, A(1,IDX(J)+1), 1, CS(J*2-1),
     $         CS(J*2) )
 20      CONTINUE
      END IF
C
      RETURN
C *** Last line of APPLROT ***
      END
C     
      LOGICAL FUNCTION INFDEF( N, POS, TNORM, T, LDT )
C     
C     Returns true if A(POS,POS) satisfies the deflation criterium for
C     an infinite eigenvalue.
C     Auxiliary function.
C     
C     .. Parameters ..
      DOUBLE PRECISION  ZERO
      PARAMETER         ( ZERO = 0.0D0 )
C     .. Configuration parameters ..
      LOGICAL           PARA_INFALL, PARA_INFBLOCK
      INTEGER           PARA_INFDEF, PARA_INFDIR, PARA_INFNZ,
     $   PARA_INFWDW
      DOUBLE PRECISION  PARA_INFTOL
      COMMON /PARA_INF/ PARA_INFALL, PARA_INFBLOCK, PARA_INFDEF, 
     $   PARA_INFDIR, PARA_INFNZ, PARA_INFWDW, PARA_INFTOL
C     .. Precision paramters ..
      DOUBLE PRECISION  SMLNUM
      COMMON /PREC/     SMLNUM
C     .. Scalar Arguments ..
      INTEGER           LDT, N, POS
      DOUBLE PRECISION  TNORM
C     .. Array Arguments ..
      DOUBLE PRECISION  T(LDT,*)
C     .. Local Scalars ..
      DOUBLE PRECISION  TST
C     .. External Functions ..
      DOUBLE PRECISION  DLANTR
      EXTERNAL          DLANTR
C     .. Intrinsic Functions ..
      INTRINSIC         ABS, MAX
C     
C     .. Executable Statements ..
C     
      IF ( T(POS,POS).EQ.ZERO ) THEN
         INFDEF = .TRUE.
      ELSE
         IF ( PARA_INFDEF.EQ.1 ) THEN
            INFDEF = ( ABS( T(POS,POS) ).LE.MAX( PARA_INFTOL*TNORM,
     $           SMLNUM ) )
         ELSE
            TST = ZERO
            IF ( POS.GT.1 )  TST = TST + ABS( T(POS-1,POS) )
            IF ( POS.LT.N )  TST = TST + ABS( T(POS,POS+1) )
            IF ( TST.EQ.ZERO ) THEN
               TST = DLANTR( '1', 'Upper', 'NoUnitDiag', N, N, T, LDT,
     $              DWORK )
            END IF
            INFDEF = ( ABS( T(POS,POS) ).LE.MAX( PARA_INFTOL*TST,
     $           SMLNUM ) )
         END IF
      END IF
C     
      RETURN
C     *** Last line of INFDEF ***
      END
