      SUBROUTINE INVHSE( N, T, LDT, TAU, V, INFO )
      IMPLICIT NONE
C     
C     Computes a Householder transformation H such that the first
C     column of T*H is a multiple of the unit vector.
C     It is assumed that N <= 13.
C     Returns INFO = 1 if Gauss with partial pivoting caused problems.
C     
C     .. Parameters ..
      INTEGER           BULGMX
      PARAMETER         ( BULGMX = 12 )
      INTEGER           LDTMP
      PARAMETER         ( LDTMP = BULGMX+1 )
      DOUBLE PRECISION  ZERO, ONE, TWO
      PARAMETER         ( ZERO = 0.0D0, ONE = 1.0D0, TWO = 2.0D0 )

      LOGICAL           PARA_ITTRI
      INTEGER           PARA_ITASH, PARA_ITBULG, PARA_ITDACK, 
     $   PARA_ITDANB, PARA_ITEXC, PARA_ITEXT, PARA_ITKK, PARA_ITMAX, 
     $   PARA_ITNS, PARA_ITOPP
      COMMON /PARA_IT/  PARA_ITASH, PARA_ITBULG, PARA_ITDACK, 
     $   PARA_ITDANB, PARA_ITEXC, PARA_ITEXT, PARA_ITKK, PARA_ITMAX, 
     $   PARA_ITNS, PARA_ITOPP, PARA_ITTRI

C     .. Scalar Arguments ..
      INTEGER           INFO, LDT, N
      DOUBLE PRECISION  TAU
C     .. Array Arguments ..
      DOUBLE PRECISION  T(LDT,*), V(*)
C     .. Local Scalars ..
      INTEGER           I, IERR
C     .. Local Arrays ..
      DOUBLE PRECISION  NU(LDTMP), TMP(LDTMP,LDTMP)
      INTEGER           IPIV(LDTMP)
C     .. External Subroutines ..
      EXTERNAL          DGESV, DLACPY
C     .. Intrinsic Functions ..
      INTRINSIC         ABS, DBLE, INT, MAX, SQRT
C     
C     .. Executable Statements ..
C     
      INFO = 0
      IF ( N.EQ.0 )
     $   RETURN
C     
      IF ( PARA_ITOPP.EQ.1 ) THEN
C        
C        Gauss with partial pivoting is ok, but better use RQ.
C        
         CALL DLACPY( 'All', N, N, T, LDT, TMP, LDTMP )
         V(1) = ONE
         DO 10 I = 2, N
            V(I) = ZERO
 10      CONTINUE
         CALL DGESV( N, 1, TMP, LDTMP, IPIV, V, LDTMP, IERR )
         IF ( IERR.NE.0 ) THEN
            INFO = 1
            RETURN
         END IF
      ELSE
         CALL DLACPY( 'All', N, N, T, LDT, TMP, LDTMP )
         CALL DGERQF( N, N, TMP, LDTMP, NU, V, LDTMP, IERR )
         CALL DORGRQ( N, N, N, TMP, LDTMP, NU, V, LDTMP, IERR ) 
         DO 20 I = 1, N
            V(I) = TMP(1,I)
 20      CONTINUE
      END IF
C     
      CALL DLARFG( N, V(1), V(2), 1, TAU )
      V(1) = ONE

C     
C     Last line of INVHSE
C     
      END
