      SUBROUTINE QZFCOL( N, NSHIFT, SR, SI, SBETA, H, LDH, T, LDT, V,
     $                   INFO )
      IMPLICIT NONE
C
C     Forms a multiple of the first column of the shift polynomial for
C     the generalized eigenvalue problem.
C     To do: Avoidance of overflow while shifting. Documentation.
C     It is assumed that NSHIFT is even and >= 2, <= MIN(12,N-1).
C
C     .. Parameters ..
      INTEGER           BULGMX
      PARAMETER         ( BULGMX = 12 )
      INTEGER           LDTMP
      PARAMETER         ( LDTMP = BULGMX+1 )
      DOUBLE PRECISION  ZERO, ONE, TWO
      PARAMETER         ( ZERO = 0.0D0, ONE = 1.0D0, TWO = 2.0D0 )
C     .. Scalar Arguments ..
      INTEGER           INFO, LDH, LDT, N, NSHIFT
C     .. Array Arguments ..
      DOUBLE PRECISION  H(LDH,*),
     $                  SBETA(*), SI(*), SR(*), T(LDT,*), V(*)
C     .. Local Scalars ..
      INTEGER           I, J, NGIV
      DOUBLE PRECISION  ALPHA, CSE, SCAL, SNE
C     .. Local Arrays ..
      DOUBLE PRECISION  TMP(LDTMP,LDTMP), CS1(BULGMX+1),
     $                  CS2(BULGMX+1), SN1(BULGMX+1), SN2(BULGMX+1)
C     .. External Subroutines ..
      EXTERNAL          DAXPY, DLACPY, DLARTG, DROT
C     .. Intrinsic Functions ..
      INTRINSIC         ABS, DBLE, INT, MAX, SQRT
C
C     .. Executable Statements ..
C
      INFO = 0
      NGIV = 0
      DO 100  I = NSHIFT, 2, -2
C
         IF ( NGIV.GT.0 ) THEN
            CALL DLACPY( 'All', NGIV+1, NGIV+1, T, LDT, TMP, LDTMP )
            DO 10  J = NGIV, 1, -1
               CALL DROT( NGIV-J+2, TMP(J,J), LDTMP, TMP(J+1,J), LDTMP,
     $                    CS1(J), SN1(J) )
               ALPHA = TMP(J+1,J+1)
               CALL DLARTG( ALPHA, TMP(J+1,J), CS2(J), SN2(J),
     $                      TMP(J+1,J+1) )
               SN2(J) = -SN2(J)
               CALL DROT( J, TMP(1,J), 1, TMP(1,J+1), 1, CS2(J),
     $                    SN2(J) )
   10       CONTINUE
            SCAL = TMP(1,1)
         ELSE 
            SCAL = T(1,1)
         END IF
C
         CALL DLACPY( 'All', NGIV+2, NGIV+2, H, LDH, TMP, LDTMP )
         DO 20  J = 1, NGIV+2
            CALL DSCAL( MIN( J+1, NGIV+2 ), SBETA(I), TMP(1,J), 1 )
            CALL DAXPY( J, -SR(I), T(1,J), 1, TMP(1,J), 1 )
   20    CONTINUE
         DO 30  J = NGIV, 1, -1
            CALL DROT( NGIV+2, TMP(1,J), 1, TMP(1,J+1), 1, CS2(J),
     $                 SN2(J) )
   30    CONTINUE
         DO 40  J = NGIV+1, 1, -1
            ALPHA = TMP(J,1)
            CALL DLARTG( ALPHA, TMP(J+1,1), CS2(J), SN2(J), TMP(J,1) )
   40    CONTINUE
         NGIV = NGIV + 1
C
         CALL DLARTG( TMP(1,1), SI(I)*SCAL, CSE, SNE, ALPHA )
         CALL DLACPY( 'All', NGIV+1, NGIV+1, T, LDT, TMP, LDTMP )
         DO 50  J = NGIV, 1, -1
            CALL DROT( NGIV-J+2, TMP(J,J), LDTMP, TMP(J+1,J), LDTMP,
     $                 CS2(J), SN2(J) )
            ALPHA = TMP(J+1,J+1)
            CALL DLARTG( ALPHA, TMP(J+1,J), CS2(J), SN2(J),
     $                   TMP(J+1,J+1) )
            SN2(J) = -SN2(J)
            CALL DROT( J, TMP(1,J), 1, TMP(1,J+1), 1, CS2(J), SN2(J) )
   50    CONTINUE
         CALL DLARTG( CSE, -SNE*TMP(1,1), CSE, SNE, ALPHA )
         SNE = -SNE
C
         CALL DLACPY( 'All', NGIV+2, NGIV+2, H, LDH, TMP, LDTMP )
         DO 60  J = 1, NGIV+2
            CALL DSCAL( MIN( J+1, NGIV+2 ), SBETA(I-1), TMP(1,J), 1 )
            CALL DAXPY( J, -SR(I-1), T(1,J), 1, TMP(1,J), 1 )
   60    CONTINUE
         DO 70  J = NGIV, 1, -1
            CALL DROT( NGIV+2, TMP(1,J), 1, TMP(1,J+1), 1, CS2(J),
     $                 SN2(J) )
   70    CONTINUE
         CALL DLASET( 'All', NGIV+2, NGIV+1, ZERO, -SI(I-1), TMP(1,2),
     $                LDTMP )
         DO 80  J = NGIV-1, 1, -1
            CALL DROT( NGIV, TMP(1,J+1), 1, TMP(1,J+2), 1, CS1(J),
     $                 SN1(J) )
   80    CONTINUE
         CALL DROT( NGIV+2, TMP(1,1), 1, TMP(1,2), 1, CSE, SNE )
         IF ( I.GT.2 ) THEN
            DO 90  J = NGIV+1, 1, -1
               ALPHA = TMP(J,1)
               CALL DLARTG( ALPHA, TMP(J+1,1), CS1(J), SN1(J),
     $                      TMP(J,1) )
   90       CONTINUE
            NGIV = NGIV + 1
         ELSE
            CALL DCOPY( NGIV+2, TMP(1,1), 1, V, 1 )
         END IF
  100 CONTINUE

C
C     Last line of QZFCOL
C
      END
