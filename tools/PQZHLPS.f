********************************************************************
*                                                                  *
* PQZHLPS - Misc. helper routines for the PDGEQZ library           *
*                                                                  *
********************************************************************

***********************************************************************
*                                                                     *
*     PCOMPRES.f:                                                     *
*     Auxilliary routine in the package PDGGHRD.                      *
*                                                                     *
*     Contributors: Björn Adlerborn                                   *
*                   Lars Karlsson                                     *
*                                                                     *
*     Department of Computing Science and HPC2N, Umea University      *
*                                                                     *
*                                                                     *
***********************************************************************
      SUBROUTINE PCOMPRES( A, B, Q, Z, S, T, DESCA, DESCB, DESCQ, DESCZ,
     $   DESCS, DESCT, N, RESA, RESB, RESQ, RESZ, ISQTRANS )
      
      IMPLICIT NONE

*
*     PURPOSE
*     =======
*
*     PCOMPRES calculates the relative residual of S and T, and level of 
*     orthogonality of Q and Z by computing :
*
*     1. ||( Q*A*Z - S )||_f / ||( A )||_f  => stored in RESA
*     2. ||( Q*B*Z - T )||_f / ||( B )||_f  => stored in RESB
*     3. ||( I - Q*Q')||_f / (eps*N) => stored in RESQ
*     4. ||( I - Z*Z')||_f / (eps*N) => stored in RESZ
* 
*     ARGUMENTS
*     =========
*
*     A       (local input/output) DOUBLE PRECISION array,
*     dimension (LLD_A, LOCc(N)).
*     On entry, the orignal matrix to be compared with S.
*     Unchanged on exit.
*
*     B       (local input/output) DOUBLE PRECISION array,
*     dimension (LLD_B, LOCc(N)).
*     On entry, the orignal matrix to be compared with T. 
*     Unchanged on exit.
*
*     Q       (local input/output) DOUBLE PRECISION array,
*     dimension (LLD_Q, LOCc(N)).
*     On entry, the matrix containing left sided orthogonal accumulated
*     transformations. 
*     Unchanged on exit.
*
*     Z       (local input/output) DOUBLE PRECISION array,
*     dimension (LLD_Q, LOCc(N)).
*     On entry, the matrix containing right sided orthogonal accumulated
*     transformations. 
*     Unchanged on exit.
*
*     S       (local input/output) DOUBLE PRECISION array,
*     dimension (LLD_S, LOCc(N)).
*     On entry, the updated version of A with corresponding updates applied
*     to Q and Z. 
*     Unchanged on exit.
*
*     T       (local input/output) DOUBLE PRECISION array,
*     dimension (LLD_S, LOCc(N)).
*     On entry, the updated version of A with corresponding updates applied
*     to Q and Z. 
*     Unchanged on exit
*
*
*     DESCA   (global and local input) INTEGER array of dimension DLEN_.
*     The array descriptor for the distributed matrix A.
*
*     DESCB   (global and local input) INTEGER array of dimension DLEN_.
*     The array descriptor for the distributed matrix B.
*
*     DESCQ   (global and local input) INTEGER array of dimension DLEN_.
*     The array descriptor for the distributed matrix Q.
*
*     DESCZ   (global and local input) INTEGER array of dimension DLEN_.
*     The array descriptor for the distributed matrix Z.
*
*     DESCS   (global and local input) INTEGER array of dimension DLEN_.
*     The array descriptor for the distributed matrix S.
*
*     DESCT   (global and local input) INTEGER array of dimension DLEN_.
*     The array descriptor for the distributed matrix T.
*
*     N       (global input) INTEGER
*     The order of the matrices A, B, S, T, Q, and Z.  N >= 0.
*
*     RESA    (global output) DOUBRE PRECISION
*     On exit, contains the caluclated value of 
*     ||( Q*A*Z - S )||_f / ||( A )||_f.
*
*     RESA    (global output) DOUBRE PRECISION
*     On exit, contains the caluclated value of 
*     ||( Q*B*Z - T )||_f / ||( B )||_f.
*
*     RESQ    (global output) DOUBRE PRECISION
*     On exit, contains the caluclated value of 
*     ||( I - Q*Q' )||_f / eps.
*
*     RESZ    (global output) DOUBLE PRECISION
*     On exit, contains the caluclated value of 
*     ||( I - Z*Z' )||_f / eps.
*
*     ISQTRANS (global input) LOGICAL
*     Indicates whether Q is calculated as transpose or not.
*
*     Constants.
*     
      INTEGER BLOCK_CYCLIC_2D, CSRC_, CTXT_, DLEN_, DTYPE_, LLD_, MB_,
     $   M_, NB_, N_, RSRC_        
      DOUBLE PRECISION ZERO, ONE, MINUSONE, EPS
      PARAMETER ( ZERO = 0.D0, ONE = 1.D0, MINUSONE = -1.D0,
     $   BLOCK_CYCLIC_2D = 1, DLEN_ = 9, DTYPE_ = 1, CTXT_ = 2, M_ = 3,
     $   N_ = 4, MB_ = 5, NB_ = 6, RSRC_ = 7, CSRC_ = 8, LLD_ = 9 )

*
*     Scalar arguments.
*
      INTEGER N
      DOUBLE PRECISION RESA, RESB, ARES, BRES, RESQ, RESZ
      LOGICAL ISQTRANS
*
*     Array arguments.
*
      INTEGER DESCA(9), DESCB(9), DESCZ(9), DESCQ(9), DESCT(9), DESCS(9)
      DOUBLE PRECISION A(*), B(*), Q(*), Z(*), S(*), T(*)

*
*     Local scalars.
*     
      INTEGER ACOLS, INFO, MYCOL, NPROW, NPCOL, MYROW, ICTXT
      DOUBLE PRECISION RES1, RES2, DPDUM, RES3, RES4

*
*     Local arrays.
*     
      DOUBLE PRECISION, ALLOCATABLE :: TMP1(:), TMP2(:)

*
*     Externals.
*     
      DOUBLE PRECISION PDLANGE, PDLAMCH
      EXTERNAL PDLANGE, PDLAMCH
      INTEGER NUMROC
      EXTERNAL NUMROC  

*
*     EXECUTABLE STATEMENTS
*     =====================
*

*
*     Extract process mesh information.
*     
      ICTXT = DESCA( CTXT_ )
      CALL BLACS_GRIDINFO( ICTXT, NPROW, NPCOL, MYROW, MYCOL )

*
*     Extract numerical characteristics about the machine.
*     
      EPS = PDLAMCH( ICTXT, 'PRECISION' )

*
*     Number of local columns in A.
*           
      ACOLS = NUMROC( DESCA( N_ ), DESCA( NB_ ), MYCOL, DESCA( CSRC_ ),
     $   NPCOL )

*
*     Allocate workspace TMP1 of size LLD_A x ACOLS (size of local A).
*     
      ALLOCATE( TMP1( DESCA( LLD_ ) * ACOLS ), STAT = INFO )
      IF( INFO .NE. 0 ) THEN
         WRITE(*,*) 'Could not allocate TMP1 for res. comp. ERR= ', INFO
         RETURN
      END IF

*
*     Allocate workspace TMP2 of size LLD_A x ACOLS (size of local A).
*     
      ALLOCATE( TMP2( DESCA( LLD_ ) * ACOLS ), STAT = INFO )
      IF( INFO .NE. 0 ) THEN
         WRITE(*,*) 'Could not allocate TMP2 for res. comp. ERR= ', INFO
         RETURN
      END IF

*
*     Compute ||A||_F.
*     
      ARES = PDLANGE( 'F', N, N, A, 1, 1, DESCA, DPDUM )

*
*     Compute ||B||_F.
*     
      BRES = PDLANGE( 'F', N, N, B, 1, 1, DESCB, DPDUM )

*
*     Copy S to TMP2.
*     
      CALL PDLACPY( 'All', N, N, S, 1, 1, DESCS, TMP2, 1, 1, DESCA )

*
*     Compute TMP1 := A * Z.
*     
      CALL PDGEMM( 'N', 'N', N, N, N, ONE, A, 1, 1, DESCA, Z, 1, 1,
     $   DESCZ, ZERO, TMP1, 1, 1, DESCA )

*
*     Compute TMP2 := TMP2 - Q^T * TMP1 = S - Q^T * A * Z.
*     
      IF (ISQTRANS) THEN
         CALL PDGEMM( 'T', 'N', N, N, N, ONE, Q, 1, 1, DESCQ, TMP1, 1, 1
     $      , DESCA, MINUSONE, TMP2, 1, 1, DESCA )
      ELSE
         CALL PDGEMM( 'N', 'N', N, N, N, ONE, Q, 1, 1, DESCQ, TMP1, 1, 1
     $      , DESCA, MINUSONE, TMP2, 1, 1, DESCA )
      END IF
           
*
*     Compute ||TMP2||_F = ||S - Q^T * A * Z||_F.
*     
      RES1 = PDLANGE( 'F', N, N, TMP2, 1, 1, DESCA, DPDUM )

*
*     Compute the relative residual for A.
*     
      RESA = RES1 / ARES

*
*     Copy T to TMP2.
*     
      CALL PDLACPY( 'All', N, N, T, 1, 1, DESCT, TMP2, 1, 1, DESCB )

*
*     Compute TMP1 := B * Z.
*     
      CALL PDGEMM( 'N', 'N', N, N, N, ONE, B, 1, 1, DESCB, Z, 1, 1,
     $   DESCZ, ZERO, TMP1, 1, 1, DESCB )

*
*     Compute TMP2 := TMP2 - Q^T * TMP1 = T - Q^T * B * Z.
*     
      IF (ISQTRANS) THEN
         CALL PDGEMM( 'T', 'N', N, N, N, ONE, Q, 1, 1, DESCQ, TMP1, 1, 1
     $      , DESCB, MINUSONE, TMP2, 1, 1, DESCB )
      ELSE
         CALL PDGEMM( 'N', 'N', N, N, N, ONE, Q, 1, 1, DESCQ, TMP1, 1, 1
     $      , DESCB, MINUSONE, TMP2, 1, 1, DESCB )
      END IF

*
*     Compute ||TMP2||_F = ||T - Q^T * B * Z||_F.
*     
      RES2 = PDLANGE( 'F', N, N, TMP2, 1, 1, DESCB, DPDUM )

*
*     Compute the relative residual for B.
*     
      RESB = RES2 / BRES  

*
*     Initialize TMP1 to the identity matrix.
*     
      CALL PDLASET( 'A', N, N, ZERO, ONE, TMP1, 1, 1, DESCQ )

*
*     Compute TMP1 := TMP1 - Q * Q^T = I - Q * Q^T.
*     
      CALL PDGEMM( 'N', 'T', N, N, N, ONE, Q, 1, 1, DESCQ, Q, 1, 1,
     $   DESCQ, MINUSONE, TMP1, 1, 1, DESCQ )

*
*     Compute ||TMP1||_F = ||I - Q * Q^T||_F.
*     
      RES3 = PDLANGE( 'F', N, N, TMP1, 1, 1, DESCQ, DPDUM )
      RESQ = RES3
      
*
*     Compute the level of orthogonality
*
      RESQ = RES3 / (eps * N)

*
*     Initialize TMP1 to the identity matrix.
*     
      CALL PDLASET( 'A', N, N, ZERO, ONE, TMP1, 1, 1, DESCZ )

*
*     Compute TMP1 = TMP1 - Z * Z^T = I - Z * Z^T.
*     
      CALL PDGEMM( 'N', 'T', N, N, N, ONE, Z, 1, 1, DESCZ, Z, 1, 1,
     $   DESCZ, MINUSONE, TMP1, 1, 1, DESCZ )

*
*     Compute ||TMP1||_F = ||I - Z * Z^T||_F.
*     
      RES4 = PDLANGE( 'F', N, N, TMP1, 1, 1, DESCZ, DPDUM )
      RESZ = RES4
      
*
*     Compute the level of orthogonality
*
      RESZ = RES4 / (eps * N)

*
*     Clean up.
*     
      DEALLOCATE( TMP1, TMP2 )
      
      RETURN
      
      END



      DOUBLE PRECISION FUNCTION CHISQ(FREEDOM)
      INTEGER FREEDOM
      DOUBLE PRECISION XX
      INTEGER J
      CHISQ = 0
      DO J = 1, FREEDOM
         CALL RANDOM_NUMBER(XX)
         CHISQ = CHISQ + XX*XX
      END DO
      END


      LOGICAL FUNCTION INT2LG(VAL)
      INTEGER VAL

      IF ( VAL .EQ. 0 ) THEN
        INT2LG = .FALSE.
      ELSE
        INT2LG = .TRUE.
      END IF
      END 
