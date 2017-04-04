***********************************************************************
*                                                                     *
*     PDHGEQZ4.f:                                                     *
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
      SUBROUTINE PDHGEQZ4( ILSCHR, ILQ, ILZ, H, DESCH, T, DESCT,  Q,
     $   DESCQ, Z, DESCZ, N, IFIRST, ILAST, ILO, IHI, ALPHAI, ALPHAR,
     $   BETA, TMPA, TMPB, ILOQ, IHIQ, ILOZ, IHIZ, PU, PV, WR, WI, 
     $   WBETA, DWORK, LDWORK, IWORK, LIWORK, INFO )

      IMPLICIT NONE
*     ..
*     .. Scalar Arguments ..
*     ..

      INTEGER             IFIRST, ILAST, INFO, IHI, ILO, N, LIWORK, 
     $   LDWORK, ILOQ, IHIQ, ILOZ, IHIZ
      LOGICAL             ILQ, ILZ, ILSCHR
*     ..
*     .. Array Arguments ..
*     ..
      DOUBLE PRECISION    H( * ), T( * ), Q( * ), Z( * ), ALPHAI( * ),
     $   ALPHAR( * ), BETA( * ), WR( * ), WI( * ), WBETA( * ), 
     $   DWORK( * ), TMPA( * ), TMPB( * ), PU( * ), PV( * )
      INTEGER             DESCH( * ), DESCT( * ), DESCQ( * ), 
     $   DESCZ( * ), IWORK( * )

*     
*  Purpose
*  =======
*
*  PDHGEQZ4 is an auxiliary routine used to find the Schur decomposition
*  and or eigenvalues of a matrix pair already in Hessenberg, Triangular form 
*  from cols IFIRS to ILAST.  This routine requires that the 
*  active block is small enough, so that it can be solved 
*  serially.
*
*  All the inputs are assumed to be valid without checking(except workspace).
*
*  Notes
*  =====
*
*  Each global data object is described by an associated description
*  vector.  This vector stores the information required to establish
*  the mapping between an object element and its corresponding process
*  and memory location.
*
*  Let A be a generic term for any 2D block cyclicly distributed array.
*  Such a global array has an associated description vector DESCA.
*  In the following comments, the character _ should be read as
*  "of the global array".
*
*  NOTATION        STORED IN      EXPLANATION
*  --------------- -------------- --------------------------------------
*  DTYPE_A(global) DESCA( DTYPE_ )The descriptor type.  In this case,
*                                 DTYPE_A = 1.
*  CTXT_A (global) DESCA( CTXT_ ) The BLACS context handle, indicating
*                                 the BLACS process grid A is distribu-
*                                 ted over. The context itself is glo-
*                                 bal, but the handle (the integer
*                                 value) may vary.
*  M_A    (global) DESCA( M_ )    The number of rows in the global
*                                 array A.
*  N_A    (global) DESCA( N_ )    The number of columns in the global
*                                 array A.
*  MB_A   (global) DESCA( MB_ )   The blocking factor used to distribute
*                                 the rows of the array.
*  NB_A   (global) DESCA( NB_ )   The blocking factor used to distribute
*                                 the columns of the array.
*  RSRC_A (global) DESCA( RSRC_ ) The process row over which the first
*                                 row of the array A is distributed.
*  CSRC_A (global) DESCA( CSRC_ ) The process column over which the
*                                 first column of the array A is
*                                 distributed.
*  LLD_A  (local)  DESCA( LLD_ )  The leading dimension of the local
*                                 array.  LLD_A >= MAX(1,LOCr(M_A)).
*
*  Let K be the number of rows or columns of a distributed matrix,
*  and assume that its process grid has dimension p x q.
*  LOCr( K ) denotes the number of elements of K that a process
*  would receive if K were distributed over the p processes of its
*  process column.
*  Similarly, LOCc( K ) denotes the number of elements of K that a
*  process would receive if K were distributed over the q processes of
*  its process row.
*  The values of LOCr() and LOCc() may be determined via a call to the
*  ScaLAPACK tool function, NUMROC:
*          LOCr( M ) = NUMROC( M, MB_A, MYROW, RSRC_A, NPROW ),
*          LOCc( N ) = NUMROC( N, NB_A, MYCOL, CSRC_A, NPCOL ).
*  An upper bound for these quantities may be computed by:
*          LOCr( M ) <= ceil( ceil(M/MB_A)/NPROW )*MB_A
*          LOCc( N ) <= ceil( ceil(N/NB_A)/NPCOL )*NB_A
*
* 
*  Arguments 
*  ========= 
*
*  ILSCHR  (global input) LOGICAL
*          = TRUE: Apply all transformations to H and T. 
*          = FALSE: Only enough of H and T are updated to preserve the eigenvalues. 
*      
*  ILQ     (input) LOGICAL
*          = TRUE: Compute Q.
*          = FALSE: Q is not referenced.
*      
*  ILZ     (global input) LOGICAL
*          = TRUE: Compute Z.
*          = FALSE: Z is not referenced.
*
*  H       (local input/output) DOUBLE PRECISION array, dimension (LLD_H, LOCc(N)).
*          On entry, the upper Hessenberg matrix H.
*          On exit, if ILSCHR is .TRUE., H is upper quasi-triangular in
*          rows and columns IFIRST:ILAST, with any 2-by-2 or larger diagonal
*          blocks not yet in standard form. If ILSCHR is .FALSE., the
*          contents of H are unspecified on exit.
*
*  DESCH   (global and local input) INTEGER array of dimension DLEN_.
*          The array descriptor for the distributed matrix H. 
*      
*  T       (local input/output) DOUBLE PRECISION array, dimension (LLD_T, LOCc(N)).
*          On entry, the upper triangular matrix T.
*          On exit, if ILSCHR is .TRUE., T is maintained upper-triangular in
*          rows and columns IFIRST:ILAST.
*          ILSCHR is .FALSE., the
*          contents of T are unspecified on exit. 
*      
*  DESCT   (global and local input) INTEGER array of dimension DLEN_.
*          The array descriptor for the distributed matrix T. 
*
* 
*  Q       (local input/output) DOUBLE PRECISION array, dimension (LLD_Q, LOCc(N)).
*          If ILQ = FALSE:  Q is not referenced. 
*          If ILQ = TRUE:  on entry, Q must contain an orthogonal matrix.
*          On exit, Q is multiplied on the right with 
*          the product of givens transformations applied on H,T from the left.
*      
*  DESCQ   (global and local input) INTEGER array of dimension DLEN_.
*          The array descriptor for the distributed matrix Q. 
*      
*  Z       (local input/output) DOUBLE PRECISION array, dimension (LLD_Z, LOCc(N)).
*          If ILZ = FALSE:  Z is not referenced. 
*          If ILZ = TRUE:  on entry, Z must contain an orthogonal matrix.
*          On exit, Z is multiplied on the right with 
*          the product of givens transformations applied on H,T from the right. 
*      
*  DESCZ   (global and local input) INTEGER array of dimension DLEN_.
*          The array descriptor for the distributed matrix Z. 
*      
*  N       (global input) INTEGER 
*          The order of the matrices H,T,Q, and Z.  N >= 0. 
*
*  IFIRST  (global input) INTEGER
*  ILAST   (global input) INTEGER
*          It is assumed that A is already upper quasi-triangular in
*          rows and columns ILAST+1:N, and that A(IFIRST,IFIRST-1) = 0 (unless
*          IFIRST = 1). PDGHEQZ4 works primarily with the
*          submatrix pair in rows and columns IFIRST to ILAST, but applies
*          transformations to all of (H,T) if ILSCHR is .TRUE..
*          1 <= IFIRST <= max (1,ILAST); ILAST <= N.
*
*  ILO     (global input) INTERGER
*          Updates are not applied to rows(from left)/cols(from right) below this index. 
*
*  IHI     (global input) INTERGER
*          Updates are not applied to rows(from left)/cols(from right) greater than this index. 
*     
*
*  ALPHAR  (global output) DOUBLE PRECISION array, dimension N
*  ALPHAI  (global output) DOUBLE PRECISION array, dimension N
*  BETA    (global output) DOUBLE PRECISION array, dimension N
*          The real, imaginary and scale parts, respectively, of the computed
*          eigenvalues IFIRST to ILAST are stored in the corresponding
*          elements of ALPHAR, ALPHAI, and BETA.
*
*  TMPA    (local input) DOUBLE PRECISION, array, dimension see below.
*  TMPB    (local input) DOUBLE PRECISION, array, -||-
*          TMPA, TMPB are local workspace of minimum size
*          (ILAST - IFIRST + 1) * (ILAST - IFIRST + 1). 
*
*  ILOQ    (global input) INTEGER
*  IHIQ    (global input) INTEGER
*          Specify the rows of Q to which transformations must be
*          applied if ILQ is .TRUE.. 1 .LE. ILQQ .LE. IHIQ .LE. N.
*
*  ILOZ    (global input) INTEGER
*  IHIZ    (global input) INTEGER
*          Specify the rows of Z to which transformations must be
*          applied if ILZ is .TRUE.. 1 .LE. ILOZ .LE. IHIZ .LE. N.
*
*  PU      (local input) DOUBLE PRECISION, array, dimensions see below. 
*  PV      (local input) DOUBLE PRECISION, array, -||-
*          PU, PV are local workspace of minimum size
*          (ILAST - IFIRST + 1) * (ILAST - IFIRST + 1). 
*
*  WR      (local input) DOUBLE PRECISION array, dimension (ILAST-IFIRST+1).
*  WI      (local input) DOUBLE PRECISION array, dimension     -||-
*  WBETA   (local input) DOUBLE PRECISION array, dimension     -||-
*          WR, WI and WBETA are Local workspace of dimension ILAST-FIRST+1.
*          On output, the real, imaginary, and scale parts of the 
*          computed eigenvalues. 
*
*  DWORK   (local workspace/global output) DOUBLE PRECISION array,
*          dimension (LDWORK) 
*      
*  LDWORK  (global input) INTEGER 
*          The dimension of the array DWORK. 
*          If LDWORK = -1, then a workspace query is assmed, and optimal 
*          workspace is returned in DWORK(1)
*  
*  IWORK   (local workspace/global output) INTEGER array, 
*          dimension (LIWORK)  
*  
*  LIWORK  (global input) INTEGER  
*          The dimension of the array IWORK. 
*          If LIWORK = -1 then workspace query is assumed, and optimal  
*          workspace is returned in IWORK(1)
*     
*  INFO    (global output) INTEGER
*          = 0:  successful exit.
*          = -1:  Double workspace is not enough
*          = -2:  Integer workspace is not enough
*          = 1 : Error calling underlying KKQZ 
*          = 2 : Error updating of diagonal updates
*
*     .. Parameters ..
*     ..
      INTEGER            BLOCK_CYCLIC_2D, CSRC_, CTXT_, DLEN_, DT_,
     $     LLD_, MB_, M_, NB_, N_, RSRC_
      PARAMETER          ( BLOCK_CYCLIC_2D = 1, DLEN_ = 9, DT_ = 1,
     $     CTXT_ = 2, M_ = 3, N_ = 4, MB_ = 5, NB_ = 6,
     $     RSRC_ = 7, CSRC_ = 8, LLD_ = 9 )

      DOUBLE PRECISION   ZERO, ONE, MINUSONE
      PARAMETER          ( ZERO = 0.0D+0, ONE = 1.0D+0, 
     $                    MINUSONE = -1.0D+0)
      DOUBLE PRECISION    EXTFAC, EXTADD, NOSWP
      PARAMETER           ( EXTFAC = 1.5D0, EXTADD = -4.0D0,
     $     NOSWP = 1.5D-1 )
      INTEGER             NSBULG, NSMIN, NSMAX      
      PARAMETER           ( NSBULG = 2, NSMIN = 12, NSMAX = 20)
*     ..
*     .. Precision paramters ..
*     ..
      DOUBLE PRECISION  SMLNUM, ULP, TNORM, HNORM, HTOL, TTOL
      COMMON /PREC/     SMLNUM, ULP, TNORM, HNORM, HTOL, TTOL
      
*     ..
*     .. Local Scalars ..
*     ..
      INTEGER             DWORKMIN, IWORKMIN, NPROCS, IAM,
     $   LTOP, IERR, NPROW, NPCOL, MYROW, MYCOL, NW, I, J, ICTXT, RSRC,
     $   CSRC, LROW, LCOL
      DOUBLE PRECISION    TMP
      LOGICAL             LQUERY, TRANSPU
*     ..
*     .. Local Arrays
*     ..

      
*     ..
*     .. Externals ..
*     ..

      
      INFO = 0

      LQUERY = LDWORK .EQ. -1 .OR. LIWORK .EQ. -1


*     Estimate required workspace for KKQZ and following updates
      DWORKMIN = 6 * N +16
      IWORKMIN = 2 * N
   

      IF ( LQUERY ) THEN
         DWORK(1) = DWORKMIN
         IWORK(1) = IWORKMIN
         RETURN
      END IF
      IF ( LDWORK .LT. DWORKMIN ) THEN
         INFO = -1
         CALL PXERBLA( ICTXT, 'PDHGEQZ4', -INFO )
         RETURN
      END IF
      IF ( LIWORK .LT. IWORKMIN ) THEN
         INFO = -2
         CALL PXERBLA( ICTXT, 'PDHGEQZ4', -INFO )
         RETURN 
      END IF
 
*     The returned U(Q) from the kkqz call is transposed, so need to 
*     untranspose it.
*     Restored to orignal version, now works with old QZ!
      TRANSPU = .FALSE.

      IF( ILSCHR ) THEN
         LTOP = ILO 
      ELSE
         LTOP = IFIRST
      END IF

*     NW holds the current problem size
      NW = ILAST - IFIRST + 1


*     Extract current communication context and get grid parameters
      ICTXT = DESCH(CTXT_)
      CALL BLACS_PINFO(IAM, NPROCS)
      CALL BLACS_GRIDINFO( ICTXT, NPROW, NPCOL, MYROW, MYCOL )
      
*     Find out of who owns first element of the data
      CALL INFOG2L( IFIRST, IFIRST, DESCH, NPROW, NPCOL , MYROW, MYCOL,
     $   LROW, LCOL, RSRC, CSRC )

*     Copy the diagonal block of H and T to local stored at (RSRC,CSRC)
      CALL PDLACP3(NW, IFIRST, H, DESCH, TMPA, NW, RSRC, CSRC, 0)
      CALL PDLACP3(NW, IFIRST, T, DESCT, TMPB, NW, RSRC, CSRC, 0)
      IERR = 0
      IF ( MYROW .EQ. RSRC .AND. MYCOL .EQ. CSRC ) THEN
*        Call the serial QZ implementation
         CALL KKQZ( 'S', 'I', 'I', NW, 1, NW, TMPA, NW, TMPB, NW, WR,
     $      WI, WBETA, PU, NW, PV, NW, DWORK, LDWORK, IWORK, LIWORK,
     $      IERR )
*        Broadcast Send IERR to all so all can abort properly if something failed ( reordering for example ).
         CALL IGEBS2D(ICTXT, 'A', ' ', 1, 1, IERR, 1 )
         IF ( IERR .NE. 0 ) THEN
             INFO = -1
             GOTO 9999
         END IF
*        Transpose PU
         IF ( TRANSPU ) THEN
         DO J = 2, NW
             DO I = 1, J - 1
                TMP = PU( I + ( J - 1 ) * NW )
                PU( I + ( J - 1 ) * NW ) = PU( J + ( I - 1 ) * NW )
                PU( J + ( I - 1 ) * NW ) = TMP
             END DO
         END DO
         END IF
*        Broadcast send PU and PV to all. Not optimal but can be fixed later
         CALL DGEBS2D( ICTXT, 'A', ' ', NW, NW, PU, NW )
         CALL DGEBS2D( ICTXT, 'A', ' ', NW, NW, PV, NW )

      ELSE
*        Broadcast receive IERR from previous AED call.
         CALL IGEBR2D(ICTXT, 'A', ' ', 1, 1, IERR, 1, RSRC, CSRC ) 

         IF (IERR.NE.0) THEN
             INFO = -1
             GOTO 9999
         END IF

*        Broadcast receive PU and PV from 0,0.
         CALL DGEBR2D( ICTXT, 'A', ' ', NW, NW, PU, NW, RSRC, CSRC )
         CALL DGEBR2D( ICTXT, 'A', ' ', NW, NW, PV, NW, RSRC, CSRC )
      END IF


*     Restore updated parts of H and T.
      CALL PDLACP3( NW, IFIRST, H, DESCH, TMPA, NW, RSRC, CSRC, -1 )
      CALL PDLACP3( NW, IFIRST, T, DESCT, TMPB, NW, RSRC, CSRC, -1 )
      
*     Apply PU=Q and PV=Z to H, T, Q, and Z. 
      CALL PDHGEQZ6( 'A', ILSCHR, ILQ, ILZ, H, DESCH, T, DESCT, Q, DESCQ
     $   , Z, DESCZ, PU, PV, NW, LTOP, IFIRST, ILAST, N, NW, ILOQ, IHIQ, 
     $   ILOZ, IHIZ, LDWORK, DWORK, IERR )
      IF (IERR .NE. 0 ) THEN
          INFO = 2
          GOTO 9999
      END IF

*     Distributed eigenvalues
      IF ( MYROW .EQ. RSRC .AND. MYCOL .EQ. CSRC ) THEN
*        Broadcast send WR, WI, and WBETA to all.
         CALL DGEBS2D( ICTXT, 'A', ' ', NW, 1, WR, 1 )
         CALL DGEBS2D( ICTXT, 'A', ' ', NW, 1, WI, 1 )
         CALL DGEBS2D( ICTXT, 'A', ' ', NW, 1, WBETA, 1 )
      ELSE
*        Broadcast receive WR, WI, and WBETA from (RSRC,CSRC).
         CALL DGEBR2D( ICTXT, 'A', ' ', NW, 1, WR, 1, RSRC, CSRC )
         CALL DGEBR2D( ICTXT, 'A', ' ', NW, 1, WI, 1, RSRC, CSRC )
         CALL DGEBR2D( ICTXT, 'A', ' ', NW, 1, WBETA, 1, RSRC, CSRC )
      END IF
*     Extract Eigenvalues
      DO I = 1, NW
         ALPHAR( IFIRST + I - 1 ) = WR( I )
         ALPHAI( IFIRST + I - 1 ) = WI( I )
         BETA( IFIRST + I - 1 ) = WBETA( I )
      END DO
9999  CONTINUE
      


      RETURN
*     End of PDHGEQZ4

      END 
