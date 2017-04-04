***********************************************************************
*                                                                     *
*     PDHGEQZ2.f:                                                     *
*         Auxiliary routine in the package PDHGEQZ.                   *
*                                                                     *
*     Contributors: Bjorn Adlerborn                                   *
*                   Bo Kagstrom                                       *
*                   Daniel Kressner                                   *
*                                                                     *
*     Department of Computing Science and HPC2N, Umea University      *
*     MATHICSE ANCHP, EPF Lausanne                                    *
*                                                                     *
***********************************************************************
      RECURSIVE SUBROUTINE PDHGEQZ2( ILSCHR, ILQ, ILZ, ILO, N, KBOT, NW,
     $   H, DESCH, T, DESCT, Q, DESCQ, Z, DESCZ, ALPHAR, ALPHAI, BETA, 
     $   ILOQ, IHIQ, ILOZ, IHIZ, NS, ND, TMPA, TMPB, PU, PV, WR, WI, 
     $   WBETA, DWORK, LDWORK, INFO )
      IMPLICIT NONE

*     ..
*     .. Scalar Arguments ..
*     ..
      INTEGER             INFO, LDWORK, N, ILO, ILOQ, IHIQ, ILOZ, IHIZ, 
     $   NS, ND, NW, KBOT
      LOGICAL             ILQ, ILZ, ILSCHR
*     ..
*     .. Array Arguments .. 
*     ..
      DOUBLE PRECISION    H( * ), T( * ), ALPHAI( * ), ALPHAR( * ), 
     $   BETA( * ), Q( * ), DWORK( * ), Z( * ), TMPA( * ), TMPB( * ), 
     $   PU(*), PV(*), WR(*), WI(*), WBETA(*)
      INTEGER             DESCH(9), DESCT(9), DESCQ(9), DESCZ(9)


*   Purpose
*   =======
*     
*  Aggressive early deflation:
*     
*  PDHGEQZ2 accepts as input an upper Hessenberg matrix A and performs an
*  orthogonal equivalanze transformation designed to detect and deflate
*  fully converged eigenvalues from a trailing principal submatrix.  On
*  output A has been overwritten by a new Hessenberg matrix that is a
*  perturbation of an orthogonal similarity transformation of A.  It is
*  to be hoped that the final version of H has many zero subdiagonal
*  entries.
*     
*  This routine handles small deflation windows which is affordable by
*  one processor. 
*
*  All the inputs assumed to be valid without checking(except workspace).
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
*  ILO     (global input) INTERGER
*          Updates from right are not applied to rows below this index. 
*      
*  N       (global input) INTEGER 
*          The order of the matrices H,T,Q, and Z.  N >= 0. 
*      
*  KBOT    (global input) INTEGER 
*  NW      (global input) INTEGER 
*          It is assumed without a check that either
*          KBOT = N or H(KBOT+1,KBOT)=0.  KBOT and NW
*          determine an isolated block along the diagonal of the
*          Hessenberg-triangular matrix pair. However, H(KBOT-NW+1,KBOT-NW)=0 is not
*          essentially necessary if ILSCHR is .TRUE. .
*      
*  H       (local input/output) DOUBLE PRECISION array, dimension (LLD_H, LOCc(N)).
*          On input the initial N-by-N section of H stores the
*          Hessenberg matrix undergoing aggressive early deflation.
*          On output H has been transformed by an orthogonal
*          transformation, perturbed, and the returned
*          to Hessenberg form that (it is to be hoped) has some
*          zero subdiagonal entries.
*
*  DESCH   (global and local input) INTEGER array of dimension DLEN_.
*          The array descriptor for the distributed matrix H. 
*      
*  T       (local input/output) DOUBLE PRECISION array, dimension (LLD_T, LOCc(N)).
*          On input the initial N-by-N section of T stores the
*          triangular matrix undergoing aggressive early deflation.
*          On output T has been transformed by an orthogonal
*          transformation, perturbed, and the returned 
*          T is again in triangular.
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
*
*  ALPHAR  (global output) DOUBLE PRECISION array, dimension N
*  ALPHAI  (global output) DOUBLE PRECISION array, dimension N
*  BETA    (global output) DOUBLE PRECISION array, dimension N
*          The real, imaginary and scale parts of converged
*          eigenvalues are stored in ALPHAR(KBOT-ND+1) through ALHPAR(KBOT),
*          ALPHAI(KBOT-ND+1) through ALHPAI(KBOT), and 
*          BETA(KBOT-ND+1) through BETA(KBOT), respectively.
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
*  NS      (global output) INTEGER
*          The number of unconverged (i.e. approximate) eigenvalues
*          returned in ALPHAR, ALPHAI, and BETA that may be used as shifts by the
*          calling subroutine.
*
*  ND      (global output) INTEGER
*          The number of converged eigenvalues discovered by this
*          subroutine.
*
*  TMPA    (local input) DOUBLE PRECISION, array, dimension see below.
*  TMPB    (local input) DOUBLE PRECISION, array, -||-
*  PU      (local input) DOUBLE PRECISION, array, -||-
*  PV      (local input) DOUBLE PRECISION, array, -||-
*          TMPA, TMPB, PU, PV are local workspace of minimum size
*          NW * NW. 
*
*  WR      (local input) DOUBLE PRECISION array, dimension NW.
*  WI      (local input) DOUBLE PRECISION array, dimension NW.
*  WBETA   (local input) DOUBLE PRECISION array, dimension NW.
*          WR, WI and WBETA are Local workspace of dimension NW.
*          On output, the real, imaginary, and scale parts of the 
*          approximate eigenvalues that may be used for shifts are stored in 
*          WR(1) through WR(NS), 
*          WI(1) through WINS), and 
*          WBETA(1) through WBETA(NS), respectively.
*     
*  DWORK   (local workspace/global output) DOUBLE PRECISION array,
*          dimension (LDWORK)
*     
*  LDWORK  (global input) INTEGER
*          The dimension of the array DWORK.
*          If LDWORK = - 1 then a workspace query is assumed. Optimal workspace
*          is then returned in DWORK(1)
*
*  INFO     (global output) INTEGER
*          = 0 : Succesful call
*          = -1 : Insufficient workspace provided
*          = 1 : Error calling underlying QZEARLY.
*          = 2 : Error updating of diagonal updates
*
*  ================================================================
*     ..
*     .. Parameters ..
*     ..
      INTEGER            BLOCK_CYCLIC_2D, CSRC_, CTXT_, DLEN_, DT_,
     $     LLD_, MB_, M_, NB_, N_, RSRC_
      PARAMETER          ( BLOCK_CYCLIC_2D = 1, DLEN_ = 9, DT_ = 1,
     $     CTXT_ = 2, M_ = 3, N_ = 4, MB_ = 5, NB_ = 6,
     $     RSRC_ = 7, CSRC_ = 8, LLD_ = 9 )

      DOUBLE PRECISION    ZERO, ONE, MINUSONE
      PARAMETER           ( ZERO = 0.0D+0, ONE = 1.0D+0, 
     $     MINUSONE = -1.0D+0)
      
*     ..
*     .. Precision paramters ..
*     ..
      DOUBLE PRECISION  SMLNUM, ULP, TNORM, HNORM, HTOL, TTOL
      COMMON /PREC/     SMLNUM, ULP, TNORM, HNORM, HTOL, TTOL      
*     ..
*     .. Local Scalars ..
*     ..
      DOUBLE PRECISION    SCAL, TMP
      INTEGER             KWTOP, LTOP,  I, J, LDWKMIN, NPROW, NPCOL,
     $   MYROW, MYCOL, IAM, NPROCS, ICTXT, IERR, RSRC, CSRC, LROW, LCOL
      LOGICAL             LQUERY, TRANSPU
*     ..
*     .. Local Arrays ..
*     ..

*     ..
*     .. Externals 
*     ..

*     ..
*     .. Executable Statements ..
*     ..      
      

      INFO = 0
      IF( N .LE. 0 ) 
     $    RETURN
      LQUERY = LDWORK .EQ. -1

*     Extract current communication context and get grid parameters.
      ICTXT = DESCH(CTXT_)
      CALL BLACS_PINFO(IAM, NPROCS)
      CALL BLACS_GRIDINFO( ICTXT, NPROW, NPCOL, MYROW, MYCOL )

*     Estimate required workspace for QZEARLY and following updates.
      LDWKMIN = 6 * N +16
      IF ( LQUERY ) THEN
         DWORK(1) = LDWKMIN
         RETURN
      END IF
      IF ( LDWORK .LT. LDWKMIN ) THEN
         INFO = -1
         CALL PXERBLA( ICTXT, 'PDHGEQZ2', -INFO )
         RETURN
      END IF
      
*     The returned U(Q) from the aed call is transposed, so need to 
*     untranspose it.
*     Restored to original version - now works with old QZ!
      TRANSPU = .FALSE.

*     KWTOP points to the left and top most element of the AED window.
      KWTOP = KBOT - NW + 1

      IF( ILSCHR ) THEN
         LTOP = ILO 
      ELSE
         LTOP = KWTOP
      END IF
*     Find out of who owns first element of the data
      CALL INFOG2L( KWTOP, KWTOP, DESCH, NPROW, NPCOL , MYROW, MYCOL,
     $   LROW, LCOL, RSRC, CSRC )

*     Copy distributed H and T to local arrays, stored at (RSRC,CSRC).
      CALL PDLACP3( NW, KWTOP, H, DESCH, TMPA, NW, RSRC, CSRC, 0 )
      CALL PDLACP3( NW, KWTOP, T, DESCT, TMPB, NW, RSRC, CSRC, 0 )
      

*     Get the scaling factor unless the AED window resides at the top of the matrix pair.
      IF( KWTOP .EQ. ILO ) THEN
         SCAL = ZERO
      ELSE
         CALL PDELGET( 'A', ' ', SCAL, H, KWTOP, KWTOP - 1, DESCH )         
      END IF
      IERR = 0
      ND = 0
      IF ( MYROW .EQ. RSRC .AND. MYCOL .EQ. CSRC ) THEN
         CALL DLASET( 'All', NW, NW, ZERO, ONE, PU, NW )
         CALL DLASET( 'All', NW, NW, ZERO, ONE, PV, NW )
*        Set number of unconverged elements, initially to the AED window size.
         NS = NW
       
*        Procesors (0,0) performs AED and distributes results.
         CALL QZEARLY( NW, SCAL, HNORM, TMPA, NW, TMPB, NW, NS, WR, WI, 
     $      WBETA, PU, NW, PV, NW, DWORK, LDWORK, IERR )

*        Broadcast Send IERR to all so all can abort properly if something failed ( reordering for example ).
         CALL IGEBS2D(ICTXT, 'A', ' ', 1, 1, IERR, 1 )
         IF ( IERR .NE. 0 ) THEN
            INFO = 1
            GOTO 9999
         END IF
*        Transpose PU.
         IF ( TRANSPU ) THEN
            DO J = 2, NW
               DO I = 1, J - 1
                  TMP = PU( I + ( J - 1 ) * NW )
                  PU( I + ( J - 1 ) * NW ) = PU( J + ( I - 1 ) * NW )
                  PU( J + ( I - 1 ) * NW ) = TMP
               END DO
            END DO
         END IF
*        Send number of unconverged elements to all.
         CALL IGEBS2D(ICTXT, 'A', ' ', 1, 1,
     $      NS, 1 )
*        Broadcast send PU, and PV to all. Not optimal but can be fixed later.
         CALL DGEBS2D(ICTXT, 'A', ' ', NW, NW, PU, NW )
         CALL DGEBS2D(ICTXT, 'A', ' ', NW, NW, PV, NW )
      ELSE
*        Broadcast receive IERR from previous AED call.
         CALL IGEBR2D(ICTXT, 'A', ' ', 1, 1,
     $        IERR, 1, RSRC, CSRC )         
         IF ( IERR .NE. 0 ) THEN
            INFO = 1
            GOTO 9999
         END IF
*        Recive number of unconverged eigenvalues.
         CALL IGEBR2D(ICTXT, 'A', ' ', 1, 1,
     $      NS, 1, RSRC, CSRC )  
*        Broadcast receive PU and PV from 0,0.
         CALL DGEBR2D(ICTXT, 'A', ' ', NW, NW, PU, NW, RSRC, CSRC )
         CALL DGEBR2D(ICTXT, 'A', ' ', NW, NW, PV, NW, RSRC, CSRC )         
      END IF

*     If some of eigenvalues converged so ...
      IF ( NS .LT. NW ) THEN
*        Set number of converged elements .
         ND = NW - NS       

*        Restore Updated parts of distributed H and T.
         IF ( KWTOP .GT. ILO ) THEN
            CALL PDELSET( H, KWTOP, KWTOP - 1, DESCH, SCAL * PU( 1 ) )
         END IF     
         CALL PDLACP3( NW, KWTOP, H, DESCH, TMPA, NW, RSRC, CSRC, -1 )
         CALL PDLACP3( NW, KWTOP, T, DESCT, TMPB, NW, RSRC, CSRC, -1 )
         
*        Apply PU ( =Q) and PV ( =Z) to H, T, Q, Z
         INFO = 0
         CALL PDHGEQZ6( 'A', ILSCHR, ILQ, ILZ, H, DESCH, T, DESCT, Q, 
     $      DESCQ, Z, DESCZ, PU, PV, NW, LTOP, KWTOP, KBOT, N, NW,
     $      ILOQ, IHIQ, ILOZ, IHIZ, LDWORK, DWORK, IERR )
         IF (IERR .NE. 0 ) THEN
            INFO = 2
            GOTO 9999
         END IF         
      END IF
        
*     Distribute all computed eigenvalues ( to be used during the bulge chase part later on ).
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
*     Extract converged eigenvalues.
      DO I = 1, ND
         ALPHAR( KBOT - I+ 1) = WR( NW - I + 1)
         BETA( KBOT - I + 1) = WBETA( NW - I + 1)
         ALPHAI( KBOT - I + 1) = WI(NW - I + 1)
      END DO
      
9999  CONTINUE
      
      RETURN
*     
*     End of PDHGEQZ2.
*     
      END
