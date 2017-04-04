***********************************************************************
*                                                                     *
*     PDHGEQZ3.f:                                                     *
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

      RECURSIVE SUBROUTINE PDHGEQZ3( WANTHT, ILQ, ILZ, N, H, DESCH, T, 
     $   DESCT, Q, DESCQ, Z, DESCZ, ILOQ, IHIQ, ILOZ, IHIZ, tH, tDESCH, 
     $   tT, tDESCT, tQ, tDESCQ, tZ, tDESCZ, KTOP, KBOT, SR, SI, SBETA, 
     $   NW, ND, KWTOP, DWORK, LDWORK, IWORK, LIWORK, INFO, UPDTIME, 
     $   RLVL )
      IMPLICIT NONE
*     ..
*     .. Scalar Arguments ..
*     ..
      INTEGER             INFO, KBOT, LDH, LDQ, LDT, LDWORK, LDZ, N,
     $   ICTXT, NPROW, NPCOL, MYROW, MYCOL,NW,
     $   KWTOP, IROFFH, NS, KTOP, ND, LIWORK,
     $   RLVL, ILOQ, IHIQ, ILOZ, IHIZ
      LOGICAL             WANTHT, ILQ, ILZ
*     ..
*     .. Array Arguments ..
*     ..
      INTEGER             DESCQ( 9 ), DESCZ( 9 ), DESCH( 9 ), 
     $   DESCT( 9 ), DESCV( 9 ), DESCTAU( 9 ), IWORK( * ),
     $   tDESCQ(9), tDESCZ(9), tDESCH(9), tDESCT(9)
      DOUBLE PRECISION    DWORK( * ), Q( * ),H( * ),
     $   T( * ), Z( *), tQ( * ), tH( * ),
     $   tT( * ), tZ( * ), SBETA(*), SI(*), SR(*)
*     Purpose
*     =======
*     
*     Aggressive early deflation:
*     
*     This subroutine accepts as input an upper Hessenberg matrix H, an upper
*     triangular matrix T, and
*     performs orthogonal equivalence transformations designed to detect
*     and deflate fully converged eigenvalues from a trailing principal
*     submatrix of the pair (H,T).  On output H has been overwritten by a new 
*     Hessenberg matrix that is a perturbation of orthogonal transformation of H.
*     
*     T is also perturbed, but maintained upper triangular.
*     It is to be hoped that the final version of H
*     has many zero subdiagonal entries.     
*     
*     All inputs are assumed valid without checking (except workspace).
*     
*     Notes
*     =====
*     
*     Each global data object is described by an associated description
*     vector.  This vector stores the information required to establish
*     the mapping between an object element and its corresponding process
*     and memory location.
*     
*     Let A be a generic term for any 2D block cyclicly distributed array.
*     Such a global array has an associated description vector DESCA.
*     In the following comments, the character _ should be read as
*     "of the global array".
*     
*     NOTATION        STORED IN      EXPLANATION
*     --------------- -------------- --------------------------------------
*     DTYPE_A(global) DESCA( DTYPE_ )The descriptor type.  In this case,
*     DTYPE_A = 1.
*     CTXT_A (global) DESCA( CTXT_ ) The BLACS context handle, indicating
*     the BLACS process grid A is distribu-
*     ted over. The context itself is glo-
*     bal, but the handle (the integer
*     value) may vary.
*     M_A    (global) DESCA( M_ )    The number of rows in the global
*     array A.
*     N_A    (global) DESCA( N_ )    The number of columns in the global
*     array A.
*     MB_A   (global) DESCA( MB_ )   The blocking factor used to distribute
*     the rows of the array.
*     NB_A   (global) DESCA( NB_ )   The blocking factor used to distribute
*     the columns of the array.
*     RSRC_A (global) DESCA( RSRC_ ) The process row over which the first
*     row of the array A is distributed.
*     CSRC_A (global) DESCA( CSRC_ ) The process column over which the
*     first column of the array A is
*     distributed.
*     LLD_A  (local)  DESCA( LLD_ )  The leading dimension of the local
*     array.  LLD_A >= MAX(1,LOCr(M_A)).
*     
*     Let K be the number of rows or columns of a distributed matrix,
*     and assume that its process grid has dimension p x q.
*     LOCr( K ) denotes the number of elements of K that a process
*     would receive if K were distributed over the p processes of its
*     process column.
*     Similarly, LOCc( K ) denotes the number of elements of K that a
*     process would receive if K were distributed over the q processes of
*     its process row.
*     The values of LOCr() and LOCc() may be determined via a call to the
*     ScaLAPACK tool function, NUMROC:
*     LOCr( M ) = NUMROC( M, MB_A, MYROW, RSRC_A, NPROW ),
*     LOCc( N ) = NUMROC( N, NB_A, MYCOL, CSRC_A, NPCOL ).
*     An upper bound for these quantities may be computed by:
*     LOCr( M ) <= ceil( ceil(M/MB_A)/NPROW )*MB_A
*     LOCc( N ) <= ceil( ceil(N/NB_A)/NPCOL )*NB_A
*     
*     Arguments
*     =========
*     
*     WANTHT  (global input) LOGICAL
*     If .TRUE., then the matrix pair H,T are fully updated
*     so that the quasi-triangular Schur factor may be
*     computed (in cooperation with the calling subroutine).
*     If .FALSE., then only enough of H,T is updated to preserve
*     
*     ILQ      (input) LOGICAL
*     = TRUE: Compute Q.
*     = FALSE: Q is not referenced.
*     
*     ILZ     (global input) LOGICAL
*     = TRUE: Compute Z.
*     = FALSE: Z is not referenced.
*     
*     N       (input) INTEGER
*     The order of the matrices H, T, Q, and Z.  N >= 0.
*     
*     H       (local input/output) DOUBLE PRECISION array, dimension (DESCH(LLD_),*)
*     On input the initial N-by-N section of H stores the
*     Hessenberg matrix undergoing aggressive early deflation.
*     On output H has been transformed by an orthogonal
*     equivalence transformations, perturbed, and the returned
*     to Hessenberg form that (it is to be hoped) has some
*     zero subdiagonal entries.
*     
*     DESCH   (global and local input) INTEGER array of dimension DLEN_.
*     The array descriptor for the distributed matrix H.
*     
*     T       (local input/output) DOUBLE PRECISION array, dimension (tDESCT(LLD_),*) 
*     On input the initial N-by-N section of T stores the 
*     triangular matrix undergoing aggressive early deflation. 
*     On output T has been transformed by an orthogonal 
*     equivalence transformations, perturbed, and then returned 
*     to triangular form. 
*     
*     DESCT   (global and local input) INTEGER array of dimension DLEN_.
*     The array descriptor for the distributed matrix T. 
*     
*     Q       (input/output) DOUBLE PRECISION array, dimension (DESCQ(LLD_),*)
*     On entry, the leading N-by-N part of this array must
*     contain an orthogonal matrix Q.
*     On exit, the leading N-by-N part of this array contains
*     the matrix Q post-multiplied by transpose of the
*     transformations which are applied to H and T on the left.
*     
*     DESCQ   (global and local input) INTEGER array of dimension DLEN_.
*     The array descriptor for the distributed matrix Q.
*     
*     Z       (input/output) DOUBLE PRECISION array, dimension (DESCZ(LLD_),*)
*     On entry, the leading N-by-N part of this array must
*     contain an orthogonal matrix Z.
*     On exit, the leading N-by-N part of this array contains
*     the matrix Z post-multiplied by the
*     transformations which are applied to H and T on the right.
*     
*     DESCZ   (global and local input) INTEGER array of dimension DLEN_.
*     The array descriptor for the distributed matrix Z.
*     
*     ILOQ    (global input) INTEGER
*     IHIQ    (global input) INTEGER
*     Specify the rows of Q to which transformations must be
*     applied if ILQ is .TRUE.. 1 .LE. ILQQ .LE. IHIQ .LE. N.
*     
*     ILOZ    (global input) INTEGER
*     IHIZ    (global input) INTEGER
*     Specify the rows of Z to which transformations must be
*     applied if ILZ is .TRUE.. 1 .LE. ILOZ .LE. IHIZ .LE. N.
*     
*     tH      (global workspace) DOUBLE PRECISION array, 
*     dimension (tDESCH(LLD_),*)
*     
*     tDESCH  (global and local input) INTEGER array of dimension DLEN_.
*     The array descriptor for the distributed matrix tH.
*     
*     tT      (global workspace) DOUBLE PRECISION array, 
*     dimension (tDESCT(LLD_),*) 
*     
*     tDESCT  (global and local input) INTEGER array of dimension DLEN_.
*     The array descriptor for the distributed matrix tT. 
*     
*     tQ      (global workspace) DOUBLE PRECISION array, 
*     dimension (tDESCQ(LLD_),*)
*     
*     tDESCQ  (global and local input) INTEGER array of dimension DLEN_.
*     The array descriptor for the distributed matrix tQ.
*     
*     tZ      (global workspace) DOUBLE PRECISION array, 
*     dimension (tDESCZ(LLD_),*)
*     
*     tDESCZ  (global and local input) INTEGER array of dimension DLEN_.
*     The array descriptor for the distributed matrix tZ.
*     
*     SR      (global output) DOUBLE PRECISION array, dimension N
*     SI      (global output) DOUBLE PRECISION array, dimension N
*     BETA    (global output) DOUBLE PRECISION array, dimension N
*     On output, the real, imaginary and scale parts of approximate
*     eigenvalues that may be used for shifts are stored in
*     SR(KBOT-ND-NS+1) through SR(KBOT-ND),
*     SI(KBOT-ND-NS+1) through SI(KBOT-ND), and 
*     BETA(KBOT-ND+NS+1) through BETA(KBOT-ND) respectively.
*     The real, imaginary and scale parts of converged eigenvalues
*     are stored in SR(KBOT-ND+1) through SR(KBOT),
*     SI(KBOT-ND+1) through SI(KBOT), and BETA(KBOT-ND+1) 
*     through BETA(KBOT) respectively.
*     
*     KTOP    (global input) INTEGER
*     It is assumed that either KTOP = 1 or H(KTOP,KTOP-1)=0.
*     KBOT and KTOP together determine an isolated block
*     along the diagonal of the Hessenberg matrix.
*     
*     KBOT    (global input) INTEGER
*     It is assumed without a check that either
*     KBOT = N or H(KBOT+1,KBOT)=0.  KBOT and KTOP together
*     determine an isolated block along the diagonal of the
*     Hessenberg matrix.
*     
*     NW      (global input) INTEGER
*     Deflation window size.  1 .LE. NW .LE. (KBOT-KTOP+1).
*     
*     ND      (global output) INTEGER
*     The number of converged eigenvalues uncovered by this
*     subroutine.
*     KTOPW   (global input) INTEGER
*     Redundant, same as KBOT-NW+1
*     
*     DWORK   (local workspace/global output) DOUBLE PRECISION array,
*     dimension (LWORK) 
*     
*     LDWORK  (global input) INTEGER 
*     The dimension of the array WORK. 
*     If LWORK = -1, then a workspace query is assumed, and optimal 
*     workspace is returned in WORK(1)
*     
*     IWORK   (local workspace/global output) INTEGER array, 
*     dimension (LIWORK)  
*     
*     LIWORK  (global input) INTEGER  
*     The dimension of the array IWORK. 
*     If LIWORK = -1 then workspace query is assumed, and optimal  
*     workspace is returned in IWORK(1)
*     
*     INFO    (global output) INTEGER
*     = 0:  successful exit.
*     = -1:  Double workspace is not enough
*     = -2:  Integer workspace is not enough
*     > 0: Some problem occured. 1, 2, and 3, means
*     some sub routine failed.
*     
*     RVLV    (global input) INTEGER
*     Indicates what level of recursion PDHGEQZ3 currently executes in.
*     0 or 1 are accepted values. Recurssion is stopped if RVLV >= 1.
*     
*     
*     
*     ================================================================

*     ..     
*     .. Parameters ..
*     ..
      DOUBLE PRECISION    ZERO, ONE, TWO
      PARAMETER           ( ZERO = 0.0D0, ONE = 1.0D0, TWO = 2.0D0)
      INTEGER             BLOCK_CYCLIC_2D, CSRC_, CTXT_, DLEN_, DT_, 
     $   LLD_, MB_, M_, NB_, N_, RSRC_
      PARAMETER          ( BLOCK_CYCLIC_2D = 1, DLEN_ = 9, DT_ = 1,
     $   CTXT_ = 2, M_ = 3, N_ = 4, MB_ = 5, NB_ = 6,
     $   RSRC_ = 7, CSRC_ = 8, LLD_ = 9 )
      DOUBLE PRECISION    EXTFAC, EXTADD, NOSWP
      PARAMETER           ( EXTFAC = 1.5D0, EXTADD = -4.0D0,
     $   NOSWP = 1.5D-1 )
      INTEGER             NSBULG, NSMIN, NSMAX      
      PARAMETER           ( NSBULG = 2, NSMIN = 12, NSMAX = 20)      
*     ..                        
*     .. Local Scalars ..
*     ..
      INTEGER             HTZQROWS0, HTZQCOLS0, J, I, IFST, ILST, IR, 
     $   NB, IPWORK, DUMMYM, NR, NC, IPTAU, NRTAU, NCTAU, JW, NPMIN, 
     $   ICTXT_NEW, MYROW_NEW, MYCOL_NEW, LC, LR, LIWORK1,IPTH, IPTT, 
     $   IPTQ, IPTZ, RSRC, CSRC, SMALLNB, NB_NEW, OFF_C, OFF_R, LTOP, 
     $   KLN, NPROCS, IAM, NMIN1, NMIN3, IDUM, IERR, IPT, LILST0, LILST,
     $   ltQROWS, ltQCOLS, IPTQ2
      DOUBLE PRECISION    ALPHA, SCAL, ELEM, TMPQ1, TMPQ2, T1, T2, TST1, 
     $   RTDET, LDWORK1

      LOGICAL             BULGE, DFLATE, DEBUG, RUNDHGEQZ1
      CHARACTER           CILQ, CILZ, CILSCHR
*     ..
*     .. Local Arrays ..
*     ..
      INTEGER             DESC_TMP( 9 ), PARA( 6 ),  PMAP( 64 * 64 ),
     $   DESCHTZQ( 9 )
      DOUBLE PRECISION    TMPH( 2 , 2 ), UPDTIME( 40 )
      DOUBLE PRECISION, ALLOCATABLE :: tQ2( : )
      
*     ..
*     .. Precision paramters ..
*     ..
      DOUBLE PRECISION    SMLNUM, ULP, TNORM, HNORM, HTOL, TTOL
      COMMON /PREC/       SMLNUM, ULP, TNORM, HNORM, HTOL, TTOL
      
*     ..      
*     .. Externals ..
*     ..
      EXTERNAL            DLAMCH, INDXG2L, INDXG2P,NUMROC, ICEIL, 
     $   MPI_WTIME, BLACS_PNUM, PILAENVX
      DOUBLE PRECISION    DLAMCH, MPI_WTIME      
      INTEGER             INDXG2L, INDXG2P,NUMROC, ICEIL, BLACS_PNUM, 
     $   PILAENVX
      
*     .. 
*     .. Intrinsic Functions ..
*     ..
      INTRINSIC           ABS, DBLE, INT, MAX, SQRT

      INFO = 0
      CILQ = 'N'
      CILZ = 'N'
      CILSCHR = 'E'
      IF (WANTHT) CILSCHR = 'S'
      IF (ILQ) CILQ = 'V'
      IF (ILZ) CILZ = 'V'

*     
*     Set number of deflated to none initially
*     
      ND = 0
      IF ( N .LE. 0 ) THEN
         RETURN
      END IF 

*     
*     Extract current communication context and get grid parameters
*     
      ICTXT = DESCH(CTXT_)
      CALL BLACS_PINFO(IAM, NPROCS)
      CALL BLACS_GRIDINFO( ICTXT, NPROW, NPCOL, MYROW, MYCOL )

*     Set DEBUG to true/false for various info during the run
      DEBUG = .FALSE.
*     Print debuginfo only for root.
      DEBUG = DEBUG .AND. IAM .EQ. 0
      

*     
*     NB_ .EQ. MB_ is assumed, we use NB as blocking factor.
*     
      NB = DESCH(NB_)




*     
*     NMIN1 defines the minimum problem size to be solved by PDHGEQZ0. 
*     If problem size is less than this value, PDHGEQZ1 is called instead
*     
      NMIN1 = PILAENVX(ICTXT, 50, 'PDHGEQZ', '', IDUM, IDUM, IDUM, IDUM)
      IF ( DEBUG ) write(*,*)'% PQZ3 : NMIN1 = ', NMIN1

*     
*     NMIN3 is Used for defining the maximum problem size to be solved serially
*     
      NMIN3 = PILAENVX(ICTXT, 55, 'PDHGEQZ', '', IDUM, IDUM, IDUM, IDUM)
      IF ( DEBUG ) write(*,*)'% PQZ3 : NMIN3 = ', NMIN3

      LDH = DESCH(LLD_)
      LDT = DESCT(LLD_)
      LDQ = DESCQ(LLD_)
      LDZ = DESCZ(LLD_)

*     
*     Global number of rows and columns
*     
      NR = DESCH( M_ )
      NC = DESCH( N_ )
      
      IROFFH = MOD( KWTOP - 1, NB )
      JW = MIN ( NW, KBOT - KTOP + 1 )
      NS = JW

*     
*     Determine which version of // QZ to use
*     For now, we only use PDHGEQZ1
      RUNDHGEQZ1 = .FALSE.
      RUNDHGEQZ1 = ( JW + IROFFH ) .LE. NMIN1
      RUNDHGEQZ1 = .TRUE.       


*     
*     Extract parameters for the // reordering
*     
*     Maximum number of independent computational windows
      PARA( 1 ) = PILAENVX(ICTXT, 80, 'PDHGEQZ', '', JW, NB, 
     $   IDUM, IDUM)
*     Number of eigenvalues in each window
      PARA( 2 ) = PILAENVX(ICTXT, 81, 'PDHGEQZ', '', IDUM, NB, 
     $   IDUM, IDUM)
*     Computational window size
      PARA( 3 ) = PILAENVX(ICTXT, 82, 'PDHGEQZ', '', IDUM, NB, 
     $   IDUM, IDUM)
*     Minimal percentage of flops required for
*     performing matrix-matrix multiplications instead of
*     pipelined orthogonal transformations
      PARA( 4 ) = PILAENVX(ICTXT, 83, 'PDHGEQZ', '', IDUM, IDUM, 
     $   IDUM, IDUM)
*     Width of block column slabs for row-wise
*     application of pipelined orthogonal transformations in
*     their factorized form
      PARA( 5 ) = PILAENVX(ICTXT, 84, 'PDHGEQZ', '', IDUM, NB, 
     $   IDUM, IDUM)
*     Maximum number of eigenvalues to bring over
*     the block border 
      PARA( 6 ) = PILAENVX(ICTXT, 85, 'PDHGEQZ', '', IDUM, NB, 
     $   IDUM, IDUM)

      IF (DEBUG) write(*,*)'% PQZ3 : Reorder PARA=', PARA(1:6)

*     Set parameters for the parallel hessenberg-triangular reduction
*     


*     
*     Check if workspace query
*     If so, check with the major called sub routines for workspace needed
*     
      IF ( LDWORK .EQ. -1 .OR. LIWORK .EQ. -1 ) THEN      
         LIWORK1 = 0
         LDWORK1 = 0
         IF ( DEBUG ) WRITE(*,*)'% PQZ3 -> PQZ1 WSPACEQ'
         CALL PDHGEQZ1( WANTHT, ILQ, ILZ,
     $      JW + IROFFH, 1 + IROFFH,
     $      JW + IROFFH, H, DESCH, T, DESCT,
     $      SR, SI, SBETA,
     $      Q, DESCQ, Z, DESCZ,
     $      DWORK, -1, IWORK, -1, IERR )
         IF ( DEBUG ) WRITE(*,*)'% PQZ3 -> PQZ1 WSPACEQ Done', 
     $      INT( DWORK( 1 )), INT( IWORK( 1 ) ) 
         LIWORK1 = IWORK( 1 )
         LDWORK1 = DWORK( 1 )

         IF ( DEBUG ) WRITE(*,*)'% PQZ3 -> PHRD WSPACEQ'
         CALL PDGGHRD( CILQ, CILZ, JW, 1, JW, tH, DESCH, tT,
     $      DESCT, tQ, DESCQ, tZ, DESCZ, DWORK, -1, IERR )
         IF ( DEBUG ) WRITE(*,*)'% PQZ3 -> PHRD WSPACEQ Done', 
     $      INT( DWORK( 1 ) )
         LDWORK1 = MAX( LDWORK1, DWORK( 1 ) )

         IF ( DEBUG ) WRITE(*,*)'% PQZ3 -> REORDER WSPACEQ'
         DUMMYM = 0
         CALL PDTGORD( ILQ, ILZ, IWORK, PARA, JW, 
     $      tH, DESCH, tT, DESCT, tQ, DESCQ,  
     $      tZ, DESCZ, 
     $      SR, SI, SBETA,
     $      DUMMYM, DWORK, -1, IWORK, -1, IERR )
         IF ( DEBUG ) WRITE(*,*)'% PQZ3 -> REORDER WSPACEQ Done', 
     $      INT( DWORK( 1 ) ), INT( IWORK( 1 ) )
         LIWORK1 = MAX( LIWORK1, IWORK( 1 ) )
         
         IWORK( 1 ) = LIWORK1
         LDWORK1 = MAX( LDWORK1, DWORK( 1 ) )
         DWORK( 1 ) = LDWORK1
         IF ( DEBUG ) THEN
            WRITE(*,*)'% PQZ3 : LDWORK, LIWORK:', INT(LDWORK1), LIWORK1
         END IF
         RETURN
      END IF
      
      
      IF ( DEBUG ) THEN
         WRITE(*,*)'% PQZ3: Enter', KTOP, KBOT, KWTOP, NW, RLVL
         CALL FLUSH(6)
      END IF
      DUMMYM = 0
      
*     
*     Get scaling factor
*     
      SCAL = ZERO
      IF ( KWTOP .GT. 1 ) THEN
         CALL PDELGET( 'A', ' ', SCAL, H, KWTOP, KWTOP - 1,DESCH )
      END IF
      
      T1 = MPI_WTIME()  
*     
*     Clean TmpQ               
*     
      CALL PDLASET( 'A', JW + IROFFH, JW + IROFFH, ZERO, ONE, 
     $   tQ, 1, 1, tDESCQ )   
*     
*     Clean TmpZ     
*     
      CALL PDLASET( 'A', JW + IROFFH, JW + IROFFH, ZERO, ONE, 
     $   tZ, 1, 1, tDESCZ )      

*     
*     Copy H to tmpH
*     
      IF (DEBUG) THEN 
         WRITE(*,*)'% PQZ3 : Copy H to TMPH'
         CALL FLUSH(6)
      END IF
      CALL PDLASET( 'A', IROFFH, JW + IROFFH, ZERO, ONE, 
     $   tH, 1, 1, tDESCH )
      CALL PDLASET( 'A', IROFFH, JW - IROFFH, ZERO, ZERO, 
     $   tH, 1, 1 + IROFFH, tDESCH )
      CALL PDLASET( 'A', JW, IROFFH, ZERO, ZERO, 
     $   tH, 1 + IROFFH, 1, tDESCH )
      CALL PDLACPY( 'A', 1, JW, H, KWTOP, KWTOP, DESCH, 
     $   tH, 1 + IROFFH, 1 + IROFFH, tDESCH )
      CALL PDLACPY( 'U', JW - 1, JW - 1, H, KWTOP + 1, KWTOP, DESCH, 
     $   tH, 1 + IROFFH + 1, 1 + IROFFH, tDESCH )

      IF (JW.GT.2)
     $   CALL PDLASET( 'L', JW - 2, JW - 2, ZERO, ZERO, 
     $   tH, 1 + IROFFH + 2, 1 + IROFFH, tDESCH )
      CALL PDLACPY( 'A', JW - 1, 1, H, KWTOP + 1, KWTOP + JW - 1, DESCH, 
     $   tH, 1 + IROFFH + 1, 1 + IROFFH + JW - 1, tDESCH )

*     
*     Copy T to tmpT
*     
      IF (DEBUG) THEN
         WRITE(*,*)'% PQZ3 : Copy T to TMPT'
         CALL FLUSH(6)
      ENDIF
      CALL PDLASET( 'A', JW + IROFFH, JW + IROFFH, ZERO, ONE, 
     $   tT, 1, 1, tDESCT )
      CALL PDLACPY( 'U', JW, JW, T, KWTOP, KWTOP, DESCT, 
     $   tT, 1 + IROFFH, 1 + IROFFH, tDESCT )


      T2 = MPI_WTIME()  
      UPDTIME(3) = UPDTIME(3) + ( T2 - T1 )
      T1 = MPI_WTIME()        

*     
*     If the problem is small, a smaller NB might be benificial.
*     In this case, when N <= 2000, we use an NB of 60 when 
*     redistributing 
*     
      SMALLNB = 60
      
      IF ( JW+IROFFH .GT. 2000) THEN 
         NPMIN = ICEIL( JW, ICEIL( 384, NB ) * NB )
      ELSE
         NPMIN = ICEIL( JW, ICEIL( 384, SMALLNB ) * SMALLNB )
      END IF


*     Setting a high value of NPMIN will prevent redistribution
*     Uncomment for skipping redistribution
c     NPMIN = 500

*     
*     Apply the QZ algorithm on the window
*     Chooose routine depending upon problem size, and decide whether to 
*     redistribute data before continuing 
*     
      IERR = 0
      IF( MIN( NPROW, NPCOL ) .LE. NPMIN + 1 .OR.
     $   RLVL .GT. 0 .OR. JW + IROFFH .LE. NMIN3 ) THEN
         IF ( DEBUG ) THEN 
            WRITE(*,*)'% PQZ3 : Calling PQ1', RUNDHGEQZ1,
     $         1+IROFFH, JW+IROFFH
            CALL FLUSH(6)
         ENDIF
         
         IF ( .NOT. RUNDHGEQZ1 ) THEN
            CALL PDHGEQZ0( WANTHT, ILQ, ILZ,
     $         JW + IROFFH, 1 + IROFFH, JW + IROFFH,
     $         tH, tDESCH, tT, tDESCT,
     $         SR( KWTOP - IROFFH ), SI( KWTOP - IROFFH ), 
     $         SBETA( KWTOP - IROFFH ),
     $         tQ, tDESCQ, tZ, tDESCZ, 
     $         DWORK, LDWORK,
     $         IWORK, LIWORK, IERR, RLVL + 1 )
         ELSE                  
            CALL PDHGEQZ1( WANTHT, ILQ, ILZ,
     $         JW + IROFFH, 1 + IROFFH, JW + IROFFH, tH, tDESCH,
     $         tT, tDESCT, SR( KWTOP - IROFFH ), 
     $         SI( KWTOP - IROFFH ),
     $         SBETA( KWTOP - IROFFH ), tQ, tDESCQ, tZ, tDESCZ, 
     $         DWORK, LDWORK, IWORK, LIWORK, IERR )            
         END IF

         T2 = MPI_WTIME()  
         UPDTIME(4) = UPDTIME(4) + ( T2 - T1)

         IF (DEBUG) THEN
            WRITE(*,*)'% PQZ3 : After QZ', IERR, JW+IROFFH
            CALL FLUSH(6)
         END IF

         IF ( IERR.NE.0 ) THEN
            IF ( IAM .EQ. 0 ) THEN
               WRITE(*,*)'PQZ3,0:PDHGEQZ1 failed with IERR=', IERR
            END IF
            INFO = 1
            NS = JW
            ND = 0
            GOTO 9999
         END IF
         
         IF ( SCAL .EQ. ZERO )  THEN
            NS = JW
            ND = NS
            GOTO 9999
         END IF
         T1 = MPI_WTIME()  

         IF (DEBUG) THEN
            WRITE(*,*)'% PQZ3 : Reorder', JW + IROFFH
            CALL FLUSH(6)
         END IF

         
*        
*        Clean up the array SELECT for reordering.
*        
         DO I = 1, NS + IROFFH
            IWORK( I ) = 0
         END DO
         
*        
*        Outer deflation detection loop (label 18).
*        In this loop a bunch of undeflatable eigenvalues
*        are moved simultaneously.
*        
         
*        
*        Lock the first entries ..
*        
         DO I = 1, IROFFH
            IWORK( I ) = 1
         END DO
         
         NS = JW
         ILST = 1+IROFFH
         
         IF (ILST.GT.1) THEN
            CALL PDELGET( 'A', ' ', TMPH, tH, ILST, ILST - 1, tDESCH )
            BULGE = TMPH( 1, 1 ) .NE. ZERO
            IF ( BULGE ) ILST = ILST + 1
         END IF
 18      CONTINUE

         IF ( ILST .LE. NS + IROFFH ) THEN
*           
*           Find the top-left corner of the local window.
*           
            LILST = MAX( ILST, NS + IROFFH - NB + 1 )
            IF( LILST .GT. 1 ) THEN
               CALL PDELGET( 'All', ' ', TMPH, tH, 
     $            LILST, LILST - 1,
     $            tDESCH )
               BULGE = TMPH( 1, 1 ) .NE. ZERO
               IF( BULGE ) LILST = LILST + 1
            END IF
*           
*           Lock all eigenvalues outside the local window.
*           
            DO I = IROFFH+1, LILST-1
               IWORK( I ) = 1
            END DO
            LILST0 = LILST
            
*           
*           Inner deflation detection loop (label 80).
*           In this loop, the undeflatable eigenvalues are moved to the
*           top-left corner of the local window.
*           
 80         CONTINUE
            IF( LILST .LE. NS + IROFFH ) THEN                       
               CALL PDELGET( 'All', ' ', TMPQ1, tQ, 1 + IROFFH,
     $            NS + IROFFH, tDESCQ )
               IF ( NS .EQ. 1 ) THEN
                  BULGE = .FALSE.
                  CALL PDELGET( 'All', ' ', TMPH( 2, 2 ), tH, 
     $               NS + IROFFH,
     $               NS + IROFFH, tDESCH )
               ELSE
                  CALL PDLACP3( 2, NS - 1 + IROFFH, tH, tDESCH, TMPH, 
     $               2, -1, -1, 0 )
                  BULGE = TMPH( 2, 1 ).NE.ZERO
                  IF ( BULGE ) THEN
                     CALL PDELGET( 'All', ' ', TMPQ2, tQ, 
     $                  1 + IROFFH,
     $                  NS - 1 + IROFFH, tDESCQ )
                  END IF
               END IF
               IF ( .NOT. BULGE ) THEN
                  TST1 = ABS( TMPH( 2, 2 ) )
                  IF ( TST1 .EQ. ZERO ) TST1 = ABS( SCAL )            
                  DFLATE = ABS( SCAL * TMPQ1 ) .LE. 
     $               MAX( ULP * TST1, SMLNUM )
               ELSE
                  IF ( ( TMPH( 1, 1 ) .EQ.ZERO )  .AND.
     $               ( TMPH( 2, 1 ) .EQ.ZERO ) ) THEN
                     RTDET = ZERO
                  ELSE IF ( ABS( TMPH( 1, 1 ) ).GE.
     $                  ABS( TMPH( 2, 1 ) ) ) THEN
                     RTDET = SQRT( ABS( TMPH( 1, 1 ) ) ) *
     $                  SQRT( ABS( TMPH( 2, 2 ) - ( TMPH( 2, 1 )/
     $                  TMPH( 1, 1) )*TMPH( 1, 2 ) ) )
                  ELSE
                     RTDET = SQRT( ABS( TMPH( 2, 1 ) ) ) *
     $                  SQRT( ABS( TMPH( 1, 2 ) - ( TMPH( 1, 1 )/
     $                  TMPH( 2, 1 ) )*TMPH( 2, 2 ) ) )
                  END IF
                  TST1 = ABS(SCAL) + RTDET
                  DFLATE = MAX( ABS( SCAL * TMPQ1 ), 
     $               ABS( SCAL * TMPQ2 ) ).LE.
     $               MAX( ULP * TST1, SMLNUM ) 
               END IF
               IF ( DFLATE ) THEN
*                 
*                 Deflation.
*                 
                  IF ( BULGE ) THEN
                     NS = NS - 2
                  ELSE
                     NS = NS - 1
                  END IF
               ELSE
*                 
*                 Move undeflatable eigenvalue up.
*                 
                  IFST = NS

*                 
*                 Clear select vector
*                 
                  DO I = LILST, JW + IROFFH
                     IWORK( I ) = 0
                  END DO
                  IWORK( IFST + IROFFH ) = 1
                  IF ( BULGE ) THEN
                     IWORK( IFST + IROFFH - 1 ) = 1
                  END IF
                  CALL PDTGORD( ILQ, ILZ, IWORK, PARA, 
     $               JW + IROFFH, tH, tDESCH,
     $               tT, tDESCT, tQ, tDESCQ, 
     $               tZ, tDESCZ, 
     $               DWORK,
     $               DWORK(JW + IROFFH + 1),
     $               DWORK(2 * (JW + IROFFH ) + 1 ),            
     $               DUMMYM, DWORK( 3 * (JW + IROFFH) + 1 ), 
     $               LDWORK, IWORK( JW + IROFFH + 1 ), 
     $               LIWORK, IERR )
                  IF ( IERR .NE. 0 ) THEN
                     IF ( IAM .EQ. 0 ) THEN
                        WRITE(*,*)
     $                     '% PQZ3 , 0:Error swaping 1 with IERR=',
     $                     IERR, JW - NS
                     END IF
                     GOTO 20
                  END IF
                  IWORK( IFST + IROFFH ) = 0
                  IF ( BULGE ) IWORK( IFST + IROFFH - 1 ) = 0
                  IWORK( LILST ) = 1
                  IF ( BULGE ) IWORK( LILST + 1 ) = 1               
                  
                  IF ( BULGE ) THEN
                     LILST = LILST + 2
                  ELSE
                     LILST = LILST + 1
                  END IF
               END IF
*              
*              End of inner deflation detection loop.
*              
               GO TO 80            
            END IF
*           
*           Unlock the eigenvalues outside the local window.
*           Then undeflatable eigenvalues are moved to the proper position.
*           
            DO I = ILST, LILST0 - 1
               IWORK( I ) = 0
            END DO         
            
            CALL PDTGORD( ILQ, ILZ, IWORK, PARA, 
     $         JW + IROFFH, tH, tDESCH,
     $         tT, tDESCT, tQ, tDESCQ, 
     $         tZ, tDESCZ, 
     $         DWORK,
     $         DWORK( JW + IROFFH + 1 ),
     $         DWORK( 2 * ( JW + IROFFH ) + 1 ),            
     $         DUMMYM, DWORK( 3 *( JW + IROFFH ) + 1 ), 
     $         LDWORK, IWORK( JW + IROFFH + 1 ), 
     $         LIWORK, IERR ) 
            
            IF ( IERR .NE. 0 ) THEN
               IF ( IAM .EQ. 0 ) THEN
                  WRITE(*,*)'% PQZ3,0:Error swaping 2 with IERR=',
     $               IERR, JW - NS
               END IF
               GOTO 20
            END IF                      
            ILST = DUMMYM + 1
            
*           
*           End of outer deflation detection loop.
*           
            GO TO 18
         END IF
 20      CONTINUE
         IF (DEBUG) THEN
            WRITE(*,*)'% PQZ3 : After reordering', NS, JW, IERR
            CALL FLUSH(6)
         END IF
*        
*        Post-reordering step: copy output eigenvalues to output.
*        
         CALL DCOPY(JW, DWORK( 1 + IROFFH ), 1, SR( KWTOP ), 1)
         CALL DCOPY(JW, DWORK( JW + 2 * IROFFH + 1 ), 1, SI( KWTOP ), 1)
         CALL DCOPY(JW, DWORK( JW * 2 + 3 * IROFFH + 1 ), 1, 
     $      SBETA( KWTOP ), 1)
         T2 = MPI_WTIME()  
         UPDTIME(5) = UPDTIME(5) + (T2-T1)
         IF( NS .EQ. 0 ) SCAL = ZERO

         IF ( DEBUG ) THEN
            WRITE(*,*)'% PQZ3 : Fix spike'
            CALL FLUSH(6)
         END IF
*        
*        Check if have we got deflated eigenvalues
*        
         IF ( NS .LT. JW .OR. SCAL .EQ. 0 ) THEN
            IF ( NS .GT. 1 .AND. SCAL .NE. ZERO ) THEN
               T1 = MPI_WTIME()
               NR = NUMROC( NS + IROFFH, NB, MYROW, tDESCQ( RSRC_ ), 
     $            NPROW )
               NC = NUMROC(1, 1, MYCOL, tDESCQ( CSRC_ ), NPCOL )      
               CALL DESCINIT( DESCV, NS + IROFFH, 1, NB, 1,
     $            tDESCQ( RSRC_ ), tDESCQ( CSRC_ ), ICTXT, 
     $            MAX( 1, NR ), IERR )
               NRTAU = NUMROC( 1, 1, MYROW, tDESCQ( RSRC_ ), NPROW )
               NCTAU = NUMROC( JW + IROFFH, NB, MYCOL, tDESCQ( CSRC_ ), 
     $            NPCOL )
               CALL DESCINIT( DESCTAU, 1, JW + IROFFH, 1, NB,
     $            tDESCQ( RSRC_ ), tDESCQ( CSRC_ ), ICTXT, 
     $            MAX( 1,NRTAU ), IERR )
               IR = 1
               IPTAU =  IR + DESCV( LLD_ ) * NC
               IPT = IPTAU + NRTAU * NCTAU
               IPWORK = IPT + tDESCT( MB_ ) * tDESCT( MB_ )

               CALL PDLASET( 'A', NS + IROFFH, 1, ZERO, ZERO, 
     $            DWORK( IR ), 1, 1, DESCV )


               CALL PDLASET( 'A', 1, JW + IROFFH, ZERO, ZERO, 
     $            DWORK( IPTAU ), 1, 1, DESCTAU )
*              
*              Copy first column from Q to calculate  Householder reflectors
*              
               CALL PDCOPY( NS, 
     $            tQ, 1 + IROFFH, 1 + IROFFH, tDESCQ, 
     $            tDESCQ( M_ ),
     $            DWORK( IR ), 1 + IROFFH, 1, DESCV, 
     $            1 )
*              
*              Calculate the householder reflectors and TAU      
*              
               CALL PDLARFG( NS, ALPHA, 1 + IROFFH, 1, DWORK(IR),
     $            2 + IROFFH, 1, DESCV, 1, DWORK(IPTAU) )

               CALL PDELSET( DWORK( IR ),1 + IROFFH, 1, DESCV, ONE ) 
               CALL PDLASET( 'L', JW - 2, JW - 2, ZERO, ZERO, 
     $            tH, 3 + IROFFH, 1 + IROFFH, tDESCH )
*              
*              Apply the householder reflectors to H     
*              
               CALL PDLARF( 'L', NS, JW, 
     $            DWORK( IR ), 1 + IROFFH, 1, DESCV, 
     $            1, DWORK(IPTAU), 
     $            tH, 1 + IROFFH, 1 + IROFFH, tDESCH, 
     $            DWORK( IPWORK ) )
*              
*              Apply the householder reflectors to Q
*              

               CALL PDLARF( 'R', JW, NS, DWORK( IR ), 1 + IROFFH, 1, 
     $            DESCV, 1, DWORK( IPTAU ), 
     $            tQ, 1 + IROFFH, 1 + IROFFH, tDESCQ, 
     $            DWORK( IPWORK ) )
               IF ( NS .GT. 1) THEN
*                 
*                 Apply the householder reflectors to T, causing fill-in
*                 
                  CALL PDLARF( 'L', NS, JW, 
     $               DWORK( IR ), 1 + IROFFH, 1, DESCV, 
     $               1, DWORK( IPTAU ), 
     $               tT, 1 + IROFFH, 1 + IROFFH, tDESCT, 
     $               DWORK( IPWORK ) )
*                 
*                 RQ factor T                              
*                 
                  CALL PDGERQF( NS, NS, 
     $               tT, 1 + IROFFH, 1 + IROFFH, tDESCT, 
     $               DWORK( IR ), DWORK( IPWORK ), LDWORK, IERR )  
*                 
*                 apply reflections on A(H) from above RQ call  
*                 
                  CALL PDORMRQ( 'R', 'T', NS, NS, NS, 
     $               tT, 1 + IROFFH, 1 + IROFFH, tDESCT, 
     $               DWORK,
     $               tH, 1 + IROFFH, 1 + IROFFH, tDESCH,
     $               DWORK( IPWORK ), LDWORK, IERR )
*                 
*                 apply reflections on Z from above RQ call             
*                 
                  CALL PDORMRQ( 'R', 'T', JW, NS, NS, 
     $               tT, 1 + IROFFH, 1 + IROFFH, tDESCT,
     $               DWORK,
     $               tZ, 1 + IROFFH, 1 + IROFFH, tDESCZ,
     $               DWORK( IPWORK ), LDWORK, IERR )
*                 
*                 Remove elements below diagonal in T               
*                 
                  CALL PDLASET( 'L', NS - 1, NS - 1, ZERO, ZERO, 
     $               tT, 2 + IROFFH, 1 + IROFFH, tDESCT )
               END IF
               T2 = MPI_WTIME()  
               UPDTIME(6) = UPDTIME(6) + (T2-T1)
*              
*              Hessenberg-triangular reduction.      
*              
               T1 = MPI_WTIME() 
               
               IF ( DEBUG ) THEN
                  WRITE(*,*)'% PQZ3 : Before HT', JW + IROFFH, 
     $               1 + IROFFH, 
     $               NS + IROFFH, JW - NS, KWTOP, JW
               END IF
*              
*              Allocate space for transpose of Q
*              
               ltQROWS = NUMROC( JW + IROFFH, NB, MYROW, 
     $            tDESCQ( RSRC_ ), NPROW )
               ltQCOLS = NUMROC( JW + IROFFH, NB, MYCOL, 
     $            tDESCQ( CSRC_ ), NPCOL )
               

               ALLOCATE( tQ2( MAX(ltQROWS * ltQCOLS, 1 ) ) )
               
*              
*              Transpose Q before call to hessenberg-triangular reduction
               CALL PDGEADD('T', JW + IROFFH, JW + IROFFH, ONE, 
     $            tQ, 1, 1, tDESCQ, ZERO, tQ2, 1, 1, tDESCQ)
*              Reduce sub to HT form again
               CALL PDGGHRD( CILQ, CILZ, JW + IROFFH, 1 + IROFFH, 
     $            NS + IROFFH, tH, tDESCH, tT, tDESCT, tQ2, tDESCQ, tZ, 
     $            tDESCZ, DWORK, LDWORK, IERR )
*              
*              Retranspose Q after HT reduction
               CALL PDGEADD('T', JW + IROFFH, JW + IROFFH, ONE, 
     $            tQ2, 1, 1, tDESCQ, ZERO, tQ, 1, 1, tDESCQ)

*              Deallocate tQ2 since not needed anymore
               DEALLOCATE( tQ2 )

               IF ( DEBUG ) THEN
                  WRITE(*,*)'% PQZ3 : After PDGGHRD', IERR
                  CALL FLUSH(6)
               END IF
               
               T2 = MPI_WTIME()  
               UPDTIME(7) = UPDTIME(7) + ( T2 - T1 )
               IF ( IERR .NE. 0 ) THEN
                  IF (IAM.EQ.0) THEN
                     WRITE(*,*)'PQZ3,0:PDGGHRD failed with IERR=',
     $                  IERR
                  END IF
                  INFO = 4
                  NS = JW
                  ND = 0 
                  GOTO 9999 
               END IF
               IPWORK = 1
            END IF      
         END IF
*        
*        .. Redistribute data before continuing
*        
      ELSE  
         ICTXT_NEW = 0 
         UPDTIME(21) = UPDTIME(21) + 1
         T1 = MPI_WTIME()
         CALL BLACS_GET(ICTXT, 10, ICTXT_NEW)
         DO I = 0, NPMIN - 1
            DO J = 0, NPMIN - 1
               PMAP( J + 1 + I * NPMIN ) = BLACS_PNUM( ICTXT, I, J )
            END DO
         END DO
         IF ( DEBUG ) THEN
            WRITE(*,*)'% PQZ3: Redistr.', 
     $         NPROW*NPCOL,'->', NPMIN*NPMIN, JW, ICTXT, ICTXT_NEW
            CALL FLUSH(6)
         END IF
         
         
         CALL BLACS_GRIDMAP( ICTXT_NEW, PMAP, NPMIN, NPMIN, NPMIN )
         CALL BLACS_GRIDINFO( ICTXT_NEW, NPMIN, NPMIN, MYROW_NEW, 
     $      MYCOL_NEW )
         
         IF( MYROW.GE.NPMIN .OR. MYCOL.GE.NPMIN ) ICTXT_NEW = -1
         IF( ICTXT_NEW .NE. -1 ) THEN
            IF ( JW .GT. 2000 ) THEN 
               NB_NEW = NB
            ELSE
               NB_NEW = SMALLNB
            END IF
            HTZQROWS0 = NUMROC( JW, NB_NEW, MYROW_NEW, 0, NPMIN )
            HTZQCOLS0 = NUMROC( JW, NB_NEW, MYCOL_NEW, 0, NPMIN )
            
            CALL DESCINIT( DESCHTZQ, JW, JW, NB_NEW, NB_NEW, 0,
     $         0, ICTXT_NEW, MAX(1,HTZQROWS0), IERR )
            
            IPTH = 1
            IPTT = IPTH + MAX(1,HTZQROWS0)*MAX(1,HTZQCOLS0)
            IPTQ = IPTT + MAX(1,HTZQROWS0)*MAX(1,HTZQCOLS0)
            IPTQ2 = IPTQ + MAX(1,HTZQROWS0)*MAX(1,HTZQCOLS0)
            IPTZ = IPTQ2 + MAX(1,HTZQROWS0)*MAX(1,HTZQCOLS0)
            IPWORK = IPTZ + MAX(1,HTZQROWS0)*MAX(1,HTZQCOLS0) 
         ELSE
            IPTH = 1
            IPTT = 2
            IPTQ = 3
            IPTQ2 = 4
            IPTZ = 5
            IPWORK =6 
            DESCHTZQ(CTXT_ ) = -1
         END IF
         CALL PDGEMR2D( JW, JW, tH, 1+IROFFH, 1+IROFFH, tDESCH,
     $      DWORK(IPTH), 1, 1, DESCHTZQ, ICTXT )
         CALL PDGEMR2D( JW, JW, tT, 1+IROFFH, 1+IROFFH, tDESCT,
     $      DWORK(IPTT), 1, 1, DESCHTZQ, ICTXT )
         
         IF (DEBUG) THEN
            WRITE(*,*)'% PQZ3 : Redistr. finished',
     $         NPMIN, JW, ',NB=',NB_NEW
            CALL FLUSH(6)
         END IF   
         
         
         
         T2 = MPI_WTIME()  
         UPDTIME(20) = UPDTIME(20) + (T2-T1)
         
         IF( ICTXT_NEW .NE. -1 ) THEN
            T1 = MPI_WTIME()
            CALL PDLASET( 'All', JW, JW, ZERO, ONE, DWORK(IPTQ), 
     $         1, 1, DESCHTZQ )
            CALL PDLASET( 'All', JW, JW, ZERO, ONE, DWORK(IPTZ), 
     $         1, 1, DESCHTZQ )  

            IF (DEBUG) THEN 
               WRITE(*,*)'% PQZ3 : QZ', JW, RUNDHGEQZ1 
               CALL FLUSH(6)
            END IF

            IF (RUNDHGEQZ1) THEN
               CALL PDHGEQZ1( WANTHT, ILQ, ILZ, JW, 1, JW,
     $            DWORK( IPTH ), DESCHTZQ, 
     $            DWORK( IPTT ), DESCHTZQ,SR( KWTOP ), SI( KWTOP ),
     $            SBETA( KWTOP ), 
     $            DWORK( IPTQ ), DESCHTZQ, DWORK( IPTZ ), DESCHTZQ, 
     $            DWORK( IPWORK ), LDWORK, IWORK, LIWORK, IERR )
            ELSE
               CALL PDHGEQZ0( WANTHT, ILQ, ILZ, JW, 1, JW,
     $            DWORK( IPTH ), DESCHTZQ, DWORK( IPTT ), DESCHTZQ,
     $            SR( KWTOP ), SI( KWTOP ), SBETA( KWTOP ),
     $            DWORK( IPTQ ), DESCHTZQ, DWORK( IPTZ ), DESCHTZQ, 
     $            DWORK( IPWORK ), LDWORK,
     $            IWORK, LIWORK, IERR, RLVL + 1 ) 
            END IF
            T2 = MPI_WTIME()  
            UPDTIME(4) = UPDTIME(4) + (T2-T1)
            
            IF (DEBUG) THEN
               WRITE(*,*)'% PQZ3 : After PDHGEQZ', NPMIN, JW, IERR
               CALL FLUSH(6)
            END IF
            IF ( IERR.NE.0 ) THEN
               IF (IAM.EQ.0) THEN
                  WRITE(*,*)'PQZ3,1:PDHGEQZ failed with IERR=',IERR,
     $               LDWORK, DWORK( IPWORK )
               END IF 
               INFO = 5
               GOTO 50
            END IF

            IF ( SCAL .EQ. ZERO )  THEN               
               NS = 0
               ND = JW
               GOTO 100
            END IF

            T1 = MPI_WTIME()             
*           
*           Clean up the array SELECT for reordering.
*           
            DO I = 1, NS
               IWORK( I ) = 0
            END DO


*           Maximum number of independent computational windows 
            PARA( 1 ) = PILAENVX( ICTXT_NEW, 80, 'PDHGEQZ', '', JW, 
     $         NB_NEW, IDUM, IDUM)
*           Number of eigenvalues in each window
            PARA( 2 ) = PILAENVX( ICTXT_NEW, 81, 'PDHGEQZ', '', IDUM, 
     $         NB_NEW, IDUM, IDUM )
*           Computational window size
            PARA( 3 ) = PILAENVX( ICTXT_NEW, 82, 'PDHGEQZ', '', IDUM, 
     $         NB_NEW, IDUM, IDUM )
*           Width of block column slabs for row-wise
*           application of pipelined orthogonal transformations in
*           their factorized form
            PARA( 5 ) = PILAENVX( ICTXT_NEW, 84, 'PDHGEQZ', '', IDUM, 
     $         NB_NEW, IDUM, IDUM )
*           Maximum number of eigenvalues to bring over
*           the block border 
            PARA( 6 ) = PILAENVX( ICTXT_NEW, 85, 'PDHGEQZ', '', IDUM, 
     $         NB_NEW, IDUM, IDUM )
            IF ( DEBUG ) WRITE(*,*)'%PQZ3 : Reorder Para=', PARA(1:6), 
     $         JW, NB_NEW

            NS = JW
            ILST = 1              
            
*           
*           Outer deflation detection loop (label 118).
*           In this loop a bunch of undeflatable eigenvalues
*           are moved simultaneously.
*           

 118        CONTINUE
c           IF (DEBUG) WRITE(*,*)'% PQZ3 : Reorder , outer:', ILST, NS
            IF ( ILST .LE. NS ) THEN
               
*              
*              Find the top-left corner of the local window.
*              
               LILST = MAX( ILST, NS - NB_NEW + 1 )
               IF( LILST.GT.1 ) THEN
                  CALL PDELGET( 'All', ' ', TMPH, DWORK( IPTH ), 
     $               LILST, LILST - 1,
     $               DESCHTZQ )
                  BULGE = TMPH(1,1).NE.ZERO
                  IF( BULGE ) LILST = LILST+1
               END IF
*              
*              Lock all eigenvalues outside the local window.
*              
               DO I = 1, LILST - 1
                  IWORK( I ) = 1
               END DO
               LILST0 = LILST
               
*              
*              Inner deflation detection loop (label 180).
*              In this loop, the undeflatable eigenvalues are moved to the
*              top-left corner of the local window.
*              
 180           CONTINUE
c              IF (DEBUG) WRITE(*,*) '% PQZ3 : Reorder , inner:', 
c              $            LILST, NS
               IF( LILST .LE. NS ) THEN 
                  
                  CALL PDELGET( 'All', ' ', TMPQ1, 
     $               DWORK(IPTQ), 1, NS, DESCHTZQ )
                  
                  IF ( NS.EQ.1 ) THEN
                     BULGE = .FALSE.
                     CALL PDELGET( 'A', ' ', TMPH( 2, 2 ), 
     $                  DWORK( IPTH ), NS, NS, DESCHTZQ )
                  ELSE
                     CALL PDLACP3( 2, NS - 1, DWORK( IPTH ), DESCHTZQ,
     $                  TMPH, 2, -1, -1, 0 )
                     BULGE = TMPH(2, 1) .NE. ZERO
                     IF (BULGE) THEN
                        CALL PDELGET( 'A', ' ', TMPQ2, 
     $                     DWORK(IPTQ), 1, NS - 1, DESCHTZQ )
                     END IF
                  END IF
                  IF ( .NOT.BULGE ) THEN
                     
                     TST1 = ABS( TMPH( 2, 2 ) )
                     IF ( TST1 .EQ. ZERO ) TST1 = ABS( SCAL )
                     DFLATE = ABS( SCAL * TMPQ1) .LE.
     $                  MAX( ULP * TST1, SMLNUM )
                  ELSE
                     
                     IF ( ( TMPH(1,1).EQ.ZERO )  .AND.
     $                  ( TMPH(2,1).EQ.ZERO ) ) THEN
                        RTDET = ZERO
                     ELSE IF ( ABS( TMPH(1,1) ).GE.
     $                     ABS( TMPH(2,1) ) ) THEN
                        RTDET = SQRT( ABS( TMPH(1,1) ) )*
     $                     SQRT( ABS( TMPH(2,2) - ( TMPH(2,1)/
     $                     TMPH(1,1) )*TMPH(1,2) ) )
                     ELSE
                        RTDET = SQRT( ABS( TMPH(2,1) ) )*
     $                     SQRT( ABS( TMPH(1,2) - ( TMPH(1,1)/
     $                     TMPH(2,1) )*TMPH(2,2) ) )
                     END IF
                     TST1 = ABS(SCAL) + RTDET
                     DFLATE = MAX( ABS( SCAL*TMPQ1 ), 
     $                  ABS( SCAL*TMPQ2 ) ).LE.
     $                  MAX( ULP*TST1, SMLNUM ) 
                  END IF
                  IF ( DFLATE ) THEN
*                    
*                    Deflation
*                    
                     IF ( BULGE ) THEN
                        NS = NS - 2
                     ELSE
                        NS = NS - 1
                     END IF
                  ELSE
                     
*                    
*                    Move undeflatable eigenvalue up.
*                    
                     IFST = NS
*                    
*                    Clear select vector
*                    
                     DO I = LILST, JW
                        IWORK( I ) = 0
                     END DO
                     IWORK( IFST ) = 1
                     IF ( BULGE ) THEN
                        IWORK( IFST - 1 ) = 1
                     END IF                      
                     IERR = 0

                     CALL PDTGORD( ILQ, ILZ, IWORK, PARA, 
     $                  JW, DWORK( IPTH ), DESCHTZQ,
     $                  DWORK( IPTT ), DESCHTZQ, 
     $                  DWORK( IPTQ ), DESCHTZQ, 
     $                  DWORK( IPTZ ), DESCHTZQ, 
     $                  DWORK( IPWORK ),
     $                  DWORK( IPWORK + JW + 1 ),
     $                  DWORK( IPWORK + 2 * JW + 1 ),            
     $                  DUMMYM, DWORK( IPWORK + 3 * JW + 1 ), 
     $                  LDWORK, IWORK( JW + 1 ), 
     $                  LIWORK, IERR )    

                     IF ( IERR .NE. 0 ) THEN
                        IF (IAM.EQ.0) THEN
                           WRITE(*,*)
     $                        '% PQZ3,1:Error swap 1 with IERR=',
     $                        IERR, JW - NS
                        END IF   
                        GOTO 200
                     END IF
                     IWORK(IFST) = 0
                     IF ( BULGE ) IWORK(IFST-1) = 0
                     IWORK(LILST) = 1
                     IF ( BULGE ) IWORK(LILST+1) = 1               
                     IF ( BULGE ) THEN
                        LILST = LILST + 2
                     ELSE
                        LILST = LILST + 1
                     END  IF                              
                  END IF 
*                 
*                 End of inner deflation detection loop.
*                 
                  GO TO 180            
               END IF
*              
*              Unlock the eigenvalues outside the local window.
*              Then undeflatable eigenvalues are moved to the proper position.
*              
               DO I = ILST, LILST0 - 1
                  IWORK( I ) = 0
               END DO         
               IERR = 2
               IF ( DEBUG ) WRITE(*,*)'% PQZ3: calling PDTGORD', ILST, 
     $            LILST0, ILST, NS, JW
               CALL PDTGORD( ILQ, ILZ, IWORK, PARA, 
     $            JW, DWORK( IPTH ), DESCHTZQ,
     $            DWORK( IPTT ), DESCHTZQ, 
     $            DWORK( IPTQ ), DESCHTZQ, 
     $            DWORK( IPTZ ), DESCHTZQ, 
     $            DWORK( IPWORK ),
     $            DWORK( IPWORK + JW + 1 ),
     $            DWORK( IPWORK + 2 * JW + 1 ),            
     $            DUMMYM, DWORK( IPWORK + 3 * JW + 1 ), 
     $            LDWORK, IWORK( JW + 1 ), 
     $            LIWORK, IERR )    
               IF ( IERR.NE.0 ) THEN
                  IF (IAM.EQ.0) THEN
                     WRITE(*,*)
     $                  '% PQZ3,1:Error swaping 2 with IERR=',
     $                  IERR, JW - NS
                  END IF   
                  GOTO 200
               END IF
               ILST = DUMMYM + 1
               IF ( DEBUG ) WRITE(*,*)'% PQZ3: after PDTGORD'
*              
*              End of outer deflation detection loop.
*              
               GO TO 118
            END IF
            
            
*           
*           Post-reordering step: copy output eigenvalues to output.
*           
 200        CONTINUE
            
            IF ( DEBUG ) THEN
               WRITE(*,*)'% PQZ3 :  Reorder complete', NS, JW, INFO, 
     $            IERR
               CALL FLUSH(6)
            END IF              
            CALL DCOPY( JW, DWORK( IPWORK ), 1, SR( KWTOP ), 1 )
            CALL DCOPY( JW, DWORK( IPWORK + JW + 1 ), 1, 
     $         SI( KWTOP ), 1 )
            CALL DCOPY( JW, DWORK( IPWORK + JW * 2 + 1 ), 1, 
     $         SBETA(KWTOP), 1)                           
            
            T2 = MPI_WTIME()  
            UPDTIME(5) = UPDTIME(5) + (T2-T1)
            ND = JW - NS            
            
            IF( NS .EQ. 0 ) SCAL = ZERO
*           
*           Check if we have got deflated eigenvalues
*           
            IF ( ND .GT. 0 .OR. SCAL .EQ. 0 ) THEN
               IF ( NS .GT. 1 .AND. SCAL .NE. ZERO ) THEN
                  T1 = MPI_WTIME()
                  NR = NUMROC( NS , NB_NEW, MYROW_NEW, 
     $               DESCHTZQ(RSRC_), NPMIN )
                  NC = NUMROC( 1, 1, MYCOL_NEW, DESCHTZQ( CSRC_ ), 
     $               NPMIN )
                  CALL DESCINIT( DESCV, NS, 1, NB_NEW, 1,
     $               DESCHTZQ(RSRC_), DESCHTZQ(CSRC_), ICTXT_NEW, 
     $               MAX( 1, NR ), IERR )
                  NRTAU = NUMROC( 1, 1, MYROW_NEW, DESCHTZQ(RSRC_), 
     $               NPMIN )
                  NCTAU = NUMROC( JW, NB_NEW, MYCOL_NEW, 
     $               DESCHTZQ( CSRC_ ), NPMIN )
                  CALL DESCINIT( DESCTAU, 1, JW, 1, NB_NEW,
     $               DESCHTZQ(RSRC_), DESCHTZQ( CSRC_ ), ICTXT_NEW, 
     $               MAX( 1, NRTAU ), IERR )
                  IR = IPWORK
                  IPTAU =  IR + DESCV(LLD_ ) * NC
                  IPT = IPTAU + NRTAU * NCTAU
                  IPWORK = IPT + DESCHTZQ( MB_ ) * DESCHTZQ( MB_ )
                  
                  CALL PDLASET('A', NS, 1, ZERO, ZERO, 
     $               DWORK(IR), 1, 1, DESCV)
                  CALL PDLASET('A', 1, JW, ZERO, ZERO, 
     $               DWORK(IPTAU), 1, 1, DESCTAU)
*                 
*                 Copy first column from Q to calculate  Householder reflectors
*                 
                  CALL PDCOPY( NS, DWORK(IPTQ), 1, 1, DESCHTZQ, 
     $               DESCHTZQ(M_), DWORK(IR), 1, 1, DESCV, 1)
*                 
*                 Calculate the householder reflectors and TAU       
*                 
                  CALL PDLARFG( NS, ALPHA, 1, 1, 
     $               DWORK(IR), 2, 1, DESCV, 1, DWORK(IPTAU) )

                  CALL PDELSET(DWORK(IR),1,1, DESCV, ONE)      
                  CALL PDLASET('Lower',JW-2, JW-2, ZERO, ZERO, 
     $               DWORK(IPTH), 3, 1, DESCHTZQ)
*                 
*                 Apply the householder reflectors to H     
*                 
                  CALL PDLARF('L', NS, JW, DWORK(IR), 1, 1, 
     $               DESCV, 1, DWORK(IPTAU),
     $               DWORK(IPTH), 1, 1, DESCHTZQ, DWORK(IPWORK))
*                 
*                 Apply the householder reflectors to Q
*                 
                  CALL PDLARF('R', JW, NS, DWORK(IR), 1, 1, 
     $               DESCV, 1, DWORK(IPTAU), 
     $               DWORK(IPTQ), 1, 1, DESCHTZQ, DWORK(IPWORK))
                  IF (NS.GT.1) THEN
*                    
*                    Apply the householder reflectors to T, causing fill-in
*                    
                     CALL PDLARF('L', NS, JW, DWORK(IR), 1 , 1, 
     $                  DESCV, 1, DWORK(IPTAU), 
     $                  DWORK(IPTT), 1, 1, DESCHTZQ, DWORK(IPWORK))
*                    
*                    RQ factor T                              
*                    
                     CALL PDGERQF( NS, NS, DWORK(IPTT), 1, 1, 
     $                  DESCHTZQ, DWORK(IR), DWORK(IPWORK), LDWORK, 
     $                  IERR)   

*                    
*                    Apply reflections on A(H) from above RQ call  
*                    
                     CALL PDORMRQ( 'R', 'T', NS, NS, NS, 
     $                  DWORK(IPTT), 1, 1, DESCHTZQ, 
     $                  DWORK(IR),
     $                  DWORK(IPTH), 1, 1, DESCHTZQ,
     $                  DWORK(IPWORK), LDWORK, IERR )

*                    
*                    Apply reflections on Z from above RQ call                
*                    
                     CALL PDORMRQ( 'R', 'T', JW, NS, NS, 
     $                  DWORK(IPTT), 1, 1, DESCHTZQ, 
     $                  DWORK(IR),
     $                  DWORK(IPTZ), 1, 1, DESCHTZQ,
     $                  DWORK(IPWORK), LDWORK, IERR )
*                    
*                    Remove elements below diagonal in T               
*                    
                     CALL PDLASET( 'L', NS-1, NS-1, ZERO, ZERO, 
     $                  DWORK(IPTT),  2, 1, DESCHTZQ )
                  END IF
                  T2 = MPI_WTIME()  
                  UPDTIME(6) = UPDTIME(6) + (T2-T1)
                  
*                 
*                 Hessenberg-triangular reduction. 
*                 
                  T1 = MPI_WTIME()
                  IF (DEBUG) THEN
                     WRITE(*,*)'% PQZ3 : Before PDGGHRD',
     $                  IERR, JW, 1, NS
                     CALL FLUSH(6)
                  END IF

*                 Transpose Q before call to hessenberg-triangular reduction
                  CALL PDGEADD( 'T', JW, JW, ONE, 
     $               DWORK( IPTQ ), 1, 1, DESCHTZQ, 
     $               ZERO, DWORK( IPTQ2 ), 1, 1, DESCHTZQ )

                  CALL PDGGHRD( CILQ, CILZ, JW, 1, NS,
     $               DWORK( IPTH ), DESCHTZQ, DWORK( IPTT ), 
     $               DESCHTZQ, DWORK( IPTQ2 ), DESCHTZQ, DWORK( IPTZ ), 
     $               DESCHTZQ, DWORK( IPWORK ), LDWORK, IERR )

                  CALL PDGEADD( 'T', JW, JW, ONE, 
     $               DWORK( IPTQ2 ), 1, 1, DESCHTZQ, 
     $               ZERO, DWORK( IPTQ ), 1, 1, DESCHTZQ )

                  T2 = MPI_WTIME()  
                  IF ( DEBUG ) THEN
                     WRITE(*,*)'% PQZ3 : After PDGGHRD', 
     $                  IERR, JW, NS, MYROW, MYCOL
                     CALL FLUSH(6)
                  END IF
                  
                  UPDTIME(7) = UPDTIME(7) + (T2-T1)
                  
                  IF ( IERR .NE. 0 ) THEN
                     IF (IAM .EQ. 0 ) THEN
                        WRITE(*,*)
     $                     '% PQZ3,1:PDGGHRD failed with IERR=',
     $                     IERR
                     END IF   
                     INFO = 8
                     GOTO 50
                  END IF                  
               END IF      
            END IF            
         END IF
         

         
         IERR = 0
         GOTO 100
         
 50      CONTINUE
         ND = 0
         NS = JW         
         
         
 100     CONTINUE
         IPWORK = 1
         IF ( IAM .NE. 0 ) THEN
            NS = 0
            ND = 0
            IERR = 0 
            INFO = 0
         END IF
         CALL IGSUM2D( ICTXT, 'All', ' ', 1, 1, INFO, 1, -1, -1 )
         CALL IGSUM2D( ICTXT, 'All', ' ', 1, 1, IERR, 1, -1, -1 )
         CALL IGSUM2D( ICTXT, 'All', ' ', 1, 1, NS, 1, -1, -1 )
         CALL IGSUM2D( ICTXT, 'All', ' ', 1, 1, ND, 1, -1, -1 )
         
*        
*        Restore to original distribution
*        
         IF ( IERR .EQ. 0 .AND. (NS .LT. JW .OR. SCAL .EQ. 0 ) ) THEN 
            
            IF (DEBUG) THEN
               WRITE(*,*)'% PQZ3 : Restore redist.', MYROW, MYCOL
               CALL FLUSH(6)
            END IF

            T1 = MPI_WTIME()
            CALL PDGEMR2D( JW, JW, DWORK( IPTH ), 1, 1, DESCHTZQ,
     $         tH, 1 + IROFFH, 1 + IROFFH, tDESCH, ICTXT )
            CALL PDGEMR2D( JW, JW, DWORK( IPTT ), 1, 1, DESCHTZQ,
     $         tT, 1 + IROFFH, 1 + IROFFH, tDESCT, ICTXT ) 
            CALL PDGEMR2D( JW, JW, DWORK( IPTQ ), 1, 1, DESCHTZQ,
     $         tQ, 1 + IROFFH, 1 + IROFFH, tDESCQ, ICTXT )
            CALL PDGEMR2D( JW, JW, DWORK( IPTZ ), 1, 1, DESCHTZQ,
     $         tZ, 1 + IROFFH, 1 + IROFFH, tDESCZ, ICTXT )               
            T2 = MPI_WTIME()
            UPDTIME(20) = UPDTIME(20) + (T2-T1) 
            
         END IF

         IF( ICTXT_NEW .NE. - 1 ) CALL BLACS_GRIDEXIT( ICTXT_NEW )
         IF( IAM .NE. 0 ) THEN
            DO J = 0, JW - 1
               SR( KWTOP + J ) = ZERO
               SI( KWTOP + J ) = ZERO
               SBETA( KWTOP + J ) = ZERO
            END DO           
         END IF
         CALL DGSUM2D( ICTXT, 'All', ' ', JW, 1, SR( KWTOP ), JW, 
     $      -1, -1 )
         CALL DGSUM2D( ICTXT, 'All', ' ', JW, 1, SI( KWTOP ), JW, 
     $      -1, -1 )
         CALL DGSUM2D( ICTXT, 'All', ' ', JW, 1, SBETA( KWTOP ), JW, 
     $      -1, -1 )  
         
         IF ( IERR .NE. 0 ) GOTO 9999
      END IF

      IPWORK = 1

      IF ( NS .LT. JW .OR. SCAL .EQ. 0 ) THEN

         CALL BLACS_BARRIER( ICTXT, 'A' )
         IF (DEBUG) THEN
            WRITE(*,*)'% PQZ3 : Restore data', 
     $         JW, NS, IROFFH, MYROW, MYCOL
            CALL FLUSH(6)
         END IF
         
*        
*        Restore data
*        
         T1 = MPI_WTIME()  

*        H
*        First the upper diagonal H(KWTOP+1, KWTOP)
         CALL PDLACPY('U', JW - 1, JW - 1, 
     $      tH, 1 + IROFFH + 1, 1 + IROFFH, tDESCH,
     $      H, KWTOP + 1, KWTOP, DESCH)
*        
*        First row of H(KWTOP, KWTOP)
         CALL PDLACPY('A', 1, JW, 
     $      tH, 1 + IROFFH, 1 + IROFFH, tDESCH, 
     $      H, KWTOP, KWTOP, DESCH)
*        
*        Last column of H(KWTOP, KWTOP)
         CALL PDLACPY('A', JW - 1, 1, 
     $      tH, 1 + IROFFH + 1, 1 + IROFFH + JW - 1, tDESCH, 
     $      H, KWTOP + 1, KWTOP + JW - 1, DESCH )
         
*        
*        T
         CALL PDLACPY('U', JW, JW, 
     $      tT, 1 + IROFFH, 1 + IROFFH, tDESCT, 
     $      T, KWTOP, KWTOP, DESCT )

*        
*        Apply top element of Q to H
         IF ( KWTOP .GT. 1 ) THEN
            CALL PDELGET( 'A', ' ', ELEM, 
     $         tQ, 1 + IROFFH, 1 + IROFFH, tDESCQ )   

            CALL PDELSET( H, KWTOP, KWTOP - 1, DESCH,  SCAL * ELEM )
         END IF 
         T2 = MPI_WTIME()  
         UPDTIME(8) = UPDTIME(8) + (T2-T1)         
         
         CALL BLACS_BARRIER( ICTXT, 'A' )

         IF (DEBUG) THEN
            WRITE(*,*)'% PQZ3 : Apply off-diagonal updates'
            CALL FLUSH(6)
         END IF

         LTOP = 1         
         T1 = MPI_WTIME()  
         
         IF ( .NOT. WANTHT ) LTOP = KTOP
         
         KLN = MAX( 0, KWTOP - LTOP )
         OFF_R = MOD( LTOP - 1, NB)
         OFF_C = MOD( KWTOP -1, NB)
         RSRC = INDXG2P( LTOP, NB, MYROW, DESCH( RSRC_ ), NPROW )
         CSRC = INDXG2P( KWTOP, NB, MYCOL, DESCH( CSRC_ ), NPCOL )
         LR = NUMROC( KLN + OFF_R, NB, MYROW, RSRC, NPROW )
         LC = NUMROC( JW + OFF_C, NB, MYCOL, CSRC, NPCOL )

         CALL BLACS_BARRIER( ICTXT, 'A' )
         IF (DEBUG) THEN
            WRITE(*,*)'% PQZ3 : tZ on H,T Update info', JW, KTOP, LTOP, 
     $         KWTOP, KLN, OFF_R, OFF_C, KLN*JW
            CALL FLUSH(6)
         END IF
*        Apply tZ from right on H and T
*        
         CALL DESCINIT( DESC_TMP , KLN + OFF_R, JW + OFF_C, NB, NB,
     $      RSRC, CSRC, ICTXT, MAX(1, LR), IERR )
         CALL PDGEMM( 'N', 'N', KLN, JW, JW, ONE, H, LTOP,
     $      KWTOP, DESCH, tZ, 1 + IROFFH, 1 + IROFFH, tDESCZ, ZERO,
     $      DWORK, 1 + OFF_R, 1 + OFF_C, DESC_TMP )
         CALL PDLACPY( 'A', KLN, JW, DWORK, 1 + OFF_R, 1 + OFF_C,
     $      DESC_TMP, H, LTOP, KWTOP, DESCH )
         
         CALL PDGEMM( 'N', 'N', KLN, JW, JW, ONE, T, LTOP,
     $      KWTOP, DESCT, tZ, 1 + IROFFH, 1 + IROFFH, tDESCZ, ZERO,
     $      DWORK, 1 + OFF_R, 1 + OFF_C, DESC_TMP )
         CALL PDLACPY( 'A', KLN, JW, DWORK, 1 + OFF_R, 1 + OFF_C,
     $      DESC_TMP, T, LTOP, KWTOP, DESCT )
         

*        Apply tQ transp from left on H and T
*        
         KLN = N - KBOT
         OFF_R = MOD( KWTOP - 1, NB )
         OFF_C = MOD( KBOT, NB )

         CALL BLACS_BARRIER( ICTXT, 'A' )
         IF (DEBUG) THEN
            WRITE(*,*)'% PQZ3 : tQ on H,T Update info', JW, KTOP, LTOP, 
     $         KWTOP, KLN, OFF_R, OFF_C, KLN*JW
            CALL FLUSH(6)
         END IF

         RSRC = INDXG2P( KWTOP, NB, MYROW, DESCH( RSRC_ ), NPROW )
         CSRC = INDXG2P( KBOT + 1, NB, MYCOL, DESCH( CSRC_ ), NPCOL )
         LR = NUMROC( JW + OFF_R, NB, MYROW, RSRC, NPROW )
         LC = NUMROC( KLN + OFF_C, NB, MYCOL, CSRC, NPCOL )
         CALL DESCINIT( DESC_TMP, JW + OFF_R, KLN + OFF_C, NB, NB,
     $      RSRC, CSRC, ICTXT, MAX(1, LR), IERR )
         CALL PDGEMM( 'T', 'N', JW, KLN, JW, ONE, tQ,
     $      1 + IROFFH, 1 + IROFFH, tDESCQ, H, KWTOP, KBOT + 1,
     $      DESCH, ZERO, DWORK, 1 + OFF_R, 1 + OFF_C, DESC_TMP )
         CALL PDLACPY( 'A', JW, KLN, DWORK, 1 + OFF_R, 1 + OFF_C,
     $      DESC_TMP, H, KWTOP, KBOT + 1, DESCH )
         
         CALL PDGEMM( 'T', 'N', JW, KLN, JW, ONE, tQ,
     $      1 + IROFFH, 1 + IROFFH, tDESCQ, T, KWTOP, KBOT + 1,
     $      DESCT, ZERO, DWORK, 1 + OFF_R, 1 + OFF_C, DESC_TMP )
         CALL PDLACPY( 'A', JW, KLN, DWORK, 1 + OFF_R, 1 + OFF_C,
     $      DESC_TMP, T, KWTOP, KBOT + 1, DESCT )

*        Apply tZ from the right on Z
*        
         IF ( ILZ ) THEN
            KLN = IHIZ - ILOZ + 1
            OFF_R = MOD( ILOZ - 1, NB )
            OFF_C = MOD( KWTOP - 1, NB )

            CALL BLACS_BARRIER( ICTXT, 'A' )
            IF (DEBUG) THEN
               WRITE(*,*)'% PQZ3 : tZ on Z Update info', JW, KTOP, LTOP, 
     $            KWTOP, KLN, OFF_R, OFF_C, KLN*JW
               CALL FLUSH(6)
            END IF
            RSRC = INDXG2P( ILOZ, NB, MYROW, DESCZ( RSRC_), NPROW )
            CSRC = INDXG2P( KWTOP, NB, MYCOL, DESCZ( CSRC_), NPCOL )
            LR = NUMROC( KLN + OFF_R, NB, MYROW, RSRC, NPROW )
            LC = NUMROC( JW + OFF_C, NB, MYCOL, CSRC, NPCOL )
            CALL DESCINIT( DESC_TMP, KLN + OFF_R, JW + OFF_C, NB, NB,
     $         RSRC, CSRC, ICTXT, MAX( 1, LR ), IERR )
            CALL PDGEMM( 'N', 'N', KLN, JW, JW, ONE, Z, ILOZ,
     $         KWTOP, DESCZ, tZ, 1 + IROFFH, 1 + IROFFH, tDESCZ,
     $         ZERO, DWORK, 1 + OFF_R, 1 + OFF_C, DESC_TMP )
            CALL PDLACPY( 'A', KLN, JW, DWORK, 1 + OFF_R, 1 + OFF_C,
     $         DESC_TMP, Z, ILOZ, KWTOP, DESCZ )                  
         END IF

*        Apply tQ from the right on Q
*        
         IF ( ILQ ) THEN
            KLN = IHIQ - ILOQ + 1
            OFF_R = MOD( ILOQ - 1, NB )
            OFF_C = MOD( KWTOP - 1, NB )

            CALL BLACS_BARRIER( ICTXT, 'A' )
            IF (DEBUG) THEN
               WRITE(*,*)'% PQZ3 : tQ on Q Update info', JW, KTOP, LTOP, 
     $            KWTOP, KLN, OFF_R, OFF_C, KLN*JW
               CALL FLUSH(6)
            END IF
            RSRC = INDXG2P( ILOQ, NB, MYROW, DESCQ( RSRC_ ), NPROW )
            CSRC = INDXG2P( KWTOP, NB, MYCOL, DESCQ( CSRC_ ), NPCOL )
            LR = NUMROC( KLN + OFF_R, NB, MYROW, RSRC, NPROW )
            LC = NUMROC( JW + OFF_C, NB, MYCOL, CSRC, NPCOL )
            CALL DESCINIT( DESC_TMP, KLN + OFF_R, JW + OFF_C, NB, NB,
     $         RSRC, CSRC, ICTXT, MAX( 1, LR ), IERR )
            CALL PDGEMM( 'N', 'N', KLN, JW, JW, ONE, Q, ILOQ,
     $         KWTOP, DESCQ, tQ, 1 + IROFFH, 1 + IROFFH, tDESCQ,
     $         ZERO, DWORK, 1 + OFF_R, 1 + OFF_C, DESC_TMP )
            CALL PDLACPY( 'A', KLN, JW, DWORK, 1 + OFF_R, 1 + OFF_C,
     $         DESC_TMP, Q, ILOQ, KWTOP, DESCQ )                  
         END IF
         T2 = MPI_WTIME()  
         UPDTIME(9) = UPDTIME(9) + (T2-T1)
      END IF      
      ND = JW - NS
      IF (DEBUG) THEN
         WRITE(*,*)'% PQZ3 : Done. Deflation=', ND,
     $      INFO, IERR
         CALL FLUSH(6)
      END IF


 9999 CONTINUE
      RETURN
      END

