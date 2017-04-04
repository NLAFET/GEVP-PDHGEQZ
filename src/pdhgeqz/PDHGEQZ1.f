***********************************************************************
*                                                                     *
*     PDHGEQZ1.f:                                                     *
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
      SUBROUTINE PDHGEQZ1(ILSCHR, ILQ, ILZ, N, ILO, IHI, H, DESCH, T,
     $   DESCT, ALPHAR, ALPHAI, BETA, Q, DESCQ, Z, DESCZ, DWORK, LDWORK,
     $   IWORK, LIWORK, INFO )
      IMPLICIT NONE
*     ..
*     .. Scalar Arguments ..
*     ..
      INTEGER             IHI, ILO, INFO, LDWORK, N, LIWORK
      LOGICAL             ILQ, ILSCHR, ILZ
*     ..
*     .. Array Arguments .. 
*     ..
      DOUBLE PRECISION    T( * ), H( * ), Q( * ), Z( * ), DWORK( * ),
     $   BETA( * ), ALPHAI( * ), ALPHAR( * )
      INTEGER             DESCH( * ), DESCT( * ), DESCQ( * ), DESCZ( * )
      INTEGER             IWORK( * )


*     Purpose
*     =======
*     
*     PDHGEQZ1 : implements a parallel single-/double-shift version of the QZ method for
*     finding the generalized eigenvalues - suited for smaller problems.
*     
*     w(j)=(ALPHAR(j) + i*ALPHAI(j))/BETAR(j)   of the equation
*     
*     det( H - w(i) T ) = 0
*     
*     In addition, the pair H,T may be reduced to generalized Schur form:
*     T is upper triangular, and H is block upper triangular, where the
*     diagonal blocks are either 1-by-1 or 2-by-2, the 2-by-2 blocks having
*     complex generalized eigenvalues (see the description of the argument
*     ILSHUR.)
*     
*     If ILSCHR=TRUE, then the pair (H,T) is simultaneously reduced to Schur
*     form by applying one orthogonal tranformation (usually called Q) on
*     the left and another (usually called Z) on the right.  The 2-by-2
*     upper-triangular diagonal blocks of T corresponding to 2-by-2 blocks
*     of H will be reduced to positive diagonal matrices.  (I.e.,
*     if H(j+1,j) is non-zero, then T(j+1,j)=T(j,j+1)=0 and T(j,j) and
*     T(j+1,j+1) will be positive.)
*     
*     If ILSCHR=FALSE, then at each iteration, the same transformations
*     are computed, but they are only applied to those parts of H and T
*     which are needed to compute ALPHAR, ALPHAI, and BETAR.
*     
*     If ILSCHR=ILQ=ILZ=TRUE, then the orthogonal
*     transformations used to reduce (A,B) are accumulated into the arrays
*     Q and Z s.t.:
*     
*     Q(in) H(in) Z(in)* = Q(out) H(out) Z(out)*
*     Q(in) T(in) Z(in)* = Q(out) T(out) Z(out)*
*     
*     Ref: C.B. Moler & G.W. Stewart, "An Algorithm for Generalized Matrix
*     Eigenvalue Problems", SIAM J. Numer. Anal., 10(1973),
*     pp. 241--256.
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
*     
*     Arguments 
*     ========= 
*     
*     ILSCHR  (global input) LOGICAL
*     = TRUE: Transform (H,T) to generalized Schur form
*     = FALSE: H,T are undefined on output
*     
*     ILQ      (input) LOGICAL
*     = TRUE: Compute Q.
*     = FALSE: Q is not referenced.
*     
*     ILZ     (global input) LOGICAL
*     = TRUE: Compute Z.
*     = FALSE: Z is not referenced.
*     
*     N       (global input) INTEGER 
*     The order of the matrices H,T,Q, and Z.  N >= 0. 
*     
*     ILO     (global input) INTEGER 
*     IHI     (global input) INTEGER 
*     It is assumed that A is already upper triangular in rows and
*     columns 1:ILO-1 and IHI+1:N.  ILO and IHI are normally set
*     by a previous call to DGGBAL; otherwise they should be set
*     to 1 and N respectively. 
*     1 <= ILO <= IHI <= N, if N > 0; ILO=1 and IHI=0, if N=0.
*     
*     H       (local input/output) DOUBLE PRECISION array, dimension (LLD_H, LOCc(N)).
*     On entry, H contains the upper Hessenberg matrix H.
*     On exit, if ILSCHR=TRUE, H is quasi-upper triangular in rows and columns 
*     (ILO:IHI), with 1 x 1 and 2 x 2 blocks on the diagonal where 
*     the 2 x 2 blocks corresponds to complex conjugated pairs of 
*     eigenvalues. If ILSCHR=FALSE, H is unspecified on exit.      
*     
*     DESCH   (global and local input) INTEGER array of dimension DLEN_.
*     The array descriptor for the distributed matrix H. 
*     
*     T       (local input/output) DOUBLE PRECISION array, dimension (LLD_T, LOCc(N)).
*     On entry, the N-by-N upper triangular matrix T. 
*     On exit, if ILSHUR the upper triangular matrix the updated tringular 
*     matrix T.
*     
*     DESCT   (global and local input) INTEGER array of dimension DLEN_.
*     The array descriptor for the distributed matrix T. 
*     
*     ALPHAR  (global output) DOUBLE PRECISION array, dimension (N)
*     ALPHAR(1:N) will be set to real parts of the diagonal 
*     elements of H that would result from reducing H and T to
*     Schur form and then further reducing them both to triangular
*     form using unitary transformations s.t. the diagonal of T
*     was non-negative real.  Thus, if H(j,j) is in a 1-by-1 block
*     (i.e., H(j+1,j)=H(j,j+1)=0), then ALPHAR(j)=H(j,j). 
*     Note that the (real or complex) values 
*     (ALPHAR(j) + i*ALPHAI(j))/BETA(j), j=1,...,N, are the 
*     generalized eigenvalues of the matrix pencil H - wT. 
*     
*     ALPHAI  (global output) DOUBLE PRECISION array, dimension (N)
*     ALPHAI(1:N) will be set to imaginary parts of the diagonal
*     elements of H that would result from reducing H and T to
*     Schur form and then further reducing them both to triangular
*     form using unitary transformations s.t. the diagonal of T
*     was non-negative real.  Thus, if H(j,j) is in a 1-by-1 block
*     (i.e., H(j+1,j)=J(j,j+1)=0), then ALPHAI(j)=0. 
*     Note that the (real or complex) values 
*     (ALPHAR(j) + i*ALPHAI(j))/BETA(j), j=1,...,N, are the 
*     generalized eigenvalues of the matrix pencil H - wT. 
*     
*     BETA    (global output) DOUBLE PRECISION array, dimension (N)
*     BETA(1:N) will be set to the (real) diagonal elements of T
*     that would result from reducing H and T to Schur form and
*     then further reducing them both to triangular form using
*     unitary transformations s.t. the diagonal of T was 
*     non-negative real.  Thus, if H(j,j) is in a 1-by-1 block
*     (i.e., H(j+1,j)=H(j,j+1)=0), then BETA(j)=T(j,j). 
*     Note that the (real or complex) values 
*     (ALPHAR(j) + i*ALPHAI(j))/BETA(j), j=1,...,N, are the 
*     generalized eigenvalues of the matrix pencil H - wT. 
*     (Note that BETA(1:N) will always be non-negative, and no
*     BETAI is necessary.)
*     
*     Q       (local input/output) DOUBLE PRECISION array, dimension (LLD_Q, LOCc(N)).
*     If ILQ = FALSE:  Q is not referenced. 
*     If ILQ = TRUE:  on entry, Q must contain an orthogonal matrix.
*     On exit, Q is multiplied on the right with 
*     the product of givens transformations applied on H,T from the left.
*     
*     DESCQ   (global and local input) INTEGER array of dimension DLEN_.
*     The array descriptor for the distributed matrix Q. 
*     
*     Z       (local input/output) DOUBLE PRECISION array, dimension (LLD_Z, LOCc(N)).
*     If ILZ = FALSE:  Z is not referenced. 
*     If ILZ = TRUE:  on entry, Z must contain an orthogonal matrix.
*     On exit, Z is multiplied on the right with 
*     the product of givens transformations applied on H,T from the right. 
*     
*     DESCZ   (global and local input) INTEGER array of dimension DLEN_.
*     The array descriptor for the distributed matrix Z. 
*     
*     DWORK   (local workspace/global output) DOUBLE PRECISION array,
*     dimension (LDWORK) 
*     
*     LDWORK  (global input) INTEGER 
*     The dimension of the array DWORK. 
*     If LDWORK = -1, then a workspace query is assmed, and optimal
*     workspace is returned in DWORK(1)
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
*     < 0:  if INFO = -i, the i-th argument had an illegal value.
*     = 1..5 one subroutine failed
*     =====================================================================



*     ..
*     .. Parameters ..
*     ..
      DOUBLE PRECISION    ZERO, ONE, TWO
      PARAMETER           ( ZERO = 0.0D+0, ONE = 1.0D+0, TWO = 1.0D+0)
      DOUBLE PRECISION    CONST
      PARAMETER           ( CONST = 1.50D+0 )
      
      INTEGER             BLOCK_CYCLIC_2D, CSRC_, CTXT_, DLEN_, DT_, 
     $   LLD_, MB_, M_, NB_, N_, RSRC_
      PARAMETER           ( BLOCK_CYCLIC_2D = 1, DLEN_ = 9, DT_ = 1,
     $   CTXT_ = 2, M_ = 3, N_ = 4, MB_ = 5, NB_ = 6,
     $   RSRC_ = 7, CSRC_ = 8, LLD_ = 9 )



*     ..
*     .. Precision paramters ..
*     ..
      DOUBLE PRECISION    SMLNUM, ULP, TFNORM, T1NORM, HFNORM, HTOL,
     $   TTOL
      COMMON /PREC/       SMLNUM, ULP, TFNORM, T1NORM, HFNORM, HTOL, 
     $   TTOL
      
*     ..
*     .. Local Scalars ..
*     ..
      DOUBLE PRECISION    SWAP, CSL, SNL, CSR, SNR, T1
      INTEGER             IFIRST, ILAST, NIBBLE, I, IERR, MAXIT, ICTXT,
     $   NPROW, NPCOL, MYROW, MYCOL, IAM, NB, NUMWIN, IPWORK, NRPROCS,
     $   KBOT, KTOP, NS, ND, NU, NH, KS,OLDKTOP, DBLK,
     $   LDWORKMIN, LIWORKMIN, SHIFTSPERWND, OLDILAST, OLDIFIRST,
     $   NDINF, TOTNDINF, USEDSHIFTS, IDUM, NMIN3, ITER0,
     $   ITER2, ITER1, ITER3, ILOQ, IHIQ, ILOZ, IHIZ
      LOGICAL             DEBUG, LQUERY

*     ..
*     .. Local Arrays ..
*     ..    
      DOUBLE PRECISION, ALLOCATABLE :: TMP1( :, : ), TMP2( :, : ), 
     $   TMP3( :, : ), TMP4( :, : ), TAR( : ), TAI( : ), TBETA( : )
      DOUBLE PRECISION    UPDTIME(40)

*     ..
*     .. Externals ..
*     ..
      EXTERNAL            ICEIL, NUMROC, INDXG2L, INDXG2P, ILCM,
     $   PDLAMCH, PDLANHS, PDLANTR, PILAENVX, MPI_WTIME
      INTEGER             ICEIL, NUMROC, INDXG2L, INDXG2P, ILCM,
     $   PILAENVX
      DOUBLE PRECISION    PDLAMCH, PDLANHS, PDLANTR, MPI_WTIME
*     ..
*     .. Intrinsic Functions ..
*     ..
      INTRINSIC           ABS, DBLE, MAX, MIN, SQRT, MOD, LOG, NINT

      
      
      
*     ..
*     .. Executable Statements ..
*     ..
      INFO = 0

      IF( N.LE.0 .OR. IHI .LT. ILO ) THEN
         RETURN
      END IF


*     Check for workspace query
      LQUERY = LDWORK .EQ. -1 .OR. LIWORK .EQ. -1


*     Extract current communication context and get grid parameters
      ICTXT = DESCH(CTXT_)
      CALL BLACS_PINFO(IAM, NRPROCS)
      CALL BLACS_GRIDINFO( ICTXT, NPROW, NPCOL, MYROW, MYCOL )

*     Set DEBUG to true/false for various information printed during the run
      DEBUG = .FALSE.
*     Debugprint for root only
      DEBUG = DEBUG .AND. IAM .EQ. 0

*     Define current problem size
      NH = IHI - ILO + 1
      ILOQ = 1
      ILOZ = 1
      IHIQ = N
      IHIZ = N

      IF ( DEBUG ) THEN
         IF (LQUERY) THEN
            WRITE(*,*)'% Entering PQZ1 for WSPACE query', NH, ILO,
     $         IHI
         ELSE
            WRITE(*,*)'% Entering PQZ1', NH, ILO, IHI            
         END IF
         CALL FLUSH( 6 )
      END IF
      
*     NMIN3 is Used for defining the maximum problem size to be solved serially
      NMIN3 = PILAENVX( ICTXT, 55, 'PDHGEQZ1', '', IDUM, IDUM, IDUM,
     $   IDUM)
      IF ( DEBUG ) THEN 
         WRITE(*,*)'% PQZ1 : NMIN3=', NMIN3
         CALL FLUSH( 6 ) 
      END IF

*     NIBBLE for when to skip QZ sweep 
      NIBBLE = PILAENVX( ICTXT, 52, 'PDHGEQZ1', '', IDUM, IDUM, IDUM,
     $   IDUM)
      IF ( DEBUG) THEN 
         WRITE(*,*)'% PQZ1 : NIBBLE', NIBBLE
         CALL FLUSH( 6 ) 
      END IF

*     NB_ .EQ. MB_ is assumed, we use NB as blocking factor.
      NB = DESCH( NB_ )

*     Calculate double precision workspace
      LDWORKMIN = N * NB

*     Integer workspace
      LIWORKMIN = 2 * N

      IF ( LQUERY ) THEN
         DWORK( 1 ) = LDWORKMIN
         IWORK( 1 ) = LIWORKMIN
         RETURN
      END IF
      IF ( LDWORK .LT. LDWORKMIN ) THEN        
         INFO = -1
         RETURN
      ELSE IF ( LIWORK .LT. LIWORKMIN ) THEN        
         INFO = -2
         RETURN
      END IF

      
      IPWORK = 1
      
*     Allocated temporary storage, used for AED and solving problems serially 
*     when remaining problem is small enough
      ALLOCATE( TMP1( NMIN3, NMIN3 ) )
      ALLOCATE( TMP2( NMIN3, NMIN3 ) )
      ALLOCATE( TMP3( NMIN3, NMIN3 ) )
      ALLOCATE( TMP4( NMIN3, NMIN3 ) )
      ALLOCATE( TAR( NMIN3 ) )
      ALLOCATE( TAI( NMIN3 ) )
      ALLOCATE( TBETA( NMIN3 ) )
      CALL BLACS_BARRIER(ICTXT, 'A')
      IF ( DEBUG ) THEN 
         WRITE(*,*)'% PQZ1 : Alloc done', NMIN3*3 + NMIN3*NMIN3*4
         CALL FLUSH( 6 ) 
      END IF
      ILAST = IHI
      IFIRST = ILO
      
      OLDILAST = ILAST
      OLDIFIRST = IFIRST


*     
*     Counters
*     
      USEDSHIFTS = 0
      TOTNDINF = 0 
*     Iter1 holds number of AED calls
      ITER1 = 0
*     ITER2 holds number of QZ sweeps
      ITER2 = 0
*     ITER3 holds number of small problems solved
      ITER3 = 0
*     ITER0 holds total number of iterations, see 360 below


*     If the diagonal block is small enough, and is not 
*     shared by more than 2 diagonal processors, 
*     copy it to local memory and
*     solve the problem serially.
      IF( NH .LT. NMIN3 .AND. NH .GT. 2 .AND. MOD( IFIRST - 1, NB ) + NH
     $   .LE. NB * 2 ) THEN   
         
         IF ( DEBUG ) THEN
            WRITE(*,*)'% PQZ1 : PQZ4', NH, ILAST, IFIRST, ILO, IHI,
     $         NMIN3
            CALL FLUSH(6)
         END IF
         CALL PDHGEQZ4( ILSCHR, ILQ, ILZ, H, DESCH, T, DESCT, Q, DESCQ,
     $      Z, DESCZ, N, IFIRST, ILAST, ILO, IHI, ALPHAI, ALPHAR, BETA,
     $      TMP1, TMP2, ILOQ, IHIQ, ILOZ, IHIZ, TMP3, TMP4, TAR, TAI,
     $      TBETA, DWORK, LDWORK, IWORK, LIWORK, IERR )
         IF ( IERR .NE. 0 ) THEN
            IF ( IAM .EQ. 0 ) THEN
               WRITE(*,*)'% PQZ1 : PDHGEQZ4 failed. INFO=', IERR
               CALL FLUSH( 6 )
            END IF
            INFO = 1
         END IF            
         GOTO 420
      END IF
      
      
*     The current problem is defined by KTOP:KBOT, updated during the 360 loop.
      KBOT = ILAST
      KTOP = IFIRST
      MAXIT = 30*( KBOT - KTOP + 1 )

      
      
      IF ( DEBUG ) THEN
         WRITE(*,*)'% PQZ1 init', IFIRST, ILAST
         CALL FLUSH(6)
      END IF


*     MAXIT = 1
      DO 360 ITER0 = 1, MAXIT
         IF ( KBOT .LT. IFIRST ) GOTO 380         
*        Update current problem size
         NH = KBOT - IFIRST + 1
*        If we fit the whole remaning problem, skip to search for deflations
         IF ( NH .LE. NMIN3 ) THEN
            KTOP = IFIRST
         ELSE          
*           Search for deflations, bottom to top.
            CALL PDLASMSUB( H, DESCH, KBOT, IFIRST, KTOP, SMLNUM, DWORK,
     $         LDWORK )
         END IF
         IF( KTOP .GT. IFIRST ) THEN
*           A( KTOP, KTOP - 1) is negligible
            CALL PDELSET( H, KTOP, KTOP - 1, DESCH, ZERO )
         END IF

*        Search for and deflate infinte eigenvalues, NDINF will hold the amount.
         NDINF = 0
         OLDKTOP = KTOP
         IF (DEBUG) WRITE(*,*)'% PQZ1 : PQZ7', KTOP, KBOT
         CALL PDHGEQZ7( ILSCHR, ILQ, ILZ, KTOP, KBOT, ILO, IHI, H, DESCH
     $      , T, DESCT, ALPHAR, ALPHAI, BETA, Q, DESCQ, Z, DESCZ, ILOQ,
     $      IHIQ, ILOZ, IHIZ, MIN(NPROW, NPCOL, 40), NDINF, DWORK,
     $      LDWORK, IWORK, IERR )               

         IF ( IERR .NE. 0 ) THEN
            IF ( IAM .EQ. 0 ) THEN
               WRITE(*,*)'% PQZ1 : PDHGEQZ7 failed. INFO=',
     $            IERR
               CALL FLUSH( 6 )
            END IF
            INFO = 2
            GOTO 420
         END IF 

*        If deflation were made at top, adjust where to start next iteration.
         IF (KTOP.NE.OLDKTOP) THEN           
            IF (OLDKTOP.EQ.IFIRST) THEN
               IFIRST = KTOP
            END IF
         END IF

*        Try another standard deflation if inf. eigvs have been found.
         IF ( NDINF .GT. 0 ) THEN
            TOTNDINF = TOTNDINF + NDINF
            IF ( DEBUG ) THEN
               WRITE(*,*)'% PQZ1 : Intermediate deflation ', NDINF, ' oo
     $            eigenvalues.', KBOT, KTOP
               CALL FLUSH(6)
               
            ENDIF          
            
            IF (KTOP.GE.KBOT) THEN
               IF (KBOT.EQ.KTOP) THEN
                  CALL PDELGET('A', ' ', ALPHAR( KBOT ), H, KBOT, KBOT,
     $               DESCH )
                  CALL PDELGET('A', ' ', BETA( KBOT ), T, KBOT, KBOT,
     $               DESCT )
                  ALPHAI( KBOT ) = ZERO
                  
               END IF
               KBOT = OLDKTOP - 1
               KTOP = IFIRST
            END IF            
            GOTO 350
         END IF
         
*        Now work on KTOP:KBOT.
         NH = KBOT - KTOP + 1
         
*        Set DBLK to the AED window size.
         DBLK = NB         
         DBLK = MIN( NH, NMIN3, DBLK )
*        Exit from loop if a small submatrix has split off.
         IF ( NH .LT. NMIN3 ) THEN
            IF ( ( MOD( KTOP - 1, NB ) + NH ) . LE. ( NB * 2 ) ) THEN
               GO TO 300
            END IF
         END IF          

*        ---------------------------
*        Time for AED
*        ---------------------------

*        Set number of deflated elements initially to none.
         ND = 0
         INFO = 0
         ITER1 = ITER1 + 1
         IF ( DEBUG ) WRITE(*,*)'% PQZ1 : AED', KBOT, DBLK
         T1 = MPI_WTIME()
         CALL PDHGEQZ2( ILSCHR, ILQ, ILZ, ILO, N, KBOT, DBLK, H, DESCH,
     $      T, DESCT, Q, DESCQ, Z, DESCZ, ALPHAR, ALPHAI, BETA, ILOQ,
     $      IHIQ, ILOZ, IHIZ, NU, ND, TMP1, TMP2, TMP3, TMP4, TAR, TAI,
     $      TBETA, DWORK, LDWORK, IERR )
         IF ( IERR .NE. 0 ) THEN
            IF ( IAM .EQ. 0 ) THEN
               WRITE(*,*)'% PQZ1 : PDHGEQZ2 failed. INFO=',
     $            IERR
               CALL FLUSH( 6 )
            END IF
            INFO = 3
            GOTO 420
         END IF 
         
         IF ( DEBUG ) THEN
            WRITE(*,'(A, I6, I6, I6, I6, I6)') '% PQZ1 : IIter=', 
     $         ITER0, KTOP,
     $         KBOT, DBLK, ND 
            CALL FLUSH(6)
         END IF
         UPDTIME(1) = UPDTIME(1) + (MPI_WTIME() - T1)
         
*        Adjust problem size for made deflations.
         KBOT = KBOT - ND
         IF ( KBOT .LT. IFIRST ) GOTO 380                           
         IF ( KBOT .LT. KTOP ) GOTO 350
         
*        Skip QZ sweep if number of converged eigenvalues are high enough.
         IF( 100 * ND .GT. NIBBLE * DBLK )  THEN
            IF (DEBUG) THEN
               WRITE(*,'(A, I6, I6)')'% PQZ1 : Retry', 100 * ND, 
     $            NIBBLE * DBLK
               CALL FLUSH(6)
            END IF
         ELSE
            
*           ---------------------------
*           Time for a QZ sweep
*           --------------------------
*           
            ITER2 = ITER2 + 1
*           Set NS to the number of computed eigenvalues.
            NS = DBLK
            
            
            
*           Shuffle shifts into pairs of real shifts
*           and pairs of complex conjugate shifts
*           assuming complex conjugate shifts are
*           already adjacent to one another. (Yes,
*           they are.)           
            DO I = NS, 3, -2               
               IF( TAI( I ) .EQ. 0 .AND. TAI( I - 1 ) .NE. 0 ) THEN
                  SWAP = TAR(I)
                  TAR(I) = TAR(I-1)
                  TAR(I-1) = TAR(I-2)
                  TAR(I-2) = SWAP
                  
                  SWAP = TAI(I)
                  TAI(I) = TAI(I-1)
                  TAI(I-1) = TAI(I-2)
                  TAI(I-2) = SWAP
                  
                  SWAP = TBETA(I)
                  TBETA(I) = TBETA(I-1)
                  TBETA(I-1) = TBETA(I-2)
                  TBETA(I-2) = SWAP                    
               END IF
            END DO
            
*           Use only undeflated eigenvalues as shift, make sure even number.    
            NS = NS - ND
            NS = NS - MOD(NS, 2)

*           Use 3/4 of the provided shifts, from bottom. 
            KS = NS / 4
            KS = KS - ( 1 - MOD( KS, 2 ) )
*           Depening on the amount of available shifts and current grid size, 
*           set number for shifts to use per window.
            SHIFTSPERWND = MAX( MIN( NS - KS, NB / 3 - 1 ), 2 )
            SHIFTSPERWND = SHIFTSPERWND - MOD( SHIFTSPERWND, 2 )
            NUMWIN = MAX(MIN(40, NPROW, NPCOL, ( ( NS - KS ) /
     $         SHIFTSPERWND ) ), 1 )

            IF ( DEBUG ) THEN
               WRITE(*,'(A, I4, I4, I4, I4, I6, I6)') '% PQZ1 : Sweep,
     $            Using shifts/wnd=',SHIFTSPERWND, NUMWIN, DBLK, NS, KS,
     $NS-KS
               CALL FLUSH(6)
            END IF

*           Introduce and chase the shifts in H and T at IFIRST.
            CALL PDHGEQZ5( ILSCHR, ILQ, ILZ, NS, KTOP, KBOT, N, ILO, IHI
     $         , H, DESCH, T, DESCT, TAR, TAI, TBETA, Q, DESCQ, Z, DESCZ
     $         , ILOQ, IHIQ, ILOZ, IHIZ, NUMWIN, USEDSHIFTS, HTOL, DWORK
     $         , LDWORK, IERR, UPDTIME, SHIFTSPERWND )
            
            IF ( IERR .NE. 0 ) THEN
               IF ( IAM .EQ. 0 ) THEN
                  WRITE(*,*)'% PQZ1 : PDHGEQZ5 failed. INFO=', IERR
                  CALL FLUSH( 6 )
               END IF
               INFO = 4
               GOTO 420
            END IF 

         END IF
         GOTO 350
 300     CONTINUE
*        Small sub matrix split
*        We have three cases: The sub matrix is 1 element, 2 elements, or larger.

         ITER3 = ITER3 + 1

*        Find the eigenvalues in A(IFIRST:ILAST,IFIRST:ILAST), IFIRST < ILAST-1
         IF ( DEBUG ) THEN
            WRITE(*,9998)' % PQZ1 : Small Sub, Jiter=',
     $         ITER0,KBOT-KTOP+1, 0
            CALL FLUSH(6)
         END IF

         IF( KTOP .EQ. KBOT ) THEN
*           A(ILAST,ILAST-1) is negligible: one (real) eigenvalue has converged.
            CALL PDELGET('A', ' ', ALPHAR( KBOT ), H, KBOT, KBOT, 
     $         DESCH )
            CALL PDELGET('A', ' ', BETA( KBOT ), T, KBOT, KBOT, 
     $         DESCT )
            
            ALPHAI( KBOT ) = ZERO
         ELSE IF( KTOP .EQ. KBOT - 1 ) THEN
*           A(ILAST-1,ILAST-2) is negligible: a pair of eigenvalues have converged.
            CALL PDLACP3( 2, KBOT - 1, H, DESCH, TMP1, NMIN3, -1, -1, 
     $         0 )
            CALL PDLACP3( 2, KBOT - 1, T, DESCT, TMP2, NMIN3, -1, -1, 
     $         0 )
            
            CALL DLAGV2( TMP1, NMIN3, TMP2, NMIN3, ALPHAR( KBOT - 1 ),
     $         ALPHAI( KBOT - 1 ), BETA( KBOT -1 ), CSL, SNL, CSR, SNR )

            CALL PDLACP3( 2, KBOT - 1, H, DESCH, TMP1, NMIN3, 0, 0, -1 )
            CALL PDLACP3( 2, KBOT - 1, T, DESCT, TMP2, NMIN3, 0, 0, -1 )
*           Update cols of A and B              
            CALL PDROT( KBOT - 2, H, 1, KBOT-1, DESCH, 1, H, 1, KBOT,
     $         DESCH, 1, CSR, SNR, DWORK, LDWORK, INFO )
            
            CALL PDROT( KBOT - 2, T, 1, KBOT-1, DESCT, 1, T, 1, KBOT,
     $         DESCT, 1, CSR, SNR, DWORK, LDWORK, INFO )

*           Update rows of A and B
            IF ( IHI .GT. KBOT ) THEN
               CALL PDROT( IHI - KBOT, H, KBOT - 1, KBOT + 1, DESCH, N,
     $            H, KBOT, KBOT + 1, DESCH, N, CSL, SNL, DWORK, LDWORK,
     $            INFO)
               
               CALL PDROT( IHI - KBOT, T, KBOT - 1, KBOT + 1, DESCT, N,
     $            T, KBOT, KBOT + 1, DESCT, N, CSL, SNL, DWORK, LDWORK,
     $            INFO )
            END IF

*           Update cols of Q
            IF (ILQ) THEN
               CALL PDROT(N, Q, 1, KBOT - 1, DESCQ, 1, Q, 1, KBOT,
     $            DESCQ, 1, CSL, SNL, DWORK, LDWORK, INFO)
            END IF
*           Update cols of Z
            IF (ILZ) THEN
               CALL PDROT(N, Z, 1, KBOT - 1, DESCZ, 1, Z, 1, KBOT, 
     $            DESCZ, 1, CSR, SNR, DWORK, LDWORK, INFO )
            END IF

         ELSE     
*           A small submatrix has been identied. Solve the problem serially and 
*           perform updates in parallel
            NH = KBOT - KTOP + 1      
            IF( NH .LT. NMIN3 ) THEN
               IF ( DEBUG ) THEN
                  WRITE(*,9999)'% Calling QZ4', KTOP, KBOT, NMIN3, 
     $               ILO, IHI
                  CALL FLUSH(6)                  
               END IF

               CALL PDHGEQZ4( ILSCHR, ILQ, ILZ, H, DESCH, T, DESCT, Q,
     $            DESCQ, Z, DESCZ, N, KTOP, KBOT, ILO, IHI, ALPHAI,
     $            ALPHAR, BETA, TMP1, TMP2, ILOQ, IHIQ, ILOZ, IHIZ, 
     $            TMP3, TMP4, TAR, TAI, TBETA, DWORK, LDWORK, IWORK, 
     $            LIWORK, IERR )

               IF ( IERR .NE. 0 ) THEN
                  IF ( IAM .EQ. 0 ) THEN
                     WRITE(*,*)'% PQZ1:PDHGEQZ5 failed. INFO=', IERR
                     CALL FLUSH( 6 )
                  END IF
                  INFO = 5
                  GOTO 420
               END IF 
               
            END IF    
         END IF
         
         KBOT = KTOP - 1
         IF ( KBOT .LT. IFIRST ) GOTO 380
         KTOP = IFIRST
         GOTO 350              
 350     CONTINUE    

*        End of iteration loop ..      
 360  CONTINUE
      
*     Drop-through = non-convergence
      INFO = ILAST
      GO TO 420      

*     Successful completion of all QZ steps  ..
 380  CONTINUE
      
*     Normal Termination
      INFO = 0

*     Exit (perhaps with INFO<>0)
 420  CONTINUE

      DEALLOCATE( TMP1, TMP2, TMP3, TMP4, TAR, TAI, TBETA )

 9998 FORMAT( A, I6, I6, I6 )
 9999 FORMAT( A, I6, I6, I6, I6, I6 )

      RETURN

      IF ( DEBUG ) THEN
         WRITE(*,*)'Leaving PQZ1'   
      END IF 
*     
*     End of PDHGEQZ1
*  
      END
