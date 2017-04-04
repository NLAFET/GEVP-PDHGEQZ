***********************************************************************
*                                                                     *
*     PDHGEQZ0.f:                                                     *
*         Auxillary routine in the package PDHGEQZ.                   *
*                                                                     *
*     Contributors: Bjorn Adlerborn                                   *
*                   Bo Kagstrom                                       *
*                   Daniel Kressner                                   *
*                                                                     *
*     Department of Computing Science and HPC2N, Umea University      *
*     MATHICSE ANCHP, EPF Lausanne                                    *
*	                                                              * 
***********************************************************************
      RECURSIVE SUBROUTINE PDHGEQZ0(ILSCHR, ILQ, ILZ, N, ILO, IHI, H, DE
     $   SCH, T, DESCT, ALPHAR, ALPHAI, BETA, Q, DESCQ, Z, DESCZ, WORK, 
     $   LWORK, IWORK, LIWORK, INFO, RLVL)

      IMPLICIT NONE
*     ..
*     .. Scalar Arguments ..
*     ..
      LOGICAL            ILSCHR, ILQ, ILZ
      INTEGER            IHI, ILO, INFO, LWORK, N, LIWORK, RLVL
*     ..
*     .. Array Arguments .. 
*     ..
      DOUBLE PRECISION   H( * ), T( * ), Q( * ), Z( * ), ALPHAI( * ), 
     $   ALPHAR( * ), BETA( * ), WORK( * )
      INTEGER            IWORK( * ), DESCH( 9 ), DESCT( 9 ), DESCQ( 9 ), 
     $   DESCZ( 9 )
      
      

*     Purpose
*     =======
*     
*     PDHGEQZ : implements a parallel single-/double-shift version of the QZ method for
*     finding the generalized eigenvalues
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
*     On exit, if JOB = 'S', H is quasi-upper triangular in rows and columns 
*     (ILO:IHI), with 1 x 1 and 2 x 2 blocks on the diagonal where 
*     the 2 x 2 blocks corresponds to complex conjugated pairs of 
*     eigenvalues. If JOB='E', H is unspecified on exit.      
*     
*     DESCH   (global and local input) INTEGER array of dimension DLEN_.
*     The array descriptor for the distributed matrix H.
*     
*     T       (local input/output) DOUBLE PRECISION array, dimension (LLD_T, LOCc(N)).
*     On entry, the N-by-N upper triangular matrix B.
*     On exit, the upper triangular matrix T = Q' B Z.  The
*     elements below the diagonal are set to zero.
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
*     elements of A that would result from reducing H and T to
*     Schur form and then further reducing them both to triangular
*     form using unitary transformations s.t. the diagonal of T
*     was non-negative real.  Thus, if H(j,j) is in a 1-by-1 block
*     (i.e., H(j+1,j)=H(j,j+1)=0), then ALPHAI(j)=0.
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
*     If COMPQ='N':  Q is not referenced.
*     If COMPQ='I':  on entry, Q need not be set, and on exit it
*     contains the orthogonal matrix Q, where Q'
*     is the product of the Givens transformations
*     which are applied to A and B on the left.
*     If COMPQ='V':  on entry, Q must contain an orthogonal matrix
*     Q1, and on exit this is overwritten by Q1*Q.
*     
*     DESCQ   (global and local input) INTEGER array of dimension DLEN_.
*     The array descriptor for the distributed matrix Q.
*     
*     Z       (local input/output) DOUBLE PRECISION array, dimension (LLD_Z, LOCc(N)).
*     If COMPZ='N':  Z is not referenced.
*     If COMPZ='I':  on entry, Z need not be set, and on exit it
*     contains the orthogonal matrix Z, which is
*     the product of the Givens transformations
*     which are applied to A and B on the right.
*     If COMPZ='V':  on entry, Z must contain an orthogonal matrix
*     Z1, and on exit this is overwritten by Z1*Z.
*     
*     DESCZ   (global and local input) INTEGER array of dimension DLEN_.
*     The array descriptor for the distributed matrix Z.
*     
*     WORK    (local workspace/global output) DOUBLE PRECISION array,
*     dimension (LWORK) 
*     
*     LWORK   (global input) INTEGER 
*     The dimension of the array WORK. 
*     If LWORK = -1, then a workspace query is assmed, and optimal 
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
*     > 0: The problem did not converge as expected. 1, 2, and 3, means
*     some sub routine failed.
*     RVLV    (global input) INTEGER
*     Indicates what level of recursion PDHGEQZ currently executes in.
*     0 or 1 are accepted values, and should initially be set to 0.
*     
*     =====================================================================



*     ..
*     .. Parameters ..
*     ..
      DOUBLE PRECISION   ZERO, ONE, TWO
      PARAMETER          ( ZERO = 0.0D+0, ONE = 1.0D+0, TWO = 2.0D+0 )
*     .. 
*     .. Parameters with common block ..
*     ..
      DOUBLE PRECISION    SMLNUM, ULP, TFNORM, T1NORM, HFNORM, HTOL,
     $   TTOL
      COMMON /PREC/       SMLNUM, ULP, TFNORM, T1NORM, HFNORM, HTOL,
     $   TTOL

      INTEGER            BLOCK_CYCLIC_2D, CSRC_, CTXT_, DLEN_, DT_, 
     $   LLD_, MB_, M_, NB_, N_, RSRC_
      PARAMETER          ( BLOCK_CYCLIC_2D = 1, DLEN_ = 9, DT_ = 1,
     $   CTXT_ = 2, M_ = 3, N_ = 4, MB_ = 5, NB_ = 6,
     $   RSRC_ = 7, CSRC_ = 8, LLD_ = 9 )

      
*     ..
*     .. Local Scalars ..
*     ..      
      LOGICAL             SORTED, LQUERY, DEBUG
      INTEGER             NMIN2, ITER1, ITER2, ITER3, 
     $   ITER4, J, JITER, MAXIT, ICTXT, NPROW, NPCOL, MYROW, 
     $   MYCOL, II, KLN, LWKOPT, LIWKOPT, KBOT, KTOP, NH, NS, ND, KWTOP, 
     $   IPTMPZ, IPTMPQ, IPTMPA, IPTMPB, I, IERR, NB, IPWORK, 
     $   USEDSHIFTS, IFIRST, IROFFH, ACSRC, ARSRC, NROWS, NCOLS, IAM, 
     $   NRPROCS, NIBBLE, KS, JW, LS, NSR, NWR, NW, K, S
     $   HIFTPERWND, NDINFTOT, ICOFFH, NUMWIN, 
     $   LLDTMP, IDUM, OLDKTOP, NDINF, ILOQ, IHIQ, ILOZ, IHIZ
      DOUBLE PRECISION    D1, D2, SWAP, T1, TMPT1, TMPT2,
     $   CSL, SNL, CSR, SNR, STARTTIME


*     .. Local Arrays ..
*     ..
      DOUBLE PRECISION    UPDTIME( 40 ), TMPA( 2, 2 ), TMPB( 2, 2 )
      INTEGER             DESCTMPQ( 9 ), DESCTMPZ( 9 ), DESCTMPA( 9 ),
     $   DESCTMPB( 9 ), DESCW( 9 )
      
*     ..
*     .. Externals ..
*     ..      
      EXTERNAL            LSAME, PDLAMCH, PDLANHS, PDLANTR, ICEIL, 
     $   NUMROC, INDXG2L, INDXG2P, PDLANGE, PILAENVX, MPI_WTIME
      DOUBLE PRECISION    PDLAMCH, PDLANHS, PDLANTR, PDLANGE, MPI_WTIME
      INTEGER             INDXG2L, INDXG2P, ICEIL, NUMROC, PILAENVX
      LOGICAL             LSAME

*     ..
*     .. Intrinsic Functions ..
*     ..
      INTRINSIC           ABS, DBLE, MAX, MIN, SQRT, MOD, LOG, NINT
      
*     ..
*     .. Executable Statements ..
*     ..


      STARTTIME = MPI_WTIME()

      IF( N .LE. 0 .OR. IHI .LT. ILO ) THEN
         RETURN
      END IF
      
      
*     
*     Extract current communication context and get grid parameters
*     
      ICTXT = DESCH(CTXT_)
      CALL BLACS_PINFO(IAM, NRPROCS)
      CALL BLACS_GRIDINFO( ICTXT, NPROW, NPCOL, MYROW, MYCOL )
*     Set debug to false/true for various info printed during the run      
      DEBUG = RLVL .EQ. 0
      DEBUG = .FALSE.
      DEBUG = DEBUG .AND. IAM .EQ. 0

*     
*     NH holds the current problem size, for now the whole matrix
*     
      NH = IHI - ILO + 1

      IPWORK = 1
      KBOT = IHI
      KTOP = ILO  
      ILOQ = 1
      IHIQ = N
      ILOZ = 1
      IHIZ = N
*     
*     NMIN2 defines the minium problem size to be solved by PDHGEQZ0
*     If remaming problem is less than this value, execution is passed 
*     to PDHGEQZ1
*     
      NMIN2 = PILAENVX(ICTXT, 54, 'PDHGEQZ', '', IDUM, IDUM, IDUM, IDUM)
      IF ( DEBUG ) WRITE(*,*)'% PQZ0: NMIN2 = ', NMIN2

*     
*     NB_ .EQ. MB_ is assumed, we use NB as blocking factor.
*     
      NB = DESCH(NB_)


      
      T1 = 0
      DO II = 1, 40
         UPDTIME(II) = 0
      END DO
      
      
*     
*     Check for Workspace query
*     
      LQUERY = LWORK .EQ. -1 . OR. LIWORK .EQ. -1
      IF ( DEBUG ) WRITE(*,*)'% PQZ0 : LQUERY = ', LQUERY
      LWKOPT = 0
      LIWKOPT = 0
*     
*     Workspace query call to PDHGEQZ3
*     
      CALL PDHGEQZ3( ILSCHR, ILQ, ILZ, N, 
     $   H, DESCH, T, DESCT, Q, DESCQ, Z, DESCZ, 
     $   ILOQ, IHIQ, ILOZ, IHIZ,
     $   WORK, DESCTMPA,
     $   WORK, DESCTMPB,
     $   WORK, DESCTMPQ,
     $   WORK, DESCTMPZ,
     $   KTOP, KBOT, ALPHAR, ALPHAI, BETA,
     $   N, ND,
     $   N, WORK, -1, IWORK, -1 , IERR,
     $   UPDTIME, RLVL )
      IF ( DEBUG ) WRITE(*,*)'% PQZ0 -> PQZ3 WQUERY: ', 
     $   INT( WORK( 1 ) ), INT( IWORK( 1 ) )
      LWKOPT = LWKOPT + INT( WORK( 1 ) )
      LIWKOPT = IWORK(1) + N
*     
*     Workspace query call to PDHGEQZ5
*     
      CALL PDHGEQZ5(ILSCHR, ILQ, ILZ, 
     $   NW, ILO, IHI, N, ILO, IHI,
     $   H, DESCH, T, DESCT,
     $   ALPHAR, ALPHAI, BETA, Q, DESCQ, Z, DESCZ,
     $   ILOQ, IHIQ, ILOZ, IHIZ,
     $   1, USEDSHIFTS, HTOL, WORK,
     $   -1, IERR, UPDTIME, 1 )

      IF ( DEBUG ) WRITE(*,*)'% PQZ0 -> PQZ5 WQUERY:', INT( WORK( 1 ) )
      LWKOPT = MAX( LWKOPT, INT( WORK( 1 ) ) )
*     
*     Quick return if workspace query
*     
      IF ( LQUERY ) THEN
         WORK(1) = LWKOPT
         IWORK( 1 ) = LIWKOPT
         RETURN
      END IF

*     
*     Check Argument Values
*     
      INFO = 0
      IF( .NOT. LQUERY .AND. LWORK .LT. MAX ( 1, LWKOPT ) ) THEN
         INFO = -1
      ELSE IF ( .NOT. LQUERY .AND. LIWORK .LT. MAX ( 1, LIWKOPT ) ) THEN
         INFO = -2
      END IF
      IF( INFO .NE. 0 ) THEN
         CALL PXERBLA( ICTXT, 'PDHGEQZ0', -INFO )
         RETURN
      END IF


      
*     
*     NIBBLE for when to skip QZ sweep 
*     
      NIBBLE = PILAENVX(ICTXT, 52, 'PDHGEQZ0', '', 
     $   N, ILO, IHI, IDUM)
      IF ( DEBUG ) WRITE(*,*)'% PQZ0 : NIBBLE = ', NIBBLE  
      
*     
*     Set number of shifts to use as a function of problem size 
*     
      NS = PILAENVX( ICTXT, 53, 'PDHGEQZ0', '', N, ILO, IHI,
     $   NB )
      NS = MIN( NS, IHI - ILO )
      NS = MAX( 2, NS - MOD( NS, 2 ) )
      IF ( DEBUG ) WRITE(*,*)'% PQZ0 : NS = ', NS


*     
*     Set deflation window size, store in NW
*     
      IF( NH .LE. 500 ) THEN
         NWR = NS
      ELSE
         NWR = 3 * NS / 2
      END IF        
      NWR = MAX( 2, NWR )
      NWR = MIN( IHI - ILO + 1, NWR )
      NW = NWR
      
*     
*     NSR = recommended number of simultaneous shifts
*     
      NSR = NS
      MAXIT = 30*( IHI-ILO+1 )
      KBOT = IHI
      KTOP = ILO        

      
      NDINF = 0
*     
*     Some Counters
*     ITER1 : Number of times smaller problems have been encounterd
*     ITER2 : Number of QZ sweeps 
*     ITER3 : Number of AED calls 
*     ITER4 : Number of times more shifts are retreived 
*     ( AED returned to few ) 
*     JITER : Holds the number of total iterations
*     USEDSHIFTS : Number of shifts applied
*     NDINFTOT : Total number for deflated infinite eigenvalues
*     

      ITER1 = 0
      ITER2 = 0            
      ITER3 = 0
      ITER4 = 0   
      USEDSHIFTS = 0
      NDINFTOT = 0
      

      IFIRST = KTOP      
      
      IF ( DEBUG ) THEN
         WRITE(*,9996)'% PQZ0 init:', KTOP, KBOT, NW, NSR, 0, 
     $      MPI_WTIME() - STARTTIME

         CALL FLUSH(6)
      END IF
*     Set MAXIT to limit the number of iteration, for debugging purpose
      DO 360 JITER = 1, MAXIT
         
         IF (KBOT .LT. IFIRST) GOTO 380
         
*        
*        Locate active block.
*        
         IPWORK = 1
         CALL PDLASMSUB( H, DESCH, KBOT, IFIRST,  J , SMLNUM, 
     $      WORK( IPWORK ), LWORK)          
         KTOP = J  

         IF( KTOP .GT. IFIRST ) THEN
*           
*           H( KTOP, KTOP - 1) is negligible
*           
            CALL PDELSET( H, KTOP, KTOP - 1, DESCH, ZERO )
         END IF

         IF ( DEBUG ) THEN
            WRITE(*,9996)'% PQZ0 : Call PQZ7:', JITER, KTOP, 
     $         KBOT, 0, 0,
     $         MPI_WTIME()-STARTTIME

            CALL FLUSH(6)
         END IF

         NDINF = 0
         OLDKTOP = KTOP
         TMPT1 = MPI_WTIME() 
*        
*        Check for and deflate infinte eigenvalues
*        NDINF will be updated to reflect number of found infinte 
*        eigenvalues - KTOP and KBOT might also change.
*        
         CALL PDHGEQZ7( ILSCHR, ILQ, ILZ, 
     $      KTOP, KBOT, ILO, IHI, 
     $      H, DESCH, T, DESCT,
     $      ALPHAR, ALPHAI, BETA, 
     $      Q, DESCQ, Z, DESCZ,
     $      ILOQ, IHIQ, ILOZ, IHIZ,
     $      MIN(NPROW, NPCOL, 40), NDINF,
     $      WORK, LWORK, IWORK, IERR )
         IF (IERR.NE.0) THEN
            IF ( IAM.EQ.0 ) THEN
               WRITE(*,*)'PQZ0 : PDHGEQZ7
     $            (Small sub) failed with INFO=', IERR
            END IF
            INFO = 1
            RETURN
         END IF 

         TMPT2 = MPI_WTIME()
         UPDTIME(19) = UPDTIME(19) + (TMPT2-TMPT1)

*        
*        If deflation made at top, adjust where to start next iteration
*        
         IF ( KTOP .NE. OLDKTOP ) THEN
            IF ( OLDKTOP .EQ. IFIRST ) THEN
               IFIRST = KTOP
            END IF
         END IF
*        
*        Try another standard deflation if inf. eigvs found
*        
         IF ( NDINF .GT. 0 ) THEN
            NDINFTOT = NDINFTOT + NDINF
            IF (DEBUG) THEN
               WRITE(*,*)'% PQZ0 : // Intermediate deflation ', NDINF,
     $            ' oo eigenvalues.', KBOT, KTOP, 
     $            MPI_WTIME()-STARTTIME
            ENDIF
            IF (KTOP.GE.KBOT) THEN
               IF (KBOT.EQ.KTOP) THEN
                  CALL PDELGET('A', ' ', ALPHAR( KBOT ), H, KBOT, KBOT,
     $               DESCH)
                  CALL PDELGET('A', ' ', BETA( KBOT ), T, KBOT, KBOT,
     $               DESCT)
                  ALPHAI( KBOT ) = ZERO
                  
               END IF
               KBOT = OLDKTOP - 1
               KTOP = IFIRST
            END IF
            GOTO 565 
         END IF




         NH = KBOT - KTOP + 1         
         IF ( NH .LE. NMIN2 ) THEN
            ITER1 = ITER1 + 1
            JW = NH
            IROFFH = MOD(KTOP - 1, NB)
            ICOFFH = IROFFH
            
            IF ( DEBUG ) THEN
               WRITE(*,9996)'% PQZ0 : Small Sub Iter=', JITER, JW, KTOP,
     $            KBOT, 0,
     $            MPI_WTIME()-STARTTIME
               CALL FLUSH(6)
            END IF              
            IF( JW.EQ.1 ) THEN
               TMPT1 = MPI_WTIME()
*              
*              H(ILAST,ILAST-1) is negligible: one eigenvalue has converged.
*              
               CALL PDELGET('A', ' ', ALPHAR( KBOT ), H, KBOT, KBOT,
     $            DESCH)
               CALL PDELGET('A', ' ', BETA( KBOT ), T, KBOT, KBOT,
     $            DESCT)
               ALPHAI( KBOT ) = ZERO
               TMPT2 = MPI_WTIME()
               UPDTIME(13) = UPDTIME(13) + (TMPT2-TMPT1)

            ELSE IF( JW.EQ.2 ) THEN
*              
*              H(ILAST-1,ILAST-2) is negligible: a pair of eigenvalues have converged.
*              
               TMPT1 = MPI_WTIME()
               CALL PDLACP3(2, KBOT-1, H, DESCH, TMPA, 2, -1, -1, 0)
               CALL PDLACP3(2, KBOT-1, T, DESCT, TMPB, 2, -1, -1, 0)
               
               CALL DLAGV2( TMPA, 2, TMPB, 2, 
     $            ALPHAR( KBOT-1), ALPHAI( KBOT-1), BETA( KBOT-1),
     $            CSL, SNL, CSR, SNR )
               CALL PDLACP3(2, KBOT-1, H, DESCH, TMPA, 2, 0, 0, -1)
               CALL PDLACP3(2, KBOT-1, T, DESCT, TMPB, 2, 0, 0, -1)
               TMPT2 = MPI_WTIME()
               UPDTIME(13) = UPDTIME(13) + (TMPT2-TMPT1)
               TMPT1 = MPI_WTIME()
*              
*              Update cols of H and T
*              
               CALL PDROT(KBOT-2, H, 1, KBOT-1, DESCH, 1,
     $            H, 1, KBOT, DESCH, 1,
     $            CSR, SNR, WORK(IPWORK), LWORK, INFO)
               
               CALL PDROT(KBOT-2, T, 1, KBOT-1, DESCT, 1,
     $            T, 1, KBOT, DESCT, 1,
     $            CSR, SNR, WORK(IPWORK), LWORK, INFO)

*              
*              Update rows of H and T
*              
               IF (IHI.GT.KBOT) THEN
                  CALL PDROT(IHI-KBOT, H, KBOT-1, KBOT+1, DESCH, N,
     $               H, KBOT, KBOT+1, DESCH, N,
     $               CSL, SNL, WORK(IPWORK), LWORK, INFO)
                  
                  CALL PDROT(IHI-KBOT, T, KBOT-1, KBOT+1, DESCT, N,
     $               T, KBOT, KBOT+1, DESCT, N,
     $               CSL, SNL, WORK(IPWORK), LWORK, INFO)
               END IF
*              
*              Update cols of Q and Z
*              
               IF (ILQ) THEN
                  CALL PDROT(N, Q, 1, KBOT-1, DESCQ, 1,
     $               Q, 1, KBOT, DESCQ, 1,
     $               CSL, SNL, WORK(IPWORK), LWORK, INFO)
               END IF
               
               IF (ILZ) THEN
                  CALL PDROT(N, Z, 1, KBOT-1, DESCZ, 1,
     $               Z, 1, KBOT, DESCZ, 1,
     $               CSR, SNR, WORK(IPWORK), LWORK, INFO)
               END IF

               TMPT2 = MPI_WTIME()
               UPDTIME(14) = UPDTIME(14) + (TMPT2-TMPT1)

            ELSE                            
*              
*              A smaller sub problem has been identified
*              
               TMPT1 = MPI_WTIME()
               ARSRC = INDXG2P( KTOP, NB, MYROW, DESCH(RSRC_), 
     $            NPROW )
               ACSRC = INDXG2P( KTOP, NB, MYCOL, DESCH(CSRC_), 
     $            NPCOL )
               
               NCOLS = NUMROC( JW, NB, MYCOL, ACSRC, NPCOL )
               NROWS = NUMROC( JW, NB, MYROW, ARSRC, NPROW )

               CALL DESCINIT( DESCTMPQ, JW, JW, 
     $            NB, NB, ARSRC, ACSRC, ICTXT,  
     $            MAX( 1, NROWS ), INFO )
               
               CALL DESCINIT( DESCTMPZ, JW, JW, 
     $            NB, NB, ARSRC, ACSRC, ICTXT,
     $            MAX( 1, NROWS ), INFO )     
               
               CALL DESCINIT( DESCTMPA, JW, JW, 
     $            NB, NB, ARSRC, ACSRC, ICTXT,  
     $            MAX( 1, NROWS ), INFO )
               
               CALL DESCINIT( DESCTMPB, JW, JW, 
     $            NB, NB, ARSRC, ACSRC, ICTXT,
     $            MAX( 1, NROWS ), INFO )     

               IPTMPA = 1
               IPTMPB = IPTMPA + DESCTMPA( LLD_ )*NCOLS
               IPTMPZ = IPTMPB + DESCTMPB( LLD_ )*NCOLS
               IPTMPQ = IPTMPZ + DESCTMPZ( LLD_ )*NCOLS
               IPWORK = IPTMPQ + DESCTMPQ( LLD_ )*NCOLS    
               
               CALL PDLASET('A', JW, JW, ZERO, ONE, 
     $            WORK(IPTMPQ), 1, 1, DESCTMPQ)                
               CALL PDLASET('A', JW, JW, ZERO, ONE, 
     $            WORK(IPTMPZ), 1, 1, DESCTMPZ)
               
               IF (KTOP.EQ.1) THEN
                  CALL PDLACPY('A', JW, JW, H, 1, 1, DESCH,
     $               WORK(IPTMPA), 1, 1, DESCTMPA )
                  CALL PDLACPY('A', JW, JW, T, 1, 1, DESCT,
     $               WORK(IPTMPB), 1, 1, DESCTMPB )
               ELSE
                  CALL PDGEMR2D( JW, JW, H, KTOP, KTOP, DESCH,
     $               WORK(IPTMPA), 1, 1, DESCTMPA, 
     $               ICTXT )              
                  CALL PDGEMR2D( JW, JW, T, KTOP, KTOP, DESCT,
     $               WORK(IPTMPB), 1, 1, DESCTMPB, 
     $               ICTXT )              
               END IF

               CALL PDHGEQZ1( .TRUE., .TRUE., .TRUE., JW, 
     $            1, JW, 
     $            WORK(IPTMPA), DESCTMPA, WORK(IPTMPB), DESCTMPB,
     $            ALPHAR(KTOP), ALPHAI(KTOP), 
     $            BETA(KTOP),
     $            WORK(IPTMPQ), DESCTMPQ, WORK(IPTMPZ), DESCTMPZ,
     $            WORK(IPWORK), LWORK, IWORK, LIWORK, IERR )
               IF (IERR.NE.0) THEN
                  IF ( IAM.EQ.0 ) THEN
                     WRITE(*,*)'PQZ0 : PDHGEQZ1
     $                  (Small sub) failed with INFO=', 
     $IERR
                  END IF
                  INFO = 2
                  RETURN
               END IF               

               IF (KTOP.EQ.1) THEN
                  CALL PDLACPY( 'A', JW, JW, 
     $               WORK(IPTMPA), 1, 1, DESCTMPA,
     $               H, 1, 1, DESCH )
                  CALL PDLACPY( 'A', JW, JW, 
     $               WORK(IPTMPB), 1, 1, DESCTMPB,
     $               T, 1, 1, DESCT )
               ELSE                                       
                  CALL PDGEMR2D( JW, JW, WORK(IPTMPA), 1, 
     $               1, DESCTMPA, H, KTOP, KTOP, DESCT,
     $               ICTXT )              
                  CALL PDGEMR2D( JW, JW, WORK(IPTMPB), 1, 
     $               1, DESCTMPB, T, KTOP, KTOP, DESCT, 
     $               ICTXT )    
               END IF
               
               
               TMPT2 = MPI_WTIME()
               UPDTIME(13) = UPDTIME(13) + (TMPT2-TMPT1)
               TMPT1 = MPI_WTIME()
               
               IROFFH = 0
               ICOFFH = 0
*              
*              Update horizontal slab in A and B with TMPQ 
*              
               KLN = N - (KTOP + JW) + 1
               IF (KLN.GT.0) THEN
                  
                  LLDTMP = NUMROC( JW, NB, MYROW, ARSRC, NPROW )
                  LLDTMP = MAX( 1, LLDTMP )
                  
                  CALL DESCINIT( DESCW, JW, KLN, 
     $               NB, NB, 
     $               ARSRC, ACSRC, ICTXT, LLDTMP, INFO )
                  
                  CALL PDGEMM( 'T', 'N', JW, KLN, JW, ONE, 
     $               WORK(IPTMPQ), 1+IROFFH, 1+ICOFFH, DESCTMPQ, 
     $               H, KTOP, KTOP+JW, DESCH, 
     $               ZERO, WORK(IPWORK), 1, 1, DESCW )
                  
                  CALL PDGEMR2D( JW, KLN,  WORK(IPWORK), 1, 1, 
     $               DESCW, H, KTOP, KTOP+JW, DESCH, ICTXT )  
                  
                  CALL PDGEMM( 'T', 'N', JW, KLN, JW, ONE, 
     $               WORK(IPTMPQ), 1+IROFFH, 1+ICOFFH, DESCTMPQ, 
     $               T, KTOP, KTOP+JW, DESCT, 
     $               ZERO, WORK(IPWORK), 1, 1, DESCW )
                  
                  CALL PDGEMR2D( JW, KLN,  WORK(IPWORK), 1, 1,
     $               DESCW, T, KTOP, KTOP+JW, DESCT, ICTXT )  

               END IF
*              
*              Update vertical slab in A and B with TMPZ..
*              
               KLN = KTOP - ILO
               IF (KLN.GT.0) THEN
                  
                  LLDTMP = NUMROC( KLN, NB, MYROW, ARSRC, NPROW )
                  LLDTMP = MAX( 1, LLDTMP )              
                  CALL DESCINIT( DESCW, KLN, JW, NB, NB, 
     $               ARSRC, ACSRC, ICTXT, LLDTMP, INFO )
                  
                  CALL PDGEMM( 'N', 'N', KLN, JW, JW, 
     $               ONE, H, ILO, KTOP, DESCH, 
     $               WORK(IPTMPZ), 1+IROFFH, 1+ICOFFH, DESCTMPZ,  
     $               ZERO, WORK(IPWORK), 1, 1, DESCW ) 
                  
                  CALL PDGEMR2D( KLN, JW, WORK(IPWORK), 1, 1, DESCW,
     $               H, ILO, KTOP, DESCH, ICTXT )  

                  
                  CALL PDGEMM( 'N', 'N', KLN, JW, JW, 
     $               ONE, T, ILO, KTOP, DESCT, 
     $               WORK(IPTMPZ), 1+IROFFH, 1+ICOFFH, DESCTMPZ,  
     $               ZERO, WORK(IPWORK), 1, 1, DESCW ) 

                  CALL PDGEMR2D( KLN, JW, WORK(IPWORK), 1, 1, DESCW,
     $               T, ILO, KTOP, DESCT, ICTXT )  
               END IF
*              
*              Update vertical slab in Q with TMPQ ..
*              
               IF( ILQ ) THEN
                  KLN = IHI - ILO + 1
                  LLDTMP = NUMROC( KLN, NB, MYROW, ARSRC, NPROW )
                  LLDTMP = MAX( 1, LLDTMP )
                  CALL DESCINIT( DESCW, KLN, JW, NB, NB,
     $               ARSRC, ACSRC, ICTXT, LLDTMP, INFO )
                  CALL PDGEMM( 'N', 'N', KLN, JW, JW, 
     $               ONE, Q, ILO, KTOP, DESCQ, 
     $               WORK(IPTMPQ), 1+IROFFH, 1+ICOFFH, DESCTMPQ, 
     $               ZERO, WORK(IPWORK), 1, 1, DESCW )   
                  CALL PDGEMR2D( KLN, JW, WORK(IPWORK), 1, 1, DESCW,
     $               Q, ILO, KTOP, DESCQ, ICTXT )              
               END IF
               
*              
*              Update vertical slab in Z with TMPZ ..
*              
               IF( ILZ ) THEN
                  KLN = IHI - ILO + 1
                  LLDTMP = NUMROC( KLN, NB, MYROW, ARSRC, NPROW )
                  LLDTMP = MAX( 1, LLDTMP )
                  CALL DESCINIT( DESCW, KLN, JW, NB, NB,
     $               ARSRC, ACSRC, ICTXT, LLDTMP, INFO )
                  CALL PDGEMM( 'N', 'N', KLN, JW, JW, 
     $               ONE, Z, ILO, KTOP, DESCZ, 
     $               WORK(IPTMPZ), 1+IROFFH, 1+ICOFFH, DESCTMPZ, 
     $               ZERO, WORK(IPWORK), 1, 1, DESCW )   
                  CALL PDGEMR2D( KLN, JW, WORK(IPWORK), 1, 1, DESCW,
     $               Z, ILO, KTOP, DESCZ, ICTXT )
     $               
               ENDIF  
               TMPT2 = MPI_WTIME()
               UPDTIME(14) = UPDTIME(14) + (TMPT2-TMPT1)

            ENDIF
            
            KBOT = KTOP - 1
            KTOP = IFIRST
            IPWORK = 1
            GOTO 350
            
         ENDIF
*        
*        Prepare for AED by selecting AED window size JW
*        
         JW = MIN( NW, NH )
         KWTOP = KBOT - JW + 1
         IROFFH = MOD(KWTOP-1, NB)
         ICOFFH = IROFFH         
         
         ARSRC = INDXG2P( KWTOP, NB, MYROW, DESCH(RSRC_), NPROW )
         ACSRC = INDXG2P( KWTOP, NB, MYCOL, DESCH(CSRC_), NPCOL )
         
         
         NCOLS = NUMROC( JW+IROFFH, NB, MYCOL, ACSRC, NPCOL )   
         NROWS = NUMROC( JW+IROFFH, NB, MYROW, ARSRC, NPROW )   

         CALL DESCINIT( DESCTMPQ, JW+IROFFH, JW+IROFFH, 
     $      NB, NB, ARSRC, ACSRC, ICTXT,  
     $      MAX( 1, NROWS ), INFO )
         
         CALL DESCINIT( DESCTMPZ, JW+IROFFH, JW+IROFFH, 
     $      NB, NB, ARSRC, ACSRC, ICTXT,
     $      MAX( 1, NROWS ), INFO )     
         
         CALL DESCINIT( DESCTMPA, JW+IROFFH, JW+IROFFH, 
     $      NB, NB, ARSRC, ACSRC, ICTXT,
     $      MAX( 1, NROWS ), INFO )

         CALL DESCINIT( DESCTMPB, JW+IROFFH, JW+IROFFH, 
     $      NB, NB, ARSRC, ACSRC, ICTXT,
     $      MAX( 1, NROWS ), INFO )  
         
         

         IPTMPA = 1
         IPTMPB = IPTMPA + DESCTMPA( LLD_ )*NCOLS
         IPTMPZ = IPTMPB + DESCTMPB( LLD_ )*NCOLS
         IPTMPQ = IPTMPZ + DESCTMPZ( LLD_ )*NCOLS
         IPWORK = IPTMPQ + DESCTMPQ( LLD_ )*NCOLS  
         
         
         TMPT1 = MPI_WTIME()

         
         ND = 0
*        
*        Time for AED, ND will hold number of deflated eigenvalues
*        
         IF ( DEBUG ) THEN
            WRITE(*,9996)'% PQZ0 : Call PQZ3', JITER, KTOP, KWTOP, 
     $         KBOT, 0,
     $         MPI_WTIME()-STARTTIME

            CALL FLUSH(6)
         END IF
         IERR = 0

         CALL PDHGEQZ3( ILSCHR, ILQ, ILZ, N, 
     $      H, DESCH, T, DESCT, Q, DESCQ, Z, DESCZ, 
     $      ILOQ, IHIQ, ILOZ, IHIZ, 
     $      WORK( IPTMPA ), DESCTMPA, 
     $      WORK( IPTMPB ), DESCTMPB, 
     $      WORK( IPTMPQ ), DESCTMPQ,
     $      WORK( IPTMPZ ), DESCTMPZ,
     $      KTOP, KBOT, ALPHAR, ALPHAI, BETA, 
     $      JW, ND, KWTOP, 
     $      WORK( IPWORK ), LWORK, IWORK, LIWORK, IERR,
     $      UPDTIME, RLVL )

         IF ( IERR .NE. 0 ) THEN
            IF ( DEBUG ) THEN
               WRITE(*,*)'% PQZ0 : PQZ3 failed with INFO=', IERR
            END IF
            ND = 0
         ENDIF

         ITER3 = ITER3 + 1

         IF ( DEBUG ) THEN
            WRITE(*,9996)'% PQZ0 : Iter=', JITER, KTOP, KBOT, JW, ND, 
     $         MPI_WTIME() - STARTTIME
            CALL FLUSH(6)
         ENDIF            

         IPWORK  = 1

         TMPT2 = MPI_WTIME()
         
         UPDTIME(1) = UPDTIME(1) + ( TMPT2 - TMPT1 )
*        
*        Adjust KBOT accounting for new deflations.
*        
         KBOT = KBOT - ND
         IF (KBOT .LT. IFIRST) GOTO 380                           
         IF (KBOT .LT. KTOP) THEN
            KTOP = IFIRST
            GOTO 350          
         ENDIF
*        
*        Number of unconverged elements          
*        
         LS  = JW - ND
*        
*        KS points to the shifts
*        
         KS = KBOT - LS + 1
         
*        
*        Skip an expensive QZ sweep if there is a (partly
*        heuristic) reason to expect that many eigenvalues
*        will deflate without it.  Here, the QZ sweep is
*        skipped if many eigenvalues have just been deflated
*        or if the remaining active block is small.
*        
         IF( ( ND.EQ.0 ) .OR. ( ( 100*ND.LE.JW*NIBBLE ) .AND. 
     $      ( KBOT - KTOP + 1 .GT. NMIN2 ) ) ) THEN
            ITER2 = ITER2 + 1
*           
*           NS = nominal number of simultaneous shifts
*           
            NS = MIN( NSR, MAX( 2, KBOT - KTOP - 1 ) )
            NS = NS - MOD( NS, 2 )
            IF ( DEBUG ) THEN
               WRITE(*,9993)'% PQZ0 : Chase ', JW, NS, KS, 0,0,0,0
               CALL FLUSH(6)
            END IF  
*           
*           If needed, get more shifts from the trailing submatrix
*           This part of the code works, but might not be needed so commented out for now
*           
c           IF ( KBOT - KS + 1 .LE. (NS / 2 ) ) THEN
c           TMPT1 = MPI_WTIME()
c           IF (DEBUG) THEN
c           WRITE(*,9993)'% // Fetching more shifts ', 
c           $                 KBOT - NS + 1, MOD(KBOT - NS + 1, NB), 0,0,0,0,0
c           CALL FLUSH(6)
c           END IF
c           KS = KBOT - NS + 1
c           IROFFH = MOD( KS - 1, NB )
c           ICOFFH = IROFFH                  
c           HRSRC = INDXG2P( KS, NB, MYROW, DESCH(RSRC_),
c           $              NPROW )
c           HCSRC = INDXG2P( KS, NB, MYROW, DESCH(CSRC_),
c           $              NPCOL )
c           NROWS = NUMROC( NS+IROFFH, NB, MYROW, HRSRC,
c           $              NPROW )
c           NCOLS = NUMROC( NS+ICOFFH, NB, MYCOL, HCSRC,
c           $              NPCOL )
c           CALL DESCINIT( DESCTMPA, NS+IROFFH, NS+ICOFFH, NB,
c           $              NB, HRSRC, HCSRC, ICTXT, MAX(1, NROWS), INFO )
c           CALL DESCINIT( DESCTMPB, NS+IROFFH, NS+ICOFFH, NB,
c           $              NB, HRSRC, HCSRC, ICTXT, MAX(1, NROWS), INFO )
c           
c           IPTMPA = 1
c           IPTMPB = IPTMPA + DESCTMPA(LLD_) * NCOLS  
c           IPWORK = IPTMPB + DESCTMPB(LLD_) * NCOLS  
c           CALL PDLASET('A', NS+IROFFH, NS+ICOFFH, ZERO, ZERO, 
c           $              WORK(IPTMPA), 1, 1, DESCTMPA)
c           CALL PDLASET('A', NS+IROFFH, NS+ICOFFH, ZERO, ZERO, 
c           $              WORK(IPTMPB), 1, 1, DESCTMPB)
c           
c           CALL PDLACPY( 'All', NS, NS, H, KS, KS, DESCH,
c           $              WORK(IPTMPA), 1+IROFFH, 1+ICOFFH, DESCTMPA )
c           CALL PDLACPY( 'All', NS, NS, T, KS, KS, DESCT,
c           $              WORK(IPTMPB), 1+IROFFH, 1+ICOFFH, DESCTMPB )
c           INFO = 0
c           
c           
c           ITER4 = ITER4 + 1
c           CALL PDHGEQZ('E', 'N', 'N', 
c           $              IROFFH+NS, 1+IROFFH, IROFFH+NS,
c           $              WORK(IPTMPA), DESCTMPA, WORK(IPTMPB),DESCTMPB,
c           $              ALPHAR(KS-IROFFH), ALPHAI(KS-IROFFH), 
c           $              BETA(KS-IROFFH),
c           $              WORK(IPWORK), DESCQ, WORK(IPWORK), DESCZ, 
c           $              1, NS, 1, NS,
c           $              WORK(IPWORK), LWORK,
c           $              IWORK, LIWORK, IERR, RLVL+1)
c           
c           TMPT2 = MPI_WTIME()
c           UPDTIME(12) = UPDTIME(12) + (TMPT2-TMPT1)
c           
c           IF ( IERR .NE. 0 ) THEN
c           IF (IAM.EQ.0) THEN
c           WRITE(*,*)'PDHGEQZ:PDHGEQZ(1)(more eig) failed with INFO=', 
c           $                    IERR
c           END IF
c           INFO = 2
c           RETURN
c           END IF                        
c           IPWORK = 1
c           END IF
            
            IF( KBOT - KS + 1 .GT. NS ) THEN
*              
*              Sort the shifts (helps a little)
*              Bubble sort keeps complex conjugate
*              pairs together.
*              
               SORTED = .TRUE.
               DO K = KBOT, KS + 1, -1
                  IF( SORTED ) GO TO 90
                  SORTED = .TRUE.
                  DO I = KS, K - 1
                     D1 = ABS( ALPHAR( I ) )+ABS( ALPHAI( I ) )
                     D2 = ABS( ALPHAR( I+1 ) )+ABS( ALPHAI( I+1 ) )
                     D1 = ABS(D1 / BETA (I))
                     D2 = ABS(D2 / BETA (I + 1))
                     
                     IF( D1 .LT. D2 ) THEN
                        SORTED = .FALSE.
*                       
                        SWAP = ALPHAR( I )
                        ALPHAR( I ) = ALPHAR( I+1 )
                        ALPHAR( I+1 ) = SWAP
*                       
                        SWAP = ALPHAI( I )
                        ALPHAI( I ) = ALPHAI( I+1 )
                        ALPHAI( I+1 ) = SWAP
*                       
                        SWAP = BETA( I )
                        BETA( I ) = BETA( I+1 )
                        BETA( I+1 ) = SWAP      
                     END IF
                  END DO
               END DO
 90            CONTINUE
            END IF

*           
*           Shuffle shifts into pairs of real shifts
*           and pairs of complex conjugate shifts
*           assuming complex conjugate shifts are
*           already adjacent to one another. (Yes,
*           they are.)
*           
            DO I = KBOT, KS + 2, -2               
               IF(ALPHAI( I ).EQ.ZERO.AND.ALPHAI( I-1 ).NE.ZERO) THEN
                  SWAP = ALPHAR( I )
                  ALPHAR( I ) = ALPHAR( I-1 )
                  ALPHAR( I-1 ) = ALPHAR( I-2 )
                  ALPHAR( I-2 ) = SWAP
                  SWAP = ALPHAI( I )
                  ALPHAI( I ) = ALPHAI( I-1 )
                  ALPHAI( I-1 ) = ALPHAI( I-2 )
                  ALPHAI( I-2 ) = SWAP                     
                  SWAP = BETA( I )
                  BETA( I ) = BETA( I-1 )
                  BETA( I-1 ) = BETA( I-2 )
                  BETA( I-2 ) = SWAP                     
               END IF
            END DO
*           
*           Use up to NS of the the smallest magnatiude
*           shifts.  If there aren't NS shifts available,
*           then use them all, possibly dropping one to
*           make the number of shifts even.
*           
            NS = MIN( NS, KBOT - KS + 1 )
            NS = NS - MOD( NS, 2 )
            KS = KBOT - NS + 1
*           
*           QZ step
*           Introduce and chase the shifts in H and T at KTOP
*           
 100        CONTINUE

*           
*           Calculate number of shifts to use per window
*           
            SHIFTPERWND = NS
            SHIFTPERWND = MIN( NB / 3 - 1, SHIFTPERWND )
            SHIFTPERWND = SHIFTPERWND - MOD( SHIFTPERWND, 2 )
            SHIFTPERWND = MAX( MIN ( SHIFTPERWND, NS ), 2 )
*           
*           Calculate number of concurrent windows to use
*           
            NUMWIN = MAX( MIN( 40, NPROW, NPCOL, NS / SHIFTPERWND ), 1 )
            IF ( DEBUG ) THEN
               WRITE(*,9993)'% PQZ0 : B4 Chase ', NS, NUMWIN, 
     $            SHIFTPERWND, KTOP, KBOT, MYROW, MYCOL
               CALL FLUSH(6)
            END IF            

            IERR = 0

            CALL PDHGEQZ5( ILSCHR, ILQ, ILZ, 
     $         NS, KTOP, KBOT, N, ILO, IHI, 
     $         H, DESCH, T, DESCT,
     $         ALPHAR(KS), ALPHAI(KS), BETA(KS),  
     $         Q, DESCQ, Z, DESCZ,
     $         ILOQ, IHIQ, ILOZ, IHIZ,
     $         NUMWIN, USEDSHIFTS, HTOL, WORK( IPWORK ),
     $         LWORK, IERR, UPDTIME, SHIFTPERWND )
            
            IF ( IERR .NE. 0 ) THEN
               IF ( IAM.EQ.0 ) THEN
                  WRITE(*,*)'PQZ0 : PDHGEQZ5
     $               failed with INFO=', IERR
               END IF
               INFO = 3
               RETURN
            END IF 

*           
*           Update the number of shifts available ( since some are used now )
*           
            NS = NS - NUMWIN * SHIFTPERWND
*           
*           We do not setup another set of windows if the number of available
*           shifts are too few
*           
            IF ( NS .GT. SHIFTPERWND ) GOTO 100
            
         ELSE
            IF ( DEBUG ) THEN
               WRITE(*,9998)'% PQZ0 : AED again ,NIB=', NIBBLE, JW, 
     $            100*ND
               CALL FLUSH(6)
            END IF
         END IF
 565     CONTINUE

*        
*        End of iteration loop ..
*        
 350     CONTINUE         
 360  CONTINUE
*     
*     Drop-through = non-convergence
*     

      INFO = KBOT
      GO TO 420     
*     
*     Successful completion of all QZ steps
*     
 380  CONTINUE
*     
*     Normal Termination
*     
      INFO = 0
*     
*     Exit (other than argument error) -- return optimal workspace size
*     
 420  CONTINUE


 9996 FORMAT(A,I8,I8,I8,I8,I8, F13.1)  
 9993 FORMAT(A,I8,I8,I8,I8,I8,I8,I8)       
 9998 FORMAT(A,I8,I8,I8)      


      RETURN
*     
*     End of PDHGEQZ0
*     
      END SUBROUTINE
