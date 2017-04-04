      SUBROUTINE PDTGORD( WANTQ, WANTZ, SELECT, PARA, N, S,  
     $   DESCS, T, DESCT, Q, DESCQ, Z,  
     $   DESCZ, ALPHAR, ALPHAI, BETA, M, DWORK, LDWORK, 
     $   IWORK, LIWORK, INFO )
C     
C     -- ScaLAPACK-style routine --
C     Preliminary version.
C     Dept. Computing Science and HPC2N, Univ. of Umeå, Sweden
C     written by Robert Granat, June 7, 2007.
C     Modified by Bjorn Adlerborn, Dec 2013. 
      IMPLICIT NONE
C     
C     .. Scalar Arguments ..
      LOGICAL            WANTQ, WANTZ
      INTEGER            INFO, LIWORK, LDWORK, M, N 
C     ..
C     .. Array Arguments ..
      INTEGER            SELECT( * )
      INTEGER            PARA( 6 ), DESCS( 9 ), DESCT( 9 ), DESCQ( 9 ),
     $   DESCZ( 9 ), IWORK( * )
      DOUBLE PRECISION   Q( * ), Z( * ), S( * ), T( * ), BETA( * ),
     $   ALPHAI( * ), DWORK( * ), ALPHAR( * )
C     ..
C     
C     Purpose
C     =======
C     
C     PBDTGSEN reorders the generalized real Schur decomposition of a real
C     matrix pair (S, T) (in terms of an orthonormal equivalence trans-
C     formation Q' * (S, T) * Z), so that a selected cluster of eigenvalues
C     appears in the leading diagonal blocks of the upper quasi-triangular
C     matrix S and the upper triangular T. The leading columns of Q and
C     Z form orthonormal bases of the corresponding left and right eigen-
C     spaces (deflating subspaces). (S, T) must be in generalized real
C     Schur canonical form, i.e. S is block upper triangular with 1-by-1 and 
C     2-by-2 diagonal blocks. T is upper triangular.
C     
C     PBDTGSEN also computes the generalized eigenvalues
C     
C     w(j) = (ALPHAR(j) + i*ALPHAI(j))/BETA(j)
C     
C     of the reordered matrix pair (S, T).
C     
C     Optionally, the routine computes the estimates of reciprocal condition
C     numbers for eigenvalues and eigenspaces. These are Difu[(S11,T11),
C     (S22,T22)] and Difl[(S11,T11), (S22,T22)], i.e. the separation(s)
C     between the matrix pairs (S11, T11) and (S22,T22) that correspond to
C     the selected cluster and the eigenvalues outside the cluster, resp.,
C     and norms of "projections" onto left and right eigenspaces w.r.t.
C     the selected cluster in the (1,1)-block.
C     
C     This subroutine uses a delay and accumulate procedure for performing
C     the off-diagonal updates (see references for details).
C     
C     Notes
C     =====
C     
C     Each global data object is described by an associated description
C     vector.  This vector stores the information required to establish
C     the mapping between an object element and its corresponding process
C     and memory location.
C     
C     Let A be a generic term for any 2D block cyclicly distributed array.
C     Such a global array has an associated description vector DESCA.
C     In the following comments, the character _ should be read as
C     "of the global array".
C     
C     NOTATION        STORED IN      EXPLANATION
C     --------------- -------------- --------------------------------------
C     DTYPE_A(global) DESCA( DTYPE_ )The descriptor type.  In this case,
C     DTYPE_A = 1.
C     CTXT_A (global) DESCA( CTXT_ ) The BLACS context handle, indicating
C     the BLACS process grid A is distribu-
C     ted over. The context itself is glo-
C     bal, but the handle (the integer
C     value) may vary.
C     M_A    (global) DESCA( M_ )    The number of rows in the global
C     array A.
C     N_A    (global) DESCA( N_ )    The number of columns in the global
C     array A.
C     MB_A   (global) DESCA( MB_ )   The blocking factor used to distribute
C     the rows of the array.
C     NB_A   (global) DESCA( NB_ )   The blocking factor used to distribute
C     the columns of the array.
C     RSRC_A (global) DESCA( RSRC_ ) The process row over which the first
C     row of the array A is distributed.
C     CSRC_A (global) DESCA( CSRC_ ) The process column over which the
C     first column of the array A is
C     distributed.
C     LLD_A  (local)  DESCA( LLD_ )  The leading dimension of the local
C     array.  LLD_A >= MAX(1,LOCr(M_A)).
C     
C     Let K be the number of rows or columns of a distributed matrix,
C     and assume that its process grid has dimension p x q.
C     LOCr( K ) denotes the number of elements of K that a process
C     would receive if K were distributed over the p processes of its
C     process column.
C     Similarly, LOCc( K ) denotes the number of elements of K that a
C     process would receive if K were distributed over the q processes of
C     its process row.
C     The values of LOCr() and LOCc() may be determined via a call to the
C     ScaLAPACK tool function, NUMROC:
C     LOCr( M ) = NUMROC( M, MB_A, MYROW, RSRC_A, NPROW ),
C     LOCc( N ) = NUMROC( N, NB_A, MYCOL, CSRC_A, NPCOL ).
C     An upper bound for these quantities may be computed by:
C     LOCr( M ) <= ceil( ceil(M/MB_A)/NPROW )*MB_A
C     LOCc( N ) <= ceil( ceil(N/NB_A)/NPCOL )*NB_A
C     
C     Arguments
C     =========
C     
C     IJOB    (global input) INTEGER
C     Specifies whether condition numbers are required for the
C     cluster of eigenvalues (PL and PR) or the deflating subspaces
C     (Difu and Difl):
C     =0: Only reorder w.r.t. SELECT. No extras.
C     =1: Reciprocal of norms of "projections" onto left and right
C     eigenspaces w.r.t. the selected cluster (PL and PR).
C     =2: Upper bounds on Difu and Difl. F-norm-based estimate
C     (DIF(1:2)).
C     =3: Estimate of Difu and Difl. 1-norm-based estimate
C     (DIF(1:2)).
C     About 5 times as expensive as IJOB = 2.
C     =4: Compute PL, PR and DIF (i.e. 0, 1 and 2 above): Economic
C     version to get it all.
C     =5: Compute PL, PR and DIF (i.e. 0, 1 and 3 above)
C     
C     WANTQ   (global input) LOGICAL
C     .TRUE. : update the left transformation matrix Q;
C     .FALSE.: do not update Q.
C     
C     WANTZ   (global input) LOGICAL
C     .TRUE. : update the right transformation matrix Z;
C     .FALSE.: do not update Z.
C     
C     SELECT  (global input/output) INTEGER/LOGICAL  array, dimension (N)
C     SELECT specifies the eigenvalues in the selected cluster. To
C     select a real eigenvalue w(j), SELECT(j) must be set to
C     .TRUE. (1). To select a complex conjugate pair of eigenvalues
C     w(j) and w(j+1), corresponding to a 2-by-2 diagonal block,
C     either SELECT(j) or SELECT(j+1) or both must be set to
C     .TRUE. (1); a complex conjugate pair of eigenvalues must be
C     either both included in the cluster or both excluded.
C     On output, the (partial) reordering is displayed.
C     
C     PARA    (global input) INTEGER*6
C     Block parameters (some should be replaced by calls to
C     PILAENV and others by meaningful default values):
C     PARA(1) = maximum number of concurrent computational windows
C     allowed in the algorithm; 
C     0 < PARA(1) <= min(NPROW,NPCOL) must hold;
C     PARA(2) = number of eigenvalues in each window; 
C     0 < PARA(2) < PARA(3) must hold; 
C     PARA(3) = window size; PARA(2) < PARA(3) < DESCT(MB_) 
C     must hold;
C     PARA(4) = minimal percentage of flops required for
C     performing matrix-matrix multiplications instead
C     of pipelined orthogonal transformations;
C     0 <= PARA(4) <= 100 must hold;
C     PARA(5) = width of block column slabs for row-wise
C     application of pipelined orthogonal
C     transformations in their factorized form; 
C     0 < PARA(5) <= DESCT(MB_) must hold.
C     PARA(6) = the maximum number of eigenvalues moved together
C     over a process border; in practice, this will be
C     approximately half of the cross border window size
C     0 < PARA(6) <= PARA(2) must hold; 
C     
C     N       (global input) INTEGER
C     The order of the globally distributed matrices S and T. N >= 0.
C     
C     S       (local input/output) DOUBLE PRECISION array, 
C     dimension (LLD_S,LOCc(N)).
C     On entry, the local pieces of the global distributed 
C     upper quasi-triangular matrix S, in Schur form. On exit, S is 
C     overwritten by the local pieces of the reordered matrix S, 
C     again in Schur form, with the selected eigenvalues in the 
C     globally leading diagonal blocks.
C     
C     IS      (global input) INTEGER
C     JS      (global input) INTEGER
C     The row and column index in the global array T indicating the
C     first column of sub( S ). IS = JS = 1 must hold. 
C     
C     DESCS   (global and local input) INTEGER array of dimension DLEN_.
C     The array descriptor for the global distributed matrix T.
C     
C     
C     T       (local input/output) DOUBLE PRECISION array, 
C     dimension (LLD_T,LOCc(N)).
C     On entry, the local pieces of the global distributed 
C     upper triangular matrix T, in Schur form. On exit, T is 
C     overwritten by the local pieces of the reordered matrix T, 
C     again in Schur form, with the selected eigenvalues in the 
C     globally leading diagonal blocks.
C     
C     IT      (global input) INTEGER
C     JT      (global input) INTEGER
C     The row and column index in the global array T indicating the
C     first column of sub( T ). IT = JT = 1 must hold.
C     
C     DESCT   (global and local input) INTEGER array of dimension DLEN_.
C     The array descriptor for the global distributed matrix T.
C     
C     Q       (local input/output) DOUBLE PRECISION array, 
C     dimension (LLD_Q,LOCc(N)).
C     On entry, if WANTQ = .TRUE., the local pieces of the global
C     distributed matrix Q of generalized Schur vectors.
C     On exit, WANTQ = .TRUE., Q has been postmultiplied by the
C     global orthogonal transformation matrix which reorders (S, T); 
C     the leading M columns of Q and Z form orthonormal bases for the
C     specified deflating subspaces.
C     If WANTQ = .FALSE., Q is not referenced.
C     
C     IQ      (global input) INTEGER
C     JQ      (global input) INTEGER
C     The column index in the global array Q indicating the
C     first column of sub( Q ). IQ = JQ = 1 must hold.
C     
C     DESCQ   (global and local input) INTEGER array of dimension DLEN_.
C     The array descriptor for the global distributed matrix Q.
C     
C     Z       (local input/output) DOUBLE PRECISION array, 
C     dimension (LLD_Z,LOCc(N)).
C     On entry, if WANTZ = .TRUE., the local pieces of the global
C     distributed matrix Z of generalized Schur vectors.
C     On exit, WANTZ = .TRUE., Z has been postmultiplied by the
C     global orthogonal transformation matrix which reorders (S, T); 
C     the leading M columns of Q and Z form orthonormal bases for the
C     specified deflating subspaces.
C     If WANTQ = .FALSE., Z is not referenced.
C     
C     IZ      (global input) INTEGER
C     JZ      (global input) INTEGER
C     The column index in the global array Z indicating the
C     first column of sub( Z ). IZ = JZ = 1 must hold.
C     
C     DESCZ   (global and local input) INTEGER array of dimension DLEN_.
C     The array descriptor for the global distributed matrix Z.
C     
C     ALPHAR  (global output) DOUBLE PRECISION array, dimension (N)
C     ALPHAI  (global output) DOUBLE PRECISION array, dimension (N)
C     BETA    (global output) DOUBLE PRECISION array, dimension (N)
C     On exit, (ALPHAR(j) + ALPHAI(j)*i)/BETA(j), j=1,...,N, will
C     be the generalized eigenvalues.  ALPHAR(j) + ALPHAI(j)*i
C     and BETA(j),j=1,...,N  are the diagonals of the complex Schur
C     form (S,T) that would result if the 2-by-2 diagonal blocks of
C     the real generalized Schur form of (A,B) were further reduced
C     to triangular form using complex unitary transformations.
C     If ALPHAI(j) is zero, then the j-th eigenvalue is real; if
C     positive, then the j-th and (j+1)-st eigenvalues are a
C     complex conjugate pair, with ALPHAI(j+1) negative.
C     Note also that if a complex eigenvalue is sufficiently 
C     ill-conditioned, then its value may differ significantly 
C     from its value before reordering. 
C     
C     M       (global output) INTEGER
C     The dimension of the specified pair of left and right eigen-
C     spaces (deflating subspaces). 0 <= M <= N.
C     
C     PL, PR  (global output) DOUBLE PRECISION
C     If IJOB = 1, 4 or 5, PL, PR are lower bounds on the
C     reciprocal of the norm of "projections" onto left and right
C     eigenspaces with respect to the selected cluster.
C     0 < PL, PR <= 1.
C     If M = 0 or M = N, PL = PR  = 1.
C     If IJOB = 0, 2 or 3, PL and PR are not referenced.
C     
C     DIF     (global output) DOUBLE PRECISION array, dimension (2).
C     If IJOB >= 2, DIF(1:2) store the estimates of Difu and Difl.
C     If IJOB = 2 or 4, DIF(1:2) are F-norm-based upper bounds on
C     Difu and Difl. If IJOB = 3 or 5, DIF(1:2) are 1-norm-based
C     estimates of Difu and Difl.
C     If M = 0 or N, DIF(1:2) = F-norm([S, T]).
C     If IJOB = 0 or 1, DIF is not referenced.
C     
C     DWORK   (local workspace/output) DOUBLE PRECISION array, 
C     dimension (LDWORK) 
C     On exit, if INFO = 0, DWORK(1) returns the optimal LDWORK.
C     
C     LDWORK  (local input) INTEGER
C     The dimension of the array DWORK.
C     
C     If LDWORK = -1, then a workspace query is assumed; the routine
C     only calculates the optimal size of the DWORK array, returns
C     this value as the first entry of the DWORK array, and no error
C     message related to LDWORK is issued by PXERBLA.
C     
C     IWORK   (local workspace/output) INTEGER array, dimension (LIWORK)
C     
C     LIWORK  (local input) INTEGER
C     The dimension of the array IWORK.
C     
C     If LIWORK = -1, then a workspace query is assumed; the
C     routine only calculates the optimal size of the IWORK array,
C     returns this value as the first entry of the IWORK array, and
C     no error message related to LIWORK is issued by PXERBLA.
C     
C     INFO    (global output) INTEGER
C     = 0: successful exit
C     < 0: if INFO = -i, the i-th argument had an illegal value.
C     If the i-th argument is an array and the j-entry had
C     an illegal value, then INFO = -(i*1000+j), if the i-th
C     argument is a scalar and had an illegal value, then INFO = -i.
C     > 0: here we have several possibilites
C     *) Reordering of (S,T) failed because some eigenvalues are too
C     close to separate (the problem is very ill-conditioned);
C     (S,T) may have been partially reordered, and the output
C     eigenvalues are in the same order as in (S,T); PL, PR and
C     DIF (if requested) are set to zero. The process that 
C     failed in the reordering will return INFO = {the index 
C     of T where the swap failed}; all others will return 
C     INFO = 1.
C     *) A 2-by-2 block to be reordered split into two 1-by-1
C     blocks and the second block failed to swap with an
C     adjacent block. The process that failed in the 
C     reordering will return INFO = {the index of T where the 
C     swap failed}; all others will return INFO = 1.
C     *) If INFO = 2, the routines used in the calculation of the
C     condition numbers raised a positive warning flag (see
C     the documentation for PGEGCSYD and PGCSYCON of the
C     SCASY library).
C     *) If INFO = 3, there is no valid BLACS context (see the
C     BLACS documentation for details).
C     *) If INFO = 333, PGEGCSYD raised an input error flag;
C     please report this bug to the authors (see below).
C     If INFO = 444, PGCSYCON raised an input error flag;
C     please report this bug to the authors (see below).
C     In a future release this subroutine may distinguish between
C     the case 1 and 2 above.  
C     
C     Method
C     ======
C     
C     This routine performs parallel eigenvalue reordering in gen. Schur
C     form by parallelizing the approach proposed in [3]. The condition
C     number estimation part is performed by using techniques and code
C     from SCASY, see http://www.cs.umu.se/research/parallel/scasy.
C     
C     Additional requirements
C     =======================
C     
C     The following alignment requirements must hold:
C     (a) DESC[S,T]( MB_ ) = DESC[S,T]( NB_ ) = DESC[Q,Z]( MB_ ) = 
C     DESC[Q,Z]( NB_ )
C     (b) DESC[S,T]( RSRC_ ) = DESC[Q,Z]( RSRC_ )
C     (c) DESC[S,T]( CSRC_ ) = DESC[Q,Z]( CSRC_ ) 
C     
C     All matrices must be blocked by a block factor larger than or
C     equal to two (3). This to simplify reordering across processor
C     borders in the presence of 2-by-2 blocks.
C     
C     Limitations
C     ===========
C     
C     This algorithm cannot work on submatrices of S, T, Q and Z i.e.,
C     IS = JS = IT = JT = IQ = JQ = IZ = JZ = 1 must hold.
C     
C     References
C     ==========
C     
C     [1] B. Kagstrom; A Direct Method for Reordering Eigenvalues in the
C     Generalized Real Schur Form of a Regular Matrix Pair (A, B), in
C     M.S. Moonen et al (eds), Linear Algebra for Large Scale and
C     Real-Time Applications, Kluwer Academic Publ. 1993, pp 195-218.
C     
C     [2] B. Kagstrom and P. Poromaa; Computing Eigenspaces with Specified
C     Eigenvalues of a Regular Matrix Pair (A, B) and Condition
C     Estimation: Theory, Algorithms and Software,
C     Report UMINF - 94.04, Department of Computing Science, Umea
C     University, S-901 87 Umea, Sweden, 1994. Also as LAPACK Working
C     Note 87. To appear in Numerical Algorithms, 1996.
C     
C     [3] B. Kagstrom and P. Poromaa, LAPACK-Style Algorithms and Software
C     for Solving the Generalized Sylvester Equation and Estimating the
C     Separation between Regular Matrix Pairs, Report UMINF - 93.23,
C     Department of Computing Science, Umea University, S-901 87 Umea,
C     Sweden, December 1993, Revised April 1994, Also as LAPACK Working
C     Note 75. To appear in ACM Trans. on Math. Software, Vol 22, No 1,
C     1996.
C     
C     [4] D. Kressner; Block algorithms for reordering standard and
C     generalized Schur forms, 2005. LAPACK Working Note ???.
C     
C     [5] R. Granat, B. Kågström, D. Kressner; Parallel eigenvalue 
C     reordering in real Schur form, 2006. LAPACK Working Note ???.
C     
C     Parallel execution recommendations
C     ==================================
C     
C     Use a square grid, if possible, for maximum performance. The block
C     parameters in PARA should be kept well below the data distribution
C     block size. In particular, see [4,5] for recommended settings for
C     these parameters.
C     
C     In general, the parallel algorithm strives to perform as much work 
C     as possible without crossing the block borders on the main block 
C     diagonal.
C     
C     Contributors
C     ============
C     
C     Implemented by Robert Granat, Umea University and HPC2N, June. 2007.
C     
C     Revisions
C     =========
C     
C     Please send bug-reports to granat@cs.umu.se
C     
C     Keywords
C     ========
C     
C     Generalized Schur form, eigenvalue reordering, Generalized coupled 
C     Sylvester matrix equation
C     
C     =====================================================================
C     ..
C     .. Parameters ..
      INTEGER            BLOCK_CYCLIC_2D, CSRC_, CTXT_, DLEN_, DTYPE_,
     $   LLD_, MB_, M_, NB_, N_, RSRC_
      DOUBLE PRECISION   ZERO, ONE
      PARAMETER          ( BLOCK_CYCLIC_2D = 1, DLEN_ = 9, DTYPE_ = 1,
     $   CTXT_ = 2, M_ = 3, N_ = 4, MB_ = 5, NB_ = 6,
     $   RSRC_ = 7, CSRC_ = 8, LLD_ = 9,
     $   ZERO = 0.0D+0, ONE = 1.0D+0 ) 
C     ..
C     .. Local Scalars ..
      LOGICAL            LQUERY, PAIR, SWAP, 
     $   ISHH, FIRST, SKIP1CR, BORDER,
     $   LASTWAIT
      INTEGER            NPROW, NPCOL, MYROW, MYCOL, NB, NPROCS, IERR,
     $   DIM1, INDX, IDUM1(2), IDUM2(2), LLDT, TRSRC, TCSRC, ILOC1,
     $   JLOC1, MYIERR, ICTXT, RSRC1, CSRC1, ILOC3, JLOC3, ILOC, JLOC,
     $   FLOPS, I, ILO, IHI, J, K, KK, KKS, KS, LIWMIN, LWMIN, MMULT,
     $   NCB, NDTRAF, NITRAF, NWIN, NUMWIN, PDTRAF, PITRAF, WINEIG,
     $   WINSIZ, LLDQ, RSRC, CSRC, ILILO, ILIHI, ILSEL, IRSRC, ICSRC,
     $   IPIW, IPW1, IPW2, IPW3, TIHI, TILO, LIHI, WINDOW, LILO, LSEL,
     $   BUFFER, NMWIN2, BUFFLEN, LROWS, LCOLS, ILOC2, JLOC2, WNEICR,
     $   WINDOW0, RSRC4, CSRC4, LIHI4, RSRC3, CSRC3, RSRC2, CSRC2, LIHIC
     $   , LIHI1, ILEN4, SELI4, ILEN1, DIM4, IPW4, QZROWS, STROWS,
     $   STCOLS, IPW5, IPW6, IPW7, IPW8, JLOC4, EAST, WEST, ILOC4, SOUTH
     $   , NORTH, INDXS, ITT, JTT, WRK1, IWRK1, ILEN, DLEN, INDXE, ILOS,
     $   DIR, TLIHI, TLILO, TLSEL, ROUND, LAST, WIN0S, WIN0E, WINE, LLDS
     $   , LLDZ, PDU, PDV, IAM

      DOUBLE PRECISION   ELEM, EPS, SMLNUM, BETA2
C     ..
C     .. Local Arrays ..
      INTEGER            IBUFF( 8 )
C     ..
C     .. External Functions ..
      LOGICAL            LSAME,INT2LG
      INTEGER            NUMROC, INDXG2P, INDXG2L
      DOUBLE PRECISION   PDLANGE, MPI_WTIME
      EXTERNAL           LSAME, PDLANGE, NUMROC, INDXG2P, INDXG2L,
     $   MPI_WTIME,INT2LG
C     ..
C     .. External Subroutines ..
      EXTERNAL           PDLACPY, PXERBLA, DGEMM, DLACPY, ILACPY,
     $   CHK1MAT, PCHK2MAT, PCHK1MAT,  INFOG2L, DGSUM2D, DGESD2D,
     $   DGERV2D, DGEBS2D, DGEBR2D, IGSUM2D, BLACS_GRIDINFO, IGEBS2D,
     $   IGEBR2D, IGAMX2D, IGAMN2D
C     ..
C     .. Intrinsic Functions ..
      INTRINSIC          ABS, MAX, SQRT, MIN
C     ..
C     .. Local Functions .. 
      INTEGER            ICEIL


      LOGICAL            DEBUG
      INTEGER            RRCNT
C     .. 
C     .. Executable Statements ..
C     
C     Get grid parameters 
C     
      RRCNT = 1
      ICTXT = DESCT( CTXT_ )
      CALL BLACS_GRIDINFO( ICTXT, NPROW, NPCOL, MYROW, MYCOL )
      CALL BLACS_PINFO(IAM, NPROCS)
      NB = DESCS( MB_ )
      DEBUG = .FALSE.

C     
C     Test if grid is O.K., i.e., the context is valid
C     
      ISHH = .TRUE.
      INFO = 0



      LQUERY = LDWORK.EQ.-1.OR.LIWORK.EQ.-1
      SMLNUM = 0
      EPS = 0

      NDTRAF = 0
      PITRAF = 0
      PDTRAF = 0
      IERR = 0
      IPW3 = 0
      SKIP1CR = .FALSE.
C     
C     Decode and test the input parameters
C     
C     
      IF( INFO .EQ. 0 ) THEN
C        

C        
C        Extract local leading dimension
C        
         LLDS = DESCS( LLD_ )
         LLDT = DESCT( LLD_ )
         LLDQ = DESCQ( LLD_ )
         LLDZ = DESCZ( LLD_ )
C        
C        Check the SELECT vector for consistency and set M to the dimension 
C        of the specified invariant subspace.
C        
         M = 0
         IF ( .NOT. LQUERY) THEN
            DO 10 K = 1, N
               IF( K.LT.N ) THEN
                  CALL INFOG2L( K+1, K, DESCS, NPROW, NPCOL, 
     $               MYROW, MYCOL, ITT, JTT, TRSRC, TCSRC )
                  IF( MYROW.EQ.TRSRC .AND. MYCOL.EQ.TCSRC ) THEN
                     IF( S( (JTT-1)*LLDS + ITT ).NE.ZERO ) THEN
                        IF( INT2LG(SELECT(K)).AND. .NOT.INT2LG(SELECT(K
     $                     +1)) ) THEN
                           SELECT(K+1) = 1
                        ELSEIF( .NOT.INT2LG(SELECT(K)).AND.
     $                        INT2LG(SELECT(K+1)) ) THEN
                           SELECT(K) = 1
                        END IF
                     END IF
                  END IF
               END IF
               IF( SELECT(K).EQ.1 ) M = M + 1
 10         CONTINUE
            
            IF( NPROCS.GT.1 ) THEN
               CALL IGAMX2D( ICTXT, 'All', ' ', N, 1, 
     $            SELECT, 1, -1, -1, 
     $            -1, -1, -1 )  
               CALL IGAMX2D( ICTXT, 'All', ' ', 1, 1, M, 1, 
     $            -1, -1, -1, 
     $            -1, -1 )           
            END IF
         END IF


C        Compute needed workspace
C        
C        
         STROWS = NUMROC( N, NB, MYROW, DESCS(RSRC_), NPROW )
         STCOLS = NUMROC( N, NB, MYCOL, DESCS(CSRC_), NPCOL )
         WRK1 = 2*N + 14*NB**2 + 4*STROWS*PARA( 3 ) + 
     $      2*STCOLS*PARA( 3 ) + MAX( 2*STROWS*PARA( 3 ), 
     $      2*STCOLS*PARA( 3 ) )
         call flush(6)
         IWRK1 = 5*PARA( 1 ) + PARA(2)*PARA(3) - 
     $      PARA(2) * (PARA(2) + 1 ) / 2
C        
         LWMIN = MAX( 1, WRK1  )
         LIWMIN = MAX( 1, IWRK1 )
C        
         IF( LDWORK.LT.LWMIN .AND. .NOT.LQUERY ) THEN
            INFO = -31
         ELSE IF( LIWORK.LT.LIWMIN .AND. .NOT.LQUERY ) THEN
            INFO = -33
         END IF
         
      END IF
C     
C     Global maximum on info
C     
      IF( NPROCS.GT.1 )
     $   CALL IGAMX2D( ICTXT, 'All', ' ', 1, 1, INFO, 1, -1, -1, -1, 
     $   -1, -1 )
C     
C     Return if some argument is incorrect
C     
      IF( INFO .NE. 0 ) THEN
         CALL PXERBLA( ICTXT, 'PDTGORD', -INFO )
         RETURN
      ELSEIF( LQUERY ) THEN
         DWORK( 1 ) = DBLE(LWMIN)
         IWORK( 1 ) = LIWMIN
         RETURN
      END IF
C     
C     Quick return if possible.
C     
      IF( M.EQ.N .OR. M.EQ.0 ) THEN
         GO TO 545
      END IF
C     
C     Set parameters
C     

      NUMWIN = PARA( 1 )
      WINEIG = MAX( PARA( 2 ), 2 )
      WINSIZ = MIN( MAX( PARA( 3 ), PARA( 2 )*2 ), NB )
      MMULT  = PARA( 4 )
      NCB    = PARA( 5 )
      WNEICR = PARA( 6 )
      IF ( DEBUG ) WRITE(*,*)'Reorder:Start', NUMWIN, WINEIG, MMULT, 
     $   NCB, WNEICR
C     
C     Insert some pointers into INTEGER workspace
C     
      ILILO = 1
      ILIHI = ILILO + NUMWIN
      ILSEL = ILIHI + NUMWIN
      IRSRC = ILSEL + NUMWIN
      ICSRC = IRSRC + NUMWIN
      IPIW  = ICSRC + NUMWIN
C     
C     Insert some pointers into DOUBLE PRECISION workspace - for now we
C     only need two pointers
C     
      IPW1 = 1
      IPW2 = IPW1 + 4*N+16
C     
C     Collect the selected blocks at the top-left corner of T.
C     
C     Globally: ignore eigenvalues that are already in order.
C     ILO is a global variable and is kept updated to be consistent
C     throughout the process mesh.
C     
      ILO = 0
 40   CONTINUE
      ILO = ILO + 1
      IF( ILO.LE.N ) THEN
         IF( INT2LG(SELECT(ILO)) ) GO TO 40
      END IF
C     
C     Globally: start the collection at the top of the matrix. Here,
C     IHI is a global variable and is kept updated to be consistent
C     throughout the process mesh.
C     
      IHI = N
C     
C     Globally:  While ( ILO <= M ) do
 50   CONTINUE
      IF ( DEBUG ) WRITE(*,*)'Reorder: loop', ILO, M
      IF( ILO.LE.M ) THEN
C        
C        Depending on the value of ILO, find the diagonal block index J,
C        such that (S,T)(1+(J-1)*NB:1+J*NB,1+(J-1)*NB:1+J*NB) contains the
C        first unsorted eigenvalue. Check that J does not point to a block
C        with only one selected eigenvalue in the last position which
C        belongs to a splitted 2-by-2 block.
C        
         ILOS = ILO - 1
 52      CONTINUE
         ILOS = ILOS + 1
         IF( .NOT. INT2LG(SELECT(ILOS)) ) GO TO 52
         IF( ILOS.LT.N ) THEN
            IF( INT2LG(SELECT( ILOS+1 )) .AND. MOD(ILOS,NB).EQ.0 ) THEN
               CALL PDELGET( 'All', ' ', ELEM, S, ILOS+1, ILOS, DESCS )
               IF( ELEM.NE.ZERO ) GO TO 52
            END IF
         END IF
         J = ICEIL(ILOS,NB)
C        
C        Globally: Set startvalues of LILO and LIHI for all processes. 
C        Choose also the number of selected eigenvalues at top of each
C        diagonal block such that the number of eigenvalues which remain 
C        to be reordered is an integer multiple of WINEIG. 
C        
C        All the information is saved into the INTEGER workspace such that 
C        all processors are aware of each others operations.  
C        
C        Compute the number of concurrent windows
C        
         NMWIN2 = (ICEIL(IHI,NB)*NB - (ILO-MOD(ILO,NB)+1)+1) / NB
         NMWIN2 = MIN( MIN( NUMWIN, NMWIN2 ), ICEIL(N,NB) - J + 1 )
C        
C        For all windows, set LSEL = 0 and find a proper start value of
C        LILO such that LILO points at the first non-selected entry in
C        the corresponding diagonal block of (S,T)
C        
         DO 80 K = 1, NMWIN2
            IWORK( ILSEL+K-1) = 0
            IWORK( ILILO+K-1) = MAX( ILO, (J-1)*NB+(K-1)*NB+1 )
            LILO = IWORK( ILILO+K-1 )
 82         CONTINUE
            IF( INT2LG(SELECT(LILO)) .AND. LILO.LT.(J+K-1)*NB ) THEN
               LILO = LILO + 1
               IF( LILO.LE.N ) GO TO 82
            END IF
            IWORK( ILILO+K-1 ) = LILO
C           
C           Fix each LILO to avoid that no 2-by-2 block is cut in top of the 
C           submatrices (LILO:LIHI,LILO:LIHI)
C           
            LILO = IWORK(ILILO+K-1)

            IF( LILO.GT.NB ) THEN
               CALL PDELGET( 'All', ' ', ELEM, S, LILO, LILO-1, DESCS ) 
               IF( ELEM.NE.ZERO ) THEN
                  IF( LILO.LT.(J+K-1)*NB ) THEN
                     IWORK(ILILO+K-1) = IWORK(ILILO+K-1) + 1
                  ELSE
                     IWORK(ILILO+K-1) = IWORK(ILILO+K-1) - 1
                  END IF 
               END IF
            END IF
C           
C           Set a proper LIHI value for each window. Also find the processors
C           corresponding to the corresponding windows.
C           
            IWORK( ILIHI+K-1 ) =  IWORK( ILILO+K-1 )

            IWORK( IRSRC+K-1 ) = INDXG2P( IWORK(ILILO+K-1), NB, MYROW, 
     $         DESCS( RSRC_ ), NPROW )
            IWORK( ICSRC+K-1 ) = INDXG2P( IWORK(ILILO+K-1), NB, MYCOL, 
     $         DESCS( CSRC_ ), NPCOL )
            TILO = IWORK(ILILO+K-1)
            TIHI = MIN( N, ICEIL( TILO, NB ) * NB )
            DO 90 KK = TIHI, TILO, -1
               IF( INT2LG(SELECT( KK )) ) THEN
                  IWORK(ILIHI+K-1) = MAX(IWORK(ILIHI+K-1) , KK )

                  IWORK(ILSEL+K-1) = IWORK(ILSEL+K-1) + 1
                  IF( IWORK(ILSEL+K-1).GT.WINEIG ) THEN
                     IWORK(ILIHI+K-1) = KK
                     IWORK(ILSEL+K-1) = 1
                  END IF
               END IF
 90         CONTINUE
C           
C           Fix each LIHI to avoid that bottom of window cuts 2-by-2 block. 
C           We exclude such a block if located on block (process) border and 
C           on window border or if an inclusion would cause violation on the 
C           maximum number of eigenvalues to reorder inside each window. If 
C           only on window border, we include it. The excluded block is 
C           included automatically later when a subcluster is reordered into 
C           the block from South-East.
C           
            LIHI = IWORK(ILIHI+K-1)
            IF( LIHI.LT.N ) THEN
               CALL PDELGET( 'All', ' ', ELEM, S, LIHI+1, LIHI, DESCS )
               IF( ELEM.NE.ZERO ) THEN
                  IF( ICEIL( LIHI, NB ) .NE. ICEIL( LIHI+1, NB ) .OR. 
     $               IWORK( ILSEL+K-1 ).EQ.WINEIG ) THEN
                     IWORK( ILIHI+K-1 ) = IWORK( ILIHI+K-1 ) - 1

                     IF( IWORK( ILSEL+K-1 ).GT.2 ) IWORK( ILSEL+K-1 ) =
     $                  IWORK( ILSEL+K-1 ) - 1
                  ELSE
                     IWORK( ILIHI+K-1 ) = IWORK( ILIHI+K-1 ) + 1
                     IF( INT2LG(SELECT(LIHI+1)) ) IWORK( ILSEL+K-1 ) =
     $                  IWORK( ILSEL+K-1 ) + 1
                  END IF
               END IF
            END IF
 80      CONTINUE
C        
C        Fix the special cases of LSEL = 0 and LILO = LIHI for each window 
C        by assuring that the stop-condition for local reordering is 
C        fulfilled directly. Do this by setting LIHI = startposition for the 
C        corresponding block and LILO = LIHI + 1.
C        
         DO 85 K = 1, NMWIN2
            LILO = IWORK( ILILO + K - 1 )
            LIHI = IWORK( ILIHI + K - 1 )
            LSEL = IWORK( ILSEL + K - 1 )
            IF( LSEL.EQ.0 .OR. LILO.EQ.LIHI ) THEN
               LIHI = IWORK( ILIHI + K - 1 )
               IWORK( ILIHI + K - 1 ) = (ICEIL(LIHI,NB)-1)*NB + 1
               IWORK( ILILO + K - 1 ) = IWORK( ILIHI + K - 1 ) + 1
            END IF
 85      CONTINUE


C        

C        
C        Associate all processors with the first computational window 
C        that should be activated, if possible. 
C        
         LILO = IHI
         LIHI = ILO
         LSEL = M
         FIRST = .TRUE.
         DO 95 WINDOW = 1, NMWIN2
            RSRC = IWORK(IRSRC+WINDOW-1)
            CSRC = IWORK(ICSRC+WINDOW-1)
            IF( MYROW.EQ.RSRC .OR. MYCOL.EQ.CSRC ) THEN
               TLILO = IWORK( ILILO + WINDOW - 1 )
               TLIHI = IWORK( ILIHI + WINDOW - 1 )
               TLSEL = IWORK( ILSEL + WINDOW - 1 )

               IF( (.NOT. ( LIHI .GE. LILO + LSEL ) ) .AND.
     $            ( (TLIHI .GE. TLILO + TLSEL) .OR. FIRST ) ) THEN
                  IF( FIRST ) FIRST = .FALSE.
                  LILO = TLILO
                  LIHI = TLIHI
                  LSEL = TLSEL                  
                  GO TO 97
               END IF
            END IF
 95      CONTINUE
 97      CONTINUE

         
C        
C        Exclude all processors that are not involved in any computational
C        window right now
C        
         IF( LILO.EQ.IHI .AND. LIHI.EQ.ILO .AND. LSEL.EQ.M ) THEN

            GO TO 114
         END IF
C        
C        Make sure all processors associated with a compuational window
C        enter the local reordering the first time 
C        
         FIRST = .TRUE.

C        Globally for all computational windows:
C        While ( LIHI >= LILO + LSEL ) do
         ROUND = 1
 130     CONTINUE
         IF( FIRST .OR. ( LIHI .GE. LILO + LSEL ) ) THEN
            IF( ROUND.GT.20 ) THEN
               STOP
            END IF
C           
C           Perform computations in parallel: loop thorugh all compuational
C           windows, do local reordering and accumulate transformations,
C           broadcast them in the corresponding block row and columns and
C           compute the corresponding updates
C           

            DO 110 WINDOW = 1, NMWIN2
               RSRC = IWORK(IRSRC+WINDOW-1)
               CSRC = IWORK(ICSRC+WINDOW-1)
C              
C              The process on the block diagonal computes the reordering 
C              
               IF( MYROW.EQ.RSRC .AND. MYCOL.EQ.CSRC ) THEN
                  LILO = IWORK(ILILO+WINDOW-1)
                  LIHI = IWORK(ILIHI+WINDOW-1)
                  LSEL = IWORK(ILSEL+WINDOW-1)
C                 
C                 Compute the local value of I
C                 
                  I = MAX( LILO, LIHI - WINSIZ + 1 )
C                 
C                 Fix my I to avoid that top of window cuts 2-by-2 block.
C                 
                  IF( I.GT.LILO ) THEN
                     CALL INFOG2L( I, I-1, DESCS, NPROW, NPCOL, MYROW,
     $                  MYCOL, ILOC, JLOC, RSRC, CSRC )
                     IF( S( LLDS*(JLOC-1) + ILOC ).NE.ZERO )   
     $                  I = I + 1
                  END IF
C                 
C                 Compute local indicies for submatrix to operate on
C                 
                  CALL INFOG2L( I, I, DESCS, NPROW, NPCOL, 
     $               MYROW, MYCOL, ILOC1, JLOC1, RSRC, 
     $               CSRC )
C                 
C                 The active window is ( I:LIHI, I:LIHI ). Reorder eigenvalues
C                 within this window and pipeline transformations.
C                 

                  NWIN = LIHI - I + 1
                  KS = 0
                  PITRAF = IPIW
                  PDTRAF = IPW2
C                 
                  PAIR = .FALSE.
                  DO 140 K = I, LIHI
                     IF( PAIR ) THEN
                        PAIR = .FALSE.
                     ELSE
                        SWAP = INT2LG(SELECT( K ))
                        IF( K.LT.LIHI ) THEN
                           CALL INFOG2L( K+1, K, DESCS, NPROW, NPCOL,
     $                        MYROW, MYCOL, ILOC, JLOC, RSRC, CSRC )
                           IF( S( LLDS*(JLOC-1) + ILOC ).NE.ZERO )
     $                        PAIR = .TRUE.
                        END IF
                        IF( SWAP ) THEN
                           KS = KS + 1
C                          
C                          Swap the K-th block to position I+KS-1.
C                          
                           
                           KK  = K - I + 1
                           KKS = KS
                           IF( KK.NE.KS ) THEN
                              NITRAF = LIWORK - PITRAF + 1
                              NDTRAF = LDWORK - PDTRAF + 1

                              IERR = 0
                              CALL BDTGEXC( NWIN, S(LLDS*(JLOC1-1) +
     $                           ILOC1), LLDS, T(LLDT*(JLOC1-1) + ILOC1)
     $                           , LLDT, KK, KKS, NITRAF, IWORK( PITRAF
     $                           ),NDTRAF, DWORK( PDTRAF ), DWORK(
     $                           PDTRAF + 18*NWIN*NWIN ), DWORK( IPW1 ),
     $                           4*N+16, IERR )

                              PITRAF = PITRAF + NITRAF
                              PDTRAF = PDTRAF + NDTRAF
C                             
C                             Update array SELECT.
C                             
                              IF ( PAIR ) THEN
                                 DO 150 J = I+KK-1, I+KKS, -1
                                    SELECT(J+1) = SELECT(J-1)
 150                             CONTINUE
                                 SELECT(I+KKS-1) = 1
                                 SELECT(I+KKS) = 1
                              ELSE
                                 DO 160 J = I+KK-1, I+KKS, -1
                                    SELECT(J) = SELECT(J-1)
 160                             CONTINUE
                                 SELECT(I+KKS-1) = 1
                              END IF
C                             
                              IF ( IERR.EQ.1 .OR. IERR.EQ.2 ) THEN

C                                
C                                Some blocks are too close to swap: prepare to leave in a clean 
C                                fashion. If IERR.EQ.2, we must update SELECT to account for the 
C                                fact that the 2 by 2 block to be reordered did split and the first 
C                                part of this block is already reordered.
C                                
                                 IF ( IERR.EQ.2 ) THEN
                                    SELECT( I+KKS-3 ) = 1
                                    SELECT( I+KKS-1 ) = 0
                                    KKS = KKS + 1
                                 END IF
C                                
C                                Update off-diagonal blocks immediately.
C                                
                                 GO TO 170
                              END IF
                              KS = KKS     
                           END IF
                           IF( PAIR )
     $                        KS = KS + 1
                        END IF
                     END IF
 140              CONTINUE
               END IF
 110        CONTINUE
 170        CONTINUE
C           
C           The on-diagonal processes save their information from the local
C           reordering in the integer buffer. This buffer is broadcasted to
C           updating processors, see below.
C           
            DO 175 WINDOW = 1, NMWIN2
               RSRC = IWORK(IRSRC+WINDOW-1)
               CSRC = IWORK(ICSRC+WINDOW-1)
               IF( MYROW.EQ.RSRC .AND. MYCOL.EQ.CSRC ) THEN
                  IBUFF( 1 ) = I
                  IBUFF( 2 ) = NWIN
                  IBUFF( 3 ) = PITRAF
                  IBUFF( 4 ) = KS
                  IBUFF( 5 ) = PDTRAF
                  IBUFF( 6 ) = NDTRAF
                  ILEN = PITRAF - IPIW 
                  DLEN = PDTRAF - IPW2
                  IBUFF( 7 ) = ILEN
                  IBUFF( 8 ) = DLEN
               END IF
 175        CONTINUE
C           
C           For the updates with respect to the local reordering, we organize
C           the updates in two phases where the update "direction" (controlled
C           by the DIR variable below) is first chosen to be the corresponding
C           rows, then the corresponding columns
C           
            DO 1111 DIR = 1, 2

C              
C              Broadcast information about the reordering and the accumulated
C              transformations: I, NWIN, PITRAF, NITRAF, PDTRAF, NDTRAF. If no
C              broadcast is performed, use a artificial value of KS to prevent
C              updating indicies for windows already finished (use KS = -1)
C              
               DO 111 WINDOW = 1, NMWIN2
                  RSRC = IWORK(IRSRC+WINDOW-1)
                  CSRC = IWORK(ICSRC+WINDOW-1)
                  IF( MYROW.EQ.RSRC .OR. MYCOL.EQ.CSRC ) THEN
                     LILO = IWORK(ILILO+WINDOW-1)
                     LIHI = IWORK(ILIHI+WINDOW-1)
                     LSEL = IWORK(ILSEL+WINDOW-1)
                  END IF

                  IF( MYROW.EQ.RSRC .AND. MYCOL.EQ.CSRC ) THEN
                     IF( NPCOL.GT.1 .AND. DIR.EQ.1 )
     $                  CALL IGEBS2D( ICTXT, 'Row', ' ', 8, 1, IBUFF, 
     $                  8)
                     IF( NPROW.GT.1 .AND. DIR.EQ.2 )
     $                  CALL IGEBS2D( ICTXT, 'Col', ' ', 8, 1, IBUFF, 
     $                  8)
                  ELSEIF( MYROW.EQ.RSRC .OR. MYCOL.EQ.CSRC ) THEN
                     IF( NPCOL.GT.1 .AND. DIR.EQ.1 .AND. MYROW.EQ.RSRC )
     $                  THEN
                        IF( FIRST .OR. (LIHI .GE. LILO + LSEL) ) THEN 

                           CALL IGEBR2D( ICTXT, 'Row', ' ', 8, 1,IBUFF,
     $                        8, RSRC, CSRC )

                           I = IBUFF( 1 ) 
                           NWIN = IBUFF( 2 ) 
                           PITRAF = IBUFF( 3 ) 
                           KS = IBUFF( 4 ) 
                           PDTRAF = IBUFF( 5 ) 
                           NDTRAF = IBUFF( 6 ) 
                           ILEN = IBUFF( 7 )
                           DLEN = IBUFF( 8 )
                        ELSE
                           ILEN = 0
                           DLEN = 0
                           KS = -1
                        END IF
                     END IF
                     IF( NPROW.GT.1 .AND. DIR.EQ.2 .AND. MYCOL.EQ.CSRC )
     $                  THEN
                        IF( FIRST .OR. (LIHI .GE. LILO + LSEL) ) THEN 

                           CALL IGEBR2D( ICTXT, 'Col', ' ', 8, 1, IBUFF,
     $                        8, RSRC, CSRC )

                           I = IBUFF( 1 ) 
                           NWIN = IBUFF( 2 ) 
                           PITRAF = IBUFF( 3 ) 
                           KS = IBUFF( 4 ) 
                           PDTRAF = IBUFF( 5 ) 
                           NDTRAF = IBUFF( 6 )
                           ILEN = IBUFF( 7 )
                           DLEN = IBUFF( 8 ) 
                        ELSE
                           ILEN = 0
                           DLEN = 0
                           KS = -1
                        END IF
                     END IF
                  END IF

C                 
C                 Broadcast the accumulated transformations - copy all information
C                 from IWORK(IPIW:PITRAF-1) and DWORK(IPW2:PDTRAF-1) and 
C                 DWORK(IPW2+18*NWIN*NWIN:PDTRAF+18*NWIN*NWIN-1) to a buffer and 
C                 broadcast this buffer in the corresponding block row and column.
C                 On arrival, copy the information back to the correct part of the
C                 workspace. This step is avoided if no computations were performed
C                 at the diagonal processor, i.e., BUFFLEN = 0.
C                 
                  IF( MYROW.EQ.RSRC .AND. MYCOL.EQ.CSRC ) THEN
                     BUFFER = PDTRAF+18*NWIN*NWIN
                     BUFFLEN = 2*DLEN + ILEN
                     IF( BUFFLEN.NE.0 ) THEN

                        DO 180 INDX = 1, ILEN
                           DWORK( BUFFER+INDX-1 ) = DBLE( IWORK(IPIW
     $                        +INDX-1) )
 180                    CONTINUE
                        CALL DLACPY( 'All', DLEN, 1, DWORK( IPW2 ), DLEN
     $                     , DWORK(BUFFER+ILEN), DLEN )
                        CALL DLACPY( 'All', DLEN, 1, DWORK( IPW2+18*NWIN
     $                     *NWIN ), DLEN, DWORK(BUFFER+ILEN+DLEN), DLEN
     $                     )

                        IF( NPCOL.GT.1 .AND. DIR.EQ.1 ) THEN
                           CALL DGEBS2D( ICTXT, 'Row', ' ', BUFFLEN, 1,
     $                        DWORK(BUFFER), BUFFLEN )
                        END IF
                        IF( NPROW.GT.1 .AND. DIR.EQ.2 ) THEN
                           CALL DGEBS2D( ICTXT, 'Col', ' ', BUFFLEN, 1,
     $                        DWORK(BUFFER), BUFFLEN )
                        END IF
                     END IF
                  ELSEIF( MYROW.EQ.RSRC .OR. MYCOL.EQ.CSRC ) THEN
                     IF( NPCOL.GT.1 .AND. DIR.EQ.1 .AND. MYROW.EQ.RSRC ) 
     $                  THEN
                        BUFFER = PDTRAF+18*NWIN*NWIN
                        BUFFLEN = 2*DLEN + ILEN
                        IF( BUFFLEN.NE.0 ) THEN 
                           CALL DGEBR2D( ICTXT, 'Row', ' ', BUFFLEN, 1,
     $                        DWORK(BUFFER), BUFFLEN, RSRC, CSRC )
                        END IF
                     END IF
                     IF( NPROW.GT.1 .AND. DIR.EQ.2 .AND. MYCOL.EQ.CSRC ) 
     $                  THEN
                        BUFFER = PDTRAF+18*NWIN*NWIN
                        BUFFLEN = 2*DLEN + ILEN
                        IF( BUFFLEN.NE.0 ) THEN
                           CALL DGEBR2D( ICTXT, 'Col', ' ', BUFFLEN, 1,
     $                        DWORK(BUFFER), BUFFLEN, RSRC, CSRC )
                        END IF
                     END IF
                     IF((NPCOL.GT.1.AND.DIR.EQ.1.AND.MYROW.EQ.RSRC).OR.
     $                  (NPROW.GT.1.AND.DIR.EQ.2.AND.MYCOL.EQ.CSRC ) ) 
     $                  THEN
                        IF( BUFFLEN.NE.0 ) THEN
                           DO 190 INDX = 1, ILEN
                              IWORK(IPIW+INDX-1) = INT(DWORK( BUFFER
     $                           +INDX-1 ))
 190                       CONTINUE
                           CALL DLACPY( 'All', DLEN, 1, DWORK( BUFFER
     $                        +ILEN ), DLEN, DWORK( IPW2 ), DLEN )
                           CALL DLACPY( 'All', DLEN, 1, DWORK( BUFFER
     $                        +ILEN+DLEN ), DLEN, DWORK( IPW2+18*NWIN
     $                        *NWIN ), DLEN )
                        END IF
                     END IF
                  END IF
 111           CONTINUE
C              
C              Now really perform the updates by applying the orthogonal 
C              transformations to the out-of-window parts of (S,T) and 
C              (Q,Z). This step is avoided if no reordering was performed 
C              by the on-diagonal processor from the beginning, i.e., 
C              BUFFLEN = 0.
C              
C              Count number of operations to decide whether to use matrix-matrix 
C              multiplications for updating off-diagonal parts or not. 
C              
C              ISHH = .TRUE.
               DO 112 WINDOW = 1, NMWIN2
                  RSRC = IWORK(IRSRC+WINDOW-1)
                  CSRC = IWORK(ICSRC+WINDOW-1)
C                 
                  IF( (MYROW.EQ.RSRC .AND. DIR.EQ.1 ).OR.
     $               (MYCOL.EQ.CSRC .AND. DIR.EQ.2 ) ) THEN
                     LILO = IWORK(ILILO+WINDOW-1)
                     LIHI = IWORK(ILIHI+WINDOW-1)
                     LSEL = IWORK(ILSEL+WINDOW-1)
C                    
C                    Skip update part for current WINDOW if BUFFLEN = 0
C                    
                     IF( BUFFLEN.EQ.0 ) GO TO 295
C                    
                     NITRAF = PITRAF - IPIW
C                    ISHH = .FALSE.
                     FLOPS = 0
                     DO 200 K = 1, NITRAF
                        IF( IWORK( IPIW + K - 1 ).LE.NWIN ) THEN
                           FLOPS = FLOPS + 6
                        ELSE
                           FLOPS = FLOPS + 11
c                          ISHH = .TRUE.
                        END IF
 200                 CONTINUE
c                    ISHH = .TRUE.
C                    
C                    Compute amount of work space necessary for performing matrix-
C                    matrix multiplications.
C                    
                     PDU = BUFFER
                     PDV = PDU + NWIN*NWIN 
                     IPW3 = PDV + NWIN*NWIN
                  ELSE
                     FLOPS = 0
                  END IF
C                 
                  IF( FLOPS.NE.0 .AND.
     $               ( FLOPS*100 ) / ( 2*NWIN*NWIN ) .GE. MMULT ) 
     $               THEN

                     
C                    
C                    Update off-diagonal blocks of (S,T) and Q and Z using matrix-
C                    matrix multiplications; if there are no Householder reflectors
C                    it is preferable to take the triangular block structure
C                    of the transformation matrix into account.
C                    
                     CALL DLASET( 'All', NWIN, NWIN, ZERO, ONE, DWORK(
     $                  PDU ), NWIN )
                     CALL BDLAGPP( 1, NWIN, NWIN, NCB, DWORK( PDU ),
     $                  NWIN, NITRAF, IWORK(IPIW), DWORK( IPW2 ),
     $                  DWORK(IPW3))
                     CALL DLASET( 'All', NWIN, NWIN, ZERO, ONE, DWORK(
     $                  PDV ), NWIN )
                     CALL BDLAGPP( 1, NWIN, NWIN, NCB, DWORK( PDV ),
     $                  NWIN, NITRAF, IWORK(IPIW), DWORK( IPW2+18*NWIN
     $                  *NWIN ), DWORK(IPW3) )
C                    
                     IF( ISHH .OR. DIM1.NE.KS .OR. DIM4.NE.KS ) THEN
C                       
C                       Loop through the local blocks of the distributed matrices T and 
C                       Q and update them according to the performed reordering
C                       
C                       Update the columns of (S,T), Q and Z affected by the reordering
C                       
                        IF( DIR.EQ.2 ) THEN
                           DO 210 INDX = 1, I-1, NB
                              CALL INFOG2L( INDX, I, DESCS, NPROW, NPCOL
     $                           , MYROW, MYCOL, ILOC, JLOC, RSRC1,
     $                           CSRC1 )
                              IF( MYROW.EQ.RSRC1 .AND. MYCOL.EQ.CSRC1 )
     $                           THEN
                                 LROWS = MIN(NB,I-INDX)
                                 CALL DGEMM( 'N', 'N', LROWS, NWIN, NWIN
     $                              , ONE, S((JLOC-1)*LLDS+ILOC), LLDS,
     $                              DWORK( PDV ), NWIN, ZERO,
     $                              DWORK(IPW3), LROWS )
                                 CALL DLACPY( 'All', LROWS, NWIN,
     $                              DWORK(IPW3), LROWS, S((JLOC-1)*LLDS
     $                              +ILOC), LLDS )

                              END IF
                              CALL INFOG2L( INDX, I, DESCT, NPROW, NPCOL
     $                           , MYROW, MYCOL, ILOC, JLOC, RSRC1,
     $                           CSRC1 )
                              IF( MYROW.EQ.RSRC1 .AND. MYCOL.EQ.CSRC1 ) 
     $                           THEN
                                 LROWS = MIN(NB,I-INDX)
                                 CALL DGEMM( 'N', 'N', LROWS, NWIN, NWIN
     $                              , ONE, T((JLOC-1)*LLDT+ILOC), LLDT,
     $                              DWORK( PDV ), NWIN, ZERO,
     $                              DWORK(IPW3), LROWS )
                                 CALL DLACPY( 'All', LROWS, NWIN,
     $                              DWORK(IPW3), LROWS, T((JLOC-1)*LLDT
     $                              +ILOC), LLDT )
                              END IF
 210                       CONTINUE
                           IF( WANTQ ) THEN
                              DO 220 INDX = 1, N, NB
                                 CALL INFOG2L( INDX, I, DESCQ, NPROW,
     $                              NPCOL, MYROW, MYCOL, ILOC, JLOC,
     $                              RSRC1, CSRC1)
                                 IF( MYROW.EQ.RSRC1.AND.MYCOL.EQ.CSRC1 ) 
     $                              THEN
                                    LROWS = MIN(NB,N-INDX+1)
                                    CALL DGEMM( 'N', 'N', LROWS, NWIN,
     $                                 NWIN, ONE, Q((JLOC-1)*LLDQ+ILOC),
     $                                 LLDQ, DWORK( PDU ), NWIN, ZERO,
     $                                 DWORK(IPW3), LROWS )
                                    CALL DLACPY( 'All', LROWS, NWIN,
     $                                 DWORK(IPW3), LROWS, Q((JLOC-1)
     $                                 *LLDQ+ILOC), LLDQ )
                                 END IF
 220                          CONTINUE
                           END IF
                           IF( WANTZ ) THEN
                              DO 225 INDX = 1, N, NB
                                 CALL INFOG2L( INDX, I, DESCZ, NPROW,
     $                              NPCOL, MYROW, MYCOL, ILOC, JLOC,
     $                              RSRC1, CSRC1)
                                 IF( MYROW.EQ.RSRC1.AND.MYCOL.EQ.CSRC1 ) 
     $                              THEN
                                    LROWS = MIN(NB,N-INDX+1)
                                    CALL DGEMM( 'N', 'N', LROWS, NWIN,
     $                                 NWIN, ONE, Z((JLOC-1)*LLDZ+ILOC),
     $                                 LLDZ, DWORK( PDV ), NWIN, ZERO,
     $                                 DWORK(IPW3), LROWS )
                                    CALL DLACPY( 'All', LROWS, NWIN,
     $                                 DWORK(IPW3), LROWS, Z((JLOC-1)
     $                                 *LLDZ+ILOC), LLDZ )
                                 END IF
 225                          CONTINUE
                           END IF
                        END IF
C                       
C                       Update the rows of (S,T) affected by the reordering
C                       
                        IF( DIR.EQ.1 ) THEN
                           IF( LIHI.LT.N ) THEN
                              IF( MOD(LIHI,NB).GT.0 ) THEN
                                 INDX = LIHI + 1
                                 CALL INFOG2L( I, INDX, DESCS, NPROW,
     $                              NPCOL, MYROW, MYCOL, ILOC, JLOC,
     $                              RSRC1, CSRC1)
                                 IF( MYROW.EQ.RSRC1.AND.MYCOL.EQ.CSRC1 ) 
     $                              THEN
                                    LCOLS = MOD( MIN( NB-MOD(LIHI,NB), 
     $                                 N-LIHI ), NB )
                                    CALL DGEMM( 'T', 'No T', NWIN, LCOLS
     $                                 , NWIN, ONE, DWORK( PDU ), NWIN,
     $                                 S((JLOC-1)*LLDS+ILOC), LLDS, ZERO
     $                                 , DWORK(IPW3), NWIN )
                                    CALL DLACPY( 'All', NWIN, LCOLS,
     $                                 DWORK(IPW3), NWIN, S((JLOC-1)
     $                                 *LLDS+ILOC), LLDS )
                                 END IF
                                 CALL INFOG2L( I, INDX, DESCT, NPROW,
     $                              NPCOL, MYROW, MYCOL, ILOC, JLOC,
     $                              RSRC1, CSRC1)
                                 IF( MYROW.EQ.RSRC1.AND.MYCOL.EQ.CSRC1 ) 
     $                              THEN
                                    LCOLS = MOD( MIN( NB-MOD(LIHI,NB), 
     $                                 N-LIHI ), NB )
                                    CALL DGEMM( 'T', 'No T', NWIN, LCOLS
     $                                 , NWIN, ONE, DWORK( PDU ), NWIN,
     $                                 T((JLOC-1)*LLDT+ILOC), LLDT, ZERO
     $                                 , DWORK(IPW3), NWIN )
                                    CALL DLACPY( 'All', NWIN, LCOLS,
     $                                 DWORK(IPW3), NWIN, T((JLOC-1)
     $                                 *LLDT+ILOC), LLDT )
                                 END IF
                              END IF
                              INDXS = ICEIL(LIHI,NB)*NB + 1
                              DO 230 INDX = INDXS, N, NB 
                                 CALL INFOG2L( I, INDX, DESCS, NPROW,
     $                              NPCOL, MYROW, MYCOL, ILOC, JLOC,
     $                              RSRC1, CSRC1)
                                 IF( MYROW.EQ.RSRC1.AND.MYCOL.EQ.CSRC1 ) 
     $                              THEN
                                    LCOLS = MIN( NB, N-INDX+1 )
                                    CALL DGEMM( 'T', 'No T', NWIN, LCOLS
     $                                 , NWIN, ONE, DWORK( PDU ), NWIN,
     $                                 S((JLOC-1)*LLDS+ILOC), LLDS, ZERO
     $                                 , DWORK(IPW3), NWIN )
                                    CALL DLACPY( 'All', NWIN, LCOLS,
     $                                 DWORK(IPW3), NWIN, S((JLOC-1)
     $                                 *LLDS+ILOC), LLDS )
                                 END IF
 230                          CONTINUE
                              DO 235 INDX = INDXS, N, NB 
                                 CALL INFOG2L( I, INDX, DESCT, NPROW,
     $                              NPCOL, MYROW, MYCOL, ILOC, JLOC,
     $                              RSRC1, CSRC1)
                                 IF( MYROW.EQ.RSRC1.AND.MYCOL.EQ.CSRC1 ) 
     $                              THEN
                                    LCOLS = MIN( NB, N-INDX+1 )
                                    CALL DGEMM( 'T', 'No T', NWIN, LCOLS
     $                                 , NWIN, ONE, DWORK( PDU ), NWIN,
     $                                 T((JLOC-1)*LLDT+ILOC), LLDT, ZERO
     $                                 , DWORK(IPW3), NWIN )
                                    CALL DLACPY( 'All', NWIN, LCOLS,
     $                                 DWORK(IPW3), NWIN, T((JLOC-1)
     $                                 *LLDT+ILOC), LLDT )
                                 END IF
 235                          CONTINUE
                           END IF
                        END IF
                     ELSE
                        IF (IAM.EQ.0) THEN
                           WRITE(*,*)'Reorder1:Not supposed to happen'
                        END IF
                        CALL BLACS_BARRIER(ICTXT, 'A')
                        STOP

                     END IF
                  ELSEIF( FLOPS.NE.0 ) THEN

C                    
C                    Update off-diagonal blocks of (S,T), Q and Z using the pipelined 
C                    elementary transformations.
C                    
                     IF( DIR.EQ.2 ) THEN
                        DO 270 INDX = 1, I-1, NB
                           CALL INFOG2L( INDX, I, DESCS, NPROW, NPCOL,
     $                        MYROW, MYCOL, ILOC, JLOC, RSRC1, CSRC1 )
                           IF( MYROW.EQ.RSRC1 .AND. MYCOL.EQ.CSRC1 )
     $                        THEN
                              LROWS = MIN(NB,I-INDX)
                              CALL BDLAGPP( 1, LROWS, NWIN, NCB, S((JLOC
     $                           -1)*LLDS+ILOC ), LLDS, NITRAF,
     $                           IWORK(IPIW), DWORK( IPW2+18*NWIN*NWIN )
     $                           , DWORK(IPW3) )
                           END IF
                           CALL INFOG2L( INDX, I, DESCT, NPROW, NPCOL,
     $                        MYROW, MYCOL, ILOC, JLOC, RSRC1, CSRC1 )
                           IF( MYROW.EQ.RSRC1 .AND. MYCOL.EQ.CSRC1 ) 
     $                        THEN
                              LROWS = MIN(NB,I-INDX)
                              CALL BDLAGPP( 1, LROWS, NWIN, NCB, T((JLOC
     $                           -1)*LLDT+ILOC ), LLDT, NITRAF,
     $                           IWORK(IPIW), DWORK( IPW2+18*NWIN*NWIN )
     $                           , DWORK(IPW3) )
                           END IF
 270                    CONTINUE
                        IF( WANTQ ) THEN
                           DO 280 INDX = 1, N, NB
                              CALL INFOG2L( INDX, I, DESCQ, NPROW, NPCOL
     $                           , MYROW, MYCOL, ILOC, JLOC, RSRC1,
     $                           CSRC1 )
                              IF( MYROW.EQ.RSRC1 .AND. MYCOL.EQ.CSRC1 ) 
     $                           THEN
                                 LROWS = MIN(NB,N-INDX+1)
                                 CALL BDLAGPP( 1, LROWS, NWIN, NCB,
     $                              Q((JLOC-1)*LLDQ+ILOC), LLDQ, NITRAF,
     $                              IWORK(IPIW), DWORK( IPW2 ),
     $                              DWORK(IPW3) )
                              END IF
 280                       CONTINUE
                        END IF 
                        IF( WANTZ ) THEN
                           DO 285 INDX = 1, N, NB
                              CALL INFOG2L( INDX, I, DESCZ, NPROW, NPCOL
     $                           , MYROW, MYCOL, ILOC, JLOC, RSRC1,
     $                           CSRC1 )
                              IF( MYROW.EQ.RSRC1 .AND. MYCOL.EQ.CSRC1 ) 
     $                           THEN
                                 LROWS = MIN(NB,N-INDX+1)
                                 CALL BDLAGPP( 1, LROWS, NWIN, NCB,
     $                              Z((JLOC-1)*LLDZ+ILOC), LLDZ, NITRAF,
     $                              IWORK(IPIW), DWORK(IPW2+18*NWIN*NWIN
     $                              ), DWORK(IPW3) )
                              END IF
 285                       CONTINUE
                        END IF
                     END IF
                     IF( DIR.EQ.1 ) THEN
                        IF( LIHI.LT.N ) THEN
                           IF( MOD(LIHI,NB).GT.0 ) THEN
                              INDX = LIHI + 1
                              CALL INFOG2L( I, INDX, DESCS, NPROW, NPCOL
     $                           , MYROW, MYCOL, ILOC, JLOC, RSRC1,
     $                           CSRC1 )
                              IF( MYROW.EQ.RSRC1 .AND. MYCOL.EQ.CSRC1 )
     $                           THEN
                                 LCOLS = MOD( MIN( NB-MOD(LIHI,NB), N
     $                              -LIHI ), NB )
                                 CALL BDLAGPP( 0, NWIN, LCOLS, NCB,
     $                              S((JLOC-1)*LLDS+ILOC), LLDS, NITRAF,
     $                              IWORK(IPIW), DWORK( IPW2 ),
     $                              DWORK(IPW3) )
                              END IF
                           END IF
                           INDXS = ICEIL(LIHI,NB)*NB + 1
                           DO 290 INDX = INDXS, N, NB 
                              CALL INFOG2L( I, INDX, DESCS, NPROW, NPCOL
     $                           , MYROW, MYCOL, ILOC, JLOC, RSRC1,
     $                           CSRC1 )
                              IF( MYROW.EQ.RSRC1 .AND. MYCOL.EQ.CSRC1 ) 
     $                           THEN
                                 LCOLS = MIN( NB, N-INDX+1 )
                                 CALL BDLAGPP( 0, NWIN, LCOLS, NCB,
     $                              S((JLOC-1)*LLDS+ILOC), LLDS, NITRAF,
     $                              IWORK(IPIW), DWORK( IPW2 ),
     $                              DWORK(IPW3) )
                              END IF
 290                       CONTINUE
                           IF( MOD(LIHI,NB).GT.0 ) THEN
                              INDX = LIHI + 1
                              CALL INFOG2L( I, INDX, DESCT, NPROW, NPCOL
     $                           , MYROW, MYCOL, ILOC, JLOC, RSRC1,
     $                           CSRC1 )
                              IF( MYROW.EQ.RSRC1 .AND. MYCOL.EQ.CSRC1 ) 
     $                           THEN
                                 LCOLS = MOD( MIN( NB-MOD(LIHI,NB), N
     $                              -LIHI ), NB )
                                 CALL BDLAGPP( 0, NWIN, LCOLS, NCB,
     $                              T((JLOC-1)*LLDT+ILOC), LLDT, NITRAF,
     $                              IWORK(IPIW), DWORK( IPW2 ),
     $                              DWORK(IPW3) )
                              END IF
                           END IF
                           INDXS = ICEIL(LIHI,NB)*NB + 1
                           DO 292 INDX = INDXS, N, NB 
                              CALL INFOG2L( I, INDX, DESCT, NPROW, NPCOL
     $                           , MYROW, MYCOL, ILOC, JLOC, RSRC1,
     $                           CSRC1 )
                              IF( MYROW.EQ.RSRC1 .AND. MYCOL.EQ.CSRC1 ) 
     $                           THEN
                                 LCOLS = MIN( NB, N-INDX+1 )
                                 CALL BDLAGPP( 0, NWIN, LCOLS, NCB,
     $                              T((JLOC-1)*LLDT+ILOC), LLDT, NITRAF,
     $                              IWORK(IPIW), DWORK( IPW2 ),
     $                              DWORK(IPW3) )
                              END IF
 292                       CONTINUE
                        END IF
                     END IF
                  END IF
C                 
C                 If I was not involved in the updates for the current window or
C                 the window was fully processed, I go here and try again for the
C                 next window
C                 
 295              CONTINUE
C                 
C                 Update LILO and LIHI depending on the number of eigenvalues really 
C                 moved - for on-diagonal processes we do this update only once since
C                 each on-diagonal process is only involved with one window at one
C                 time. The indicies are updated in three cases:
C                 1) When some reordering was really performed - indicated by
C                 BUFFLEN > 0.
C                 2) When no selected eigenvalues was found in the current window
C                 - indicated by KS = 0.
C                 3) When some selected eigenvalues was found in the current window 
C                 but no one of them was moved (KS > 0 and BUFFLEN = 0)
C                 False index updating is avoided by sometimes setting KS = -1. This
C                 will affect processors involved in more than one window and where
C                 the first one ends up with KS = 0 and for the second one is done
C                 already. 
C                 
                  
                  IF( MYROW.EQ.RSRC.AND.MYCOL.EQ.CSRC ) THEN
                     IF( DIR.EQ.2 ) THEN
                        IF( BUFFLEN.NE.0 .OR. KS.EQ.0 .OR. (
     $                     BUFFLEN.EQ.0 .AND. KS.GT.0 ) ) LIHI = I + KS
     $                     - 1


                        IWORK( ILIHI+WINDOW-1 ) = LIHI
                        IF( .NOT. LIHI.GE.LILO+LSEL ) THEN
                           LILO = LILO + LSEL
                           IWORK( ILILO+WINDOW-1 ) = LILO
                        END IF
                     END IF
                  ELSEIF( MYROW.EQ.RSRC .AND. DIR.EQ.1 ) THEN
                     IF( BUFFLEN.NE.0 .OR. KS.EQ.0 .OR. ( BUFFLEN.EQ.0
     $                  .AND. KS.GT.0 ) ) LIHI = I + KS - 1
                     IWORK( ILIHI+WINDOW-1 ) = LIHI
                     IF( .NOT. LIHI.GE.LILO+LSEL ) THEN
                        LILO = LILO + LSEL
                        IWORK( ILILO+WINDOW-1 ) = LILO
                     END IF
                  ELSEIF( MYCOL.EQ.CSRC .AND. DIR.EQ.2 ) THEN
                     IF( BUFFLEN.NE.0 .OR. KS.EQ.0 .OR. ( BUFFLEN.EQ.0
     $                  .AND. KS.GT.0 ) ) LIHI = I + KS - 1
                     IWORK( ILIHI+WINDOW-1 ) = LIHI
                     IF( .NOT. LIHI.GE.LILO+LSEL ) THEN
                        LILO = LILO + LSEL
                        IWORK( ILILO+WINDOW-1 ) = LILO
                     END IF
                  END IF
C                 
 112           CONTINUE
C              
C              End of direction loop for updates with respect to local reordering
C              
 1111       CONTINUE
C           
C           Associate each process with one of the corresponding computational 
C           windows such that the test for another round of local reordering
C           is carried out properly. Since the column updates were computed 
C           after the row updates, it is sufficient to test for changing the
C           association to the window in the corresponding process row
C           
            DO 113 WINDOW = 1, NMWIN2
               RSRC = IWORK( IRSRC + WINDOW - 1 )
               IF( MYROW.EQ.RSRC .AND. (.NOT. LIHI.GE.LILO+LSEL )) THEN
                  LILO = IWORK( ILILO + WINDOW - 1 )
                  LIHI = IWORK( ILIHI + WINDOW - 1 )
                  LSEL = IWORK( ILSEL + WINDOW - 1 )
               END IF
 113        CONTINUE
C           
C           End While ( LIHI >= LILO + LSEL )
            ROUND = ROUND + 1
            IF( FIRST ) FIRST = .FALSE.
            GO TO 130
         END IF
C        
C        All processors excluded from the local reordering go here
C        
 114     CONTINUE

C        
C        Barrier to collect the processes before proceeding
C        

         
         CALL BLACS_BARRIER( ICTXT, 'All' )
         

C        
C        Compute global maximum of IERR so that we know if some process
C        experienced a failure in the reordering
C        
         MYIERR = IERR
         IF( NPROCS.GT.1 )
     $      CALL IGAMX2D( ICTXT, 'All', ' ', 1, 1, IERR, 1, -1, 
     $      -1, -1, -1, -1 )
C        
         IF( IERR.EQ.1 .OR. IERR.EQ.2 ) THEN
C           
C           When calling BDTREXC, the block at position I+KKS-1 failed to swap.
C           
            IF( MYIERR.EQ.1 .OR. MYIERR.EQ.2 ) THEN
               INFO = I+KKS-1
            ELSE 
               INFO = IERR
            END IF
            GO TO 300
         END IF
C        
C        
C        
C        Now, for each compuational window, move the selected eigenvalues
C        across the process border. Do this by forming the processors into
C        groups of four working together to bring the window over the 
C        border. The processes are numbered as follows
C        
C        1 | 2
C        -----
C        3 | 4
C        
C        where '|' and '-' denotes the process (and block) borders.
C        This implies that the cluster to be reordered over the border
C        is held by process 4, process 1 will receive the cluster after
C        the reordering, process 3 holds the local (2,1)th element of a 
C        2-by-2 diagonal block located on the block border and process 2
C        holds the closest off-diagonal part of the windo that is affected 
C        by the cross-border reordering.
C        
C        The active window is now ( I : LIHI[4], I : LIHI[4] ), where 
C        I = MAX( ILO, LIHI - 2*MOD(LIHI,NB) ). If this active window is to 
C        large compared to the value of PARA( 6 ), it will be truncated in 
C        both ends such that a maximum of PARA( 6 ) eigenvalues is reordered 
C        across the border this time.
C        
C        The active window will be collected and built in workspace at 
C        process 1 and 4, which both compute the reordering and return the 
C        updated parts to the corresponding processes 2-3. Next, the 
C        accumulated transformations are broadcasted for updates in the 
C        block rows and column that corresponds to the process rows and 
C        columns where process 1 and 4 reside.  
C        
C        The off-diagonal blocks are updated by the processes receiving
C        from the broadcasts of the orthogonal transformations. Since the
C        active window is split over the process borders, the updates
C        of (S,T), Q and Z requires that stripes of block rows of columns 
C        are exchanged between neighboring processes in the corresponding
C        process rows and columns.
C        
C        First, form each group of processors involved in the crossborder
C        reordering. Do this in two (or three) phases: 
C        1) Reorder each odd window over the border.
C        2) Reorder each even window over the border.
C        3) Reorder the last odd window over the border, if it was not
C        processed in the first phase.
C        
C        When reordering the odd windows over the border, we must make sure
C        that no process row or column is involved in both the first and 
C        the last window at the same time. This happens when the total 
C        number of windows is odd, greater than one and equal to the 
C        minumum process mesh dimension. Therefore the last odd window may 
C        be reordered over the border at last.
C        
         LASTWAIT = NMWIN2.GT.1 .AND. MOD(NMWIN2,2).EQ.1 .AND. 
     $      NMWIN2.EQ.MIN(NPROW,NPCOL)

         LAST = 0
 308     CONTINUE
         IF( LASTWAIT ) THEN
            IF( LAST.EQ.0 ) THEN
               WIN0S = 1
               WIN0E = 2
               WINE = NMWIN2 - 1
            ELSE
               WIN0S = NMWIN2
               WIN0E = NMWIN2
               WINE = NMWIN2
            END IF
         ELSE
            WIN0S = 1
            WIN0E = 2
            WINE = NMWIN2
         END IF

         DO 310 WINDOW0 = WIN0S, WIN0E
            DO 320 WINDOW = WINDOW0, WINE, 2
C              
C              Define the process holding the down-right part of the window
C              
               RSRC4 = IWORK(IRSRC+WINDOW-1)
               CSRC4 = IWORK(ICSRC+WINDOW-1)
C              
C              Define the other processes in the group of four
C              
               RSRC3 = RSRC4
               CSRC3 = MOD( CSRC4 - 1 + NPCOL, NPCOL )
               RSRC2 = MOD( RSRC4 - 1 + NPROW, NPROW )
               CSRC2 = CSRC4
               RSRC1 = RSRC2
               CSRC1 = CSRC3
               IF( ( MYROW.EQ.RSRC1 .AND. MYCOL.EQ.CSRC1 ) .OR. (
     $            MYROW.EQ.RSRC2 .AND. MYCOL.EQ.CSRC2 ) .OR. (
     $            MYROW.EQ.RSRC3 .AND. MYCOL.EQ.CSRC3 ) .OR. (
     $            MYROW.EQ.RSRC4 .AND. MYCOL.EQ.CSRC4 ) ) THEN
C                 
C                 Compute the correct active window - for reordering into a block
C                 that has not been active at all before, we try to reorder as
C                 many of our eigenvalues over the border as possible without
C                 knowing of the situation on the other side - this may cause
C                 very few eigenvalues to be reordered over the border this time
C                 (perhaps not any) but this should be an initial problem. 
C                 Anyway, the bottom-right position of the block will be at
C                 position LIHIC.
C                 
                  IF( MYROW.EQ.RSRC4 .AND. MYCOL.EQ.CSRC4 ) THEN
                     LIHI4 = ( IWORK( ILILO + WINDOW - 1 ) + IWORK(
     $                  ILIHI + WINDOW - 1 ) ) / 2 

                     LIHIC = MIN(LIHI4,(ICEIL(LIHI4,NB)-1)*NB+WNEICR)

C                    
C                    Fix LIHIC to avoid that bottom of window cuts 2-by-2 block and
C                    make sure all processors in the group knows about the correct
C                    value
C                    
                     IF( (.NOT. LIHIC.LE.NB) .AND. LIHIC.LT.N ) THEN
                        ILOC = INDXG2L( LIHIC+1, NB, MYROW, DESCS( RSRC_
     $                     ), NPROW )
                        JLOC = INDXG2L( LIHIC, NB, MYCOL, DESCS( CSRC_ )
     $                     , NPCOL )
                        IF( S( (JLOC-1)*LLDS+ILOC ).NE.ZERO ) THEN
                           IF( MOD( LIHIC, NB ).EQ.1 .OR. ( MOD( LIHIC,
     $                        NB ).EQ.2 .AND. .NOT. INT2LG(SELECT(LIHIC
     $                        -2)) ) ) THEN
                              LIHIC = LIHIC + 1
                           ELSE
                              LIHIC = LIHIC - 1
                           END IF
                        END IF
                     END IF

                     IF( RSRC4.NE.RSRC1 .OR. CSRC4.NE.CSRC1 ) CALL
     $                  IGESD2D( ICTXT, 1, 1, LIHIC, 1, RSRC1, CSRC1 )
                     IF( RSRC4.NE.RSRC2 .OR. CSRC4.NE.CSRC2 ) CALL
     $                  IGESD2D( ICTXT, 1, 1, LIHIC, 1, RSRC2, CSRC2 )
                     IF( RSRC4.NE.RSRC3 .OR. CSRC4.NE.CSRC3 ) CALL
     $                  IGESD2D( ICTXT, 1, 1, LIHIC, 1, RSRC3, CSRC3 ) 
                  END IF
                  IF( MYROW.EQ.RSRC1 .AND. MYCOL.EQ.CSRC1 ) THEN
                     IF( RSRC4.NE.RSRC1 .OR. CSRC4.NE.CSRC1 ) THEN
                        CALL IGERV2D( ICTXT, 1, 1, LIHIC, 1, RSRC4,
     $                     CSRC4 )
                     END IF
                  END IF
                  IF( MYROW.EQ.RSRC2 .AND. MYCOL.EQ.CSRC2 ) THEN
                     IF( RSRC4.NE.RSRC2 .OR. CSRC4.NE.CSRC2 ) THEN
                        CALL IGERV2D( ICTXT, 1, 1, LIHIC, 1, RSRC4,
     $                     CSRC4 )
                     END IF
                  END IF
                  IF( MYROW.EQ.RSRC3 .AND. MYCOL.EQ.CSRC3 ) THEN
                     IF( RSRC4.NE.RSRC3 .OR. CSRC4.NE.CSRC3 ) THEN
                        CALL IGERV2D( ICTXT, 1, 1, LIHIC, 1, RSRC4,
     $                     CSRC4 )
                     END IF
                  END IF
C                 
C                 Avoid going over the border with the first window if it resides
C                 in the block where the last global position S(ILO,ILO) is or
C                 ILO has been updated to point to a position right of 
C                 S(LIHIC,LIHIC)
C                 
                  SKIP1CR = WINDOW.EQ.1 .AND. 
     $               ICEIL(LIHIC,NB).LE.ICEIL(ILO,NB) 
C                 
C                 Decide I, where to put top of window, such that top of window 
C                 does not cut 2-by-2 block. Make sure that we do not end up
C                 in a situation where a 2-by-2 block splitted on the border
C                 is left in its original place - this can cause infinite loops.
C                 Remedy: make sure that the part of the window that resides left
C                 to the border is at least of dimension two (2) in case we have
C                 2-by-2 blocks in top of the cross border window.
C                 
C                 Also make sure all processors in the group knows about the correct 
C                 value of I. When skipping the crossborder reordering, just set
C                 I = LIHIC.
C                 
                  IF( .NOT. SKIP1CR ) THEN
                     IF( MYROW.EQ.RSRC1 .AND. MYCOL.EQ.CSRC1 ) THEN
                        IF( WINDOW.EQ.1 ) THEN
                           LIHI1 = ILO 
                        ELSE
                           LIHI1 = IWORK( ILIHI + WINDOW - 2 )
                        END IF

                        I = MAX( LIHI1, 
     $                     MIN( LIHIC-2*MOD(LIHIC,NB) + 1, 
     $                     (ICEIL(LIHIC,NB)-1)*NB - 1  ) )

                        ILOC = INDXG2L( I, NB, MYROW, DESCS( RSRC_ ),
     $                     NPROW )
                        JLOC = INDXG2L( I-1, NB, MYCOL, DESCS( CSRC_ ),
     $                     NPCOL )
                        IF( S( (JLOC-1)*LLDS+ILOC ).NE.ZERO ) THEN
                           I = I - 1
                        END IF

                        IF( RSRC1.NE.RSRC4 .OR. CSRC1.NE.CSRC4 ) CALL
     $                     IGESD2D( ICTXT, 1, 1, I, 1, RSRC4, CSRC4 )
                        IF( RSRC1.NE.RSRC2 .OR. CSRC1.NE.CSRC2 ) CALL
     $                     IGESD2D( ICTXT, 1, 1, I, 1, RSRC2, CSRC2 )
                        IF( RSRC1.NE.RSRC3 .OR. CSRC1.NE.CSRC3 ) CALL
     $                     IGESD2D( ICTXT, 1, 1, I, 1, RSRC3, CSRC3 ) 
                     END IF
                     IF( MYROW.EQ.RSRC2 .AND. MYCOL.EQ.CSRC2 ) THEN
                        IF( RSRC1.NE.RSRC2 .OR. CSRC1.NE.CSRC2 ) THEN 
                           CALL IGERV2D( ICTXT, 1, 1, I, 1, RSRC1, CSRC1
     $                        )
                        END IF
                     END IF
                     IF( MYROW.EQ.RSRC3 .AND. MYCOL.EQ.CSRC3 ) THEN
                        IF( RSRC1.NE.RSRC3 .OR. CSRC1.NE.CSRC3 ) THEN
                           CALL IGERV2D( ICTXT, 1, 1, I, 1, RSRC1, CSRC1
     $                        )
                        END IF
                     END IF
                     IF( MYROW.EQ.RSRC4 .AND. MYCOL.EQ.CSRC4 ) THEN
                        IF( RSRC1.NE.RSRC4 .OR. CSRC1.NE.CSRC4 ) THEN
                           CALL IGERV2D( ICTXT, 1, 1, I, 1, RSRC1, CSRC1
     $                        )

                        ENDIF
                     END IF            
                  ELSE
                     I = LIHIC
                  END IF
C                 
C                 Finalize computation of window size: active window is now
C                 (I:LIHIC,I:LIHIC)
C                 
                  NWIN = LIHIC - I + 1
                  KS = 0
C                 
C                 Skip rest of this part if appropriate
C                 
                  IF( SKIP1CR ) GO TO 360
C                 
C                 Divide workspace - put active window in 
C                 DWORK(IPW2:IPW2+2*NWIN**2-1) and orthogonal transformations 
C                 in DWORK(IPW3:...). The active window is stored as follows
C                 WIN = [S(I:LIHIC,I:LIHIC),T(I:LIHIC,I:LIHIC)]
C                 
                  CALL DLASET( 'All', NWIN, NWIN, ZERO, ZERO, DWORK(
     $               IPW2 ), NWIN )
                  CALL DLASET( 'All', NWIN, NWIN, ZERO, ZERO, DWORK(
     $               IPW2+NWIN*NWIN ), NWIN )
C                 
                  PITRAF = IPIW
                  IPW3 = IPW2 + 2*NWIN*NWIN
                  PDTRAF = IPW3
C                 
C                 Exchange the current view of SELECT for the active window between 
C                 process 1 and 4 to make sure that exactly the same job is performed
C                 for both processes.
C                 
                  IF( RSRC1.NE.RSRC4 .OR. CSRC1.NE.CSRC4 ) THEN
                     ILEN4 = MOD(LIHIC,NB)
                     SELI4 = ICEIL(I,NB)*NB+1
                     ILEN1 = NWIN - ILEN4 
                     IF( MYROW.EQ.RSRC1 .AND. MYCOL.EQ.CSRC1 ) THEN
                        CALL IGESD2D( ICTXT, ILEN1, 1, SELECT(I), ILEN1,
     $                     RSRC4, CSRC4 )

                        CALL IGERV2D( ICTXT, ILEN4, 1, SELECT(SELI4),
     $                     ILEN4, RSRC4, CSRC4 )

                     END IF
                     IF( MYROW.EQ.RSRC4 .AND. MYCOL.EQ.CSRC4 ) THEN
                        CALL IGESD2D( ICTXT, ILEN4, 1, SELECT(SELI4),
     $                     ILEN4, RSRC1, CSRC1 )

                        CALL IGERV2D( ICTXT, ILEN1, 1, SELECT(I), ILEN1,
     $                     RSRC1, CSRC1 )

                     END IF
                  END IF
C                 
C                 Form the active window by a series of point-to-point sends
C                 and receives.  
C                 
                  DIM1 = NB - MOD(I-1,NB) 
                  DIM4 = NWIN - DIM1
                  IF( MYROW.EQ.RSRC1 .AND. MYCOL.EQ.CSRC1 ) THEN
                     ILOC = INDXG2L( I, NB, MYROW, DESCS( RSRC_ ), NPROW
     $                  )
                     JLOC = INDXG2L( I, NB, MYCOL, DESCS( CSRC_ ), NPCOL
     $                  ) 
                     CALL DLACPY( 'All', DIM1, DIM1, S((JLOC-1)*LLDS
     $                  +ILOC), LLDS, DWORK(IPW2), NWIN )
                     CALL DLACPY( 'All', DIM1, DIM1, T((JLOC-1)*LLDT
     $                  +ILOC), LLDT, DWORK(IPW2+NWIN*NWIN), NWIN )
                     IF( RSRC1.NE.RSRC4 .OR. CSRC1.NE.CSRC4 ) THEN
                        CALL DGESD2D( ICTXT, DIM1, DIM1, DWORK(IPW2),
     $                     NWIN, RSRC4, CSRC4 )

                        CALL DGERV2D( ICTXT, DIM4, DIM4, DWORK(IPW2+DIM1
     $                     *NWIN+DIM1), NWIN, RSRC4, CSRC4 )

                        CALL DGESD2D( ICTXT, DIM1, DIM1, DWORK(IPW2+NWIN
     $                     *NWIN), NWIN, RSRC4, CSRC4)

                        CALL DGERV2D( ICTXT, DIM4, DIM4, DWORK(IPW2+NWIN
     $                     *NWIN+DIM1*NWIN+DIM1),NWIN, RSRC4, CSRC4 )

                     END IF
                  END IF
                  IF( MYROW.EQ.RSRC4 .AND. MYCOL.EQ.CSRC4 ) THEN
                     ILOC = INDXG2L( I+DIM1, NB, MYROW, DESCS( RSRC_ ),
     $                  NPROW )
                     JLOC = INDXG2L( I+DIM1, NB, MYCOL, DESCS( CSRC_ ),
     $                  NPCOL ) 
                     CALL DLACPY( 'All', DIM4, DIM4, S((JLOC-1)*LLDS
     $                  +ILOC), LLDS, DWORK(IPW2+DIM1*NWIN+DIM1), NWIN )
                     CALL DLACPY( 'All', DIM4, DIM4, T((JLOC-1)*LLDT
     $                  +ILOC), LLDT, DWORK(IPW2+NWIN*NWIN+DIM1*NWIN
     $                  +DIM1), NWIN )
                     IF( RSRC4.NE.RSRC1 .OR. CSRC4.NE.CSRC1 ) THEN
                        CALL DGESD2D( ICTXT, DIM4, DIM4, DWORK(IPW2+DIM1
     $                     *NWIN+DIM1), NWIN, RSRC1, CSRC1 ) 

                        CALL DGERV2D( ICTXT, DIM1, DIM1, DWORK(IPW2),
     $                     NWIN, RSRC1, CSRC1 )

                        CALL DGESD2D( ICTXT, DIM4, DIM4, DWORK(IPW2+NWIN
     $                     *NWIN+DIM1*NWIN+DIM1),NWIN, RSRC1, CSRC1 ) 

                        CALL DGERV2D( ICTXT, DIM1, DIM1, DWORK(IPW2+NWIN
     $                     *NWIN), NWIN, RSRC1, CSRC1)

                     END IF
                  END IF
                  IF( MYROW.EQ.RSRC2 .AND. MYCOL.EQ.CSRC2 ) THEN
                     ILOC = INDXG2L( I, NB, MYROW, DESCS( RSRC_ ), 
     $                  NPROW )
                     JLOC = INDXG2L( I+DIM1, NB, MYCOL, DESCS( CSRC_ ),
     $                  NPCOL )
                     CALL DLACPY( 'All', DIM1, DIM4, S((JLOC-1)*LLDS
     $                  +ILOC), LLDS, DWORK(IPW2+DIM1*NWIN), NWIN )
                     CALL DLACPY( 'All', DIM1, DIM4, T((JLOC-1)*LLDT
     $                  +ILOC), LLDT, DWORK(IPW2+NWIN*NWIN+DIM1*NWIN),
     $                  NWIN )
                     IF( RSRC2.NE.RSRC1 .OR. CSRC2.NE.CSRC1 ) THEN
                        CALL DGESD2D( ICTXT, DIM1, DIM4, DWORK(IPW2+DIM1
     $                     *NWIN), NWIN, RSRC1, CSRC1)
                        CALL DGESD2D( ICTXT, DIM1, DIM4, DWORK(IPW2+NWIN
     $                     *NWIN+DIM1*NWIN), NWIN, RSRC1, CSRC1 )
                     END IF
                  END IF
                  IF( MYROW.EQ.RSRC2 .AND. MYCOL.EQ.CSRC2 ) THEN 
                     IF( RSRC2.NE.RSRC4 .OR. CSRC2.NE.CSRC4 ) THEN
                        CALL DGESD2D( ICTXT, DIM1, DIM4, DWORK(IPW2+DIM1
     $                     *NWIN), NWIN, RSRC4, CSRC4)
                        CALL DGESD2D( ICTXT, DIM1, DIM4, DWORK(IPW2+NWIN
     $                     *NWIN+DIM1*NWIN), NWIN, RSRC4, CSRC4 )
                     END IF
                  END IF
                  IF( MYROW.EQ.RSRC3 .AND. MYCOL.EQ.CSRC3 ) THEN
                     ILOC = INDXG2L( I+DIM1, NB, MYROW, DESCS( RSRC_ ),
     $                  NPROW )
                     JLOC = INDXG2L( I+DIM1-1, NB, MYCOL, 
     $                  DESCS( CSRC_ ), NPCOL )
                     CALL DLACPY( 'All', 1, 1, S((JLOC-1)*LLDS+ILOC),
     $                  LLDS, DWORK(IPW2+(DIM1-1)*NWIN+DIM1), NWIN )
                     IF( RSRC3.NE.RSRC1 .OR. CSRC3.NE.CSRC1 ) THEN 
                        CALL DGESD2D( ICTXT, 1, 1, DWORK(IPW2+(DIM1-1)
     $                     *NWIN+DIM1), NWIN, RSRC1, CSRC1 )
                     END IF
                  END IF
                  IF( MYROW.EQ.RSRC3 .AND. MYCOL.EQ.CSRC3 ) THEN
                     IF( RSRC3.NE.RSRC4 .OR. CSRC3.NE.CSRC4 ) THEN 
                        CALL DGESD2D( ICTXT, 1, 1, DWORK(IPW2+(DIM1-1)
     $                     *NWIN+DIM1), NWIN, RSRC4, CSRC4 )
                     END IF
                  END IF
                  IF( MYROW.EQ.RSRC1 .AND. MYCOL.EQ.CSRC1 ) THEN
                     IF( RSRC1.NE.RSRC2 .OR. CSRC1.NE.CSRC2 ) THEN

                        CALL DGERV2D( ICTXT, DIM1, DIM4, DWORK(IPW2+DIM1
     $                     *NWIN), NWIN, RSRC2, CSRC2 )
                        CALL DGERV2D( ICTXT, DIM1, DIM4, DWORK(IPW2+NWIN
     $                     *NWIN+DIM1*NWIN), NWIN, RSRC2, CSRC2 )

                     END IF
                     IF( RSRC1.NE.RSRC3 .OR. CSRC1.NE.CSRC3 ) THEN

                        CALL DGERV2D( ICTXT, 1, 1, DWORK(IPW2+(DIM1-1)
     $                     *NWIN+DIM1), NWIN, RSRC3, CSRC3 )

                     END IF
                  END IF
                  IF( MYROW.EQ.RSRC4 .AND. MYCOL.EQ.CSRC4 ) THEN
                     IF( RSRC4.NE.RSRC2 .OR. CSRC4.NE.CSRC2 ) THEN

                        CALL DGERV2D( ICTXT, DIM1, DIM4, DWORK(IPW2+DIM1
     $                     *NWIN), NWIN, RSRC2, CSRC2 )
                        CALL DGERV2D( ICTXT, DIM1, DIM4, DWORK(IPW2+NWIN
     $                     *NWIN+DIM1*NWIN), NWIN, RSRC2, CSRC2 )

                     END IF
                     IF( RSRC4.NE.RSRC3 .OR. CSRC4.NE.CSRC3 ) THEN
                        
                        CALL DGERV2D( ICTXT, 1, 1, DWORK(IPW2+(DIM1-1)
     $                     *NWIN+DIM1), NWIN, RSRC3, CSRC3 )

                     END IF
                  END IF
C                 
C                 Compute the reordering (just as in the total local case) and
C                 accumulate the transformations (ONLY ON-DIAGONAL PROCESSES)
C                 
                  IF( ( MYROW.EQ.RSRC1 .AND. MYCOL.EQ.CSRC1 ) .OR.
     $               ( MYROW.EQ.RSRC4 .AND. MYCOL.EQ.CSRC4 ) ) THEN 
                     PAIR = .FALSE.
                     DO 330 K = I, LIHIC
                        IF( PAIR ) THEN
                           PAIR = .FALSE.
                        ELSE
                           SWAP = INT2LG(SELECT( K ))
                           IF( K.LT.LIHIC ) THEN
                              ELEM = DWORK(IPW2+(K-I)*NWIN+K-I+1)
                              IF( ELEM.NE.ZERO ) PAIR = .TRUE.
                           END IF
                           IF(SWAP ) THEN
                              KS = KS + 1
C                             
C                             Swap the K-th block to position I+KS-1.
C                             
                              
                              KK  = K - I + 1
                              KKS = KS
                              IF( KK.NE.KS ) THEN
                                 NITRAF = LIWORK - PITRAF + 1
                                 NDTRAF = LDWORK - PDTRAF + 1

                                 IERR = 0
                                 CALL BDTGEXC( NWIN, DWORK(IPW2), NWIN,
     $                              DWORK(IPW2+NWIN*NWIN), NWIN, KK, KKS
     $                              , NITRAF, IWORK( PITRAF ), NDTRAF,
     $                              DWORK( PDTRAF ), DWORK( PDTRAF + 18
     $                              *NWIN*NWIN ), DWORK( IPW1 ), 4*N+16,
     $                              IERR )

                                 PITRAF = PITRAF + NITRAF
                                 PDTRAF = PDTRAF + NDTRAF

C                                
C                                Update array SELECT.
C                                
                                 IF ( PAIR ) THEN
                                    DO 340 J = I+KK-1, I+KKS, -1
                                       SELECT(J+1) = SELECT(J-1)
 340                                CONTINUE
                                    SELECT(I+KKS-1) = 1
                                    SELECT(I+KKS) = 1
                                 ELSE
                                    DO 350 J = I+KK-1, I+KKS, -1
                                       SELECT(J) = SELECT(J-1)
 350                                CONTINUE
                                    SELECT(I+KKS-1) = 1
                                 END IF
C                                
                                 IF ( IERR.EQ.1 .OR. IERR.EQ.2 ) THEN
C                                   
                                    IF ( IERR.EQ.2 ) THEN
                                       SELECT( I+KKS-3 ) = 1
                                       SELECT( I+KKS-1 ) = 0
                                       KKS = KKS + 1
                                    END IF
                                    GO TO 360
                                 END IF
                                 KS = KKS     
                              END IF
                              IF( PAIR )
     $                           KS = KS + 1
                           END IF
                        END IF
 330                 CONTINUE
                  END IF
 360              CONTINUE
C                 
C                 Save information about the reordering 
C                 
                  IF( ( MYROW.EQ.RSRC1 .AND. MYCOL.EQ.CSRC1 ) .OR. (
     $               MYROW.EQ.RSRC4 .AND. MYCOL.EQ.CSRC4 ) ) THEN 
                     IBUFF( 1 ) = I
                     IBUFF( 2 ) = NWIN
                     IBUFF( 3 ) = PITRAF
                     IBUFF( 4 ) = KS
                     IBUFF( 5 ) = PDTRAF
                     IBUFF( 6 ) = NDTRAF 
                     ILEN = PITRAF - IPIW + 1
                     DLEN = PDTRAF - IPW3 + 1
                     IBUFF( 7 ) = ILEN
                     IBUFF( 8 ) = DLEN
C                    
C                    Put reordered data back into global matrix if a reordering
C                    took place
C                    
                     IF( .NOT. SKIP1CR ) THEN 
                        IF( MYROW.EQ.RSRC1 .AND. MYCOL.EQ.CSRC1 ) THEN
                           ILOC = INDXG2L( I, NB, MYROW,DESCS( RSRC_ ),
     $                        NPROW )
                           JLOC = INDXG2L( I, NB, MYCOL,DESCS( CSRC_ ),
     $                        NPCOL ) 
                           CALL DLACPY( 'All', DIM1, DIM1, DWORK(IPW2),
     $                        NWIN, S((JLOC-1)*LLDS+ILOC), LLDS )
                           CALL DLACPY( 'All', DIM1, DIM1, DWORK(IPW2
     $                        +NWIN*NWIN), NWIN, T((JLOC-1)*LLDT+ILOC),
     $                        LLDT )
                        END IF
                        IF( MYROW.EQ.RSRC4 .AND. MYCOL.EQ.CSRC4 ) THEN
                           ILOC = INDXG2L( I+DIM1, NB, MYROW, DESCS(
     $                        RSRC_ ), NPROW )
                           JLOC = INDXG2L( I+DIM1, NB, MYCOL, DESCS(
     $                        CSRC_ ), NPCOL ) 
                           CALL DLACPY( 'All', DIM4, DIM4, DWORK(IPW2
     $                        +DIM1*NWIN+DIM1), NWIN, S((JLOC-1)*LLDS
     $                        +ILOC), LLDS )
                           CALL DLACPY( 'All', DIM4, DIM4, DWORK(IPW2
     $                        +NWIN*NWIN+DIM1*NWIN+DIM1), NWIN, T((JLOC
     $                        -1)*LLDT+ILOC), LLDT )
                        END IF
                     END IF
                  END IF
C                 
C                 Break if appropriate - IBUFF(3:8) may now contain nonsens, but
C                 that's no problem. The processors outside the cross border group
C                 only needs to know about I and NWIN to get a correct value of
C                 SKIP1CR (see below) and to skip the cross border updates if
C                 necessary.
C                 
                  IF( WINDOW.EQ.1 .AND. SKIP1CR ) GO TO 325
C                 
C                 Return reordered data to process 2 and 3
C                 
                  IF( MYROW.EQ.RSRC1 .AND. MYCOL.EQ.CSRC1 ) THEN
                     IF( RSRC1.NE.RSRC3 .OR. CSRC1.NE.CSRC3 ) THEN
                        CALL DGESD2D( ICTXT, 1, 1, DWORK( IPW2+(DIM1-1)
     $                     *NWIN+DIM1 ), NWIN, RSRC3, CSRC3 )
                     END IF
                  END IF
                  IF( MYROW.EQ.RSRC4 .AND. MYCOL.EQ.CSRC4 ) THEN
                     IF( RSRC4.NE.RSRC2 .OR. CSRC4.NE.CSRC2 ) THEN
                        CALL DGESD2D( ICTXT, DIM1, DIM4, DWORK( IPW2
     $                     +DIM1*NWIN), NWIN, RSRC2, CSRC2 )
                        CALL DGESD2D( ICTXT, DIM1, DIM4, DWORK( IPW2
     $                     +NWIN*NWIN+DIM1*NWIN), NWIN, RSRC2, CSRC2 )
                     END IF
                  END IF
                  IF( MYROW.EQ.RSRC2 .AND. MYCOL.EQ.CSRC2 ) THEN
                     ILOC = INDXG2L( I, NB, MYROW, DESCS( RSRC_ ), NPROW
     $                  )
                     JLOC = INDXG2L( I+DIM1, NB, MYCOL, DESCS( CSRC_ ),
     $                  NPCOL )
                     IF( RSRC2.NE.RSRC4 .OR. CSRC2.NE.CSRC4 ) THEN 

                        CALL DGERV2D( ICTXT, DIM1, DIM4, 
     $                     DWORK(IPW2+DIM1*NWIN), NWIN, RSRC4, CSRC4)

                     END IF
                     CALL DLACPY( 'All', DIM1, DIM4, 
     $                  DWORK( IPW2+DIM1*NWIN ), NWIN,
     $                  S((JLOC-1)*LLDS+ILOC), LLDS )
                     IF( RSRC2.NE.RSRC4 .OR. CSRC2.NE.CSRC4 ) THEN 

                        CALL DGERV2D( ICTXT, DIM1, DIM4, DWORK(IPW2+NWIN
     $                     *NWIN+DIM1*NWIN), NWIN, RSRC4, CSRC4 )

                     END IF
                     CALL DLACPY( 'All', DIM1, DIM4, DWORK( IPW2+NWIN
     $                  *NWIN+DIM1*NWIN ), NWIN, T((JLOC-1)*LLDT+ILOC),
     $                  LLDT )
                  END IF
                  IF( MYROW.EQ.RSRC3 .AND. MYCOL.EQ.CSRC3 ) THEN
                     ILOC = INDXG2L( I+DIM1, NB, MYROW, DESCS( RSRC_ ),
     $                  NPROW )
                     JLOC = INDXG2L( I+DIM1-1, NB, MYCOL, DESCS( CSRC_ )
     $                  , NPCOL )
                     IF( RSRC3.NE.RSRC1 .OR. CSRC3.NE.CSRC1 ) THEN

                        CALL DGERV2D( ICTXT, 1, 1, DWORK( IPW2+(DIM1-1)
     $                     *NWIN+DIM1 ), NWIN, RSRC1, CSRC1 )

                     END IF
                     S((JLOC-1)*LLDS+ILOC) = DWORK( IPW2+(DIM1-1)*NWIN
     $                  +DIM1 )
                  END IF
               END IF
C              
 325           CONTINUE
C              
 320        CONTINUE
C           
C           For the crossborder updates, we use the same directions as in
C           the local reordering case above
C           
            DO 2222 DIR = 1, 2
C              
C              Broadcast informations about the reordering
C              
               DO 321 WINDOW = WINDOW0, WINE, 2
                  RSRC4 = IWORK(IRSRC+WINDOW-1)
                  CSRC4 = IWORK(ICSRC+WINDOW-1)
                  RSRC1 = MOD( RSRC4 - 1 + NPROW, NPROW )
                  CSRC1 = MOD( CSRC4 - 1 + NPCOL, NPCOL )
                  IF( MYROW.EQ.RSRC1 .AND. MYCOL.EQ.CSRC1 ) THEN
                     IF( NPCOL.GT.1 .AND. DIR.EQ.1 ) CALL IGEBS2D( ICTXT
     $                  , 'Row', ' ', 8, 1, IBUFF, 8 )
                     IF( NPROW.GT.1 .AND. DIR.EQ.2 ) CALL IGEBS2D( ICTXT
     $                  , 'Col', ' ', 8, 1, IBUFF, 8 )
                     SKIP1CR = WINDOW.EQ.1 .AND. ICEIL(LIHIC,
     $                  NB).LE.ICEIL(ILO, NB)
                  ELSEIF( MYROW.EQ.RSRC1 .OR. MYCOL.EQ.CSRC1 ) THEN
                     IF( NPCOL.GT.1 .AND. DIR.EQ.1 .AND. MYROW.EQ.RSRC1
     $                  ) THEN

                        CALL IGEBR2D( ICTXT, 'Row', ' ', 8, 1, IBUFF, 8,
     $                     RSRC1, CSRC1 )

                        I = IBUFF( 1 ) 
                        NWIN = IBUFF( 2 ) 
                        PITRAF = IBUFF( 3 ) 
                        KS = IBUFF( 4 ) 
                        PDTRAF = IBUFF( 5 ) 
                        NDTRAF = IBUFF( 6 ) 
                        ILEN = IBUFF( 7 )
                        DLEN = IBUFF( 8 )
                        BUFFLEN = ILEN + 2*DLEN
                        IPW3 = IPW2 + 2*NWIN*NWIN
                        DIM1 = NB - MOD(I-1,NB) 
                        DIM4 = NWIN - DIM1
                        LIHIC = NWIN + I - 1
                        SKIP1CR = WINDOW.EQ.1 .AND. 
     $                     ICEIL(LIHIC,NB).LE.ICEIL(ILO,NB)
                     END IF
                     IF( NPROW.GT.1 .AND. DIR.EQ.2 .AND. MYCOL.EQ.CSRC1
     $                  ) THEN

                        CALL IGEBR2D( ICTXT, 'Col', ' ', 8, 1, IBUFF, 8,
     $                     RSRC1, CSRC1 ) 

                        I = IBUFF( 1 ) 
                        NWIN = IBUFF( 2 ) 
                        PITRAF = IBUFF( 3 ) 
                        KS = IBUFF( 4 ) 
                        PDTRAF = IBUFF( 5 ) 
                        NDTRAF = IBUFF( 6 )
                        ILEN = IBUFF( 7 )
                        DLEN = IBUFF( 8 )
                        BUFFLEN = ILEN + 2*DLEN
                        IPW3 = IPW2 + 2*NWIN*NWIN
                        DIM1 = NB - MOD(I-1,NB) 
                        DIM4 = NWIN - DIM1
                        LIHIC = NWIN + I - 1
                        SKIP1CR = WINDOW.EQ.1 .AND. ICEIL(LIHIC
     $                     ,NB).LE.ICEIL(ILO,NB)
                     END IF
                  END IF
                  IF( RSRC1.NE.RSRC4 ) THEN
                     IF( MYROW.EQ.RSRC4 .AND. MYCOL.EQ.CSRC4 ) THEN
                        IF( NPCOL.GT.1 .AND. DIR.EQ.1 ) CALL IGEBS2D(
     $                     ICTXT, 'Row', ' ', 8, 1, IBUFF, 8 )
                        SKIP1CR = WINDOW.EQ.1 .AND. ICEIL(LIHIC
     $                     ,NB).LE.ICEIL(ILO,NB)
                     ELSEIF( MYROW.EQ.RSRC4 ) THEN
                        IF( NPCOL.GT.1 .AND. DIR.EQ.1 ) THEN

                           CALL IGEBR2D( ICTXT, 'Row', ' ', 8, 1, IBUFF,
     $                        8, RSRC4, CSRC4 )

                           I = IBUFF( 1 ) 
                           NWIN = IBUFF( 2 ) 
                           PITRAF = IBUFF( 3 ) 
                           KS = IBUFF( 4 ) 
                           PDTRAF = IBUFF( 5 ) 
                           NDTRAF = IBUFF( 6 ) 
                           ILEN = IBUFF( 7 )
                           DLEN = IBUFF( 8 )
                           BUFFLEN = ILEN + 2*DLEN
                           IPW3 = IPW2 + 2*NWIN*NWIN
                           DIM1 = NB - MOD(I-1,NB) 
                           DIM4 = NWIN - DIM1
                           LIHIC = NWIN + I - 1
                           SKIP1CR = WINDOW.EQ.1 .AND. ICEIL(LIHIC
     $                        ,NB).LE.ICEIL(ILO,NB)
                        END IF
                     END IF
                  END IF
                  IF( CSRC1.NE.CSRC4 ) THEN
                     IF( MYROW.EQ.RSRC4 .AND. MYCOL.EQ.CSRC4 ) THEN
                        IF( NPROW.GT.1 .AND. DIR.EQ.2 ) CALL IGEBS2D(
     $                     ICTXT, 'Col', ' ', 8, 1, IBUFF, 8 )
                        SKIP1CR = WINDOW.EQ.1 .AND. 
     $                     ICEIL(LIHIC,NB).LE.ICEIL(ILO,NB)
                     ELSEIF( MYCOL.EQ.CSRC4 ) THEN
                        IF( NPROW.GT.1 .AND. DIR.EQ.2 ) THEN

                           CALL IGEBR2D( ICTXT, 'Col', ' ', 8, 1, IBUFF,
     $                        8, RSRC4, CSRC4 )

                           I = IBUFF( 1 ) 
                           NWIN = IBUFF( 2 ) 
                           PITRAF = IBUFF( 3 ) 
                           KS = IBUFF( 4 ) 
                           PDTRAF = IBUFF( 5 ) 
                           NDTRAF = IBUFF( 6 ) 
                           ILEN = IBUFF( 7 )
                           DLEN = IBUFF( 8 )
                           BUFFLEN = ILEN + 2*DLEN
                           IPW3 = IPW2 + 2*NWIN*NWIN
                           DIM1 = NB - MOD(I-1,NB) 
                           DIM4 = NWIN - DIM1
                           LIHIC = NWIN + I - 1
                           SKIP1CR = WINDOW.EQ.1 .AND. ICEIL(LIHIC
     $                        ,NB).LE.ICEIL(ILO,NB)
                        END IF
                     END IF
                  END IF
C                 
C                 Skip rest of broadcasts and updates if appropriate
C                 
                  IF( SKIP1CR ) GO TO 326 
C                 
C                 Broadcast the orthogonal transformations
C                 
                  IF( MYROW.EQ.RSRC1 .AND. MYCOL.EQ.CSRC1 ) THEN
                     BUFFER = PDTRAF+18*NWIN*NWIN
                     BUFFLEN = 2*DLEN + ILEN
                     IF( (NPROW.GT.1 .AND. DIR.EQ.2) .OR. 
     $                  (NPCOL.GT.1 .AND. DIR.EQ.1) ) THEN
                        DO 370 INDX = 1, ILEN
                           DWORK( BUFFER+INDX-1 ) = DBLE( IWORK(IPIW
     $                        +INDX-1) )
 370                    CONTINUE
                        CALL DLACPY( 'All', DLEN, 1, DWORK( IPW3 ), 
     $                     DLEN, DWORK(BUFFER+ILEN), DLEN )
                        CALL DLACPY( 'All', DLEN, 1, DWORK( IPW3+18*NWIN
     $                     *NWIN ), DLEN, DWORK(BUFFER+ILEN+DLEN), DLEN
     $                     )
                     END IF
                     IF( NPCOL.GT.1 .AND. DIR.EQ.1 ) THEN
                        CALL DGEBS2D( ICTXT, 'Row', ' ', BUFFLEN, 1, 
     $                     DWORK(BUFFER), BUFFLEN )
                     END IF
                     IF( NPROW.GT.1 .AND. DIR.EQ.2 ) THEN
                        CALL DGEBS2D( ICTXT, 'Col', ' ', BUFFLEN, 1,
     $                     DWORK(BUFFER), BUFFLEN )
                     END IF
                  ELSEIF( MYROW.EQ.RSRC1 .OR. MYCOL.EQ.CSRC1 ) THEN
                     IF( NPCOL.GT.1 .AND. DIR.EQ.1 .AND. 
     $                  MYROW.EQ.RSRC1 ) 
     $                  THEN
                        BUFFER = PDTRAF+18*NWIN*NWIN
                        BUFFLEN = 2*DLEN + ILEN

                        CALL DGEBR2D( ICTXT, 'Row', ' ', BUFFLEN, 1,
     $                     DWORK(BUFFER), BUFFLEN, RSRC1, CSRC1 )

                     END IF
                     IF( NPROW.GT.1 .AND. DIR.EQ.2 .AND. 
     $                  MYCOL.EQ.CSRC1 ) 
     $                  THEN
                        BUFFER = PDTRAF+18*NWIN*NWIN
                        BUFFLEN = 2*DLEN + ILEN

                        CALL DGEBR2D( ICTXT, 'Col', ' ', BUFFLEN, 1,
     $                     DWORK(BUFFER), BUFFLEN, RSRC1, CSRC1 )

                     END IF
                     IF( (NPCOL.GT.1.AND.DIR.EQ.1.AND.
     $                  MYROW.EQ.RSRC1) 
     $                  .OR.
     $                  (NPROW.GT.1.AND.DIR.EQ.2.AND. 
     $                  MYCOL.EQ.CSRC1) ) 
     $                  THEN
                        DO 380 INDX = 1, ILEN
                           IWORK(IPIW+INDX-1) = INT( DWORK( BUFFER+INDX
     $                        -1 ) )
 380                    CONTINUE
                        CALL DLACPY( 'All', DLEN, 1, DWORK( BUFFER+ILEN
     $                     ), DLEN, DWORK( IPW3 ), DLEN )
                        CALL DLACPY( 'All', DLEN, 1, DWORK( BUFFER+ILEN
     $                     +DLEN ), DLEN, DWORK( IPW3+18*NWIN*NWIN ),
     $                     DLEN )
                     END IF
                  END IF
                  IF( RSRC1.NE.RSRC4 ) THEN
                     IF( MYROW.EQ.RSRC4 .AND. MYCOL.EQ.CSRC4 ) THEN
                        BUFFER = PDTRAF+18*NWIN*NWIN
                        BUFFLEN = 2*DLEN + ILEN
                        IF( NPCOL.GT.1 .AND. DIR.EQ.1 ) THEN
                           DO 390 INDX = 1, ILEN
                              DWORK( BUFFER+INDX-1 ) = DBLE( IWORK(IPIW
     $                           +INDX-1) )
 390                       CONTINUE
                           CALL DLACPY( 'All', DLEN, 1, DWORK( IPW3 ),
     $                        DLEN, DWORK(BUFFER+ILEN), DLEN )
                           CALL DLACPY( 'All', DLEN, 1, DWORK( IPW3+18
     $                        *NWIN*NWIN ), DLEN, DWORK(BUFFER+ILEN
     $                        +DLEN), DLEN )
                           CALL DGEBS2D( ICTXT, 'Row', ' ', BUFFLEN, 1,
     $                        DWORK(BUFFER), BUFFLEN )
                        END IF
                     ELSEIF( MYROW.EQ.RSRC4.AND.DIR.EQ.1.AND.
     $                     NPCOL.GT.1 ) 
     $                     THEN
                        BUFFER = PDTRAF+18*NWIN*NWIN
                        BUFFLEN = 2*DLEN + ILEN

                        CALL DGEBR2D( ICTXT, 'Row', ' ', BUFFLEN, 1,
     $                     DWORK(BUFFER), BUFFLEN, RSRC4, CSRC4 )

                        DO 400 INDX = 1, ILEN
                           IWORK(IPIW+INDX-1) = INT( DWORK( BUFFER+INDX
     $                        -1 ) )
 400                    CONTINUE
                        CALL DLACPY( 'All', DLEN, 1, DWORK( BUFFER+ILEN
     $                     ), DLEN, DWORK( IPW3 ), DLEN )
                        CALL DLACPY( 'All', DLEN, 1, DWORK( BUFFER+ILEN
     $                     +DLEN ), DLEN, DWORK( IPW3+18*NWIN*NWIN ),
     $                     DLEN )
                     END IF
                  END IF
                  IF( CSRC1.NE.CSRC4 ) THEN
                     IF( MYROW.EQ.RSRC4 .AND. MYCOL.EQ.CSRC4 ) THEN
                        BUFFER = PDTRAF+18*NWIN*NWIN
                        BUFFLEN = 2*DLEN + ILEN
                        IF( NPROW.GT.1 .AND. DIR.EQ.2 ) THEN 
                           DO 395 INDX = 1, ILEN
                              DWORK( BUFFER+INDX-1 ) = DBLE( IWORK(IPIW
     $                           +INDX-1) )
 395                       CONTINUE
                           CALL DLACPY( 'All', DLEN, 1, DWORK( IPW3 ),
     $                        DLEN, DWORK(BUFFER+ILEN), DLEN )
                           CALL DLACPY( 'All', DLEN, 1, DWORK( IPW3+18
     $                        *NWIN*NWIN ), DLEN, DWORK(BUFFER+ILEN
     $                        +DLEN), DLEN )
                           CALL DGEBS2D( ICTXT, 'Col', ' ', BUFFLEN, 1,
     $                        DWORK(BUFFER), BUFFLEN )
                        END IF
                     ELSEIF( MYCOL.EQ.CSRC4.AND.DIR.EQ.2.AND.
     $                     NPROW.GT.1 ) 
     $                     THEN
                        BUFFER = PDTRAF+18*NWIN*NWIN
                        BUFFLEN = 2*DLEN + ILEN

                        CALL DGEBR2D( ICTXT, 'Col', ' ', BUFFLEN, 1,
     $                     DWORK(BUFFER), BUFFLEN, RSRC4, CSRC4 )

                        DO 402 INDX = 1, ILEN
                           IWORK(IPIW+INDX-1) = INT( DWORK( BUFFER+INDX
     $                        -1 ) )
 402                    CONTINUE
                        CALL DLACPY( 'All', DLEN, 1, DWORK( BUFFER+ILEN
     $                     ), DLEN, DWORK( IPW3 ), DLEN )
                        CALL DLACPY( 'All', DLEN, 1, DWORK( BUFFER+ILEN
     $                     +DLEN ), DLEN, DWORK( IPW3+18*NWIN*NWIN ),
     $                     DLEN )
                     END IF
                  END IF
C                 
 326              CONTINUE
C                 
 321           CONTINUE
C              
C              Compute crossborder updates
C              
               DO 322 WINDOW = WINDOW0, WINE, 2
                  IF( WINDOW.EQ.1 .AND. SKIP1CR ) GO TO 327
                  RSRC4 = IWORK(IRSRC+WINDOW-1)
                  CSRC4 = IWORK(ICSRC+WINDOW-1)
                  RSRC1 = MOD( RSRC4 - 1 + NPROW, NPROW )
                  CSRC1 = MOD( CSRC4 - 1 + NPCOL, NPCOL )
C                 
C                 Prepare workspaces for updates: 
C                 IPW3 holds now the orthogonal transformations
C                 IPW4 holds the explicit orthogonal matrices, if formed
C                 IPW5 holds the crossborder block columns of (S,T) 
C                 IPW6 holds the crossborder block rows of (S,T)
C                 IPW7 holds the crossborder block columns of Q (if WANTQ=.TRUE.)
C                 and Z (if WANTZ=.TRUE.)
C                 IPW8 points to the leftover workspace used as lhs in matmults 
C                 
                  IF(((MYCOL.EQ.CSRC1.OR.MYCOL.EQ.CSRC4) .AND. 
     $               DIR.EQ.2)
     $               .OR.
     $               ((MYROW.EQ.RSRC1.OR.MYROW.EQ.RSRC4).AND.
     $               DIR.EQ.1)) 
     $               THEN
                     IPW4 = BUFFER 
                     IF( DIR.EQ.2 ) THEN
                        IF( WANTQ .OR. WANTZ ) THEN
                           QZROWS = NUMROC( N, NB, MYROW, 
     $                        DESCQ( RSRC_ ), 
     $                        NPROW )
                        ELSE
                           QZROWS = 0
                        END IF
                        STROWS = NUMROC( I-1, NB, MYROW, 
     $                     DESCS( RSRC_ ), NPROW )
                     ELSE
                        QZROWS = 0
                        STROWS = 0
                     END IF
                     IF( DIR.EQ.1 ) THEN
                        STCOLS = NUMROC( N - (I+DIM1-1), NB, MYCOL, 
     $                     CSRC4,
     $                     NPCOL )
                        IF( MYCOL.EQ.CSRC4 ) STCOLS = STCOLS - DIM4
                     ELSE
                        STCOLS = 0
                     END IF
                     IPW5 = IPW4 +  2 * NWIN*NWIN
                     IPW6 = IPW5 + 2 * STROWS*NWIN
                     IF( WANTQ .OR. WANTZ ) THEN
                        IPW7 = IPW6 + 2 * NWIN*STCOLS
                        IPW8 = IPW7 + 2 * QZROWS*NWIN
                     ELSE
                        IPW8 = IPW6 + 2 * NWIN*STCOLS
                     END IF
                  END IF
C                 
C                 Let each process row and column involved in the updates
C                 exchange data in (S,T), Q and Z with their neighbours 
C                 
                  IF( DIR.EQ.2 ) THEN
                     IF( MYCOL.EQ.CSRC1 .OR. MYCOL.EQ.CSRC4 ) THEN
                        DO 410 INDX = 1, NPROW
                           IF( MYCOL.EQ.CSRC1 ) THEN
                              CALL INFOG2L( 1+(INDX-1)*NB, I, DESCS,
     $                           NPROW, NPCOL, MYROW, MYCOL, ILOC, JLOC1
     $                           , RSRC, CSRC1 )
                              IF( MYROW.EQ.RSRC ) THEN
                                 CALL DLACPY( 'All', STROWS, DIM1,
     $                              S((JLOC1-1)*LLDS+ILOC), LLDS,
     $                              DWORK(IPW5), 2*STROWS )
                                 CALL DLACPY( 'All', STROWS, DIM1,
     $                              T((JLOC1-1)*LLDT+ILOC), LLDT,
     $                              DWORK(IPW5+STROWS), 2*STROWS )
                                 IF( NPCOL.GT.1 ) THEN
                                    EAST = MOD( MYCOL + 1, NPCOL )
                                    CALL DGESD2D( ICTXT, 2*STROWS, DIM1,
     $                                 DWORK(IPW5), 2*STROWS, RSRC, EAST
     $                                 )
                                    CALL DGERV2D( ICTXT, 2*STROWS, DIM4,
     $                                 DWORK(IPW5+2*STROWS*DIM1), 2
     $                                 *STROWS, RSRC, EAST )
                                 END IF
                              END IF
                           END IF
                           IF( MYCOL.EQ.CSRC4 ) THEN
                              CALL INFOG2L( 1+(INDX-1)*NB, I+DIM1, DESCS
     $                           , NPROW, NPCOL, MYROW, MYCOL, ILOC,
     $                           JLOC4, RSRC, CSRC4 )
                              IF( MYROW.EQ.RSRC ) THEN
                                 CALL DLACPY( 'All', STROWS, DIM4,
     $                              S((JLOC4-1)*LLDS+ILOC), LLDS,
     $                              DWORK(IPW5+2*STROWS*DIM1), 2*STROWS
     $                              )
                                 CALL DLACPY( 'All', STROWS, DIM4,
     $                              T((JLOC4-1)*LLDT+ILOC), LLDT,
     $                              DWORK(IPW5+2*STROWS*DIM1+STROWS), 2
     $                              *STROWS )
                                 IF( NPCOL.GT.1 ) THEN
                                    WEST = MOD( MYCOL - 1 + NPCOL, 
     $                                 NPCOL )
                                    CALL DGESD2D( ICTXT, 2*STROWS, DIM4,
     $                                 DWORK(IPW5+2*STROWS*DIM1), 2
     $                                 *STROWS, RSRC, WEST )
                                    CALL DGERV2D( ICTXT, 2*STROWS, DIM1,
     $                                 DWORK(IPW5), 2*STROWS, RSRC, WEST
     $                                 )
                                 END IF
                              END IF
                           END IF
 410                    CONTINUE
                     END IF
                  END IF
C                 
                  IF( DIR.EQ.1 ) THEN
                     IF( MYROW.EQ.RSRC1 .OR. MYROW.EQ.RSRC4 ) THEN
                        DO 420 INDX = 1, NPCOL
                           IF( MYROW.EQ.RSRC1 ) THEN
                              IF( INDX.EQ.1 ) THEN
                                 CALL INFOG2L( I, LIHIC+1, DESCS, NPROW,
     $                              NPCOL, MYROW, MYCOL, ILOC1, JLOC,
     $                              RSRC1, CSRC )
                              ELSE
                                 CALL INFOG2L( I, (ICEIL(LIHIC,NB)+(INDX
     $                              -2))*NB+1, DESCS, NPROW, NPCOL,
     $                              MYROW, MYCOL, ILOC1, JLOC, RSRC1,
     $                              CSRC )
                              END IF
                              IF( MYCOL.EQ.CSRC ) THEN
                                 CALL DLACPY( 'All', DIM1, STCOLS,
     $                              S((JLOC-1)*LLDS+ILOC1), LLDS,
     $                              DWORK(IPW6), NWIN )
                                 CALL DLACPY( 'All', DIM1, STCOLS,
     $                              T((JLOC-1)*LLDT+ILOC1), LLDT,
     $                              DWORK(IPW6+NWIN*STCOLS), NWIN )
                                 IF( NPROW.GT.1 ) THEN
                                    SOUTH = MOD( MYROW + 1, NPROW )
                                    CALL DGESD2D( ICTXT, DIM1, 2*STCOLS,
     $                                 DWORK(IPW6), NWIN, SOUTH, CSRC )
                                    CALL DGERV2D( ICTXT, DIM4, 2*STCOLS,
     $                                 DWORK(IPW6+DIM1), NWIN, SOUTH,
     $                                 CSRC )
                                 END IF
                              END IF
                           END IF
                           IF( MYROW.EQ.RSRC4 ) THEN
                              IF( INDX.EQ.1 ) THEN
                                 CALL INFOG2L( I+DIM1, LIHIC+1, DESCS,
     $                              NPROW, NPCOL, MYROW, MYCOL, ILOC4,
     $                              JLOC, RSRC4, CSRC )
                              ELSE
                                 CALL INFOG2L( I+DIM1, (ICEIL(LIHIC,NB)
     $                              +(INDX-2))*NB+1, DESCS, NPROW, NPCOL
     $                              , MYROW, MYCOL, ILOC4, JLOC, RSRC4,
     $                              CSRC )
                              END IF
                              IF( MYCOL.EQ.CSRC ) THEN
                                 CALL DLACPY( 'All', DIM4, STCOLS,
     $                              S((JLOC-1)*LLDS+ILOC4), LLDS,
     $                              DWORK(IPW6+DIM1), NWIN )
                                 CALL DLACPY( 'All', DIM4, STCOLS,
     $                              T((JLOC-1)*LLDT+ILOC4), LLDT,
     $                              DWORK(IPW6+DIM1+NWIN*STCOLS), NWIN )
                                 IF( NPROW.GT.1 ) THEN
                                    NORTH = MOD( MYROW - 1 + NPROW, 
     $                                 NPROW )
                                    CALL DGESD2D( ICTXT, DIM4, 2*STCOLS,
     $                                 DWORK(IPW6+DIM1), NWIN, NORTH,
     $                                 CSRC )
                                    CALL DGERV2D( ICTXT, DIM1, 2*STCOLS,
     $                                 DWORK(IPW6), NWIN, NORTH, CSRC )
                                 END IF
                              END IF
                           END IF
 420                    CONTINUE
                     END IF
                  END IF
C                 
                  IF( DIR.EQ.2 ) THEN
                     IF( WANTQ ) THEN
                        IF( MYCOL.EQ.CSRC1 .OR. MYCOL.EQ.CSRC4 ) THEN
                           DO 430 INDX = 1, NPROW
                              IF( MYCOL.EQ.CSRC1 ) THEN
                                 CALL INFOG2L( 1+(INDX-1)*NB, I, DESCQ,
     $                              NPROW, NPCOL, MYROW, MYCOL, ILOC,
     $                              JLOC1, RSRC, CSRC1 )
                                 IF( MYROW.EQ.RSRC ) THEN
                                    CALL DLACPY( 'All', QZROWS, DIM1,
     $                                 Q((JLOC1-1)*LLDQ+ILOC), LLDQ,
     $                                 DWORK(IPW7), 2*QZROWS )
                                    CALL DLACPY( 'All', QZROWS, DIM1,
     $                                 Z((JLOC1-1)*LLDZ+ILOC), LLDZ,
     $                                 DWORK(IPW7+QZROWS), 2*QZROWS )
                                    IF( NPCOL.GT.1 ) THEN
                                       EAST = MOD( MYCOL + 1, NPCOL )
                                       CALL DGESD2D( ICTXT, 2*QZROWS,
     $                                    DIM1, DWORK(IPW7), 2*QZROWS,
     $                                    RSRC, EAST )
                                       CALL DGERV2D( ICTXT, 2*QZROWS,
     $                                    DIM4, DWORK(IPW7+2*QZROWS
     $                                    *DIM1), 2*QZROWS, RSRC, EAST )
                                    END IF
                                 END IF
                              END IF
                              IF( MYCOL.EQ.CSRC4 ) THEN
                                 CALL INFOG2L( 1+(INDX-1)*NB, I+DIM1,
     $                              DESCQ, NPROW, NPCOL, MYROW, MYCOL,
     $                              ILOC, JLOC4, RSRC, CSRC4 )
                                 IF( MYROW.EQ.RSRC ) THEN
                                    CALL DLACPY( 'All', QZROWS, DIM4,
     $                                 Q((JLOC4-1)*LLDQ+ILOC), LLDQ,
     $                                 DWORK(IPW7+2*QZROWS*DIM1), 2
     $                                 *QZROWS)
                                    CALL DLACPY( 'All', QZROWS, DIM4,
     $                                 Z((JLOC4-1)*LLDZ+ILOC), LLDZ,
     $                                 DWORK(IPW7+2*QZROWS*DIM1+ QZROWS)
     $                                 , 2*QZROWS )
                                    IF( NPCOL.GT.1 ) THEN
                                       WEST = MOD( MYCOL - 1 + NPCOL,
     $                                    NPCOL )
                                       CALL DGESD2D( ICTXT, 2*QZROWS,
     $                                    DIM4, DWORK(IPW7+2*QZROWS
     $                                    *DIM1), 2*QZROWS, RSRC, WEST )
                                       CALL DGERV2D( ICTXT, 2*QZROWS,
     $                                    DIM1, DWORK(IPW7), 2*QZROWS,
     $                                    RSRC, WEST )
                                       
                                    END IF
                                 END IF
                              END IF
 430                       CONTINUE
                        END IF
                     END IF
                  END IF
C                 
 327              CONTINUE
C                 
 322           CONTINUE
C              
               DO 323 WINDOW = WINDOW0, WINE, 2
                  RSRC4 = IWORK(IRSRC+WINDOW-1)
                  CSRC4 = IWORK(ICSRC+WINDOW-1)
                  RSRC1 = MOD( RSRC4 - 1 + NPROW, NPROW )
                  CSRC1 = MOD( CSRC4 - 1 + NPCOL, NPCOL )
                  FLOPS = 0
                  IF(((MYCOL.EQ.CSRC1.OR.MYCOL.EQ.CSRC4) .AND. DIR.EQ.2)
     $               .OR. ((MYROW.EQ.RSRC1.OR.MYROW.EQ.RSRC4) .AND.
     $               DIR.EQ.1)) THEN
C                    
C                    Skip this part of the updates if appropriate
C                    
                     IF( WINDOW.EQ.1 .AND. SKIP1CR ) GO TO 328
C                    
C                    Count number of operations to decide whether to use matrix-matrix 
C                    multiplications for updating off-diagonal parts or not.
C                    
                     NITRAF = PITRAF - IPIW
c                    ISHH = .FALSE.
                     DO 405 K = 1, NITRAF
                        IF( IWORK( IPIW + K - 1 ).LE.NWIN ) THEN
                           FLOPS = FLOPS + 6
                        ELSE
                           FLOPS = FLOPS + 11
c                          ISHH = .TRUE.
                        END IF
 405                 CONTINUE
                     PDU = IPW4
                     PDV = IPW4 + NWIN*NWIN
C                    
C                    Perform updates in parallel
C                    
                     IF( FLOPS.NE.0 .AND. ( 2*FLOPS*100 )/( 2*NWIN*NWIN)
     $                  .GE. MMULT ) THEN
C                       
                        CALL DLASET( 'All', NWIN, NWIN, ZERO, ONE,
     $                     DWORK( PDU ), NWIN )
                        CALL BDLAGPP( 1, NWIN, NWIN, NCB, DWORK( PDU ),
     $                     NWIN, NITRAF, IWORK(IPIW), DWORK( IPW3 ),
     $                     DWORK(IPW8) )
                        CALL DLASET( 'All', NWIN, NWIN, ZERO, ONE,
     $                     DWORK( PDV ), NWIN )
                        CALL BDLAGPP( 1, NWIN, NWIN, NCB, DWORK( PDV ),
     $                     NWIN, NITRAF, IWORK(IPIW), DWORK( IPW3+18
     $                     *NWIN*NWIN ), DWORK(IPW8) )
C                       
C                       Test if sparsity structure of orthogonal matrix can be exploited
C                       (see below) 
C                       
                        IF( ISHH .OR. DIM1.NE.KS .OR. DIM4. NE. KS) THEN

C                          
C                          Update the columns of (S,T), Q and Z affected by the reordering
C                          
                           IF( DIR.EQ.2 ) THEN
                              DO 440 INDX = 1, MIN(I-1,1+(NPROW-1)*NB), 
     $                           NB
                                 IF( MYCOL.EQ.CSRC1 ) THEN
                                    CALL INFOG2L( INDX, I, DESCS, NPROW
     $                                 ,NPCOL, MYROW, MYCOL,ILOC, JLOC,
     $                                 RSRC, CSRC1 )
                                    IF( MYROW.EQ.RSRC ) THEN
                                       CALL DGEMM( 'N', 'N', 2*STROWS,
     $                                    DIM1, NWIN, ONE, DWORK( IPW5 )
     $                                    , 2*STROWS, DWORK( PDV ), NWIN
     $                                    , ZERO, DWORK(IPW8), 2*STROWS)
                                       CALL DLACPY( 'All', STROWS, DIM1
     $                                    ,DWORK(IPW8), 2*STROWS,S((JLOC
     $                                    -1)*LLDS+ILOC),LLDS )
                                       CALL DLACPY( 'All', STROWS, DIM1
     $                                    ,DWORK(IPW8+STROWS),2*STROWS
     $                                    ,T((JLOC-1)*LLDT+ILOC),LLDT )
                                    END IF
                                 END IF
                                 IF( MYCOL.EQ.CSRC4 ) THEN
                                    CALL INFOG2L( INDX, I+DIM1, DESCS,
     $                                 NPROW, NPCOL, MYROW, MYCOL, ILOC,
     $                                 JLOC, RSRC, CSRC4 )
                                    IF( MYROW.EQ.RSRC ) THEN
                                       CALL DGEMM( 'N', 'N', 2*STROWS,
     $                                    DIM4, NWIN, ONE, DWORK( IPW5 )
     $                                    , 2*STROWS, DWORK( PDV+NWIN
     $                                    *DIM1 ), NWIN, ZERO,
     $                                    DWORK(IPW8), 2*STROWS )
                                       CALL DLACPY( 'All', STROWS, DIM4
     $                                    ,DWORK(IPW8), 2*STROWS,S((JLOC
     $                                    -1)*LLDS+ILOC),LLDS )
                                       CALL DLACPY( 'All', STROWS, DIM4
     $                                    ,DWORK(IPW8+STROWS),2*STROWS
     $                                    ,T((JLOC-1)*LLDT+ILOC),LLDT ) 
                                    END IF
                                 END IF
 440                          CONTINUE
C                             
                              IF( WANTQ ) THEN
                                 DO 450 INDX = 1, MIN(N,1+(NPROW-1)*NB)
     $                              ,NB
                                    IF( MYCOL.EQ.CSRC1 ) THEN
                                       CALL INFOG2L( INDX, I, DESCQ,
     $                                    NPROW, NPCOL, MYROW, MYCOL,
     $                                    ILOC, JLOC, RSRC, CSRC1 )
                                       IF( MYROW.EQ.RSRC ) THEN
                                          CALL DGEMM( 'N', 'N', QZROWS,
     $                                       DIM1, NWIN, ONE,
     $                                       DWORK(IPW7), 2*QZROWS,
     $                                       DWORK(PDU ), NWIN, ZERO,
     $                                       DWORK(IPW8), QZROWS )
                                          CALL DLACPY( 'All', QZROWS,
     $                                       DIM1, DWORK(IPW8), QZROWS,
     $                                       Q((JLOC-1)*LLDQ+ILOC),
     $                                       LLDQ)
                                          CALL DGEMM( 'N', 'N', QZROWS,
     $                                       DIM1, NWIN, ONE, DWORK(IPW7
     $                                       +QZROWS ), 2*QZROWS,
     $                                       DWORK(PDV ), NWIN, ZERO,
     $                                       DWORK(IPW8), QZROWS )
                                          CALL DLACPY( 'All', QZROWS,
     $                                       DIM1, DWORK(IPW8), QZROWS,
     $                                       Z((JLOC-1)*LLDZ+ILOC),
     $                                       LLDZ)
                                       END IF
                                    END IF
                                    IF( MYCOL.EQ.CSRC4 ) THEN
                                       CALL INFOG2L( INDX, I+DIM1, DESCQ
     $                                    , NPROW, NPCOL, MYROW, MYCOL,
     $                                    ILOC, JLOC, RSRC, CSRC4 )
                                       IF( MYROW.EQ.RSRC ) THEN
                                          CALL DGEMM( 'N', 'N', QZROWS,
     $                                       DIM4, NWIN, ONE,
     $                                       DWORK(IPW7), 2*QZROWS,
     $                                       DWORK(PDU+NWIN *DIM1 ),
     $                                       NWIN, ZERO , DWORK(IPW8),
     $                                       QZROWS )
                                          CALL DLACPY( 'All', QZROWS,
     $                                       DIM4, DWORK(IPW8), QZROWS,
     $                                       Q((JLOC-1)*LLDQ+ILOC),
     $                                       LLDQ)
                                          CALL DGEMM( 'N', 'N', QZROWS,
     $                                       DIM4, NWIN, ONE, DWORK(IPW7
     $                                       +QZROWS ), 2*QZROWS,
     $                                       DWORK(PDV+NWIN*DIM1 ), NWIN
     $                                       , ZERO , DWORK(IPW8),
     $                                       QZROWS )
                                          CALL DLACPY( 'All', QZROWS,
     $                                       DIM4, DWORK(IPW8), QZROWS,
     $                                       Z((JLOC-1)*LLDZ+ILOC),
     $                                       LLDZ)
                                       END IF
                                    END IF
 450                             CONTINUE
                              END IF
                           END IF
C                          
C                          Update the rows of T affected by the reordering
C                          
                           IF( DIR.EQ.1 ) THEN
                              IF ( LIHIC.LT.N ) THEN
                                 IF( MYROW.EQ.RSRC1.AND.MYCOL.EQ.CSRC4
     $                              .AND.
     $                              MOD(LIHIC,NB).NE.0 ) THEN
                                    INDX = LIHIC + 1
                                    CALL INFOG2L( I, INDX, DESCS, NPROW
     $                                 ,NPCOL, MYROW, MYCOL,ILOC, JLOC
     $                                 ,RSRC1, CSRC4 )
                                    CALL DGEMM( 'T', 'N', DIM1, 2*STCOLS
     $                                 , NWIN, ONE, DWORK(PDU), NWIN,
     $                                 DWORK( IPW6 ), NWIN, ZERO,
     $                                 DWORK(IPW8), DIM1 )
                                    CALL DLACPY( 'All', DIM1, STCOLS,
     $                                 DWORK(IPW8), DIM1, S((JLOC-1)
     $                                 *LLDS+ILOC), LLDS )
                                    CALL DLACPY( 'All', DIM1, STCOLS,
     $                                 DWORK(IPW8+DIM1*STCOLS), DIM1,
     $                                 T((JLOC-1)*LLDT+ILOC), LLDT )
                                 END IF
                                 IF( MYROW.EQ.RSRC4.AND.MYCOL.EQ.CSRC4
     $                              .AND.
     $                              MOD(LIHIC,NB).NE.0 ) THEN
                                    INDX = LIHIC + 1
                                    CALL INFOG2L( I+DIM1, INDX, DESCS,
     $                                 NPROW, NPCOL, MYROW, MYCOL, ILOC,
     $                                 JLOC, RSRC4, CSRC4 )
                                    CALL DGEMM( 'T', 'N', DIM4, 2*STCOLS
     $                                 , NWIN, ONE, DWORK( PDU+DIM1
     $                                 *NWIN), NWIN, DWORK( IPW6), NWIN,
     $                                 ZERO , DWORK(IPW8), DIM4 )
                                    CALL DLACPY( 'All', DIM4, STCOLS,
     $                                 DWORK(IPW8), DIM4, S((JLOC-1)
     $                                 *LLDS+ILOC), LLDS )
                                    CALL DLACPY( 'All', DIM4, STCOLS,
     $                                 DWORK(IPW8+DIM4*STCOLS), DIM4,
     $                                 T((JLOC-1)*LLDT+ILOC), LLDT )
                                 END IF
                                 INDXS = ICEIL(LIHIC,NB)*NB + 1 
                                 INDXE = MIN(N,INDXS+(NPCOL-2)*NB)
                                 DO 460 INDX = INDXS, INDXE, NB
                                    IF( MYROW.EQ.RSRC1 ) THEN
                                       CALL INFOG2L( I, INDX, DESCS,
     $                                    NPROW, NPCOL, MYROW, MYCOL,
     $                                    ILOC, JLOC, RSRC1, CSRC )
                                       IF( MYCOL.EQ.CSRC ) THEN
                                          CALL DGEMM( 'T', 'N', DIM1, 2
     $                                       *STCOLS, NWIN, ONE,
     $                                       DWORK(PDU ), NWIN, DWORK(
     $                                       IPW6 ), NWIN, ZERO,
     $                                       DWORK(IPW8), DIM1 )
                                          CALL DLACPY( 'All', DIM1,
     $                                       STCOLS, DWORK(IPW8), DIM1,
     $                                       S((JLOC-1)*LLDS+ILOC),
     $                                       LLDS)
                                          CALL DLACPY( 'All', DIM1,
     $                                       STCOLS, DWORK(IPW8+DIM1
     $                                       *STCOLS), DIM1, T((JLOC-1)
     $                                       *LLDT+ILOC), LLDT )
                                       END IF
                                    END IF
                                    IF( MYROW.EQ.RSRC4 ) THEN
                                       CALL INFOG2L( I+DIM1, INDX, DESCS
     $                                    , NPROW, NPCOL, MYROW, MYCOL,
     $                                    ILOC, JLOC, RSRC4, CSRC )
                                       IF( MYCOL.EQ.CSRC ) THEN
                                          CALL DGEMM( 'T', 'N', DIM4, 2
     $                                       *STCOLS, NWIN, ONE,
     $                                       DWORK(PDU+NWIN*DIM1 ), NWIN
     $                                       , DWORK( IPW6 ), NWIN, ZERO
     $                                       , DWORK(IPW8), DIM4 )
                                          CALL DLACPY( 'All', DIM4,
     $                                       STCOLS, DWORK(IPW8), DIM4,
     $                                       S((JLOC-1)*LLDS+ILOC),
     $                                       LLDS)
                                          CALL DLACPY( 'All', DIM4,
     $                                       STCOLS, DWORK(IPW8+DIM4
     $                                       *STCOLS), DIM4, T((JLOC-1)
     $                                       *LLDT+ILOC), LLDT )
                                       END IF
                                    END IF
 460                             CONTINUE
                              END IF
                           END IF
                        ELSE
                           IF (MYROW.EQ.0.AND.MYCOL.EQ.0) THEN
                              WRITE(*,*)
     $                           'Reorder2:Not supposed to happen'
                           END IF
                           CALL BLACS_BARRIER(ICTXT, 'A')
                           STOP
C                          
                        END IF
                     ELSEIF( FLOPS.NE.0 ) THEN

C                       
C                       Update off-diagonal blocks of (S,T), Q and Z using the pipelined 
C                       elementary transformations. Now we have a delicate problem: how to 
C                       do this without redundant work? For now, we let the processes 
C                       involved compute the whole crossborder block rows and column saving 
C                       only the part belonging to the corresponding side of the border. To
C                       make this a realistic alternative, we have modified the ratio
C                       r_flops (see Kressner paper above) to give more favor to the
C                       ordinary matrix multiplication.
C                       
                        IF( DIR.EQ.2 ) THEN
                           INDXE =  MIN(I-1,1+(NPROW-1)*NB)
                           DO 500 INDX = 1, INDXE, NB
                              CALL INFOG2L( INDX, I, DESCS, NPROW, NPCOL
     $                           , MYROW, MYCOL, ILOC, JLOC, RSRC, CSRC
     $                           )
                              IF( MYROW.EQ.RSRC .AND. MYCOL.EQ.CSRC )
     $                           THEN
                                 CALL BDLAGPP( 1, 2*STROWS, NWIN, NCB,
     $                              DWORK(IPW5), 2*STROWS, NITRAF,
     $                              IWORK(IPIW), DWORK( IPW3+18*NWIN
     $                              *NWIN), DWORK(IPW8) )
                                 CALL DLACPY( 'All', STROWS, DIM1,
     $                              DWORK(IPW5), STROWS, S((JLOC-1)*LLDS
     $                              +ILOC ), LLDS )
                                 CALL DLACPY( 'All', STROWS, DIM1,
     $                              DWORK(IPW5+STROWS), 2*STROWS,
     $                              T((JLOC-1)*LLDT+ILOC ), LLDT )
                              END IF
                              CALL INFOG2L( INDX, I+DIM1, DESCS, NPROW,
     $                           NPCOL, MYROW, MYCOL, ILOC, JLOC, RSRC,
     $                           CSRC )
                              IF( MYROW.EQ.RSRC .AND. MYCOL.EQ.CSRC ) 
     $                           THEN
                                 IF( NPCOL.GT.1 ) CALL BDLAGPP( 1, 2
     $                              *STROWS, NWIN, NCB, DWORK(IPW5), 2
     $                              *STROWS, NITRAF, IWORK(IPIW),
     $                              DWORK(IPW3+18*NWIN*NWIN),
     $                              DWORK(IPW8) )
                                 CALL DLACPY( 'All', STROWS, DIM4,
     $                              DWORK(IPW5+2*STROWS*DIM1), 2*STROWS,
     $                              S((JLOC-1)*LLDS+ILOC ), LLDS )
                                 CALL DLACPY( 'All', STROWS, DIM4,
     $                              DWORK(IPW5+STROWS+2*STROWS*DIM1), 2
     $                              *STROWS, T((JLOC-1)*LLDT+ILOC ),
     $                              LLDT )
                              END IF
 500                       CONTINUE
                           IF( WANTQ ) THEN
                              INDXE = MIN(N,1+(NPROW-1)*NB)
                              DO 510 INDX = 1, INDXE, NB
                                 CALL INFOG2L( INDX, I, DESCQ, NPROW,
     $                              NPCOL, MYROW, MYCOL, ILOC, JLOC,
     $                              RSRC, CSRC )
                                 IF( MYROW.EQ.RSRC .AND. MYCOL.EQ.CSRC)
     $                              THEN
                                    CALL BDLAGPP( 1, QZROWS, NWIN, NCB,
     $                                 DWORK(IPW7), 2*QZROWS, NITRAF,
     $                                 IWORK(IPIW), DWORK( IPW3),
     $                                 DWORK(IPW8) )
                                    CALL DLACPY( 'All', QZROWS, DIM1,
     $                                 DWORK(IPW7), 2*QZROWS, Q((JLOC-1)
     $                                 *LLDQ+ILOC ), LLDQ )
                                 END IF
                                 CALL INFOG2L( INDX, I+DIM1, DESCQ,
     $                              NPROW, NPCOL, MYROW, MYCOL, ILOC,
     $                              JLOC, RSRC, CSRC )
                                 IF( MYROW.EQ.RSRC .AND. MYCOL.EQ.CSRC)
     $                              THEN
                                    IF( NPCOL.GT.1 ) CALL BDLAGPP( 1,
     $                                 QZROWS, NWIN, NCB, DWORK(IPW7), 2
     $                                 *QZROWS, NITRAF, IWORK(IPIW),
     $                                 DWORK( IPW3 ), DWORK(IPW8) )
                                    CALL DLACPY( 'All', QZROWS, DIM4,
     $                                 DWORK(IPW7+2*QZROWS*DIM1), 2
     $                                 *QZROWS, Q((JLOC-1)*LLDQ+ILOC ),
     $                                 LLDQ )
                                 END IF  
 510                          CONTINUE
                           END IF
                           IF( WANTZ ) THEN
                              INDXE = MIN(N,1+(NPROW-1)*NB)
                              DO 515 INDX = 1, INDXE, NB
                                 CALL INFOG2L( INDX, I, DESCZ, NPROW,
     $                              NPCOL, MYROW, MYCOL, ILOC, JLOC,
     $                              RSRC, CSRC )
                                 IF( MYROW.EQ.RSRC .AND. MYCOL.EQ.CSRC)
     $                              THEN
                                    CALL BDLAGPP( 1, QZROWS, NWIN, NCB,
     $                                 DWORK(IPW7+QZROWS), 2*QZROWS,
     $                                 NITRAF, IWORK( IPIW ), DWORK(IPW3
     $                                 +18*NWIN*NWIN ), DWORK( IPW8) )
                                    CALL DLACPY( 'All', QZROWS, DIM1,
     $                                 DWORK(IPW7+QZROWS), 2*QZROWS,
     $                                 Z((JLOC-1)*LLDZ+ILOC ), LLDZ )
                                 END IF
                                 CALL INFOG2L( INDX, I+DIM1, DESCZ,
     $                              NPROW, NPCOL, MYROW, MYCOL, ILOC,
     $                              JLOC, RSRC, CSRC )
                                 IF( MYROW.EQ.RSRC .AND. MYCOL.EQ.CSRC)
     $                              THEN
                                    IF( NPCOL.GT.1 ) CALL BDLAGPP( 1,
     $                                 QZROWS, NWIN, NCB, DWORK(IPW7
     $                                 +QZROWS), 2*QZROWS, NITRAF,
     $                                 IWORK( IPIW ), DWORK( IPW3+18
     $                                 *NWIN*NWIN ), DWORK( IPW8 ) )
                                    CALL DLACPY( 'All', QZROWS, DIM4,
     $                                 DWORK( IPW7+QZROWS+ 2*QZROWS*DIM1
     $                                 ), 2*QZROWS, Z((JLOC-1)*LLDZ+ILOC
     $                                 ), LLDZ )
                                 END IF  
 515                          CONTINUE
                           END IF
                        END IF
C                       
                        IF( DIR.EQ.1 ) THEN
                           IF( LIHIC.LT.N ) THEN
                              INDX = LIHIC + 1
                              CALL INFOG2L( I, INDX, DESCS, NPROW, NPCOL
     $                           , MYROW, MYCOL, ILOC, JLOC, RSRC, CSRC
     $                           )
                              IF( MYROW.EQ.RSRC .AND. MYCOL.EQ.CSRC.AND.
     $                           MOD(LIHIC,NB).NE.0 ) THEN
                                 CALL BDLAGPP( 0, NWIN, 2*STCOLS, NCB,
     $                              DWORK( IPW6 ), NWIN, NITRAF,
     $                              IWORK(IPIW), DWORK( IPW3 ),
     $                              DWORK(IPW8) )
                                 CALL DLACPY( 'All', DIM1, STCOLS,
     $                              DWORK( IPW6 ), NWIN, S((JLOC-1)*LLDS
     $                              +ILOC), LLDS )
                                 CALL DLACPY( 'All', DIM1, STCOLS, 
     $                              DWORK( IPW6+NWIN*STCOLS ), NWIN, 
     $                              T((JLOC-1)*LLDT+ILOC), LLDT )
                              END IF
                              CALL INFOG2L( I+DIM1, INDX, DESCS, NPROW,
     $                           NPCOL, MYROW, MYCOL, ILOC, JLOC, RSRC,
     $                           CSRC )
                              IF( MYROW.EQ.RSRC .AND. MYCOL.EQ.CSRC.AND.
     $                           MOD(LIHIC,NB).NE.0 ) THEN
                                 IF( NPROW.GT.1 ) CALL BDLAGPP( 0, NWIN,
     $                              2*STCOLS, NCB, DWORK( IPW6 ), NWIN,
     $                              NITRAF, IWORK(IPIW), DWORK( IPW3 ),
     $                              DWORK(IPW8) )
                                 CALL DLACPY( 'All', DIM4, STCOLS,
     $                              DWORK( IPW6+DIM1 ), NWIN, S((JLOC-1)
     $                              *LLDS+ILOC), LLDS )
                                 CALL DLACPY( 'All', DIM4, STCOLS,
     $                              DWORK( IPW6+NWIN*STCOLS+DIM1 ), NWIN
     $                              , T((JLOC-1)*LLDT+ILOC), LLDT )
                              END IF
                              INDXS = ICEIL(LIHIC,NB)*NB + 1
                              INDXE = MIN(N,INDXS+(NPCOL-2)*NB)
                              DO 520 INDX = INDXS, INDXE, NB
                                 CALL INFOG2L( I, INDX, DESCS, NPROW,
     $                              NPCOL, MYROW, MYCOL, ILOC, JLOC,
     $                              RSRC, CSRC )
                                 IF( MYROW.EQ.RSRC .AND. MYCOL.EQ.CSRC)
     $                              THEN
                                    CALL BDLAGPP( 0, NWIN, 2*STCOLS, NCB
     $                                 , DWORK(IPW6), NWIN, NITRAF,
     $                                 IWORK(IPIW), DWORK( IPW3 ),
     $                                 DWORK(IPW8) )
                                    CALL DLACPY( 'All', DIM1, STCOLS,
     $                                 DWORK( IPW6 ), NWIN, S((JLOC-1)
     $                                 *LLDS+ILOC), LLDS )
                                    CALL DLACPY( 'All', DIM1, STCOLS,
     $                                 DWORK( IPW6+NWIN*STCOLS ), NWIN,
     $                                 T((JLOC-1)*LLDT+ILOC), LLDT )
                                 END IF
                                 CALL INFOG2L( I+DIM1, INDX, DESCS,
     $                              NPROW, NPCOL, MYROW, MYCOL, ILOC,
     $                              JLOC, RSRC, CSRC )
                                 IF( MYROW.EQ.RSRC .AND. MYCOL.EQ.CSRC)
     $                              THEN
                                    IF( NPROW.GT.1 ) CALL BDLAGPP( 0,
     $                                 NWIN, 2*STCOLS, NCB, DWORK(IPW6),
     $                                 NWIN, NITRAF, IWORK(IPIW),
     $                                 DWORK(IPW3 ), DWORK(IPW8) )
                                    CALL DLACPY( 'All', DIM4, STCOLS,
     $                                 DWORK( IPW6+DIM1 ), NWIN, S((JLOC
     $                                 -1)*LLDS+ILOC), LLDS )
                                    CALL DLACPY( 'All', DIM4, STCOLS,
     $                                 DWORK( IPW6+NWIN*STCOLS+DIM1 ),
     $                                 NWIN, T((JLOC-1)*LLDT+ILOC),
     $                                 LLDT)
                                 END IF
 520                          CONTINUE
                           END IF
                        END IF
                     END IF
                  END IF
C                 
 328              CONTINUE

C                 
 323           CONTINUE
C              
C              End of loops over directions (DIR)
C              
 2222       CONTINUE
C           
C           End of loops over diagonal blocks for reordering over the block
C           diagonal
C           
 310     CONTINUE
         LAST = LAST + 1
         IF( LASTWAIT .AND. LAST.LT.2 ) GO TO 308
C        
C        Barrier to collect the processes before proceeding
C        
         CALL BLACS_BARRIER( ICTXT, 'All' )
C        
C        
C        Compute global maximum of IERR so that we know if some process
C        experienced a failure in the reordering
C        
         MYIERR = IERR
         IF( NPROCS.GT.1 ) CALL IGAMX2D( ICTXT, 'All', ' ', 1, 1, IERR,
     $      1, -1, -1, -1, -1, -1 )
C        
         IF( IERR.EQ.1 .OR. IERR.EQ.2 ) THEN

C           
C           When calling BDTREXC, the block at position I+KKS-1 failed to swap.
C           
            IF( MYIERR.EQ.1 .OR. MYIERR.EQ.2 ) THEN
               INFO = I+KKS-1
            ELSE 
               INFO = 1
            END IF
            IF( NPROCS.GT.1 ) CALL IGAMX2D( ICTXT, 'All', ' ', 1, 1,
     $         INFO, 1, -1, -1, -1, -1, -1 )
            GO TO 300
         END IF
C        
C        Do a global update of the SELECT vector
C        
         DO 530 K = 1, N
            RSRC = INDXG2P( K, NB, MYROW, DESCS( RSRC_ ), NPROW )
            CSRC = INDXG2P( K, NB, MYCOL, DESCS( CSRC_ ), NPCOL )
            IF( MYROW.NE.RSRC .OR. MYCOL.NE.CSRC ) SELECT( K ) = 0
 530     CONTINUE 
         IF( NPROCS.GT.1 ) CALL IGSUM2D( ICTXT, 'All', ' ', N, 1, SELECT
     $      , N, -1, -1)

C        
C        Find the global minumum of ILO and IHI
C        
         ILO = ILO - 1
 523     CONTINUE
         ILO = ILO + 1
         IF( ILO .LE. N ) THEN
            IF( INT2LG(SELECT(ILO)) ) GO TO 523
         END IF
         IHI = IHI + 1
 527     CONTINUE
         IHI = IHI - 1
         IF( IHI .GE. 1 ) THEN
            IF( .NOT. INT2LG(SELECT(IHI)) ) GO TO 527
         END IF


         GO TO 50
      END IF
C     
C     
 300  CONTINUE
C     
C     In case an error occured, do an additional global update of SELECT
C     

      IF( INFO.NE.0 ) THEN
         DO 540 K = 1, N
            RSRC = INDXG2P( K, NB, MYROW, DESCS( RSRC_ ), NPROW )
            CSRC = INDXG2P( K, NB, MYCOL, DESCS( CSRC_ ), NPCOL )
            IF( MYROW.NE.RSRC .OR. MYCOL.NE.CSRC ) 
     $         SELECT( K ) = 0
 540     CONTINUE
         IF( NPROCS.GT.1 )
     $      CALL IGSUM2D( ICTXT, 'All', ' ', N, 1, SELECT, N, -1, -1)
      END IF   
C     
 545  CONTINUE
C     
C     Store the output eigenvalues in ALPHA[R,I] and BETA: first let 
C     all the processes compute the eigenvalue inside their diagonal 
C     blocks in parallel, except for the eigenvalue located next to a 
C     block border. After that, compute all eigenvalues located next to 
C     the block borders. Finally, do a global summation over so that 
C     all processors receive the result. Notice: real eigenvalues
C     extracted from a non-canonical 2-by-2 block are not stored in
C     in any particular order
C     
      DO 550 K = 1, N
         ALPHAR( K ) = ZERO
         ALPHAI( K ) = ZERO
         BETA( K ) = ZERO
 550  CONTINUE

      PAIR = .FALSE.
      DO 560 K = 1, N
         IF( .NOT. PAIR ) THEN            
            BORDER = ( K.NE.N .AND. MOD( K, NB ).EQ.0 ) .OR.
     $         ( K.NE.1 .AND. MOD( K, NB ).EQ.1 )

            IF( .NOT. BORDER ) THEN
               CALL INFOG2L( K, K, DESCS, NPROW, NPCOL, MYROW,MYCOL,
     $            ILOC1, JLOC1, RSRC1, CSRC1 )
               IF( MYROW.EQ.RSRC1 .AND. MYCOL.EQ.CSRC1 ) THEN
                  IF( K.LT.N ) THEN
                     DWORK(2) = S((JLOC1-1)*LLDS+ILOC1+1)
                  ELSE
                     DWORK(2) = ZERO
                  END IF
                  IF( DWORK(2).NE.ZERO ) THEN
                     DWORK(1) = S((JLOC1-1)*LLDS+ILOC1)
                     DWORK(3) = S((JLOC1)*LLDS+ILOC1)
                     DWORK(4) = S((JLOC1)*LLDS+ILOC1+1)
                     DWORK(5) = T((JLOC1-1)*LLDT+ILOC1)
                     DWORK(6) = T((JLOC1-1)*LLDT+ILOC1+1)
                     DWORK(7) = T((JLOC1)*LLDT+ILOC1)
                     DWORK(8) = T((JLOC1)*LLDT+ILOC1+1)
                     CALL DLAG2( DWORK, 2, DWORK( 5 ), 2,SMLNUM*EPS,
     $                  BETA( K ), BETA( K+1 ), ALPHAR( K ), ALPHAR( K+1
     $                  ), ALPHAI( K ) )
                     ALPHAI( K+1 ) = -ALPHAI( K )
                     PAIR = .TRUE.
                  ELSE 
                     IF( K.GT.1 ) THEN
                        
                        DWORK(2) = S((JLOC1-2)*LLDS+ILOC1)
                        IF( DWORK(2).NE.ZERO ) THEN
                           DWORK(1) = S((JLOC1-2)*LLDS+ILOC1-1)
                           DWORK(3) = S((JLOC1-1)*LLDS+ILOC1-1)
                           DWORK(4) = S((JLOC1-1)*LLDS+ILOC1)
                           DWORK(5) = T((JLOC1-2)*LLDT+ILOC1-1)
                           DWORK(6) = T((JLOC1-2)*LLDT+ILOC1)
                           DWORK(7) = T((JLOC1-1)*LLDT+ILOC1-1)
                           DWORK(8) = T((JLOC1-1)*LLDT+ILOC1)
                           CALL DLAG2( DWORK, 2, DWORK( 5 ), 2, SMLNUM
     $                        *EPS, BETA( K-1 ), BETA( K ), ALPHAR( K-1
     $                        ), ALPHAR( K ), ALPHAI( K-1 ) )
                           ALPHAI( K ) = -ALPHAI( K-1 )
                        ELSE
                           ALPHAR( K ) = S((JLOC1-1)*LLDS+ILOC1)
                           ALPHAI( K ) = ZERO
                           BETA( K ) = T((JLOC1-1)*LLDT+ILOC1) 
                        END IF
                     ELSE
                        ALPHAR( K ) = S((JLOC1-1)*LLDS+ILOC1)
                        ALPHAI( K ) = ZERO
                        BETA( K ) = T((JLOC1-1)*LLDT+ILOC1) 
                     END IF
                  END IF
               END IF
            END IF
         ELSE
            PAIR = .FALSE.
         END IF
 560  CONTINUE

      DO 570 K = NB, N-1, NB
         CALL INFOG2L( K, K, DESCS, NPROW, NPCOL, MYROW, MYCOL, ILOC1,
     $      JLOC1, RSRC1, CSRC1 )
         CALL INFOG2L( K+1, K, DESCS, NPROW, NPCOL, MYROW, MYCOL, ILOC2,
     $      JLOC2, RSRC2, CSRC2 )
         CALL INFOG2L( K, K+1, DESCS, NPROW, NPCOL, MYROW, MYCOL, ILOC3,
     $      JLOC3, RSRC3, CSRC3 )
         CALL INFOG2L( K+1, K+1, DESCS, NPROW, NPCOL, MYROW, MYCOL,
     $      ILOC4, JLOC4, RSRC4, CSRC4 )
         IF( MYROW.EQ.RSRC2 .AND. MYCOL.EQ.CSRC2 ) THEN
            DWORK(9) = S((JLOC2-1)*LLDS+ILOC2)
            DWORK(10) = T((JLOC2-1)*LLDT+ILOC2)
            IF( RSRC1.NE.RSRC2 .OR. CSRC1.NE.CSRC2 ) 
     $         CALL DGESD2D( ICTXT, 2, 1, DWORK(9), 2, RSRC1, CSRC1 )
         END IF
         IF( MYROW.EQ.RSRC3 .AND. MYCOL.EQ.CSRC3 ) THEN
            DWORK(11) = S((JLOC3-1)*LLDS+ILOC3)
            DWORK(12) = T((JLOC3-1)*LLDT+ILOC3)
            IF( RSRC1.NE.RSRC3 .OR. CSRC1.NE.CSRC3 ) 
     $         CALL DGESD2D( ICTXT, 2, 1, DWORK(11), 2, RSRC1, CSRC1)
         END IF
         IF( MYROW.EQ.RSRC4 .AND. MYCOL.EQ.CSRC4 ) THEN

            DWORK(13) = S((JLOC4-1)*LLDS+ILOC4) 
            DWORK(14) = T((JLOC4-1)*LLDT+ILOC4)

            IF( K+1.LT.N ) THEN
               DWORK(15) = S((JLOC4-1)*LLDS+ILOC4+1)
            ELSE
               DWORK(15) = ZERO
            END IF
            IF( RSRC1.NE.RSRC4 .OR. CSRC1.NE.CSRC4 ) CALL DGESD2D( ICTXT
     $         , 3, 1, DWORK(13), 3, RSRC1, CSRC1 )
         END IF
         IF( MYROW.EQ.RSRC1 .AND. MYCOL.EQ.CSRC1 ) THEN
            DWORK(1) = S((JLOC1-1)*LLDS+ILOC1)
            DWORK(5) = T((JLOC1-1)*LLDT+ILOC1)
            IF( RSRC1.NE.RSRC2 .OR. CSRC1.NE.CSRC2 ) CALL DGERV2D( ICTXT
     $         , 2, 1, DWORK(9), 2, RSRC2, CSRC2 )
            DWORK(2) = DWORK(9)
            DWORK(6) = DWORK(10)
            IF( RSRC1.NE.RSRC3 .OR. CSRC1.NE.CSRC3 ) CALL DGERV2D( ICTXT
     $         , 2, 1, DWORK(11), 1, RSRC3, CSRC3)
            DWORK(3) = DWORK(11)
            DWORK(7) = DWORK(12)
            IF( RSRC1.NE.RSRC4 .OR. CSRC1.NE.CSRC4 ) CALL DGERV2D( ICTXT
     $         , 3, 1, DWORK(13), 3, RSRC4, CSRC4)
            DWORK(4) = DWORK(13)
            DWORK(8) = DWORK(14)
            DWORK(9) = DWORK(15)
            IF( DWORK(9).EQ.ZERO ) THEN
               IF( DWORK(2).NE.ZERO ) THEN
                  CALL DLAG2( DWORK, 2, DWORK( 5 ), 2, SMLNUM*EPS, BETA(
     $               K ), BETA( K+1 ), ALPHAR( K ), ALPHAR( K+1 ),
     $               ALPHAI( K ) )
                  ALPHAI( K+1 ) = -ALPHAI( K )
               ELSE
                  IF( ALPHAR( K ).EQ.ZERO .AND. 
     $               ALPHAI( K ).EQ.ZERO ) THEN
                     ALPHAR( K ) = DWORK(1)
                     ALPHAI( K ) = ZERO
                     BETA( K ) = DWORK(5) 
                  END IF
                  IF( ALPHAR( K+1 ).EQ.ZERO .AND. 
     $               ALPHAI( K+1 ).EQ.ZERO ) THEN 
                     ALPHAR( K+1 ) = DWORK(4)
                     ALPHAI( K+1 ) = ZERO
                     BETA( K+1 ) = DWORK(8)
                  END IF
               END IF
            ELSE
               IF( ALPHAR( K ).EQ.ZERO .AND. ALPHAI( K ).EQ.ZERO )THEN
                  ALPHAR( K ) = DWORK(1)
                  ALPHAI( K ) = ZERO
                  BETA( K ) = DWORK(5)                
               END IF
            END IF
         END IF
 570  CONTINUE
C     
C     Global reduction
C     
      IF( NPROCS.GT.1 ) THEN
         CALL DGSUM2D( ICTXT, 'All', ' ', N, 1, ALPHAR, N, -1, -1 )
         CALL DGSUM2D( ICTXT, 'All', ' ', N, 1, ALPHAI, N, -1, -1 )
         CALL DGSUM2D( ICTXT, 'All', ' ', N, 1, BETA, N, -1, -1 )
      END IF
C     
C     If BETA is negative for a real eigenvalue, make it positive
C     
      DO 580 K = 1, N
         CALL PDELGET( 'All', ' ', BETA2, T, K, K, DESCT ) 
         IF( ALPHAI( K ).EQ.ZERO .AND. BETA2.LT.ZERO ) THEN
            BETA( K ) = -BETA( K )
            ALPHAR( K ) = -ALPHAR( K )
            CALL PDSCAL( N, -ONE, S, K, 1, DESCS, N )
            CALL PDSCAL( N, -ONE, T, K, 1, DESCT, N )
            CALL PDSCAL( N, -ONE, Q, 1, K, DESCQ, 1 )
         END IF
 580  CONTINUE
 600  CONTINUE
C     
C     Store storage requirements in workspaces
C     
      DWORK( 1 ) = DBLE(LWMIN)
      IWORK( 1 ) = LIWMIN
C     
C     Return to calling program
C     
C     
      RETURN
C     
C     End of PBDTGORD
C     
      END
C     
