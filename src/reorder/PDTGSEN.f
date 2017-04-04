      SUBROUTINE PDTGSEN( IJOB, WANTQ, WANTZ, SELECT, PARA, N, S, IS, 
     $     JS, DESCS, T, IT, JT, DESCT, Q, IQ, JQ, DESCQ, Z, IZ, JZ, 
     $     DESCZ, ALPHAR, ALPHAI, BETA, M, PL, PR, DIF, DWORK, LDWORK, 
     $     IWORK, LIWORK, INFO )
C   
C  -- ScaLAPACK-style routine --
C     Preliminary version.
C     Dept. Computing Science and HPC2N, Univ. of Umeå, Sweden
C     June 7, 2007.
C
      IMPLICIT NONE
C
C     .. Scalar Arguments ..
      LOGICAL            WANTQ, WANTZ
      INTEGER            IJOB, INFO, LIWORK, LDWORK, M, N,
     $                   IT, JT, IQ, JQ, IS, JS, IZ, JZ
      DOUBLE PRECISION   PL, PR
C     ..
C     .. Array Arguments ..
      LOGICAL            SELECT( * )
      INTEGER            PARA( 6 ), DESCS( * ), DESCT( * ), DESCQ( * ), 
     $                   DESCZ( * ), IWORK( * )
      DOUBLE PRECISION   Q( * ), Z( * ), S( * ), T( * ), BETA( * ),
     $                   ALPHAI( * ), DWORK( * ), ALPHAR( * ), DIF( 2 )
C     ..
C
C  Purpose
C  =======
C     
C  PBDTGSEN reorders the generalized real Schur decomposition of a real
C  matrix pair (S, T) (in terms of an orthonormal equivalence trans-
C  formation Q' * (S, T) * Z), so that a selected cluster of eigenvalues
C  appears in the leading diagonal blocks of the upper quasi-triangular
C  matrix S and the upper triangular T. The leading columns of Q and
C  Z form orthonormal bases of the corresponding left and right eigen-
C  spaces (deflating subspaces). (S, T) must be in generalized real
C  Schur canonical form, i.e. S is block upper triangular with 1-by-1 and 
C  2-by-2 diagonal blocks. T is upper triangular.
C
C  PBDTGSEN also computes the generalized eigenvalues
C
C              w(j) = (ALPHAR(j) + i*ALPHAI(j))/BETA(j)
C
C  of the reordered matrix pair (S, T).
C
C  Optionally, the routine computes the estimates of reciprocal condition
C  numbers for eigenvalues and eigenspaces. These are Difu[(S11,T11),
C  (S22,T22)] and Difl[(S11,T11), (S22,T22)], i.e. the separation(s)
C  between the matrix pairs (S11, T11) and (S22,T22) that correspond to
C  the selected cluster and the eigenvalues outside the cluster, resp.,
C  and norms of "projections" onto left and right eigenspaces w.r.t.
C  the selected cluster in the (1,1)-block.
C     
C  This subroutine uses a delay and accumulate procedure for performing
C  the off-diagonal updates (see references for details).
C
C  Notes
C  =====
C
C  Each global data object is described by an associated description
C  vector.  This vector stores the information required to establish
C  the mapping between an object element and its corresponding process
C  and memory location.
C
C  Let A be a generic term for any 2D block cyclicly distributed array.
C  Such a global array has an associated description vector DESCA.
C  In the following comments, the character _ should be read as
C  "of the global array".
C
C  NOTATION        STORED IN      EXPLANATION
C  --------------- -------------- --------------------------------------
C  DTYPE_A(global) DESCA( DTYPE_ )The descriptor type.  In this case,
C                                 DTYPE_A = 1.
C  CTXT_A (global) DESCA( CTXT_ ) The BLACS context handle, indicating
C                                 the BLACS process grid A is distribu-
C                                 ted over. The context itself is glo-
C                                 bal, but the handle (the integer
C                                 value) may vary.
C  M_A    (global) DESCA( M_ )    The number of rows in the global
C                                 array A.
C  N_A    (global) DESCA( N_ )    The number of columns in the global
C                                 array A.
C  MB_A   (global) DESCA( MB_ )   The blocking factor used to distribute
C                                 the rows of the array.
C  NB_A   (global) DESCA( NB_ )   The blocking factor used to distribute
C                                 the columns of the array.
C  RSRC_A (global) DESCA( RSRC_ ) The process row over which the first
C                                 row of the array A is distributed.
C  CSRC_A (global) DESCA( CSRC_ ) The process column over which the
C                                 first column of the array A is
C                                 distributed.
C  LLD_A  (local)  DESCA( LLD_ )  The leading dimension of the local
C                                 array.  LLD_A >= MAX(1,LOCr(M_A)).
C
C  Let K be the number of rows or columns of a distributed matrix,
C  and assume that its process grid has dimension p x q.
C  LOCr( K ) denotes the number of elements of K that a process
C  would receive if K were distributed over the p processes of its
C  process column.
C  Similarly, LOCc( K ) denotes the number of elements of K that a
C  process would receive if K were distributed over the q processes of
C  its process row.
C  The values of LOCr() and LOCc() may be determined via a call to the
C  ScaLAPACK tool function, NUMROC:
C          LOCr( M ) = NUMROC( M, MB_A, MYROW, RSRC_A, NPROW ),
C          LOCc( N ) = NUMROC( N, NB_A, MYCOL, CSRC_A, NPCOL ).
C  An upper bound for these quantities may be computed by:
C          LOCr( M ) <= ceil( ceil(M/MB_A)/NPROW )*MB_A
C          LOCc( N ) <= ceil( ceil(N/NB_A)/NPCOL )*NB_A
C     
C  Arguments
C  =========
C     
C  IJOB    (global input) INTEGER
C          Specifies whether condition numbers are required for the
C          cluster of eigenvalues (PL and PR) or the deflating subspaces
C          (Difu and Difl):
C           =0: Only reorder w.r.t. SELECT. No extras.
C           =1: Reciprocal of norms of "projections" onto left and right
C               eigenspaces w.r.t. the selected cluster (PL and PR).
C           =2: Upper bounds on Difu and Difl. F-norm-based estimate
C               (DIF(1:2)).
C           =3: Estimate of Difu and Difl. 1-norm-based estimate
C               (DIF(1:2)).
C               About 5 times as expensive as IJOB = 2.
C           =4: Compute PL, PR and DIF (i.e. 0, 1 and 2 above): Economic
C               version to get it all.
C           =5: Compute PL, PR and DIF (i.e. 0, 1 and 3 above)
C     
C  WANTQ   (global input) LOGICAL
C          .TRUE. : update the left transformation matrix Q;
C          .FALSE.: do not update Q.
C
C  WANTZ   (global input) LOGICAL
C          .TRUE. : update the right transformation matrix Z;
C          .FALSE.: do not update Z.
C     
C  SELECT  (global input/output) INTEGER/LOGICAL  array, dimension (N)
C          SELECT specifies the eigenvalues in the selected cluster. To
C          select a real eigenvalue w(j), SELECT(j) must be set to
C          .TRUE. (1). To select a complex conjugate pair of eigenvalues
C          w(j) and w(j+1), corresponding to a 2-by-2 diagonal block,
C          either SELECT(j) or SELECT(j+1) or both must be set to
C          .TRUE. (1); a complex conjugate pair of eigenvalues must be
C          either both included in the cluster or both excluded.
C          On output, the (partial) reordering is displayed.
C
C  PARA    (global input) INTEGER*6
C          Block parameters (some should be replaced by calls to
C          PILAENV and others by meaningful default values):
C          PARA(1) = maximum number of concurrent computational windows
C                    allowed in the algorithm; 
C                    0 < PARA(1) <= min(NPROW,NPCOL) must hold;
C          PARA(2) = number of eigenvalues in each window; 
C                    0 < PARA(2) < PARA(3) must hold; 
C          PARA(3) = window size; PARA(2) < PARA(3) < DESCT(MB_) 
C                    must hold;
C          PARA(4) = minimal percentage of flops required for
C                    performing matrix-matrix multiplications instead
C                    of pipelined orthogonal transformations;
C                    0 <= PARA(4) <= 100 must hold;
C          PARA(5) = width of block column slabs for row-wise
C                    application of pipelined orthogonal
C                    transformations in their factorized form; 
C                    0 < PARA(5) <= DESCT(MB_) must hold.
C          PARA(6) = the maximum number of eigenvalues moved together
C                    over a process border; in practice, this will be
C                    approximately half of the cross border window size
C                    0 < PARA(6) <= PARA(2) must hold; 
C
C  N       (global input) INTEGER
C          The order of the globally distributed matrices S and T. N >= 0.
C
C  S       (local input/output) DOUBLE PRECISION array, 
C          dimension (LLD_S,LOCc(N)).
C          On entry, the local pieces of the global distributed 
C          upper quasi-triangular matrix S, in Schur form. On exit, S is 
C          overwritten by the local pieces of the reordered matrix S, 
C          again in Schur form, with the selected eigenvalues in the 
C          globally leading diagonal blocks.
C
C  IS      (global input) INTEGER
C  JS      (global input) INTEGER
C          The row and column index in the global array T indicating the
C          first column of sub( S ). IS = JS = 1 must hold. 
C     
C  DESCS   (global and local input) INTEGER array of dimension DLEN_.
C          The array descriptor for the global distributed matrix T.
C     
C    
C  T       (local input/output) DOUBLE PRECISION array, 
C          dimension (LLD_T,LOCc(N)).
C          On entry, the local pieces of the global distributed 
C          upper triangular matrix T, in Schur form. On exit, T is 
C          overwritten by the local pieces of the reordered matrix T, 
C          again in Schur form, with the selected eigenvalues in the 
C          globally leading diagonal blocks.
C
C  IT      (global input) INTEGER
C  JT      (global input) INTEGER
C          The row and column index in the global array T indicating the
C          first column of sub( T ). IT = JT = 1 must hold.
C     
C  DESCT   (global and local input) INTEGER array of dimension DLEN_.
C          The array descriptor for the global distributed matrix T.
C     
C  Q       (local input/output) DOUBLE PRECISION array, 
C          dimension (LLD_Q,LOCc(N)).
C          On entry, if WANTQ = .TRUE., the local pieces of the global
C          distributed matrix Q of generalized Schur vectors.
C          On exit, WANTQ = .TRUE., Q has been postmultiplied by the
C          global orthogonal transformation matrix which reorders (S, T); 
C          the leading M columns of Q and Z form orthonormal bases for the
C          specified deflating subspaces.
C          If WANTQ = .FALSE., Q is not referenced.
C
C  IQ      (global input) INTEGER
C  JQ      (global input) INTEGER
C          The column index in the global array Q indicating the
C          first column of sub( Q ). IQ = JQ = 1 must hold.
C     
C  DESCQ   (global and local input) INTEGER array of dimension DLEN_.
C          The array descriptor for the global distributed matrix Q.
C     
C  Z       (local input/output) DOUBLE PRECISION array, 
C          dimension (LLD_Z,LOCc(N)).
C          On entry, if WANTZ = .TRUE., the local pieces of the global
C          distributed matrix Z of generalized Schur vectors.
C          On exit, WANTZ = .TRUE., Z has been postmultiplied by the
C          global orthogonal transformation matrix which reorders (S, T); 
C          the leading M columns of Q and Z form orthonormal bases for the
C          specified deflating subspaces.
C          If WANTQ = .FALSE., Z is not referenced.
C
C  IZ      (global input) INTEGER
C  JZ      (global input) INTEGER
C          The column index in the global array Z indicating the
C          first column of sub( Z ). IZ = JZ = 1 must hold.
C     
C  DESCZ   (global and local input) INTEGER array of dimension DLEN_.
C          The array descriptor for the global distributed matrix Z.
C     
C  ALPHAR  (global output) DOUBLE PRECISION array, dimension (N)
C  ALPHAI  (global output) DOUBLE PRECISION array, dimension (N)
C  BETA    (global output) DOUBLE PRECISION array, dimension (N)
C          On exit, (ALPHAR(j) + ALPHAI(j)*i)/BETA(j), j=1,...,N, will
C          be the generalized eigenvalues.  ALPHAR(j) + ALPHAI(j)*i
C          and BETA(j),j=1,...,N  are the diagonals of the complex Schur
C          form (S,T) that would result if the 2-by-2 diagonal blocks of
C          the real generalized Schur form of (A,B) were further reduced
C          to triangular form using complex unitary transformations.
C          If ALPHAI(j) is zero, then the j-th eigenvalue is real; if
C          positive, then the j-th and (j+1)-st eigenvalues are a
C          complex conjugate pair, with ALPHAI(j+1) negative.
C          Note also that if a complex eigenvalue is sufficiently 
C          ill-conditioned, then its value may differ significantly 
C          from its value before reordering. 
C     
C  M       (global output) INTEGER
C          The dimension of the specified pair of left and right eigen-
C          spaces (deflating subspaces). 0 <= M <= N.
C  
C  PL, PR  (global output) DOUBLE PRECISION
C          If IJOB = 1, 4 or 5, PL, PR are lower bounds on the
C          reciprocal of the norm of "projections" onto left and right
C          eigenspaces with respect to the selected cluster.
C          0 < PL, PR <= 1.
C          If M = 0 or M = N, PL = PR  = 1.
C          If IJOB = 0, 2 or 3, PL and PR are not referenced.
C
C  DIF     (global output) DOUBLE PRECISION array, dimension (2).
C          If IJOB >= 2, DIF(1:2) store the estimates of Difu and Difl.
C          If IJOB = 2 or 4, DIF(1:2) are F-norm-based upper bounds on
C          Difu and Difl. If IJOB = 3 or 5, DIF(1:2) are 1-norm-based
C          estimates of Difu and Difl.
C          If M = 0 or N, DIF(1:2) = F-norm([S, T]).
C          If IJOB = 0 or 1, DIF is not referenced.
C     
C  DWORK   (local workspace/output) DOUBLE PRECISION array, 
C          dimension (LDWORK) 
C          On exit, if INFO = 0, DWORK(1) returns the optimal LDWORK.
C     
C  LDWORK  (local input) INTEGER
C          The dimension of the array DWORK.
C     
C          If LDWORK = -1, then a workspace query is assumed; the routine
C          only calculates the optimal size of the DWORK array, returns
C          this value as the first entry of the DWORK array, and no error
C          message related to LDWORK is issued by PXERBLA.
C     
C  IWORK   (local workspace/output) INTEGER array, dimension (LIWORK)
C     
C  LIWORK  (local input) INTEGER
C          The dimension of the array IWORK.
C     
C          If LIWORK = -1, then a workspace query is assumed; the
C          routine only calculates the optimal size of the IWORK array,
C          returns this value as the first entry of the IWORK array, and
C          no error message related to LIWORK is issued by PXERBLA.
C     
C  INFO    (global output) INTEGER
C          = 0: successful exit
C          < 0: if INFO = -i, the i-th argument had an illegal value.
C          If the i-th argument is an array and the j-entry had
C          an illegal value, then INFO = -(i*1000+j), if the i-th
C          argument is a scalar and had an illegal value, then INFO = -i.
C          > 0: here we have several possibilites
C            *) Reordering of (S,T) failed because some eigenvalues are too
C               close to separate (the problem is very ill-conditioned);
C               (S,T) may have been partially reordered, and the output
C               eigenvalues are in the same order as in (S,T); PL, PR and
C               DIF (if requested) are set to zero. The process that 
C               failed in the reordering will return INFO = {the index 
C               of T where the swap failed}; all others will return 
C               INFO = 1.
C            *) A 2-by-2 block to be reordered split into two 1-by-1
C               blocks and the second block failed to swap with an
C               adjacent block. The process that failed in the 
C               reordering will return INFO = {the index of T where the 
C               swap failed}; all others will return INFO = 1.
C            *) If INFO = 2, the routines used in the calculation of the
C               condition numbers raised a positive warning flag (see
C               the documentation for PGEGCSYD and PGCSYCON of the
C               SCASY library).
C            *) If INFO = 3, there is no valid BLACS context (see the
C               BLACS documentation for details).
C            *) If INFO = 333, PGEGCSYD raised an input error flag;
C               please report this bug to the authors (see below).
C               If INFO = 444, PGCSYCON raised an input error flag;
C               please report this bug to the authors (see below).
C          In a future release this subroutine may distinguish between
C          the case 1 and 2 above.  
C     
C  Method
C  ======
C
C  This routine performs parallel eigenvalue reordering in gen. Schur
C  form by parallelizing the approach proposed in [3]. The condition
C  number estimation part is performed by using techniques and code
C  from SCASY, see http://www.cs.umu.se/research/parallel/scasy.
C
C  Additional requirements
C  =======================
C
C  The following alignment requirements must hold:
C  (a) DESC[S,T]( MB_ ) = DESC[S,T]( NB_ ) = DESC[Q,Z]( MB_ ) = 
C      DESC[Q,Z]( NB_ )
C  (b) DESC[S,T]( RSRC_ ) = DESC[Q,Z]( RSRC_ )
C  (c) DESC[S,T]( CSRC_ ) = DESC[Q,Z]( CSRC_ ) 
C
C  All matrices must be blocked by a block factor larger than or
C  equal to two (3). This to simplify reordering across processor
C  borders in the presence of 2-by-2 blocks.
C
C  Limitations
C  ===========
C
C  This algorithm cannot work on submatrices of S, T, Q and Z i.e.,
C  IS = JS = IT = JT = IQ = JQ = IZ = JZ = 1 must hold.
C
C  References
C  ==========
C
C  [1] B. Kagstrom; A Direct Method for Reordering Eigenvalues in the
C      Generalized Real Schur Form of a Regular Matrix Pair (A, B), in
C      M.S. Moonen et al (eds), Linear Algebra for Large Scale and
C      Real-Time Applications, Kluwer Academic Publ. 1993, pp 195-218.
C
C  [2] B. Kagstrom and P. Poromaa; Computing Eigenspaces with Specified
C      Eigenvalues of a Regular Matrix Pair (A, B) and Condition
C      Estimation: Theory, Algorithms and Software, Numerical Algorithms
C      12, 3-4, 369-407, 1996. Also as LAPACK Working Note 87.
C
C  [3] B. Kagstrom and P. Poromaa, LAPACK-Style Algorithms and Software
C      for Solving the Generalized Sylvester Equation and Estimating the
C      Separation between Regular Matrix Pairs, ACM TOMS 22(1):78-103,
C      1996. Also as LAPACK Working Note 75. 
C
C  [4] D. Kressner; Block algorithms for reordering standard and
C      generalized Schur forms, ACM TOMS, 32(4):521-532, 2006. 
C      Also as LAPACK Working Note 171.
C
C  [5] R. Granat, B. Kågström, D. Kressner; Parallel eigenvalue 
C      reordering in real Schur form, submitted to Concurrency and
C      Computations: Practice & Experience, Septemeber, 2007. Also as 
C      LAPACK Working Note ???.
C
C  Parallel execution recommendations
C  ==================================
C 
C  Use a square grid, if possible, for maximum performance. The block
C  parameters in PARA should be kept well below the data distribution
C  block size. In particular, see [4,5] for recommended settings for
C  these parameters.
C  
C  In general, the parallel algorithm strives to perform as much work 
C  as possible without crossing the block borders on the main block 
C  diagonal.
C
C  Contributors
C  ============
C  
C  Implemented by Robert Granat, Umea University and HPC2N, June 2007,
C  in collaboration with Bo Kågström and Daniel Kressner.
C 
C  Revisions
C  =========
C
C  Please send bug-reports to granat@cs.umu.se
C
C  Keywords
C  ========
C
C  Generalized Schur form, eigenvalue reordering, Generalized coupled 
C  Sylvester matrix equation
C
C  =====================================================================
C     ..
C     .. Parameters ..
      CHARACTER          TOP
      INTEGER            BLOCK_CYCLIC_2D, CSRC_, CTXT_, DLEN_, DTYPE_,
     $                   LLD_, MB_, M_, NB_, N_, RSRC_
      LOGICAL            DBG, DBG2, DBG3
      DOUBLE PRECISION   ZERO, ONE
      PARAMETER          ( TOP = ' ',
     $                     BLOCK_CYCLIC_2D = 1, DLEN_ = 9, DTYPE_ = 1,
     $                     CTXT_ = 2, M_ = 3, N_ = 4, MB_ = 5, NB_ = 6,
     $                     RSRC_ = 7, CSRC_ = 8, LLD_ = 9,
     $                     ZERO = 0.0D+0, ONE = 1.0D+0, DBG = .FALSE.,
     $                     DBG2 = .FALSE., DBG3 = .FALSE. ) 
C     ..
C     .. Local Scalars ..
      LOGICAL            LQUERY, WANTP, WANTD1, WANTD2,
     $                   WANTD
     $                   
      INTEGER            NPROW, NPCOL, MYROW, MYCOL, NB, NPROCS, 
     $                   IERR, LLDT, TRSRC, TCSRC, JLOC1,
     $                   ICOFF12, ST12RWS, ST12CLS, SPACE, 
     $                   NOEXSY, ICTXT, I, K, KKS, 
     $                   LIWMIN, LWMIN, N1, N2,
     $                   LLDQ, RSRC, CSRC, IRSRC,
     $                   IPW1, ILOC1, LIHI, IPW4, STROWS,
     $                   STCOLS, ITT, JTT, WRK1, IWRK1, WRK2, IWRK2, 
     $                   WRK3, IWRK3, LLDS, LLDZ 
      DOUBLE PRECISION   EST, SCALE, TIME1, TIME2, TIME, TIME3, TIME4, 
     $                   DSUM, RDSCAL
C     ..
C     .. Local Arrays ..
      INTEGER            DESC12( DLEN_ ), MBNB2( 2 )
      INTEGER            IDUM1( 1 ), IDUM2( 1 )
C     ..
C     .. External Functions ..
      LOGICAL            LSAME, INT2LG
      INTEGER            NUMROC, INDXG2P, INDXG2L, ILG2NT
      DOUBLE PRECISION   PDLANGE, MPI_WTIME
      EXTERNAL           LSAME, PDLANGE, NUMROC, INDXG2P, INDXG2L,
     $                   MPI_WTIME, INT2LG, ILG2NT
C     ..
C     .. External Subroutines ..
      EXTERNAL           PDLACPY, PXERBLA,
     $                   DGEMM, DLACPY, ILACPY, CHK1MAT, PCHK1MAT, 
     $                   INFOG2L, DGSUM2D, DGESD2D, DGERV2D, DGEBS2D, 
     $                   DGEBR2D, IGSUM2D, BLACS_GRIDINFO, IGEBS2D, 
     $                   IGEBR2D, IGAMX2D, IGAMN2D, PCHK2MAT
C      EXTERNAL           PSYCTCON, PGESYCTD
C     ..
C     .. Intrinsic Functions ..
      INTRINSIC          ABS, MAX, SQRT, MIN
C     ..
C     .. Local Functions .. 
C     .. 
C     .. Executable Statements ..
C
C     Get grid parameters 
C     
      ICTXT = DESCT( CTXT_ )
      CALL BLACS_GRIDINFO( ICTXT, NPROW, NPCOL, MYROW, MYCOL )
      NPROCS = NPROW*NPCOL
      IF(DBG) WRITE(*,*)MYROW,MYCOL, 'Entering PBDTGSEN'
      IF(DBG) WRITE(*,*)MYROW,MYCOL, 'NPROW,NPCOL,NPROCS=',NPROW,NPCOL,
     $     NPROCS
C     
C     Test if grid is O.K., i.e., the context is valid
C     
      INFO = 0
      IF( NPROW.EQ.-1 ) THEN
         INFO = 3
      END IF
C
C     Check if workspace query
C
      LQUERY = LDWORK.EQ.-1 .OR. LIWORK.EQ.-1 
C     
C     Test dimensions for local sanity
C     
      IF( INFO.EQ.0 ) THEN
         CALL CHK1MAT( N, 6, N, 6, IS, JS, DESCS, 10, INFO )
         CALL CHK1MAT( N, 6, N, 6, IT, JT, DESCT, 14, INFO )
      END IF
      IF( INFO.EQ.0 ) THEN
         CALL CHK1MAT( N, 6, N, 6, IQ, JQ, DESCQ, 18, INFO )
         CALL CHK1MAT( N, 6, N, 6, IZ, JZ, DESCZ, 22, INFO )
      END IF
C     
C     Check the blocking sizes for alignment requirements
C     
      IF( INFO.EQ.0 ) THEN
         IF( DESCS( MB_ ).NE.DESCS( NB_ ) ) INFO = -(1000*10 + MB_)
         IF( DESCT( MB_ ).NE.DESCT( NB_ ) ) INFO = -(1000*14 + MB_) 
      END IF
      IF( INFO.EQ.0 ) THEN
         IF( DESCQ( MB_ ).NE.DESCQ( NB_ ) ) INFO = -(1000*18 + MB_)
         IF( DESCZ( MB_ ).NE.DESCZ( NB_ ) ) INFO = -(1000*22 + MB_)
      END IF
      IF( INFO.EQ.0 ) THEN
         IF( DESCS( MB_ ).NE.DESCT( MB_ ) ) INFO = -(1000*14 + MB_)
         IF( DESCS( MB_ ).NE.DESCQ( MB_ ) ) INFO = -(1000*18 + MB_)
         IF( DESCS( MB_ ).NE.DESCZ( MB_ ) ) INFO = -(1000*22 + MB_)
      END IF
C     
C     Check the blocking sizes for minimum sizes
C
      IF( INFO.EQ.0 ) THEN
         IF( N.NE.DESCS( MB_ ) .AND. DESCS( MB_ ).LT.3 ) 
     $        INFO = -(1000*10 + MB_)
      END IF
C     
C     Check parameters in PARA  
C
      NB = DESCS( MB_ )
      IF( INFO.EQ.0 ) THEN
         IF( PARA(1).LT.1 .OR. PARA(1).GT.MIN(NPROW,NPCOL) ) 
     $        INFO = -(1000 * 5 + 1)
         IF( PARA(2).LT.1 .OR. PARA(2).GE.PARA(3) ) 
     $        INFO = -(1000 * 5 + 2)
         IF( PARA(3).LT.1 .OR. PARA(3).GT.NB ) 
     $        INFO = -(1000 * 5 + 3)
         IF( PARA(4).LT.0 .OR. PARA(4).GT.100 ) 
     $        INFO = -(1000 * 5 + 4)
         IF( PARA(5).LT.1 .OR. PARA(5).GT.NB )
     $        INFO = -(1000 * 5 + 5)  
         IF( PARA(6).LT.1 .OR. PARA(6).GT.PARA(2) )
     $        INFO = -(1000 * 5 + 6)  
      END IF
       
      
C     
C     Check requirements on IT, JT, IQ and JQ
C     
      IF( INFO.EQ.0 ) THEN
         IF( IS.NE.1 ) INFO = -8
         IF( JS.NE.IS ) INFO = -9
         IF( IT.NE.1 ) INFO = -12
         IF( JT.NE.IT ) INFO = -13
         IF( IQ.NE.1 ) INFO = -16
         IF( JQ.NE.IQ ) INFO = -17
         IF( IZ.NE.1 ) INFO = -20
         IF( JZ.NE.IZ ) INFO = -21
      END IF
C     
C     Test input parameters for global sanity
C     
      IF( INFO.EQ.0 ) THEN
         CALL PCHK1MAT( N, 6, N, 6, IS, JS, DESCS, 10, 0, IDUM1, 
     $        IDUM2, INFO )
         CALL PCHK1MAT( N, 6, N, 6, IT, JT, DESCT, 14, 0, IDUM1, 
     $        IDUM2, INFO )
      END IF
      IF( INFO.EQ.0 ) THEN
         CALL PCHK1MAT( N, 6, N, 6, IQ, JQ, DESCQ, 18, 0, IDUM1, 
     $        IDUM2, INFO )
         CALL PCHK1MAT( N, 6, N, 6, IZ, JZ, DESCZ, 22, 0, IDUM1, 
     $        IDUM2, INFO ) 
      END IF
      IF( INFO.EQ.0 ) THEN
         CALL PCHK2MAT( N, 6, N, 6, IS, JS, DESCS, 10, N, 6, N, 6, 
     $        IT, JT, DESCT, 14, 0, IDUM1, IDUM2, INFO )
         CALL PCHK2MAT( N, 6, N, 6, IS, JS, DESCS, 10, N, 6, N, 6, 
     $        IQ, JQ, DESCQ, 18, 0, IDUM1, IDUM2, INFO )
         CALL PCHK2MAT( N, 6, N, 6, IT, JT, DESCT, 14, N, 6, N, 6, 
     $        IZ, JZ, DESCZ, 22, 0, IDUM1, IDUM2, INFO )
      END IF
C     
C     Decode and test the input parameters
C     
      WANTP = IJOB.EQ.1 .OR. IJOB.GE.4
      WANTD1 = IJOB.EQ.2 .OR. IJOB.EQ.4
      WANTD2 = IJOB.EQ.3 .OR. IJOB.EQ.5
      WANTD = WANTD1 .OR. WANTD2
C
      IF( INFO.EQ.0 .OR. LQUERY ) THEN
C     
         IF( IJOB.LT.0 .OR. IJOB.GT.5 ) THEN
            INFO = -1
         ELSEIF( N.LT.0 ) THEN
            INFO = -6
         ELSE
C     
C     Extract local leading dimension
C     
            LLDS = DESCS( LLD_ )
            LLDT = DESCT( LLD_ )
            LLDQ = DESCQ( LLD_ )
            LLDZ = DESCZ( LLD_ )
C     
C     Check the SELECT vector for consistency and set M to the dimension 
C     of the specified invariant subspace.
C     
            IF(DBG) WRITE(*,*)MYROW,MYCOL, 'Calculating M'
            M = 0
            DO 10 K = 1, N
*
*              IWORK(1:N) is an integer copy of SELECT.
*
               IF( SELECT(K) ) THEN
                  IWORK(K) = 1
               ELSE
                  IWORK(K) = 0
               END IF

               IF( K.LT.N ) THEN
                  CALL INFOG2L( K+1, K, DESCS, NPROW, NPCOL, 
     $                 MYROW, MYCOL, ITT, JTT, TRSRC, TCSRC )
                  IF( MYROW.EQ.TRSRC .AND. MYCOL.EQ.TCSRC ) THEN
                     IF( S( (JTT-1)*LLDS + ITT ).NE.ZERO ) THEN
                        IF( SELECT(K).AND..NOT.
     $                       SELECT(K+1) ) THEN
                           IF(DBG) WRITE(*,*) 'SELECT ERROR 1 K=',K
C                           INFO = -3
                           SELECT(K+1) = .TRUE.
                           IWORK(K+1) = 1
                        ELSEIF( .NOT.SELECT(K).AND.
     $                           SELECT(K+1) ) THEN
                           IF(DBG) WRITE(*,*) 'SELECT ERROR 2 K=',K
C                           INFO = -3
                           SELECT(K) = .TRUE.
                           IWORK(K) = 1
                        END IF
                     END IF
                  END IF
               END IF
               IF( SELECT(K) ) M = M + 1
               IF(DBG) WRITE(*,*)MYROW,MYCOL,'M=',M 
 10         CONTINUE
            IF(DBG) WRITE(*,*)MYROW,MYCOL, 'Done Calculating M=',M
            IF(DBG) WRITE(*,*) 'INFO=',INFO
C     
C     Set parameters for deep pipelining in parallel Sylvester solver
C
            MBNB2( 1 ) = MIN( MAX( PARA( 3 ), PARA( 2 )*2 ), NB )
            MBNB2( 2 ) = MBNB2( 1 )
C     
C     Compute needed workspace
C     
            IF(DBG) WRITE(*,*) MYROW,MYCOL,'computing workspace'
            N1 = M
            N2 = N - M
            IF( WANTP ) THEN
               IF(DBG) WRITE(*,*) MYROW,MYCOL,'pgegcsyd'
c               CALL PGEGCSYD( 'Solve', 'Notranspose','Schur', 'Schur', 
c     $              'Notranspose', 'Notranspose', -1, 'Demand', N1, N2, 
c     $              S, 1, 1, DESCS, S, N1+1, N1+1, DESCS, S, 1, N1+1,
c     $              DESCS, T, 1, 1, DESCT, T, N1+1, N1+1, DESCT, T, 1, 
c     $              N1+1, DESCT, MBNB2, DWORK, -1, IWORK, -1, NOEXSY, 
c     $              SCALE, IERR )
               WRK1 = INT(DWORK(1))
               IWRK1 = IWORK(1)
               WRK1 = 0
               IWRK1 = 0               
            ELSE
               WRK1 = 0
               IWRK1 = 0
            END IF
C
            IF( WANTD ) THEN
              IF(DBG) WRITE(*,*) MYROW,MYCOL,'pgcsycon' 
c               CALL PGCSYCON( 'Notranspose', 'Notranspose', -1, 
c     $              'Demand', N1, N2, S, 1, 1, DESCS, S, N1+1, N1+1,  
c     $              DESCS, T, 1, 1, DESCT, T, N1+1, N1+1, DESCT, MBNB2, 
c     $              DWORK, -1,  IWORK, -1, EST, ITER, IERR )
               WRK2 = INT(DWORK(1))
               IWRK2 = IWORK(1)
               WRK2 = 0
               IWRK2 = 0               
            ELSE
               WRK2 = 0
               IWRK2 = 0
            END IF
C
            STROWS = NUMROC( N, NB, MYROW, DESCS(RSRC_), NPROW )
            STCOLS = NUMROC( N, NB, MYCOL, DESCS(CSRC_), NPCOL )
            WRK3 = 2*N + 14*NB**2 + 4*STROWS*PARA( 3 ) + 
     $           2*STCOLS*PARA( 3 ) + MAX( 2*STROWS*PARA( 3 ), 
     $           2*STCOLS*PARA( 3 ) )
            IWRK3 = N + 5*PARA( 1 ) + PARA(2)*PARA(3) - 
     $           PARA(2) * (PARA(2) + 1 ) / 2
C
            IF( IJOB.EQ.1 .OR. IJOB.EQ.2 .OR. IJOB.EQ.4 ) THEN
               LWMIN = MAX( 1, MAX( WRK2, WRK3) )
               LIWMIN = MAX( 1, MAX( IWRK2, IWRK3 ) )
            ELSE IF( IJOB.EQ.3 .OR. IJOB.EQ.5 ) THEN
               LWMIN = MAX( 1, WRK3 )
               LIWMIN = IWRK3
            ELSE
               LWMIN = MAX( 1, MAX( WRK1, WRK3) )
               LIWMIN = MAX( 1, MAX( IWRK1, IWRK3 ) )
            END IF
C     
            IF( LDWORK.LT.LWMIN .AND. .NOT.LQUERY ) THEN
               INFO = -31
            ELSE IF( LIWORK.LT.LIWMIN .AND. .NOT.LQUERY ) THEN
               INFO = -33
            END IF
            IF(DBG) WRITE(*,*) MYROW,MYCOL,'done computing workspace'
         END IF
      END IF
C     
C     Global maximum on info
C     
      IF( NPROCS.GT.1 )
     $     CALL IGAMX2D( ICTXT, 'All', TOP, 1, 1, INFO, 1, -1, -1, -1, 
     $     -1, -1 )
C     
C     Return if some argument is incorrect
C     
      IF(DBG) WRITE(*,*) MYROW,MYCOL,'return if argument incorrect'
      IF( INFO.NE.0 .AND. .NOT. LQUERY ) THEN
         M = 0                 
         CALL PXERBLA( ICTXT, 'PDTGSEN', -INFO )
         RETURN
      ELSEIF( LQUERY ) THEN
         DWORK( 1 ) = DBLE(LWMIN)
         IWORK( 1 ) = LIWMIN
         RETURN
      END IF
C     
C     Quick return if possible.
C    
      IF(DBG) WRITE(*,*) MYROW,MYCOL,'quick return if possible' 
      IF( M.EQ.N .OR. M.EQ.0 ) THEN
         IF( WANTP ) THEN
            PL = ONE
            PR = ONE
         END IF
         IF( WANTD ) THEN
            SCALE = ZERO
            DSUM = ONE
            DO 20 I = 1, N
               CALL PDLASSQ( N, S, 1, I, DESCS, 1, SCALE, DSUM )
               CALL PDLASSQ( N, T, 1, I, DESCT, 1, SCALE, DSUM )
   20       CONTINUE
            DIF( 1 ) = SCALE*SQRT( DSUM )
            DIF( 2 ) = DIF( 1 )
         END IF
         GO TO 50
      END IF
      
      CALL PDTGORD( WANTQ, WANTZ, IWORK, PARA, N, 
     $     S, DESCS, T, DESCT, 
     $     Q, DESCQ, Z, DESCZ, 
     $     ALPHAR, ALPHAI, BETA, M, DWORK, LDWORK, 
     $     IWORK(N+1), LIWORK-N, INFO )
      
C     
      IF( WANTP ) THEN
         TIME = MPI_WTIME()
C     
C     Solve generalized Sylvester equation for R and L and compute PL 
C     and PR.
C     
C     Copy (S12,T12) to workspace
C     
         CALL INFOG2L( 1, N1+1, DESCS, NPROW, NPCOL, MYROW, 
     $        MYCOL, ILOC1, JLOC1, RSRC, CSRC )
         ICOFF12 = MOD( N1, NB )
         ST12RWS = NUMROC( N1, NB, MYROW, RSRC, NPROW )
         ST12CLS = NUMROC( N2+ICOFF12, NB, MYCOL, CSRC, NPCOL )
         SPACE = DESC12( LLD_ ) * ST12CLS
         CALL DESCINIT( DESC12, N1, N2+ICOFF12, NB, NB, RSRC, 
     $        CSRC, ICTXT, MAX(1,ST12RWS), IERR )
         CALL PDLACPY( 'All', N1, N2, S, 1, N1+1, DESCS, DWORK, 
     $        1, 1+ICOFF12, DESC12 )
         CALL PDLACPY( 'All', N1, N2, T, 1, N1+1, DESCT, 
     $        DWORK( 1 + SPACE ), 1, 1+ICOFF12, DESC12 )
C     
C     Solve the equation to get the solution in workspace 
C     
         IPW1 = 1 + 2*SPACE
         IERR = 0
c         CALL PGEGCSYD( 'Solve', 'Notranspose','Schur', 'Schur', 
c     $        'Notranspose', 'Notranspose', -1, 'Demand', N1, N2, 
c     $        S, 1, 1, DESCS, S, N1+1, N1+1, DESCS, DWORK, 1, 1+ICOFF12, 
c     $        DESC12, T, 1, 1, DESCT, T, N1+1, N1+1, DESCT, 
c     $        DWORK( 1 + SPACE ), 1, 1+ICOFF12, DESC12, MBNB2, 
c     $        DWORK(IPW1), LDWORK-2*SPACE+1, IWORK, LIWORK, NOEXSY, 
c     $        SCALE, IERR )
         IF( IERR.LT.0 ) THEN
            INFO = 333
         ELSE
            INFO = 2
         END IF 
C     
C     Estimate the reciprocal of norms of "projections" onto left
C     and right eigenspaces.
C     
         RDSCAL = ZERO
         DSUM = ONE
         CALL PDLASSQ( N1*N2, DWORK, 1, 1, DESC12, 1, RDSCAL, DSUM )
         PL = RDSCAL*SQRT( DSUM )
         IF( PL.EQ.ZERO ) THEN
            PL = ONE
         ELSE
            PL = SCALE / ( SQRT( SCALE*SCALE / PL+PL )*SQRT( PL ) )
         END IF
         RDSCAL = ZERO
         DSUM = ONE
         CALL PDLASSQ( N1*N2, DWORK( 1 + SPACE ), 1, 1, DESC12, 1, 
     $        RDSCAL, DSUM )
         PR = RDSCAL*SQRT( DSUM )
         IF( PR.EQ.ZERO ) THEN
            PR = ONE
         ELSE
            PR = SCALE / ( SQRT( SCALE*SCALE / PR+PR )*SQRT( PR ) )
         END IF
      END IF
C     
      IF( WANTD ) THEN
         TIME = MPI_WTIME()
C     
C     Compute estimates of Difu and Difl.
C     
         IF( WANTD1 ) THEN
            I = N1 + 1
C     
C     Frobenius norm-based Difu-estimate.
C     
            WRITE(*,*) 
     $           'No support in SCASY for Frobenius norm estimates yet!'
            DIF(1) = ZERO
C     
C     Frobenius norm-based Difl-estimate.
C     
C     
            WRITE(*,*) 
     $           'No support in SCASY for Frobenius norm estimates yet!'
            DIF(2) = ZERO
         ELSE          
C     
C     1-norm-based estimate of Difu.
C     
c            CALL PGCSYCON( 'Notranspose', 'Notranspose', -1, 
c     $           'Demand', N1, N2, S, 1, 1, DESCS, S, N1+1, N1+1,  
c     $           DESCS, T, 1, 1, DESCT, T, N1+1, N1+1, DESCT, MBNB2, 
c     $           DWORK, LDWORK, IWORK, LIWORK, DIF( 1 ), ITER, IERR )
            DIF( 1 ) = DIF( 1 ) * SQRT(DBLE(2*N1*N2))
            DIF( 1 ) = ONE / DIF( 1 )
C     
C     1-norm-based estimate of Difl.
C     
c            CALL PGCSYCON( 'Notranspose', 'Notranspose', -1, 
c     $           'Demand', N2, N1, S, N1+1, N1+1, DESCS, S, 1, 1,  
c     $           DESCS, T, N1+1, N1+1, DESCT, T, 1, 1, DESCT, MBNB2, 
c     $           DWORK, LDWORK, IWORK, LIWORK, DIF( 2 ), ITER, IERR )
            DIF( 2 ) = DIF( 2 ) * SQRT(DBLE(2*N1*N2))
            DIF( 2 ) = ONE / DIF( 2 )
C     
         END IF
         IF( IERR.LT.0 ) THEN
            INFO = 444
         ELSE
            INFO = 2
         END IF
         TIME4 = MPI_WTIME() - TIME
      END IF
C     

C     
C     Store storage requirements in workspaces
C     
      IF(DBG) WRITE(*,*) MYROW,MYCOL,'Storage requirements'
      DWORK( 1 ) = DBLE(LWMIN)
      IWORK( 1 ) = LIWMIN
C     
C     Return to calling program
C     
      DWORK(2) = TIME1
      DWORK(3) = TIME2
      DWORK(4) = TIME3
      DWORK(5) = TIME4
C     
 50   CONTINUE
      
      RETURN
C     
C     End of PBDTGSEN
C     
      END
C     
