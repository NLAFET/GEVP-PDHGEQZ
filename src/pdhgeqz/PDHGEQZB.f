***********************************************************************
*                                                                     *
*     PDHGEQZB.f:                                                     *
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
      SUBROUTINE PDHGEQZB( ILSCHR, ILQ, ILZ, UPDOFFD, N, NBULGES, IFIRST
     $   , ILAST, J, WSIZE, MOVEDIST, ILO, A, DESCA, B, DESCB, Q, DESCQ, 
     $   Z, DESCZ, ILOZ, ILOQ, IHIQ, IHIZ, ALPHAR, ALPHAI, BETA, ATOL,
     $   PU, PV, LDUV, TMPA, TMPB, LDTMPAB, UPDTIME, DWORK, LDWORK, INFO
     $   )
      IMPLICIT NONE
*     ..
*     .. Scalar Arguments ..
*     ..
      INTEGER             NBULGES, IFIRST, ILAST, J, WSIZE, ILO,
     $   MOVEDIST, N, ILOZ, ILOQ, IHIQ, IHIZ, INFO, LDUV, LDWORK,
     $   LDTMPAB
      LOGICAL             ILSCHR, ILQ, ILZ, UPDOFFD
      DOUBLE PRECISION    ATOL
*     ..
*     .. Array Arguments ..
*     ..
      INTEGER             DESCA( * ), DESCB( * ), DESCZ( * ), DESCQ( * )
      DOUBLE PRECISION    A( * ) , B( * ), Q( * ) , Z( * ), DWORK( * ),
     $   UPDTIME( 40 ), ALPHAI( * ), ALPHAR( * ), BETA( * ), PU( LDUV, *
     $   ), PV( LDUV, * ), TMPA( LDTMPAB, * ), TMPB( LDTMPAB, * )       
*     
*     Purpose
*     =======
*     
*     PDHGEQZB chases bulges, movedist steps down, within a 
*     diagonal block of A and B.
*     
*     All the inputs assumed to be valid without checking.
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
*     = TRUE: Apply all transformations to A and B. 
*     = FALSE: A and B are perturbed, but only so that
*     the eigenvalues are preserved.
*     
*     ILQ     (input) LOGICAL
*     = TRUE: Compute Q.
*     = FALSE: Q is not referenced.
*     
*     ILZ     (global input) LOGICAL
*     = TRUE: Compute Z.
*     = FALSE: Z is not referenced.
*     
*     UPDOFFD (global input) LOGICAL
*     = TRUE: Update off diagonal entries after the chase is complete.
*     = FALSE: Do not update off diagonal entries.
*     
*     N       (global input) INTEGER 
*     Order of the matrices A, B, Q, and Z.
*     
*     NBULGES (global input) INTEGER
*     Number of bulges to be chased.
*     
*     IFIRST  (global input) INTEGER
*     Minimum row/column index of where to chase bulges.
*     
*     ILAST   (global input) INTEGER
*     Maximum row/column index for where to chase bulges.
*     
*     J       (global input) INTEGER 
*     Row/column index of where to start the chase.
*     
*     WSIZE   (global input) INTEGER
*     Size of the diagonal block where bulges are to be chased.         
*     
*     MOVEDIST(global input) INTEGER
*     Number of step each bulge should be chased.
*     
*     ILO     (global input) INTEGER
*     Updates wont be applied below this row/column index.
*     
*     A       (local input/output) DOUBLE PRECISION array of size 
*     (DESCA(LLD_),*)
*     On input A contains a Hessenberg matrix.  On output A
*     is perturbed s.t. shifts 
*     (ALPHAR(1:NBULGES*2)+i*ALPHAI(1:NBULGES*2)/BETA(1:NBULGES*2) is  
*     used to create bulges  and these are and chased forward so 
*     that all bulges can fit. 
*     
*     DESCA   (global and local input) INTEGER array of dimension DLEN_.
*     The array descriptor for the distributed matrix A.
*     
*     B       (local input/output) DOUBLE PRECISION array of size  
*     (DESCB(LLD_),*)
*     On input B contains an upper triangular matrix.  On output B
*     is perturbed s.t. shifts 
*     (ALPHAR(1:NBULGES*2)+i*ALPHAI(1:NBULGES*2)/BETA(1:NBULGES*2) is 
*     used to create bulges and these are and chased forward so 
*     that all bulges can fit. 
*     
*     DESCB   (global and local input) INTEGER array of dimension DLEN_. 
*     The array descriptor for the distributed matrix B.
*     
*     Q       (local input/output) DOUBLE PRECISION array of size
*     (DESCQ(LLD_),*)
*     If ILQ = .TRUE., then the QZ Sweep left orthogonal
*     eqivalence transformation is accumulated into
*     Q(1:N,ILOQ:IHIQ) from the right.
*     If ILQ = .FALSE., then Q is unreferenced.
*     
*     DESCQ   (global and local input) INTEGER array of dimension DLEN_.
*     The array descriptor for the distributed matrix Q.
*     
*     Z       (local input/output) DOUBLE PRECISION array of size
*     (DESCZ(LLD_),*)
*     If ILZ = .TRUE., then the QZ Sweep right orthogonal
*     eqivalence transformation is accumulated into
*     Z(1:N,ILOZ:IHIZ) from the right.
*     If ILZ = .FALSE., then Z is unreferenced.
*     
*     ILOQ    (global input) INTEGER
*     IHIQ    (global input) INTEGER
*     Specify the rows of Q to which transformations must be
*     applied if ILQ is .TRUE.. 1 .LE. ILOQ .LE. IHIQ .LE. N
*     
*     ILOZ    (global input) INTEGER
*     IHIZ    (global input) INTEGER
*     Specify the rows of Z to which transformations must be
*     applied if ILZ is .TRUE.. 1 .LE. ILOZ .LE. IHIZ .LE. N
*     
*     ALPHAR  (global input) DOUBLE PRECISION array of size ( NBULGES * 2 )
*     ALPHAI  (global input) DOUBLE PRECISION array of size ( NBULGES * 2 )
*     BETA    (global input) DOUBLE PRECISION array of size ( NBULGES * 2 )
*     ALPHAR( 1 : NBULGES * 2 ) contains the real parts and 
*     ALPHAI( 1 : NBUGLES * 2 ) contains the imaginary
*     parts and BETA( 1 : NBULGES * 2 ) contains the scale parts of 
*     the shifts that early was used to created the bulges that are
*     about to chased - used for recreating bulges if vigialent 
*     deflation occurs.
*     
*     ATOL    (global input) DOUBLE PRECISION 
*     Holds number where safe to set a value of A to zero, used 
*     for vigialent deflations.
*     
*     PU      (local input) DOUBLE PRECISION, array, dimension (LDUV, *) 
*     Local workspace, minium size (WSIZE, WSIZE)
*     
*     PV      (local input) DOUBLE PRECISION, array, dimension (LDUV, *)
*     Local workspace, minium size (WSIZE, WSIZE)
*     
*     LDUV    (global input) INTEGER
*     Leading dimension for the local workpaces PU and PV.
*     
*     TMPA    (local input) DOUBLE PRECISION, array, dimension (LDTMPAB, *).
*     Local workspace, minium size (WSIZE, WSIZE)
*     
*     TMPB    (local input) DOUBLE PRECISION, array, dimension (LDTMPAB, *).
*     Local workspace, minium size (WSIZE, WSIZE)
*     
*     LDTMPAB (global input) INTEGER
*     Leading dimension for the local workspace TMPA and TMPB.
*     
*     UPDTIME (local input) INTEGER, array dimension 40.
*     Storage for misc. timers.
*     
*     DWORK   (local workspace/global input) DOUBLE PRECISION array,
*     dimension (LDWORK)
*     
*     LDWORK  (global input) INTEGER
*     The dimension of the array DWORK.
*     
*     INFO     (global output) INTEGER
*     = 0 : Succesful call
*     
*     ================================================================


*     ..
*     .. Parameters ..
*     ..
      DOUBLE PRECISION    ZERO, ONE
      PARAMETER          ( ZERO = 0.0D+0, ONE = 1.0D+0 )
      INTEGER             BLOCK_CYCLIC_2D, CSRC_, CTXT_, DLEN_, DT_, 
     $   LLD_, MB_, M_, NB_, N_, RSRC_
      PARAMETER          ( BLOCK_CYCLIC_2D = 1, DLEN_ = 9, DT_ = 1,
     $   CTXT_ = 2, M_ = 3, N_ = 4, MB_ = 5, NB_ = 6,
     $   RSRC_ = 7, CSRC_ = 8, LLD_ = 9 )
      INTEGER             NSBULGE, BULGESIZE
      PARAMETER           (NSBULGE = 2,BULGESIZE = NSBULGE + 1)

*     ..
*     .. Local Arrays ..
*     ..      


*     ..
*     .. Local Scalars ..
*     ..
      INTEGER             JJ, IB, NNB, NB, ICTXT, MYROW, MYCOL, NPROW,
     $   NPCOL, DIM1, DIM4, IAM, NPROCS, NCOLS, LTOP, LDA, LDB, RSRC1,
     $   CSRC1, RSRC2, CSRC2, RSRC3, CSRC3, RSRC4, CSRC4

      LOGICAL             OWNER, PARTOWNER

      DOUBLE PRECISION    C, S, TEMP, T1


*     ..
*     .. Intrinsic Functions ..
*     ..
      INTRINSIC           ABS, DBLE, MAX, MIN, SQRT, MOD
*     ..
*     .. Externals
*     ..
      EXTERNAL            INDXG2L, INDXG2P, MPI_WTIME    
      DOUBLE PRECISION    MPI_WTIME
      INTEGER             INDXG2L, INDXG2P


      IF ( IFIRST .GE. ILAST ) RETURN
      IF ( J .GE. ILAST ) RETURN
      IF ( IFIRST + 1 .EQ. ILAST ) RETURN
      IF ( J + 1 .EQ. ILAST ) RETURN
      
      

      IF( ILSCHR ) THEN
         LTOP = ILO
      ELSE
         LTOP = J
      END IF
      
*     Extract current communication context and get grid parameters
      ICTXT = DESCA( CTXT_ )
      CALL BLACS_PINFO( IAM, NPROCS )
      CALL BLACS_GRIDINFO( ICTXT, NPROW, NPCOL, MYROW, MYCOL )


*     NB_ .EQ. MB_ is assumed, we use NB as blocking factor.
      NB = DESCA( NB_ )
      
      LDA = DESCA( LLD_ )
      LDB = DESCB( LLD_ )
      
      T1 = MPI_WTIME()

*     Extract processors for current window,
*     as below:
*     
*     1 | 2
*     -----
*     3 | 4
      RSRC1 = INDXG2P( J, NB, -1, DESCA(RSRC_), NPROW )
      CSRC1 = INDXG2P( J, NB, -1, DESCA(CSRC_), NPCOL )
      RSRC2 = RSRC1
      CSRC2 = MOD( CSRC1 + 1, NPCOL )
      RSRC3 = MOD( RSRC1 + 1, NPROW )
      CSRC3 = CSRC1
      RSRC4 = RSRC3
      CSRC4 = CSRC2
      
      
      OWNER = ( MYROW .EQ. RSRC1 .AND. MYCOL .EQ. CSRC1 ) .OR.
     $   ( MYROW .EQ. RSRC4 .AND. MYCOL .EQ. CSRC4 )
      PARTOWNER = ( MYROW .EQ. RSRC1 .AND. MYCOL .EQ. CSRC1 ) .OR.
     $   ( MYROW .EQ. RSRC2 .AND. MYCOL .EQ. CSRC2 ) .OR.
     $   ( MYROW .EQ. RSRC3 .AND. MYCOL .EQ. CSRC3 ) .OR.
     $   ( MYROW .EQ. RSRC4 .AND. MYCOL .EQ. CSRC4 )
      
      DIM1 = NB - MOD( J - 1, NB )
      IF ( J + DIM1 - 1 .GE. ILAST ) THEN
         DIM1 = WSIZE
      END IF
      DIM4 = WSIZE - DIM1
      IF ( NPROW * NPCOL .EQ. 1 ) THEN
         DIM1 = WSIZE
         DIM4 = 0
      END IF

      IF ( PARTOWNER ) THEN
*        Initialize PU and PV to the Identity.
         CALL DLASET( 'A', WSIZE, WSIZE, ZERO, ONE, PU, LDUV )
         CALL DLASET( 'A', WSIZE, WSIZE, ZERO, ONE, PV, LDUV )

*        Gather elements to chase bulges.
         CALL PDLACP4( J, ILAST, WSIZE,
     $      A, B, DESCA,
     $      TMPA, TMPB, LDTMPAB,
     $      1, 0 )

         IF ( .NOT. OWNER) GOTO 20
         IF ( DIM4 .EQ. 0 .AND. MYROW .EQ. RSRC4 .AND. MYCOL .EQ. CSRC4
     $      .AND. NPCOL * NPROW .GT. 1 ) GOTO 20
         
         DO 30 IB = NBULGES, 1, -1
            NCOLS = 3
            JJ = J + ( IB - 1 ) * BULGESIZE
*           Skip if current bulge already moved over the edge.
            IF ( JJ .GT. ( ILAST - 2 ) ) GOTO 30
            NNB = MIN( NCOLS + MOVEDIST + 1, ILAST - JJ + 1 )
            IF (NNB .LE. 0) GOTO 30
*           Check if Last block.
            IF ( ( JJ + NNB - 1 ) .EQ. ILAST ) NCOLS = 2
*           Perform the actual chase.
            CALL DHGEQZ5( NNB, TMPA, LDTMPAB, TMPB, LDTMPAB, ( IB - 1) *
     $         BULGESIZE, NCOLS, WSIZE, ( NBULGES - IB + 1 ) * 2 - 1,
     $         ALPHAR, ALPHAI, BETA, PU, PV, LDUV, ATOL ) 

*           Check if possible to move the current bulge over the edge.
            IF ( NCOLS .EQ. 2 ) THEN   
               TEMP = TMPA( WSIZE - 1, WSIZE - 2 )
               CALL DLARTG( TEMP, TMPA( WSIZE, WSIZE - 2 ), C, S, TMPA(
     $            WSIZE - 1, WSIZE - 2 ) )
               TMPA( WSIZE, WSIZE - 2 ) = ZERO

*              Update rows of A.
               CALL DROT( 2, TMPA( WSIZE - 1, WSIZE - 1 ), LDTMPAB,
     $            TMPA( WSIZE, WSIZE - 1 ), LDTMPAB, C, S )

*              Update rows of B.
               CALL DROT(2, TMPB( WSIZE - 1, WSIZE - 1 ), LDTMPAB, TMPB(
     $            WSIZE, WSIZE - 1 ), LDTMPAB, C, S )

*              Update cols of U (= Q)
               CALL DROT( WSIZE, PU( 1, WSIZE - 1 ), 1, 
     $                PU( 1, WSIZE ), 1, C, S )               

*              Annihilate fill in B.
               TEMP = TMPB( WSIZE, WSIZE )
               CALL DLARTG( TEMP ,TMPB( WSIZE, WSIZE - 1 ), C, S, TMPB(
     $            WSIZE, WSIZE ) )

               TMPB( WSIZE, WSIZE - 1 ) = ZERO

*              Update cols of A.
               CALL DROT( WSIZE, TMPA( 1, WSIZE ), 1, TMPA( 1, WSIZE - 1
     $            ), 1, C, S )

*              Update cols of B.
               CALL DROT( WSIZE - 1, TMPB( 1, WSIZE ), 1, TMPB( 1, WSIZE
     $            - 1 ), 1, C, S )

*              Update cols of PV (= Z ).
               CALL DROT( WSIZE, PV( 1, WSIZE ), 1, PV( 1, WSIZE - 1 ),
     $            1, C, S )
            END IF
 30      CONTINUE
 20      CONTINUE
         
*        Copy updated TMPA and TMPB to A and B.
         CALL PDLACP4( J, ILAST, WSIZE, A, B, DESCA, TMPA, TMPB, LDTMPAB
     $      , 1, 1 )

      END IF

      
      UPDTIME( 10 ) = UPDTIME( 10 ) + ( MPI_WTIME() - T1 )
*     Should we distribute PU and PV and perform full update of A and B ? 
      IF ( UPDOFFD ) THEN
         T1 = MPI_WTIME()
*        Distribute U (=Q) and V(=Z)
*        U both directions(needed for Q accumulations), V only columnwize  
         IF ( MYROW .EQ. RSRC1 .AND. MYCOL .EQ. CSRC1 ) THEN
            IF ( ILQ ) THEN
*              Broadcast send PU columnwise
               CALL DGEBS2D( ICTXT, 'C', ' ', WSIZE, WSIZE,
     $            PU, LDUV )
            END IF
*           Broadcast send PV columnwise
            CALL DGEBS2D( ICTXT, 'C', ' ', WSIZE, WSIZE, PV, LDUV )
         ELSE IF ( MYCOL .EQ. CSRC1 ) THEN
            IF ( ILQ ) THEN
*              Broadcast receive PU columnwise
               CALL DGEBR2D( ICTXT, 'C', ' ', WSIZE, WSIZE, PU, LDUV,
     $            RSRC1, CSRC1 )
            END IF
*           Broadcast receive PV columnwise
            CALL DGEBR2D( ICTXT, 'C', ' ',  WSIZE, WSIZE, PV, LDUV,
     $         RSRC1, CSRC1 )                 
         END IF

         IF ( MYROW .EQ. RSRC1 .AND. MYCOL .EQ. CSRC1 ) THEN
*           Broadcast send PU rowwise
            CALL DGEBS2D( ICTXT, 'R', ' ', WSIZE, WSIZE, PU, LDUV )
         ELSE IF ( MYROW .EQ. RSRC1 ) THEN
*           Broadcast receive PU rowwise
            CALL DGEBR2D( ICTXT, 'R', ' ', WSIZE, WSIZE, PU, LDUV, RSRC1
     $         , CSRC1 ) 
         END IF

         IF (DIM4.GT.0) THEN
            IF ( MYROW .EQ. RSRC4 .AND. MYCOL .EQ. CSRC4 ) THEN
               IF ( ILQ ) THEN
*                 Broadcast send PU columnwise
                  CALL DGEBS2D(ICTXT, 'C', ' ', WSIZE, WSIZE,
     $               PU, LDUV)
               END IF
*              Broadcast send PV columnwise
               CALL DGEBS2D(ICTXT, 'C', ' ', WSIZE, WSIZE, PV, LDUV )
            ELSE IF ( MYCOL .EQ. CSRC4 ) THEN
               IF ( ILQ ) THEN
*                 Broadcast receive PU columnwise
                  CALL DGEBR2D(ICTXT, 'C', ' ', WSIZE, WSIZE, PU, LDUV,
     $               RSRC4, CSRC4)
               END IF
*              Broadcast receive PV columwise
               CALL DGEBR2D( ICTXT, 'C', ' ',  WSIZE, WSIZE, PV, LDUV,
     $            RSRC4, CSRC4 )
            END IF

            IF ( MYROW .EQ. RSRC4 .AND. MYCOL .EQ. CSRC4 ) THEN
*              Broadcast send PU rowwise
               CALL DGEBS2D( ICTXT, 'R', ' ', WSIZE, WSIZE, PU, LDUV )
            ELSE IF ( MYROW .EQ. RSRC4 ) THEN
*              Broadcast receive PU rowwise
               CALL DGEBR2D( ICTXT, 'R', ' ', WSIZE, WSIZE, PU, LDUV,
     $            RSRC4, CSRC4 )
            END IF
         END IF
         
         
*        Update remaning parts of A and B, Q, and Z with dgemm
         CALL PDHGEQZ6( 'A', ILSCHR, ILQ, ILZ, A, DESCA, B, DESCB, Q,
     $      DESCQ, Z, DESCZ, PU, PV, LDUV, 1, J, J + WSIZE - 1, N, WSIZE
     $      , ILOQ, IHIQ, ILOZ, IHIZ, LDWORK, DWORK, INFO )
         
         UPDTIME( 11 ) = UPDTIME( 11 ) + ( MPI_WTIME() - T1 )

      END IF
      RETURN
*     
*     End of PDHGEQZB
*     
      END

