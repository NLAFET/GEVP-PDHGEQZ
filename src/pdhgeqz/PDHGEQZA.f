***********************************************************************
*                                                                     *
*     PDHGEQZA.f:                                                     *
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
      SUBROUTINE PDHGEQZA( ILSCHR, ILQ, ILZ, N, NBULGES, IFIRST, ILAST,
     $   J, ILO, ILOQ, IHIQ, ILOZ, IHIZ, A, DESCA, B, DESCB, Q, DESCQ, Z
     $   , DESCZ, ALPHAR, ALPHAI, BETA, ATOL, PU, PV, LDUV, TMPA, TMPB,
     $   LDTMPAB, UPDTIME, DWORK, LDWORK, INFO )

      IMPLICIT NONE
*     ..
*     .. Scalar Arguments ..
*     ..
      INTEGER             IFIRST, ILAST, J, ILO, N, LDUV, NBULGES, INFO,
     $   LDWORK, LDTMPAB, ILOQ, IHIQ, ILOZ, IHIZ
      LOGICAL             ILQ, ILZ, ILSCHR
      DOUBLE PRECISION    ATOL
*     ..
*     .. Array Arguments ..
*     ..
      DOUBLE PRECISION    A( * ), B( * ), Q( * ), Z( * ), DWORK( * ),
     $   UPDTIME( 40 ), ALPHAI( * ), ALPHAR( * ), BETA( * ), PU( LDUV, *
     $   ) , PV( LDUV, * ), TMPA( LDTMPAB, * ), TMPB( LDTMPAB, * ) 
      INTEGER             DESCA( * ), DESCB( * ), DESCQ( * ), DESCZ( * )
*     Purpose
*     =======
*     
*     PDHGEQZA creates NBULGES and chases them down within a 
*     diagonal block of H and T. The bulges are only chased so
*     much so that all NBULGES can fit. 
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
*     N       (global input) INTEGER 
*     Order of the matrices A, B, Q, and Z.
*     
*     NBULGES (global input) INTEGER
*     Number of bulges to be created.
*     
*     IFIRST  (global input) INTEGER
*     Minimum row/column index of where to chase bulges.
*     
*     ILAST   (global input) INTEGER
*     Maxiumu row/column index for where to chase bulges.
*     
*     J       (global input) INTEGER 
*     Row/column index of where to first create bulges.
*     
*     ILO     (global input) INTEGER
*     Updates wont be applied below this row/column index.
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
*     ALPHAR  (global input) DOUBLE PRECISION array of size ( NBULGES * 2 )
*     ALPHAI  (global input) DOUBLE PRECISION array of size ( NBULGES * 2 )
*     BETA    (global input) DOUBLE PRECISION array of size ( NBULGES * 2 )
*     ALPHAR( 1 : NBULGES * 2 ) contains the real parts and 
*     ALPHAI( 1 : NBUGLES * 2 ) contains the imaginary
*     parts and BETA( 1 : NBULGES * 2 ) contains the scale parts of 
*     the shifts that define the multi-shift QZ sweep.
*     
*     ATOL    (global input) DOUBLE PRECISION 
*     Holds number where safe to set a value of A to zero, used 
*     for vigialent deflations.
*     
*     PU      (local input) DOUBLE PRECISION, array, dimension (LDUV, *) 
*     Local workspace, minium size (3 * NBULGES + 1, 3 * NBULGES + 1)
*     
*     PV      (local input) DOUBLE PRECISION, array, dimension (LDUV, *)
*     Local workspace, minium size (3 * NBULGES + 1, 3 * NBULGES + 1)
*     
*     LDUV    (global input) INTEGER
*     Leading dimension for the local workpaces PU and PV.
*     
*     TMPA    (local input) DOUBLE PRECISION, array, dimension (LDTMPAB, *).
*     Local workspace, minium size (3 * NBULGES + 1, 3 * NBULGES + 1)
*     
*     TMPB    (local input) DOUBLE PRECISION, array, dimension (LDTMPAB, *).
*     Local workspace, minium size (3 * NBULGES + 1, 3 * NBULGES + 1)
*     
*     LDTMPAB (global input) INTEGER
*     Leading dimension for the local workspace TMPA and TMPB.
*     
*     UPDTIME (local input) INTEGER, array dimension 40.
*     Storage for misc. timers.
*     
*     DWORK   (local workspace/global output) DOUBLE PRECISION array,
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
      DOUBLE PRECISION   ZERO, ONE
      PARAMETER          ( ZERO = 0.0D+0, ONE = 1.0D+0)
      INTEGER            DOWN, UP, RIGHT, LEFT
      PARAMETER          ( DOWN = 1, UP = 2, RIGHT = 3, LEFT = 4 )
      INTEGER            ROWWISE, COLUMNWISE, ALL
      PARAMETER          ( ROWWISE = 1, COLUMNWISE = 2, ALL=3)
      INTEGER            BLOCK_CYCLIC_2D, CSRC_, CTXT_, DLEN_, DT_, 
     $   LLD_, MB_, M_, NB_, N_, RSRC_
      PARAMETER          ( BLOCK_CYCLIC_2D = 1, DLEN_ = 9, DT_ = 1,
     $   CTXT_ = 2, M_ = 3, N_ = 4, MB_ = 5, NB_ = 6,
     $   RSRC_ = 7, CSRC_ = 8, LLD_ = 9 )
      
      INTEGER            NSBULGE, BULGESIZE
      PARAMETER          (NSBULGE = 2, BULGESIZE = NSBULGE + 1)
      
*     ..
*     .. Local Scalars ..
*     ..
      LOGICAL             OWNER
      INTEGER             ISHIFT, IBULGE, LTOP, ICTXT, NPROW, NPCOL,
     $   MYROW, MYCOL, IAROW1, IACOL1, IAROW2, IACOL2, WSIZE, IAM,
     $   NPROCS, JA, JJ, IB, NNB, NCOLS
      DOUBLE PRECISION    TAU, T1, C, S, TEMP
*     ..
*     .. Local Arrays
*     ..

      DOUBLE PRECISION    V( 4 )



*     ..
*     .. Intrinsic Functions ..
*     ..
      INTRINSIC           MIN, MOD
*     ..
*     .. Externals ..
*     ..
      EXTERNAL            INDXG2L, INDXG2P, MPI_WTIME
      INTEGER             INDXG2L, INDXG2P
      DOUBLE PRECISION    MPI_WTIME
      

*     
*     Extract current communication context and get grid parameters
*     
      ICTXT = DESCA( CTXT_ )
      CALL BLACS_GRIDINFO( ICTXT, NPROW, NPCOL, MYROW, MYCOL )  
      CALL BLACS_PINFO( IAM, NPROCS )

      IF( ILSCHR ) THEN
         LTOP = ILO
      ELSE
         LTOP = J
      END IF
      


*     .. Number of bulges going so far
      IBULGE = 0 
      
      
      WSIZE = BULGESIZE * NBULGES + 1
      WSIZE = MIN( WSIZE, ILAST - J + 1)
      


      T1 = MPI_WTIME()
      
      IAROW1 = INDXG2P( J, DESCA( NB_ ), DESCA( RSRC_ ), 
     $   DESCA( RSRC_ ), NPROW )
      IACOL1 = INDXG2P( J, DESCA( NB_ ), DESCA( CSRC_ ), 
     $   DESCA( CSRC_ ), NPCOL )
      IAROW2 = INDXG2P( J + WSIZE - 1, DESCA( NB_ ), DESCA( RSRC_ ), 
     $   DESCA(RSRC_), NPROW )
      IACOL2 = INDXG2P( J + WSIZE - 1, DESCA( NB_ ), DESCA( CSRC_ ), 
     $   DESCA( CSRC_ ), NPCOL )
      
      OWNER = MYROW .EQ. IAROW1 .AND. MYCOL .EQ. IACOL1
      
*     Gather elements to introduce shifts - stored at 
*     process (IARWO1, IACOL1)
      CALL PDLACP3( WSIZE, J, A, DESCA,
     $   TMPA, LDTMPAB, IAROW1, IACOL1, 0)
      CALL PDLACP3( WSIZE, J, B, DESCB,
     $   TMPB, LDTMPAB, IAROW1, IACOL1, 0)

      
      
      IF ( OWNER ) THEN
*        Initialize PU and PV to the Identity.
         CALL DLASET( 'A', WSIZE, WSIZE, ZERO, ONE, PU,
     $      LDUV )
         CALL DLASET( 'A', WSIZE, WSIZE, ZERO, ONE, PV,
     $      LDUV ) 

         
         DO 10 ISHIFT=1,NBULGES * 2, 2
            
*           Number of bulges going so far.
            IBULGE = IBULGE + 1    
            
*           Introduce a bulge.
            CALL QZFCOL( 3, 2, ALPHAR( ISHIFT ), ALPHAI( ISHIFT ), BETA(
     $         ISHIFT ), TMPA, LDTMPAB, TMPB, LDTMPAB, V, INFO )
            CALL DLARFG( 3, V( 1 ), V( 2 ), 1, TAU )
            V( 1 ) = ONE
            
*           Update rows of A.
            CALL DLARFX( 'Left', 3, WSIZE, V, TAU, TMPA, LDTMPAB, DWORK
     $         )
*           Update rows of B.
            CALL DLARFX( 'Left', 3, WSIZE, V, TAU, TMPB, LDTMPAB, DWORK
     $         )
*           Update columns of PU.
            CALL DLARFX( 'Right', WSIZE, 3, V, TAU, PU, LDUV, DWORK )

*           Zero j-th column of B.
            CALL INVHSE(3, TMPB, LDTMPAB, TAU, V, INFO)
            
*           Update columns of B.
            CALL DLARFX( 'Right', 3, 3, V, TAU,
     $         TMPB, LDTMPAB, DWORK )
            TMPB(2,1) = ZERO
            TMPB(3,1) = ZERO

*           Update columns of A.
            CALL DLARFX( 'Right', 4, 3, V, TAU, TMPA, LDTMPAB, DWORK )

*           Update columns of PV.
            CALL DLARFX( 'Right', WSIZE, 3, V, TAU, PV, LDUV, DWORK )
            
            IF ( IBULGE .LT. NBULGES ) THEN
               
*              Hunt bulges one step forward to prepare for next.
               JA = IFIRST
               DO 20 IB = IBULGE, 1, -1
                  NCOLS = 3
                  JJ = JA + ( IB - 1 ) * ( BULGESIZE )
                  NNB = MIN( BULGESIZE + 1 + NCOLS, ILAST -JJ + 1 )
                  IF ( NNB .LE. 0 ) GOTO 20

*                 Check if last block
                  IF ((JJ + NNB -1) .EQ. ILAST) NCOLS = 2
*                 Preform the actual chase.
                  CALL DHGEQZ5( NNB,
     $               TMPA, LDTMPAB,
     $               TMPB, LDTMPAB,
     $               ( IB - 1 ) * BULGESIZE, NCOLS, WSIZE,
     $               ( IBULGE - IB + 1 ) * 2 - 1,
     $               ALPHAR, ALPHAI, BETA,
     $               PU, PV, LDUV, ATOL )

*                 Check if possible to push of.
                  IF ( NCOLS .EQ. 2 ) THEN
                     TEMP = TMPA( WSIZE - 1, WSIZE - 2 )
                     CALL DLARTG( TEMP, TMPA( WSIZE, WSIZE - 2 ), C, S,
     $                  TMPA( WSIZE - 1, WSIZE - 2 ) )
                     TMPA( WSIZE, WSIZE - 2 ) = ZERO

*                    Update rows of A.                
                     CALL DROT(2, TMPA( WSIZE - 1, WSIZE - 1 ), WSIZE,
     $                  TMPA(WSIZE, WSIZE - 1 ), WSIZE, C, S )

*                    Update rows of B.
                     CALL DROT( 2, TMPB( WSIZE - 1, WSIZE - 1 ), WSIZE,
     $                  TMPB( WSIZE, WSIZE - 1 ), WSIZE, C, S )

*                    Update cols of PU (=Q)                                      
                     CALL DROT( WSIZE, PU( 1, WSIZE - 1 ), 1,
     $                  PU( 1, WSIZE ), 1, C, S )

*                    Annihilate fill in B.
                     TEMP = TMPB( WSIZE, WSIZE )
                     CALL DLARTG(TEMP, TMPB( WSIZE, WSIZE - 1 ), C, S,
     $                  TMPB( WSIZE, WSIZE )) 
                     TMPB( WSIZE, WSIZE - 1 ) = ZERO

*                    Update columns of A.
                     CALL DROT(WSIZE, TMPA( 1, WSIZE ), 1, TMPA( 1,
     $                  WSIZE - 1 ), 1, C, S)

*                    Update columns of B.
                     CALL DROT( WSIZE - 1, TMPB( 1, WSIZE ), 1, TMPB( 1,
     $                  WSIZE - 1 ), 1, C, S )
                     
*                    Update columns of PV.
                     CALL DROT(WSIZE, PV( 1, WSIZE ), 1, PV( 1, WSIZE -
     $                  1 ), 1, C, S )
                  END IF
 20            CONTINUE
            END IF     
 10      CONTINUE
      END IF
      
*     Copy updated TMPA and TMPB to A and B.
      CALL PDLACP3( WSIZE, J, A, DESCA, TMPA, LDTMPAB, IAROW1, IACOL1,
     $   -1 )
      CALL PDLACP3( WSIZE, J, B, DESCB, TMPB, LDTMPAB, IAROW1, IACOL1,
     $   -1 )

      UPDTIME(10) = UPDTIME(10) + (MPI_WTIME() - T1)
      T1 = MPI_WTIME()
*     
*     Distribute U (=Q) and V(=Z).
*     U both directions(needed for Q accumulations), V only columnwize
*     But first to right/below.
*     
      IF ( IACOL1 .NE. IACOL2 .OR. IAROW1 .NE. IAROW2 ) THEN          
         IF ( MYROW .EQ. IAROW1 .AND. MYCOL .EQ. IACOL1 ) THEN
*           Send PU and PV to right/below.
            CALL DGESD2D( ICTXT, WSIZE, WSIZE, PU, LDUV, IAROW2, IACOL2
     $         )          
            CALL DGESD2D( ICTXT, WSIZE, WSIZE, PV, LDUV, IAROW2, IACOL2
     $         )              
         ELSE IF ( MYROW .EQ. IAROW2 .AND. MYCOL .EQ. IACOL2 ) THEN
*           Receive PU and PV from left/up.
            CALL DGERV2D( ICTXT, WSIZE, WSIZE, PU, LDUV, IAROW1, IACOL1
     $         )   
            CALL DGERV2D( ICTXT, WSIZE, WSIZE, PV, LDUV, IAROW1, IACOL1
     $         )               
         END IF
      END IF
      IF ( MYROW .EQ. IAROW1 .AND. MYCOL .EQ. IACOL1 ) THEN
         IF ( ILQ ) THEN
*           Broadcast send PU columnwise
            CALL DGEBS2D( ICTXT, 'C', ' ', WSIZE, WSIZE,
     $         PU, LDUV )
         END IF
*        Broadcast send PV columnwise.
         CALL DGEBS2D( ICTXT, 'C', ' ', WSIZE, WSIZE, PV, LDUV )
      ELSE IF ( MYCOL .EQ. IACOL1) THEN
         IF ( ILQ ) THEN
*           Broadcast receive PU columnwise
            CALL DGEBR2D( ICTXT, 'C', ' ', WSIZE, WSIZE, PU, LDUV,
     $         IAROW1, IACOL1 )  
         END IF
*        Broadcast receive PV columnwise.
         CALL DGEBR2D( ICTXT, 'C', ' ', WSIZE, WSIZE, PV, LDUV,
     $      IAROW1, IACOL1 )
      END IF
      IF ( MYROW .EQ. IAROW1 .AND. MYCOL .EQ. IACOL1 ) THEN
*        Broadcast send PU rowwise.
         CALL DGEBS2D(ICTXT, 'R', ' ', WSIZE, WSIZE, PU, LDUV)
      ELSE IF ( MYROW .EQ. IAROW1 ) THEN
*        Broadcast receive PU rowwise.
         CALL DGEBR2D(ICTXT, 'R', ' ', WSIZE, WSIZE, PU, LDUV,
     $      IAROW1, IACOL1)
      END IF
*     If TMPA and TMPB is split over a 2x2 grid or more (or a 2x1/1x2 grid)
      IF ( IAROW1 .NE. IAROW2 .OR. IACOL1 .NE. IACOL2 ) THEN
         IF (MYROW .EQ. IAROW2 .AND. MYCOL .EQ. IACOL2 ) THEN
            IF ( ILQ ) THEN
*              Broadcast send PU columwise
               CALL DGEBS2D( ICTXT, 'C', ' ', WSIZE, WSIZE,
     $            PU, LDUV )
            END IF
*           Broadcast send PV columnwise.
            CALL DGEBS2D( ICTXT, 'C', ' ', WSIZE, WSIZE, PV, LDUV )
         ELSE IF ( MYCOL .EQ. IACOL2) THEN
            IF ( ILQ ) THEN
*              Broadcast receive PU columnwise
               CALL DGEBR2D( ICTXT, 'C', ' ', WSIZE, WSIZE, PU, LDUV,
     $            IAROW2, IACOL2 )
            END IF
*           Broadcast receive PV columnwise.
            CALL DGEBR2D( ICTXT, 'C', ' ', WSIZE, WSIZE, PV, LDUV,
     $         IAROW2, IACOL2 )
         END IF
         IF ( MYROW .EQ. IAROW2 .AND. MYCOL .EQ. IACOL2 ) THEN
*           Broadcast send PU rowwise.
            CALL DGEBS2D( ICTXT, 'R', ' ', WSIZE, WSIZE, PU, LDUV )
         ELSE IF ( MYROW .EQ. IAROW2 ) THEN
*           Broadcast receive PU rowwise.
            CALL DGEBR2D( ICTXT, 'R', ' ', WSIZE, WSIZE, PU, LDUV,
     $         IAROW2, IACOL2 )
         END IF
      END IF    

*     Update remaning parts of A and B, and Q, Z with dgemm	
      CALL PDHGEQZ6( 'A', ILSCHR, ILQ, ILZ, A, DESCA, B, DESCB, Q, DESCQ
     $   ,Z, DESCZ, PU, PV, LDUV, 1, J, J + WSIZE - 1, N, WSIZE, ILOQ
     $   ,IHIQ, ILOZ, IHIZ, LDWORK, DWORK, INFO )

      UPDTIME( 11 ) = UPDTIME( 11 ) + ( MPI_WTIME() - T1 )
      


      RETURN
*     End of PDHGEQZA

      END 
