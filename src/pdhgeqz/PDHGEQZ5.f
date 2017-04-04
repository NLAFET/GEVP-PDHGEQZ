***********************************************************************
*                                                                     *
*     PDHGEQZ5.f:                                                     *
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
      SUBROUTINE PDHGEQZ5( ILSCHR, ILQ, ILZ, NS, IFIRST, ILAST, N, ILO,
     $   IHI, A, DESCA, B, DESCB, ALPHAR, ALPHAI, BETA, Q, DESCQ, Z,
     $   DESCZ, ILOQ, IHIQ, ILOZ, IHIZ, NUMWIN, USEDSHIFTS, ATOL, WORK,
     $   LWORK, INFO, UPDTIME, SHIFTSPERWND )
      IMPLICIT NONE
*     .. 
*     .. Scalar Arguments ..
*     ..
      INTEGER             IHI, ILO, LWORK, N, NS,
     $   NUMWIN, USEDSHIFTS,
     $   ILOQ, IHIQ, ILOZ, IHIZ,
     $   IFIRST, ILAST, SHIFTSPERWND, INFO
      LOGICAL             ILQ, ILZ, ILSCHR
      DOUBLE PRECISION    ATOL
*     ..
*     .. Array Arguments .. 
*     ..
      DOUBLE PRECISION    A( * ), ALPHAI( * ), ALPHAR( * ), B( * ),
     $   BETA( * ),Q( * ), WORK( * ), Z( * )
      INTEGER             DESCA( 9 ), DESCB( 9 ), DESCQ( 9 ), DESCZ( 9 )
      DOUBLE PRECISION    UPDTIME( 40 )

*     Purpose
*     =======
*     
*     This auxiliary subroutine called by PDHGEQZ0 and PDHGEQZ1 and performs a
*     small-bulge multi-shift QZ sweep by chasing separated
*     groups of bulges along the main block diagonal of A and B, while keeping
*     Q and Z updated from right.
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
*     NS      (global input) INTEGER
*     NS gives the number of simultaneous shifts. 
*     
*     IFIRST  (global input) INTEGER
*     Row/column index of where to first insert shifts.
*     
*     ILAST   (global input) INTEGER
*     Row/column index for the last index to apply shifts to.
*     
*     N       (global input) INTEGER 
*     Order of the matrices A, B, Q, and Z.
*     
*     ILO     (global input) INTEGER
*     Updates wont be applied below this row/column index.
*     
*     IHI     (global input) INTEGER
*     Updates wont be applied to indexes greater than this.
*     
*     A       (local input/output) DOUBLE PRECISION array of size 
*     (DESCA(LLD_),*)
*     On input A contains a Hessenberg matrix.  On output a
*     multi-shift QZ sweep with shifts (ALPHAR(J)+i*ALPHAI(J)/BETA(J) is 
*     applied to the isolated diagonal block in rows and columns IFIRST
*     through ILAST.
*     
*     DESCA   (global and local input) INTEGER array of dimension DLEN_.
*     The array descriptor for the distributed matrix A.
*     
*     
*     B       (local input/output) DOUBLE PRECISION array of size  
*     (DESCB(LLD_),*)
*     On input B contains an upper triangular matrix.  On output a
*     multi-shift QZ sweep with shifts (ALPHAR(J)+i*ALPHAI(J)/BETA(J) is 
*     applied to the isolated diagonal block in rows and columns IFIRST 
*     through ILAST. 
*     
*     DESCB   (global and local input) INTEGER array of dimension DLEN_. 
*     The array descriptor for the distributed matrix B.
*     
*     ALPHAR  (global input) DOUBLE PRECISION array of size (N)
*     ALPHAI  (global input) DOUBLE PRECISION array of size (N)
*     BETA    (global input) DOUBLE PRECISION array of size (N)
*     ALPHAR(KBOT-NS+1:KBOT) contains the real parts and 
*     ALPHAI(KBOT-NS+1:KBOT contains the imaginary
*     parts and BETA(KBOT-NS+1:KBOT) contains the scale parts of 
*     the NS shifts of origin that define the
*     multi-shift QZ sweep.
*     
*     Q       (local input/output) DOUBLE PRECISION array of size
*     (DESCQ(LLD_),*)
*     If ILQ = .TRUE., then the QZ Sweep left orthogonal
*     equivalence transformation is accumulated into
*     Q(1:N,ILOQ:IHIQ) from the right.
*     If ILZ = .FALSE., then Q is unreferenced.
*     
*     DESCQ   (global and local input) INTEGER array of dimension DLEN_.
*     The array descriptor for the distributed matrix Q.
*     
*     Z (local input/output) DOUBLE PRECISION array of size (DESCZ(LLD_),*) 
*     If ILZ = .TRUE., then the QZ Sweep right
*     orthogonal equivalence transformation is accumulated into
*     Z(1:N,ILOZ:IHIZ) from the right.  If ILZ = .FALSE., then Z is
*     unreferenced.
*     
*     DESCZ   (global and local input) INTEGER array of dimension DLEN_.
*     The array descriptor for the distributed matrix Z.
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
*     NUMWIN  (global input) INTEGER
*     Number of concurrent windows to use
*     
*     USEDSHIFTS(global input/output) INTEGER
*     Incremented value of number of shifts applied
*     
*     ATOL    (global input) DOUBLE PRECISION 
*     Holds number where safe to set a value of A to zero, used 
*     for vigialent deflations.
*     
*     WORK    (local workspace) DOUBLE PRECISION array, dimension(DWORK)
*     
*     LWORK   (local input) INTEGER
*     The length of the workspace array WORK.
*     
*     INFO    (global output) INTEGER
*     
*     UPDTIME (global input/output) DOUBLE PRECISION array, dimension 40.
*     Holds timings for different subroutines.
*     
*     SHIFTSPERWND (global input) INTEGER
*     Number of shifts to be applied for each window.
*     
*     
*     ================================================================
*     

*     ..
*     .. Parameters ..
*     ..
      INTEGER            BLOCK_CYCLIC_2D, CSRC_, CTXT_, DLEN_, DT_, 
     $   LLD_, MB_, M_, NB_, N_, RSRC_
      PARAMETER          ( BLOCK_CYCLIC_2D = 1, DLEN_ = 9, DT_ = 1,
     $   CTXT_ = 2, M_ = 3, N_ = 4, MB_ = 5, NB_ = 6,
     $   RSRC_ = 7, CSRC_ = 8, LLD_ = 9 )
      
      INTEGER            NSBULGE, BULGESIZE
      PARAMETER         (NSBULGE = 2,BULGESIZE = NSBULGE + 1)
      
      
*     ..
*     .. Local Scalars ..
*     ..
      
      INTEGER             J, CTX, WORKMIN, MOVEDIST, NUMCHASE, SHIFTS,
     $   BULGES, NPROW, NPCOL, MYROW, MYCOL, TOBORDER, CHAINSIZE,
     $   CURRWIN, IC, WSIZE, RSRC1, CSRC1, RSRC2, CSRC2, RSRC3, CSRC3, 
     $   RSRC4, CSRC4, DIM1, DIM4, SWIN, MAXWSIZE, RPT, LDUV, IAM, 
     $   NPROCS, NB, IPWORK, IERR
      DOUBLE PRECISION    T1
      LOGICAL             LQUERY, CHASECOMPLETEALL, DEBUG
*     ..
*     ..  Local Arrays
*     ..
      INTEGER             WINPOS(40), WINSHIFTPOS(40), PREVWINPOS(40),
     $   UPOS(40), VPOS(40)
      LOGICAL             CHASECOMPLETE(40), WINXBORDER(40),
     $   PREVWINXBORDER(40)
      DOUBLE PRECISION, ALLOCATABLE :: UV( :, : ), TMPA( :, : ), 
     $   TMPB( :, : )      

*     ..
*     .. External Functions ..
*     ..
      EXTERNAL            LSAME, INDXG2L, INDXG2P, MPI_WTIME
      LOGICAL             LSAME
      INTEGER             INDXG2L, INDXG2P
      DOUBLE PRECISION    MPI_WTIME

      
*     ..
*     .. Externals ..
*     ..


      
*     ..
*     .. Intrinsic Functions ..
*     ..
      INTRINSIC           ABS, DBLE, MAX, MIN, SQRT, MOD, LOG
      
      INFO = 0
      IF( N.LE.0 ) THEN
         WORK( 1 ) = DBLE( 1 )
         RETURN
      END IF
      
*     If IHI < ILO, skip QZ sweep step.
      IF( IHI .LT. ILO )  RETURN
      J = IFIRST
      

      LQUERY = LWORK .EQ. -1

*     Extract current communication context and get grid parameters
      CTX = DESCA( CTXT_ )
      CALL BLACS_PINFO( IAM, NPROCS )
      CALL BLACS_GRIDINFO( CTX, NPROW, NPCOL, MYROW, MYCOL )


*     NB_ .EQ. MB_ is assumed, we use NB as blocking factor.
      NB = DESCA( NB_ )

*     Set DEBUG to true/false for various info printed during the run
      DEBUG = .FALSE.
*     Print debuginfo only for root.
      DEBUG = DEBUG .AND. IAM .EQ. 0

      
      
      
      SHIFTS = SHIFTSPERWND
      BULGES = SHIFTS/2      

      CHAINSIZE = BULGESIZE * BULGES
      MAXWSIZE = MAX( MIN( ILAST - IFIRST + 1, NB ), 3 )
      
*     Calulate maximum needed work, for dgemm updates only      
      WORKMIN = N * NB
      
      
      CHASECOMPLETE = .FALSE.
      
      IF ( LQUERY ) THEN
         WORK( 1 ) = WORKMIN
         RETURN
      END IF
      


      LDUV = MAXWSIZE * NUMWIN * 2
      
      ALLOCATE( UV( LDUV, MAXWSIZE ) )
      ALLOCATE( TMPA( MAXWSIZE, MAXWSIZE ) )
      ALLOCATE( TMPB( MAXWSIZE, MAXWSIZE ) )
      
      IPWORK = 1

      IERR = 0

*     Prepare windows      
      DO IC = 1, NUMWIN 
         WINXBORDER( IC ) = .FALSE.
         WINPOS( IC ) = J
         WINSHIFTPOS( IC ) = NS - IC * SHIFTS + 1
         IF (IC.EQ.1) THEN
            UPOS( IC ) = 1
            VPOS( IC ) = UPOS( IC ) + MAXWSIZE
         ELSE
            UPOS( IC ) = VPOS( IC - 1 ) + MAXWSIZE
            VPOS( IC ) = UPOS( IC ) + MAXWSIZE
         END IF
         CHASECOMPLETE(IC) = .FALSE.
      END DO
      
      IF (DEBUG) THEN
         WRITE(*,*)'% PQZ5 : Intro part. Numwin=', NUMWIN
      END IF 
      DO CURRWIN = 1, NUMWIN
*        Introduce the shifts in A and B at J
         J = WINPOS( CURRWIN )
         IF ( J .LT. ILAST - 1 ) THEN
            IF (DEBUG) THEN
               WRITE(*,*)'% PQZ5 : Intro A. Currwin=', CURRWIN, 
     $            MYROW, MYCOL, J
               CALL FLUSH(6)
            END IF 
            CALL PDHGEQZA( ILSCHR, ILQ, ILZ, N, BULGES, IFIRST, ILAST,
     $         WINPOS( CURRWIN ), ILO, ILOQ, IHIQ, ILOZ, IHIZ, A, DESCA,
     $         B , DESCB, Q, DESCQ, Z, DESCZ, ALPHAR( WINSHIFTPOS(
     $         CURRWIN ) ) , ALPHAI( WINSHIFTPOS( CURRWIN ) ), BETA(
     $         WINSHIFTPOS( CURRWIN ) ), ATOL, UV( UPOS( CURRWIN ), 1 ),
     $         UV( VPOS( CURRWIN ),1 ), LDUV, TMPA, TMPB, MAXWSIZE,
     $         UPDTIME, WORK( IPWORK ), LWORK, IERR )
            
            CALL BLACS_BARRIER( CTX, 'A' )
            IF (DEBUG) THEN
               WRITE(*,*)'% PQZ5 : Intro A Done. Currwin=', 
     $            CURRWIN, MYROW, MYCOL, J
               CALL FLUSH(6)
            END IF 
            USEDSHIFTS = USEDSHIFTS + SHIFTS
         END IF
         IF ( CURRWIN .NE. NUMWIN ) THEN
            DO RPT = 1, 2
*              Chase the created bulge(s) down the diagonal of 
*              H and T. At this stage, only so next chain can be introduced
               DO IC = 1, CURRWIN
                  J = WINPOS( IC )
                  IF ( WINXBORDER( IC ) ) THEN
                     MOVEDIST = CHAINSIZE + 1
                     WINXBORDER( IC ) = .FALSE.
                     WSIZE = MIN( ILAST - J + 1, (CHAINSIZE + 1 ) * 2 )
                  ELSE
                     TOBORDER = NB - MOD( J - 1, NB )
                     MOVEDIST = TOBORDER - CHAINSIZE - 1
                     WINXBORDER(IC) = .TRUE.    
                     IF ( MOVEDIST .LE. 0 ) THEN
                        MOVEDIST = CHAINSIZE + MOVEDIST + 1
                        WINXBORDER( IC ) = .FALSE.
                        PREVWINXBORDER( IC ) = .FALSE.
                     END IF      
                     WSIZE = MIN( ILAST - J + 1, CHAINSIZE + MOVEDIST +
     $                  1 )
                  END IF
                  
                  IF ( J .LT. ILAST - 1 ) THEN
                     IF (DEBUG) THEN
                        WRITE(*,*)'% PQZ5 : Intro B. Currwin=', IC, RPT, 
     $                  MYROW, MYCOL, J
                        CALL FLUSH(6)
                     END IF 
                     CALL PDHGEQZB( ILSCHR, ILQ, ILZ, .TRUE., N, BULGES,
     $                  IFIRST, ILAST, J, WSIZE, MOVEDIST, ILO, A, DESCA
     $                  , B, DESCB, Q, DESCQ, Z, DESCZ, ILOZ, ILOQ, IHIQ
     $                  , IHIZ, ALPHAR( WINSHIFTPOS( IC ) ), ALPHAI(
     $                  WINSHIFTPOS( IC ) ), BETA( WINSHIFTPOS( IC ) ),
     $                  ATOL, UV( UPOS( IC ), 1 ), UV( VPOS( IC ), 1 ),
     $                  LDUV, TMPA , TMPB, MAXWSIZE, UPDTIME, WORK(
     $                  IPWORK ), LWORK, IERR )

                     CALL BLACS_BARRIER( CTX, 'A' )
                     IF (DEBUG) THEN
                        WRITE(*,*)'% PQZ5 : Intro B Done. Currwin=', 
     $                     IC, RPT, MYROW, MYCOL, J
                        CALL FLUSH(6)
                     END IF 
                     CHASECOMPLETE( IC ) = .FALSE.
                     WINPOS(IC)  = J + MOVEDIST
                  ELSE
                     CHASECOMPLETE( IC ) = .TRUE.
                  END IF                  
               END DO               
            END DO
         END IF
      END DO
      CHASECOMPLETEALL = .TRUE.
      DO IC = 1, NUMWIN
         IF ( .NOT. CHASECOMPLETE( IC ) ) THEN
            CHASECOMPLETEALL = .FALSE.
         END IF
      END DO

*     Now when all bulges are introduced, continue the chase until all 
*     bulges are chased of.
*     Chase odd windows first, then even
      NUMCHASE = 0
 10   CONTINUE
      NUMCHASE = NUMCHASE + 1
      IF ( DEBUG ) THEN
         WRITE(*,*)'% PQZ5 : Chase part. Chase=', NUMCHASE, 
     $      'Chase complete =', CHASECOMPLETEALL
      END IF 
      IF ( .NOT. CHASECOMPLETEALL ) THEN
         
         DO SWIN = 1, 2
            DO IC = SWIN, NUMWIN, 2
               J = WINPOS( IC )
               PREVWINXBORDER( IC ) = WINXBORDER( IC )
               PREVWINPOS( IC ) = WINPOS( IC )

               IF ( WINXBORDER( IC ) ) THEN                  
                  WINXBORDER( IC ) = .FALSE.
                  MOVEDIST = CHAINSIZE + 1                  
                  WSIZE = MIN( ILAST - J + 1, ( CHAINSIZE + 1 ) * 2 )
               ELSE
                  TOBORDER = NB - MOD( J - 1, NB ) 
                  MOVEDIST = TOBORDER - CHAINSIZE - 1
                  WINXBORDER( IC ) = .TRUE.                  
                  IF ( MOVEDIST .LE. 0 ) THEN 
                     MOVEDIST = CHAINSIZE + MOVEDIST + 1
                     WINXBORDER( IC ) = .FALSE.   
                     PREVWINXBORDER( IC ) = .FALSE. 
                  END IF
                  WSIZE = MIN( ILAST - J + 1, CHAINSIZE + MOVEDIST + 1 )
               END IF
               IF (DEBUG) THEN
                  WRITE(*,*)'% PQZ5 : Chase SWIN=', SWIN, ',J=',J, 
     $               WSIZE, 
     $               MOVEDIST, IC, ILAST,
     $               WINXBORDER( IC ), CHASECOMPLETEALL
               END IF 

               IF ( J .LT. ILAST - 1 .AND. .NOT. CHASECOMPLETE( IC ) ) 
     $            THEN

                  CALL PDHGEQZB( ILSCHR, ILQ, ILZ, .FALSE., N, BULGES,
     $               IFIRST, ILAST, J, WSIZE, MOVEDIST, ILO, A, DESCA, B
     $               , DESCB, Q, DESCQ, Z, DESCZ, ILOZ, ILOQ, IHIQ, IHIZ
     $               , ALPHAR( WINSHIFTPOS( IC ) ), ALPHAI( WINSHIFTPOS(
     $               IC ) ), BETA( WINSHIFTPOS( IC ) ), ATOL, UV( UPOS(
     $               IC ), 1 ), UV( VPOS( IC ), 1 ), LDUV, TMPA, TMPB,
     $               MAXWSIZE, UPDTIME, WORK( IPWORK ), LWORK, IERR )
                  INFO = 0
                  CHASECOMPLETE( IC ) = .FALSE.
                  WINPOS(IC)  = J + MOVEDIST
               ELSE
                  CHASECOMPLETE( IC ) = .TRUE.
               END IF                
            END DO
         END DO

         CHASECOMPLETEALL = .TRUE.
         DO IC = 1, NUMWIN
            IF ( .NOT. CHASECOMPLETE( IC ) ) THEN
               CHASECOMPLETEALL = .FALSE.
            END IF
         END DO
         
         
         T1 = MPI_WTIME()

*        Distribute U (=Q) and V(=Z) for each window
*        Start with Columnwise broadcasts

*        Extract processors for current window,
*        as below:
*        
*        1 | 2
*        -----
*        3 | 4
         DO SWIN = 1, 2
            DO IC = SWIN, NUMWIN, 2
               J = PREVWINPOS(IC)
               
               RSRC1 = INDXG2P( J, DESCA( NB_ ), DESCA( RSRC_ ),
     $            DESCA( RSRC_ ), NPROW )
               CSRC1 = INDXG2P( J, DESCA( NB_), DESCA( CSRC_),
     $            DESCA( CSRC_ ), NPCOL )
               
               RSRC2 = RSRC1
               CSRC2 = MOD( CSRC1 + 1, NPCOL )
               RSRC3 = MOD( RSRC1 + 1, NPROW )
               CSRC3 = CSRC1
               RSRC4 = RSRC3
               CSRC4 = CSRC2    
               
               IF ( PREVWINXBORDER( IC ) ) THEN
                  MOVEDIST = CHAINSIZE + 1                  
                  WSIZE = MIN( ILAST - J + 1, ( CHAINSIZE + 1 ) * 2 )
               ELSE
                  TOBORDER = NB - MOD( J - 1, NB ) 
                  MOVEDIST = TOBORDER - CHAINSIZE - 1
                  IF (MOVEDIST.LE.0) THEN
                     MOVEDIST = CHAINSIZE + MOVEDIST + 1
                  END IF
                  WSIZE = MIN( ILAST-J+1, CHAINSIZE + MOVEDIST + 1 )
               END IF
               
               DIM1 = NB - MOD( J - 1, NB )
               IF ( J +DIM1 - 1 .GE. ILAST ) THEN
                  DIM1 = WSIZE
               END IF
               DIM4 = WSIZE - DIM1
               IF ( NPROW * NPCOL .EQ. 1 ) THEN
                  DIM1 = WSIZE
                  DIM4 = 0
               END IF              
               IF ( J .LT. ILAST - 1 ) THEN                  
                  IF ( MYROW .EQ. RSRC1 .AND. MYCOL .EQ. CSRC1 ) THEN
                     IF (ILQ) THEN
                        CALL DGEBS2D(CTX, 'C', ' ', WSIZE, WSIZE,
     $                       UV(UPOS(IC),1), LDUV)
                        CALL DGEBS2D(CTX, 'C', ' ', WSIZE, WSIZE,
     $                     UV(VPOS(IC),1), LDUV)                        
                     ELSE 
                        CALL DGEBS2D( CTX, 'C', ' ', WSIZE, WSIZE, 
     $                       UV( VPOS( IC ), 1 ), LDUV )
                     END IF                      
                  ELSE IF ( MYCOL .EQ. CSRC1 ) THEN
                     IF ( ILQ ) THEN
                        CALL DGEBR2D(CTX, 'C', ' ', WSIZE, WSIZE,
     $                       UV(UPOS(IC),1), LDUV, 
     $                       RSRC1, CSRC1)
                        CALL DGEBR2D( CTX, 'C', ' ', WSIZE, WSIZE,
     $                     UV( VPOS( IC ), 1 ), LDUV, RSRC1, CSRC1 )                        
                     ELSE 
                        CALL DGEBR2D( CTX, 'C', ' ', WSIZE, WSIZE,
     $                     UV( VPOS( IC ), 1 ), LDUV, RSRC1, CSRC1 )
                     END IF                         
                  END IF              
                  IF ( DIM4 .GT. 0 ) THEN
                     IF ( MYROW .EQ. RSRC4 .AND. MYCOL .EQ. CSRC4 ) THEN
                        IF (ILQ) THEN
                           CALL DGEBS2D(CTX, 'C', ' ', WSIZE,WSIZE,
     $                          UV(UPOS(IC),1), LDUV)
                           CALL DGEBS2D( CTX, 'C', ' ', WSIZE,WSIZE,
     $                        UV( VPOS( IC ), 1 ), LDUV )                           
                        ELSE 
                           CALL DGEBS2D( CTX, 'C', ' ', WSIZE, WSIZE, 
     $                          UV( VPOS( IC ), 1 ), LDUV )
                        END IF                            
                     ELSE IF ( MYCOL .EQ. CSRC4 ) THEN
                        IF (ILQ) THEN
                           CALL DGEBR2D(CTX, 'C', ' ', WSIZE,WSIZE,
     $                          UV(UPOS(IC),1), LDUV, 
     $                          RSRC4, CSRC4)
                           CALL DGEBR2D( CTX, 'C', ' ', WSIZE,WSIZE,
     $                        UV( VPOS( IC ), 1 ), LDUV, RSRC4, CSRC4 )
                        ELSE 
                           CALL DGEBR2D( CTX, 'C', ' ', WSIZE, WSIZE,
     $                        UV( VPOS( IC ), 1 ), LDUV, RSRC4, CSRC4 )
                        END IF                          
                     END IF
                  END IF
               END IF
            END DO
         END DO

*        Now Rowwise broadcasts          
         DO SWIN = 1, 2
            DO IC = SWIN, NUMWIN, 2
               J = PREVWINPOS( IC )
               RSRC1 = INDXG2P( J, DESCA( NB_ ), DESCA( RSRC_ ), 
     $            DESCA( RSRC_ ), NPROW )
               CSRC1 = INDXG2P( J, DESCA( NB_ ), DESCA( CSRC_ ), 
     $            DESCA( CSRC_ ), NPCOL )
               
               RSRC2 = RSRC1
               CSRC2 = MOD( CSRC1 + 1, NPCOL )
               RSRC3 = MOD( RSRC1 + 1, NPROW )
               CSRC3 = CSRC1
               RSRC4 = RSRC3
               CSRC4 = CSRC2    

               IF ( PREVWINXBORDER( IC ) ) THEN
                  MOVEDIST = CHAINSIZE + 1                  
                  WSIZE = MIN( ILAST - J + 1, ( CHAINSIZE + 1 ) * 2 )
               ELSE
                  TOBORDER = NB - MOD( J - 1, NB ) 
                  MOVEDIST = TOBORDER - CHAINSIZE - 1
                  IF ( MOVEDIST .LE. 0 ) THEN
                     MOVEDIST = CHAINSIZE + MOVEDIST + 1
                  END IF
                  WSIZE = MIN( ILAST - J + 1, CHAINSIZE + MOVEDIST + 1 )
               END IF
               
               DIM1 = NB - MOD( J - 1, NB )
               IF ( J + DIM1 - 1 .GE. ILAST ) THEN
                  DIM1 = WSIZE
               END IF
               DIM4 = WSIZE - DIM1
               IF ( NPROW * NPCOL .EQ. 1 ) THEN
                  DIM1 = WSIZE
                  DIM4 = 0
               END IF              
               IF ( J .LT. ILAST-1 ) THEN   
                  IF ( MYROW .EQ. RSRC1 .AND. MYCOL .EQ. CSRC1 ) THEN
                     CALL DGEBS2D( CTX, 'R', ' ', WSIZE, WSIZE, UV(
     $                  UPOS( IC ), 1 ), LDUV )  
                  ELSE IF ( MYROW .EQ. RSRC1 ) THEN
                     CALL DGEBR2D( CTX, 'R', ' ', WSIZE, WSIZE, UV(
     $                  UPOS( IC ), 1 ), LDUV, RSRC1, CSRC1 )
                  END IF              
                  IF ( DIM4 .GT. 0 ) THEN
                     IF ( MYROW .EQ. RSRC4 .AND. MYCOL .EQ. CSRC4 ) THEN
                        CALL DGEBS2D( CTX, 'R', ' ', WSIZE, WSIZE, UV(
     $                     UPOS( IC ), 1 ), LDUV)     
                     ELSE IF ( MYROW .EQ. RSRC4 ) THEN
                        CALL DGEBR2D( CTX, 'R', ' ', WSIZE, WSIZE, UV(
     $                     UPOS( IC ), 1 ), LDUV, RSRC4, CSRC4) 
                     END IF
                  END IF
               END IF
            END DO
         END DO


*        Now the actual update of of-diagonal elements with dgemm
*        First updates from the right
         DO SWIN = 1, 2
            DO IC = SWIN, NUMWIN, 2
               J = PREVWINPOS(IC)

               RSRC1 = INDXG2P( J, DESCA( NB_ ), DESCA( RSRC_ ), 
     $            DESCA( RSRC_ ), NPROW)
               CSRC1 = INDXG2P( J, DESCA( NB_ ), DESCA( CSRC_ ), 
     $            DESCA( CSRC_ ), NPCOL)
               
               RSRC2 = RSRC1
               CSRC2 = MOD( CSRC1+1, NPCOL )
               RSRC3 = MOD( RSRC1+1, NPROW )
               CSRC3 = CSRC1
               RSRC4 = RSRC3
               CSRC4 = CSRC2    

               IF ( PREVWINXBORDER( IC ) ) THEN
                  MOVEDIST = CHAINSIZE + 1                  
                  WSIZE = MIN( ILAST - J + 1, ( CHAINSIZE + 1 ) * 2 )
               ELSE
                  TOBORDER = NB - MOD( J - 1, NB ) 
                  MOVEDIST = TOBORDER - CHAINSIZE - 1
                  IF ( MOVEDIST .LE. 0 ) THEN
                     MOVEDIST = CHAINSIZE + MOVEDIST + 1
                  END IF
                  WSIZE = MIN( ILAST -J + 1, CHAINSIZE + MOVEDIST + 1 )
               END IF
               
               DIM1 = NB - MOD( J - 1, NB )
               IF ( J + DIM1 - 1 .GE. ILAST ) THEN
                  DIM1 = WSIZE
               END IF
               DIM4 = WSIZE - DIM1
               IF ( NPROW * NPCOL.EQ.1 ) THEN
                  DIM1 = WSIZE
                  DIM4 = 0
               END IF 
               
               IF ( J .LT. ILAST - 1 ) THEN                       
                  CALL PDHGEQZ6( 'R', ILSCHR, ILQ, ILZ, A, DESCA, B,
     $               DESCB, Q, DESCQ, Z, DESCZ, UV( UPOS( IC ), 1 ), UV(
     $               VPOS( IC ), 1 ), LDUV, 1, J, J + WSIZE - 1, N,
     $               WSIZE, ILOQ, IHIQ, ILOZ, IHIZ, LWORK, WORK( IPWORK)
     $               , IERR )
               END IF
            END DO
            
         END DO

*        Left Updates
         DO SWIN = 1, 2
            DO IC = SWIN, NUMWIN, 2
               J = PREVWINPOS( IC )

               RSRC1 = INDXG2P( J, DESCA(NB_), DESCA(RSRC_),
     $            DESCA(RSRC_), NPROW)
               CSRC1 = INDXG2P( J, DESCA(NB_), DESCA(CSRC_),
     $            DESCA(CSRC_), NPCOL)
               
               RSRC2 = RSRC1
               CSRC2 = MOD( CSRC1+1, NPCOL )
               RSRC3 = MOD( RSRC1+1, NPROW )
               CSRC3 = CSRC1
               RSRC4 = RSRC3
               CSRC4 = CSRC2    

               IF ( PREVWINXBORDER( IC ) ) THEN
                  MOVEDIST = CHAINSIZE + 1                  
                  WSIZE = MIN( ILAST - J + 1, ( CHAINSIZE + 1 ) * 2 )
               ELSE
                  TOBORDER = NB - MOD( J - 1, NB ) 
                  MOVEDIST = TOBORDER - CHAINSIZE - 1
                  IF ( MOVEDIST .LE. 0 ) THEN
                     MOVEDIST = CHAINSIZE + MOVEDIST + 1
                  END IF
                  WSIZE = MIN( ILAST - J + 1, CHAINSIZE + MOVEDIST + 1 )
               END IF
               
               DIM1 = NB - MOD( J - 1, NB )
               IF (J + DIM1 - 1 .GE. ILAST) THEN
                  DIM1 = WSIZE
               END IF
               DIM4 = WSIZE - DIM1
               IF ( NPROW * NPCOL .EQ. 1 ) THEN
                  DIM1 = WSIZE
                  DIM4 = 0
               END IF 
               
               IF (J .LT. ILAST-1) THEN                       
                  CALL PDHGEQZ6( 'L', ILSCHR, ILQ, ILZ, A, DESCA, B,
     $               DESCB, Q, DESCQ, Z, DESCZ, UV( UPOS( IC ), 1 ),
     $               UV(VPOS( IC ) ,1 ), LDUV, 1, J, J + WSIZE - 1, N,
     $               WSIZE, ILOQ, IHIQ, ILOZ, IHIZ, LWORK, WORK( IPWORK
     $               ), IERR )
               END IF
            END DO
            
         END DO
         
         UPDTIME( 11 ) = UPDTIME( 11 ) + ( MPI_WTIME() - T1 )   


         GOTO 10

      END IF 
      
      DEALLOCATE( UV, TMPA, TMPB )
      INFO = IERR
      RETURN
*     
*     End of PDHGEQZ5
*     
      END



