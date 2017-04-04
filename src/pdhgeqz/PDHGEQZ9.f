***********************************************************************
*                                                                     *
*     PDHGEQZ9.f:                                                     *
*         Auxillary routine in the package PDHGEQZ.                   *
*                                                                     *
*     Contributors: Bjorn Adlerborn                                   *
*                   Bo Kagstrom                                       *
*                   Daniel Kressner                                   *
*                                                                     *
*     Department of Computing Science and HPC2N, Umea University      *
*     MATHICSE ANCHP, EPF Lausannei                                   *
*                                                                     * 
***********************************************************************
      RECURSIVE SUBROUTINE PDHGEQZ9( ILSCHR, ILQ, ILZ, IFIRST, ILAST, 
     $   ILO, IHI, H, DESCH, T, DESCT, ALPHAR, ALPHAI, BETA, Q, DESCQ, 
     $   Z, DESCZ, ILOQ, IHIQ, ILOZ, IHIZ, GMAXWIN, ND, SELECT, DWORK, 
     $   LWORK, INFO)
      IMPLICIT NONE
*     ..
*     .. Scalar Arguments ..
*     ..
      INTEGER             IHI, ILO, LWORK, ND, GMAXWIN, IFIRST, ILAST, 
     $   INFO, ILOQ, IHIQ, ILOZ, IHIZ
      LOGICAL             ILQ, ILZ, ILSCHR


*     ..
*     .. Array Arguments .. 
*     ..
      DOUBLE PRECISION    ALPHAI( * ), ALPHAR( * ), BETA( * ), H( * ), 
     $   T( * ), Q( * ), Z( * ), DWORK( * )
      INTEGER             DESCH( * ), DESCT( * ), DESCQ( * ), 
     $   DESCZ( * ), SELECT( * )

*     Purpose
*     =======
*     PDHGEQZ9 works on H(IFIRST:ILAST, IFIRST:ILAST) 
*     and T(IFIRST:MYILAST, IFIRST:MYILAST)
*     The aim is to push 0 elements in the diagonal of T to the 
*     top of T while updating H,T,Q, and Z.
*     
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
*     = FALSE: A and B are not referenced.
*     
*     ILQ     (global input) LOGICAL
*     = TRUE: Compute Q.
*     = FALSE: Q is not referenced.
*     
*     ILZ     (global input) LOGICAL
*     = TRUE: Compute Z.
*     = FALSE: Z is not referenced.
*     
*     IFIRST  (global input) INTEGRE
*     Row/column index for where to stop the process.
*     
*     ILAST   (global input) INTEGER
*     Row/column index for where to start the process.
*     
*     ILO     (global input) INTEGER
*     Updates wont be applied below this row/column index.
*     
*     IHI     (global input) INTEGER
*     Updates wont be applied to indexes greater than this.
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
*     
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
*     GMAXWIN (global input) INTEGER
*     Total number of windows to setup and use - have affect on avialble workspace
*     since these are allocated using ALLOCATE
*     
*     ND      (global output) INTEGER
*     The number of converged eigenvalues discovered by this
*     subroutine.
*     
*     DWORK   (local workspace/global output) DOUBLE PRECISION array,
*     dimension (LDWORK)
*     
*     LDWORK  (global input) INTEGER
*     The dimension of the array DWORK.
*     If LDWORK = - 1 then a workspace query is assumed. Optimal workspace
*     is then returned in DWORK(1)
*     
*     SELECT  (global input) INTEGER, array dimension N
*     Placeholder for index of discovered infite eigenvalues. The SELECT 
*     vector is kept updated to refelect new positions.
*     
*     INFO    (global output) INTEGER
*     
*     ===================================
      
*     ..
*     .. Parameters ..
*     ..
      DOUBLE PRECISION    ZERO, ONE
      PARAMETER           ( ZERO = 0.0D+0, ONE = 1.0D+0 )      
      INTEGER            BLOCK_CYCLIC_2D, CSRC_, CTXT_, DLEN_, DT_, 
     $   LLD_, MB_, M_, NB_, N_, RSRC_
      PARAMETER          ( BLOCK_CYCLIC_2D = 1, DLEN_ = 9, DT_ = 1,
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
      LOGICAL             STARTJ, FIRSTBLK, LASTBLK, BLOCKOWNER, LQUERY, 
     $   PARTOWNER
      INTEGER             J, CTX, NB, MYILAST, MYIFIRST, NUMWIN, NPROW, 
     $   NPCOL, MYROW, MYCOL, MAXWIN, LDUV, JJ, CURRWIN, SWIN, BSIZE, 
     $   LWKOPT, MAXWSIZE, MOVECNT, RSRC1, RSRC4, CSRC1, CSRC4, XBORDER, 
     $   UPCNT, DWCNT, NUM0ZEL, I, GMOVECNT, IAM, NPROCS, WINSZ

*     ..  
*     ..  Local Arrays
*     ..
      INTEGER             UPOS( GMAXWIN ), VPOS( GMAXWIN ),
     $   STARTPOS( GMAXWIN ), BLKSZ( GMAXWIN ), MOVECNTS( GMAXWIN )

      DOUBLE PRECISION, ALLOCATABLE :: UV(:,:), TMPH(:,:), TMPT(:,:)       

      
*     ..
*     .. Externals
*     ..

      EXTERNAL            INDXG2P, ICEIL, INDXG2L
      INTEGER             INDXG2P, ICEIL, INDXG2L
      
*     ..
*     .. Intrinsic Functions ..
*     ..
      INTRINSIC           ABS, DBLE, MAX, MIN, SQRT, MOD, LOG
      
      
      
*     
*     Extract current communication context and get grid parameters
*     
      CTX = DESCT(CTXT_)
      CALL BLACS_PINFO(IAM, NPROCS)
      CALL BLACS_GRIDINFO( CTX, NPROW, NPCOL, MYROW, MYCOL )
*     
*     NB_ .EQ. MB_ is assumed, we use NB as blocking factor.
*     
      NB = DESCT(NB_)
      
      INFO = 0
      LWKOPT = 1 
      DWORK ( 1 ) = LWKOPT
      LQUERY = LWORK .EQ. -1
*     
*     Set number for deflated eigenvalues to none initially
*     
      ND = 0
*     
*     Quick return if possible
*     
      IF( IHI.LT.ILO )  RETURN
      IF (ILAST.LE.IFIRST) RETURN
      
*     
*     Set current window size. 
*     Often differs from IHI - ILO + 1!
*     
      WINSZ = ILAST - IFIRST + 1
      
*     
*     Calulate maximum needed work, for dgemm updates only      
*     
      LWKOPT = NB*NB
      IF ( LQUERY ) THEN
         DWORK(1) = LWKOPT
         RETURN
      END IF
      
      
*     
*     Window size is set to be of size that fits on one process
*     
      MAXWSIZE = MAX( MIN( ILAST - IFIRST + 1, NB ) ,3 )  
      
*     
*     Set maxiumn number for windows to use 
*     
      MAXWIN = MAX( 1, MIN( GMAXWIN, NPROW, 
     $   ICEIL( ILAST - IFIRST + 1, NB ) ) )
      
      
*     
*     Storage required for the accumulated transformations
*     
      LDUV = MAXWSIZE*MAXWIN*2
      
      ALLOCATE( UV( LDUV, MAXWSIZE ), STAT = INFO ) 
      IF ( INFO .NE. 0 ) THEN
         INFO = -10
         RETURN
      END IF     

*     
*     Storage for gathered parts of H and T
*     
      ALLOCATE( TMPH( MAXWSIZE, MAXWSIZE ), STAT = INFO )
      IF ( INFO .NE. 0) THEN
         INFO = - 10
         RETURN
      END IF   
      
      ALLOCATE( TMPT( MAXWSIZE, MAXWSIZE ), STAT = INFO )
      IF ( INFO .NE. 0) THEN
         INFO = -10
         RETURN
      END IF 
      
      
      MYILAST = ILAST
      MYIFIRST = IFIRST
      
      

*     
*     Prepare windows      
*     Start with position for U and V, these will not change during the iteration
*     
      DO CURRWIN = 1, MAXWIN 
         IF (CURRWIN.EQ.1) THEN
            UPOS(CURRWIN) = 1
            VPOS(CURRWIN) = UPOS(CURRWIN) + MAXWSIZE
         ELSE
            UPOS(CURRWIN) = VPOS(CURRWIN-1) + MAXWSIZE
            VPOS(CURRWIN) = UPOS(CURRWIN) + MAXWSIZE
         END IF
      END DO

*     
*     Blocksizes and end position. These might change during the iteration
*     
      J = MYIFIRST
      NUMWIN = 0
      FIRSTBLK = .TRUE.
      LASTBLK = .FALSE.
      DO WHILE ( J .LT. MYILAST .AND. NUMWIN .LT. MAXWIN)
         NUMWIN = NUMWIN + 1
*        
*        Number of rows and column in current block is BSIZE          
*        
         BSIZE = NB - MOD( J - 1, NB )

*        
*        Check for potential last block, adjust block size         
*        
         IF ( J + BSIZE - 1 .GE. MYILAST ) THEN             
            BSIZE = MYILAST - J + 1   
            LASTBLK = .TRUE.
         END IF
         
         STARTPOS( NUMWIN ) = J
         BLKSZ( NUMWIN ) = BSIZE
         MOVECNTS( NUMWIN ) = 0
         J = J + BSIZE
      END DO
      
      
 10   CONTINUE
      GMOVECNT = 0
      

*     
*     Repeat at least 2 times. One for none crossborder (XBORDER=0) and 
*     one for crossborder (XBORDER=1)
*     
      DO XBORDER = 0, 1
 20      CONTINUE    
         
*        
*        Each window either resides on 1 process, or is shared by 4 processes.
*        
*        Extract processes for current block as below,
*        
*        .. 
*        (RSRC1, CSRC1) | (RSCR1, CSRC4)     
*        ---------------| --------------
*        (RSCR4, CSRC1) | (RSCR4, CSRC4) 
*        
*        RSRC1 = RSCR4 and CSRC1 = CSRC4 if window is not shared
*        
*        Since the SELECT values are local, (RSRC1, CSRC1) and (RSRC4, CSRC4) 
*        first exhange values so that both are aware of zeros in the diagonal 
*        of T. (RSCR1, CSRC4) and (RSCR4, CSRC1) also recive information about
*        the number of zeros in the diagonal of T.      
*        
*        If there are zeros in the diagonal of T and the window is shared, scatter 
*        the block so that both (RSRC1, CSRC1) and (RSRC4, CSRC4) will have their 
*        own copy of the full window, stored in TMPH and TMPT.
*        TMPA and TMPB are then processed by (RSRC1, CSRC1) and (RSRC4, CSRC4) while
*        updating SELECT and keeping track of how many zeros that are moved.
*        
*        (H,T) is updated w.r.t TMPH and TMPT after the local updates.
*        

*        
*        Do odd windows first, then even. This will enable us to deal with several
*        windows at the same time if we run in parallel.
*        
         DO SWIN = 1, 2
            DO CURRWIN = SWIN, NUMWIN, 2          
               JJ = STARTPOS( CURRWIN ) 
               BSIZE = BLKSZ( CURRWIN ) 
               
               RSRC1 = INDXG2P( JJ, NB, 
     $            0, DESCH( RSRC_ ), NPROW)
               CSRC1 = INDXG2P( JJ, NB, 
     $            0, DESCH( CSRC_ ), NPCOL)               
               CSRC4 = INDXG2P( JJ + BSIZE - 1, NB, 
     $            0, DESCH( RSRC_ ), NPROW)
               RSRC4 = INDXG2P( JJ + BSIZE - 1, NB, 
     $            0, DESCH( CSRC_ ), NPCOL)
               
               PARTOWNER = (MYROW.EQ.RSRC1.AND.MYCOL.EQ.CSRC1) 
     $            .OR. (MYROW.EQ.RSRC1.AND.MYCOL.EQ.CSRC4) .OR.
     $            (MYROW.EQ.RSRC4.AND.MYCOL.EQ.CSRC1) .OR.
     $            (MYROW.EQ.RSRC4.AND.MYCOL.EQ.CSRC4)
               
               BLOCKOWNER = (MYROW.EQ.RSRC1.AND.MYCOL.EQ.CSRC1) 
     $            .OR. (MYROW.EQ.RSRC4.AND.MYCOL.EQ.CSRC4)
*              
*              If two diagonal processors are invloved, exchange select values
*              
               IF ( BLOCKOWNER .AND. 
     $            ( RSRC1.NE.RSRC4.OR.CSRC1.NE.CSRC4 ) ) THEN
                  
                  
                  DWCNT = 1 + MOD( JJ + BSIZE - 1 - 1, NB )
                  UPCNT = BSIZE - DWCNT 
                  
                  
                  
                  IF ( DWCNT .LT. BSIZE ) THEN
                     
                     
                     IF ( MYROW .EQ . RSRC1 .AND. 
     $                  MYCOL .EQ. CSRC1 ) THEN
*                       
*                       DWCNT holds number of elements to receive from process below
*                       UPCNT holds number of elements to send to process below
*                       
                        
                        CALL IGESD2D( CTX, UPCNT, 1,
     $                     SELECT( JJ ), 1, 
     $                     RSRC4, CSRC4 )
                        
                        CALL IGERV2D( CTX, DWCNT, 1,
     $                     SELECT( JJ + UPCNT ), 1, 
     $                     RSRC4, CSRC4 )
                        
                        
                     ELSE IF ( MYROW .EQ . RSRC4 .AND. 
     $                     MYCOL .EQ. CSRC4 ) THEN
*                       
*                       UPCNT holds number of elements to receive from process above
*                       DWCNT holds number of elements to send to process above
*                       
                        
                        CALL IGESD2D( CTX, DWCNT, 1,
     $                     SELECT( JJ + UPCNT ), 1, 
     $                     RSRC1, CSRC1 )
                        
                        CALL IGERV2D( CTX, UPCNT, 1,
     $                     SELECT( JJ ), 1, 
     $                     RSRC1, CSRC1 )
                        
                     END IF
                     
                  END IF
               END IF    
               
               IF ( BLOCKOWNER ) THEN
*                 
*                 Check the number of non-zeros in the select vector
*                 not counting the first and last element since those potential zeros in T
*                 cant be moved right now. 
*                 
                  NUM0ZEL = 0
                  
                  DO I = JJ + 1, JJ + BSIZE - 2 
                     IF ( SELECT (I) .EQ. 1 ) THEN
                        NUM0ZEL = NUM0ZEL + 1
                     END IF
                  END DO                          
               END IF

*              
*              Let nearest neightbors take part of the select info
*              for being able to determine wether to copy data or not.
*              
               IF ( PARTOWNER .AND. 
     $            ( RSRC1 .NE. RSRC4 .OR. CSRC1 .NE. CSRC4 ) ) THEN
*                 
*                 Send RSRC1, RSRC1 --> RSRC1, RSRC4
*                 
                  IF ( MYROW .EQ . RSRC1 .AND. 
     $               MYCOL .EQ. CSRC1 ) THEN
                     
                     CALL IGESD2D( CTX, 1, 1,
     $                  NUM0ZEL, 1, 
     $                  RSRC1, CSRC4 )
*                    
*                    Recv RRSRC1, RSRC4 <-- SRC1, RSRC1
*                    
                  ELSE IF ( MYROW .EQ . RSRC1 .AND. 
     $                  MYCOL .EQ. CSRC4 ) THEN

                     CALL IGERV2D( CTX, 1, 1,
     $                  NUM0ZEL, 1, 
     $                  RSRC1, CSRC1 )

                  END IF
*                 
*                 Send RSRC4, RSRC4 --> RSRC4, RSRC1
*                 
                  IF ( MYROW .EQ . RSRC4 .AND. 
     $               MYCOL .EQ. CSRC4 ) THEN
                     
                     CALL IGESD2D( CTX, 1, 1,
     $                  NUM0ZEL, 1, 
     $                  RSRC4, CSRC1 )
*                    
*                    Recv RSRC4, RSRC1 <-- RSRC4, RSRC4
*                    
                  ELSE IF ( MYROW .EQ . RSRC4 .AND. 
     $                  MYCOL .EQ. CSRC1 ) THEN

                     CALL IGERV2D( CTX, 1, 1,
     $                  NUM0ZEL, 1, 
     $                  RSRC4, CSRC4 )

                  END IF  
               END IF
               
*              
*              Gather window to TMPH and TMPT
*              
               IF ( PARTOWNER .AND. NUM0ZEL .GT. 0 ) THEN
                  CALL PDLACP4( JJ, JJ + BSIZE - 1, 
     $               BSIZE,  H, T, DESCT,
     $               TMPH, TMPT, MAXWSIZE,      
     $               0, 0 )    
               END IF
               
               
               
               
               IF ( BLOCKOWNER .AND. NUM0ZEL .GT. 0  ) THEN
*                 
*                 Set number of moved zeros in T to none                          
*                 DHGEQZ7 will update this number if any zeros are moved
*                 
                  MOVECNT = 0                              
                  STARTJ = JJ .EQ. MYIFIRST
                  CALL DHGEQZ7 ( 'U', BSIZE, BSIZE, 
     $               TMPH,  MAXWSIZE, TMPT, MAXWSIZE, 
     $               UV( UPOS( CURRWIN ), 1 ), 
     $               UV( VPOS( CURRWIN ), 1 ),
     $               LDUV, SELECT( JJ ), 
     $               STARTJ, MOVECNT )
*                 
*                 Save number of moved zeros in T 
*                 
                  MOVECNTS( CURRWIN ) = MOVECNT            
               END IF
               
*              
*              Restore (H,T) w.r.t TMPH and TMPT
*              
               IF ( PARTOWNER .AND. NUM0ZEL .GT. 0 ) THEN
                  CALL PDLACP4( JJ, JJ + BSIZE - 1, 
     $               BSIZE,  H, T, DESCT,
     $               TMPH, TMPT, MAXWSIZE,      
     $               0, 1 )
               END IF
               
               
            END DO
         END DO
         
         
*        
*        Everyone investigate MOVECNTS to update our global move count (GMOVECNT).
*        If no zeros has been moved during XBORDER = 0..1, we are done!
*        Note that GMOVECNT is set to 0 as the loop for XBORDER = 1 begins and all 
*        windows initially have MOVECNT = 0
*        
         CALL IGAMX2D( CTX, 'A', ' ', NUMWIN, 1, 
     $      MOVECNTS, 1,  
     $      -1, -1, -1, -1, -1 )                
         MOVECNT = 0
         DO I = 1, NUMWIN
            MOVECNT = MOVECNT + MOVECNTS( I )                 
         END DO
         GMOVECNT = GMOVECNT + MOVECNT


*        
*        The updates are stored on both (RSRC1, CSRC1) and (RSRC4, CSRC4). 
*        These need to be distributed column- and row-wise so that updates  
*        can be applied on the of-digaonal blocks (using dgemm).
*        
         
*        
*        Broadcast V columnwise. If ILQ also U
*        
         DO SWIN = 1, 2
            DO CURRWIN = SWIN, NUMWIN, 2
*              
*              Check if something to do 
*              
               MOVECNT = MOVECNTS( CURRWIN ) 
               IF ( MOVECNT .EQ. 0 ) GOTO 100
               
               
               JJ = STARTPOS( CURRWIN )
               BSIZE = BLKSZ( CURRWIN )
               
               
               RSRC1 = INDXG2P( JJ, NB, 
     $            0, DESCH( RSRC_ ), NPROW)
               CSRC1 = INDXG2P( JJ, NB, 
     $            0, DESCH( CSRC_ ), NPCOL)
               
               CSRC4 = INDXG2P( JJ + BSIZE - 1, NB, 
     $            0, DESCH( RSRC_ ), NPROW)
               RSRC4 = INDXG2P( JJ + BSIZE - 1, NB, 
     $            0, DESCH( CSRC_ ), NPCOL)
               
               IF ( MYROW .EQ. RSRC1 .AND. MYCOL .EQ. CSRC1 ) 
     $            THEN
                  IF ( ILQ ) THEN
                     CALL DGEBS2D( CTX, 'C', ' ', BSIZE, BSIZE,
     $                  UV( UPOS( CURRWIN ), 1 ), LDUV )
                     CALL DGEBS2D( CTX, 'C', ' ', BSIZE, BSIZE,
     $                  UV( VPOS( CURRWIN ), 1 ), LDUV )
                  ELSE 
                     CALL DGEBS2D( CTX, 'C', ' ', BSIZE, BSIZE,
     $                  UV( VPOS( CURRWIN ), 1 ), LDUV )
                  END IF
               ELSE IF (MYCOL .EQ. CSRC1) THEN
                  IF( ILQ ) THEN
                     CALL DGEBR2D( CTX, 'C', ' ', BSIZE, BSIZE,
     $                  UV( UPOS( CURRWIN ), 1 ), LDUV, 
     $                  RSRC1, CSRC1 )
                     CALL DGEBR2D( CTX, 'C', ' ', BSIZE, BSIZE,
     $                  UV( VPOS( CURRWIN ), 1 ), LDUV, 
     $                  RSRC1, CSRC1 )                  
                  ELSE
                     CALL DGEBR2D( CTX, 'C', ' ', BSIZE, BSIZE,
     $                  UV( VPOS( CURRWIN ), 1 ), LDUV, 
     $                  RSRC1, CSRC1 )
                  END IF                  
               END IF
               
               IF ( RSRC1 .EQ. RSRC4 .AND. CSRC1 .EQ. CSRC4 ) 
     $            GOTO 100
               
               IF ( MYROW .EQ. RSRC4 .AND. MYCOL .EQ. CSRC4 ) 
     $            THEN
                  IF ( ILQ ) THEN
                     CALL DGEBS2D( CTX, 'C', ' ', BSIZE, BSIZE,
     $                  UV( UPOS( CURRWIN ), 1 ), LDUV )
                     CALL DGEBS2D( CTX, 'C', ' ', BSIZE, BSIZE,
     $                  UV( VPOS( CURRWIN ),1 ), LDUV )
                  ELSE 
                     CALL DGEBS2D( CTX, 'C', ' ', BSIZE, BSIZE,
     $                  UV( VPOS( CURRWIN ),1 ), LDUV )
                  END IF                  
               ELSE IF ( MYCOL .EQ. CSRC4 ) THEN
                  IF ( ILQ ) THEN
                     CALL DGEBR2D( CTX, 'C', ' ', BSIZE, BSIZE,
     $                  UV( UPOS( CURRWIN ), 1 ), LDUV, 
     $                  RSRC4, CSRC4 )
                     CALL DGEBR2D( CTX, 'C', ' ', BSIZE, BSIZE,
     $                  UV( VPOS( CURRWIN ),1 ), LDUV, 
     $                  RSRC4, CSRC4 )
                  ELSE 
                     CALL DGEBR2D( CTX, 'C', ' ', BSIZE, BSIZE,
     $                  UV( VPOS( CURRWIN ),1 ), LDUV, 
     $                  RSRC4, CSRC4 )
                  END IF                   
               END IF  
 100           CONTINUE
               
               
            END DO
         END DO 
*        
*        Broadcast U rowwise
*        
         DO SWIN = 1, 2
            DO CURRWIN = SWIN, NUMWIN, 2
               
*              
*              Check if something to do 
*              
               MOVECNT = MOVECNTS( CURRWIN ) 
               IF ( MOVECNT .EQ. 0 ) GOTO 200
               
               
               JJ = STARTPOS( CURRWIN )
               BSIZE = BLKSZ( CURRWIN )
               
               RSRC1 = INDXG2P( JJ, NB, 
     $            0, DESCH( RSRC_ ), NPROW)
               CSRC1 = INDXG2P( JJ, NB, 
     $            0, DESCH( CSRC_ ), NPCOL)
               
               CSRC4 = INDXG2P( JJ + BSIZE - 1, NB, 
     $            0, DESCH( RSRC_ ), NPROW)
               RSRC4 = INDXG2P( JJ + BSIZE - 1, NB, 
     $            0, DESCH( CSRC_ ), NPCOL)
               
               
               
               IF ( MYROW .EQ. RSRC1 .AND. MYCOL .EQ. CSRC1 ) 
     $            THEN
                  CALL DGEBS2D( CTX, 'R', ' ', BSIZE, BSIZE, 
     $               UV( UPOS( CURRWIN ), 1 ), LDUV )  
               ELSE IF ( MYROW .EQ. RSRC1 ) THEN
                  CALL DGEBR2D( CTX, 'R', ' ', BSIZE, BSIZE, 
     $               UV( UPOS( CURRWIN ), 1 ), LDUV, 
     $               RSRC1, CSRC1 )
               END IF              
               
               
               IF ( RSRC1 .EQ. RSRC4 .AND. CSRC1 .EQ. CSRC4 ) 
     $            GOTO 200
               
               IF ( MYROW .EQ. RSRC4 .AND. MYCOL .EQ. CSRC4 ) 
     $            THEN
                  CALL DGEBS2D( CTX, 'R', ' ', BSIZE, BSIZE, 
     $               UV( UPOS( CURRWIN ), 1 ), LDUV )   
               ELSE IF ( MYROW .EQ. RSRC4 ) THEN
                  CALL DGEBR2D( CTX, 'R', ' ', BSIZE, BSIZE, 
     $               UV( UPOS( CURRWIN ), 1 ), LDUV, 
     $               RSRC4, CSRC4 ) 
               END IF
               
               
 200           CONTINUE
               
            END DO
         END DO 
         
*        
*        Update remaning parts of A, B, Q and Z using U and V
*        First from right              
*        
         DO SWIN = 1, 2
            DO CURRWIN = SWIN, NUMWIN, 2
*              
*              Check if something to do 
*              
               MOVECNT = MOVECNTS( CURRWIN ) 
               IF ( MOVECNT .EQ. 0 ) GOTO 300
               
               JJ = STARTPOS( CURRWIN )
               BSIZE = BLKSZ( CURRWIN )
               
               CALL PDHGEQZ6( 'R', ILSCHR, ILQ, ILZ, 
     $            H, DESCH, T, DESCT, Q, DESCQ, Z, DESCZ, 
     $            UV( UPOS( CURRWIN ), 1 ), 
     $            UV( VPOS( CURRWIN ), 1 ), LDUV, 
     $            1, JJ, JJ + BSIZE - 1, IHI - ILO + 1, 
     $            BSIZE, ILOQ, IHIQ, ILOZ, IHIZ, 
     $            LWORK, DWORK, INFO )
 300           CONTINUE                      
            END DO 
         END DO
*        
*        Now from left           
*        
         DO SWIN = 1, 2
            DO CURRWIN = SWIN, NUMWIN, 2
*              
*              Check if something to do 
*              
               MOVECNT = MOVECNTS( CURRWIN ) 
               IF ( MOVECNT .EQ. 0 ) GOTO 400
               
               JJ = STARTPOS( CURRWIN )
               BSIZE = BLKSZ( CURRWIN )
               
               CALL PDHGEQZ6( 'L', ILSCHR, ILO, ILQ, 
     $            H, DESCH, T, DESCT, Q, DESCQ, Z, DESCZ, 
     $            UV( UPOS( CURRWIN ), 1 ), 
     $            UV( VPOS( CURRWIN ), 1 ), LDUV, 
     $            1, JJ, JJ + BSIZE - 1, IHI - ILO + 1, 
     $            BSIZE, ILOQ, IHIQ, ILOZ, IHIZ, 
     $            LWORK, DWORK, INFO )
 400           CONTINUE     
            END DO 
         END DO
         
*        
*        Check for deflations in first block
*        and adjust problem size accordingly 
*        If MYIFIRST has changed, recalculate maximum number of windows
*        
         IF ( STARTPOS( 1 ) .EQ. MYIFIRST .AND. 
     $      MOVECNTS( 1 ) .GT. 0 ) THEN
*           
*           Set the BETA to zero to indicate an deflated infinite eigenvalue
*           
            MOVECNT = MOVECNTS( 1 )
            DO I = MYIFIRST, MYIFIRST + MOVECNT - 1
               BETA( I ) = ZERO
               ALPHAI( I ) = ZERO
               CALL PDELGET( 'A', ' ', ALPHAR( I ),
     $            H, I,  I, DESCH)

               ND = ND + 1
            END DO
            
            MYIFIRST = MYIFIRST + MOVECNT
*           
*           Check if done                      
*           
            IF ( MYIFIRST .GT. MYILAST ) GOTO 500
            
            MAXWIN = MAX( 1, MIN( GMAXWIN, NPROW, 
     $         ICEIL( MYILAST - MYIFIRST + 1, NB ) ) )
            
         END IF           
         

         
*        
*        Recalculate windows and repeat until we have done all blocks
*        and no zeors are moved.
*        
         IF ( .NOT. LASTBLK ) THEN
            

*           
*           Prepare new windows
*           
            
            NUMWIN = 0
            DO WHILE ( J .LE. MYILAST .AND. NUMWIN .LT. MAXWIN )
               NUMWIN = NUMWIN + 1
*              
*              Number of rows and column in current block is BSIZE          
*              
               IF ( XBORDER .EQ. 0 ) THEN
                  BSIZE = NB - MOD( J - 1, NB )
               ELSE
                  BSIZE = NB
               END IF
*              
*              Check for potential last block, adjust block size         
*              
               IF ( J + BSIZE - 1 .GE. MYILAST ) THEN             
                  BSIZE = MYILAST - J + 1   
                  LASTBLK = .TRUE.
               END IF

               MOVECNTS( NUMWIN ) = 0
               STARTPOS( NUMWIN ) = J
               BLKSZ( NUMWIN ) = BSIZE

               J = J + BSIZE

               

               
            END DO
            GOTO 20
            
         ELSE
*           
*           Check if done
*           
            IF ( XBORDER .EQ. 1 .AND. GMOVECNT .EQ. 0 ) GOTO 500

            LASTBLK = .FALSE.            
            J = MYIFIRST
*           
*           Prepare new windows
*           

*           
*           First block is of size BSIZE
*           
            BSIZE = NB - MOD ( J - 1, NB )

*           
*           If next loop is a crossborder loop:
*           
            IF ( XBORDER .EQ. 0 ) THEN
*              
*              Find middle of first+1 block 
*              
               J = J + BSIZE
               J = MIN( J + NB / 2 - 1, MYILAST )
*              
*              Try to subtract NB to new STARTPOS I                  
*              
               I = MAX( J - NB + 1, MYIFIRST )
               BSIZE = J - I + 1
               J = I
            END IF

            NUMWIN = 0
            DO WHILE ( J .LE. MYILAST .AND. NUMWIN .LT. MAXWIN )
               NUMWIN = NUMWIN + 1
*              
*              Check for potential last block and adjust block size
*              
               IF ( J + BSIZE - 1 .GE. MYILAST ) THEN             
                  BSIZE = MYILAST - J + 1   
                  LASTBLK = .TRUE.
               END IF
               MOVECNTS( NUMWIN ) = 0
               STARTPOS( NUMWIN ) = J
               BLKSZ( NUMWIN ) = BSIZE
               
               J = J + BSIZE
*              
*              Number of rows and column in next block is BSIZE          
*              
               IF ( XBORDER .EQ. 0 ) THEN
                  BSIZE = NB
               ELSE
                  BSIZE = NB - MOD ( J - 1, NB )                  
               END IF 
               

            END DO   
*           
*           Restart outer loop if we just performed crossborder loop
*           
            IF ( XBORDER .EQ. 1 ) GOTO 10

         END IF
      END DO    
      
      
 500  CONTINUE          
      
      
      IFIRST = MYIFIRST
      
      DEALLOCATE( UV )
      DEALLOCATE( TMPH )
      DEALLOCATE( TMPT )
      
      RETURN
*     
*     End of PDHGEQZ7
*     
      END
