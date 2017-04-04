***********************************************************************
*                                                                     *
*     PDHGEQZ7.f:                                                     *
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
      SUBROUTINE PDHGEQZ7( ILSCHR, ILQ, ILZ, IFIRST, ILAST, ILO, IHI, 
     $   H, DESCH, T, DESCT, ALPHAR, ALPHAI, BETA, Q, DESCQ, Z, DESCZ,
     $   ILOQ, IHIQ, ILOZ, IHIZ, GMAXWIN, ND, DWORK, LWORK, SELECT, 
     $   INFO )
      IMPLICIT NONE   
*     ..
*     .. Scalar Arguments ..
*     ..
      INTEGER             IHI, ILO, LWORK, ND, ILOQ, IHIQ, ILOZ, IHIZ,
     $   GMAXWIN, IFIRST, ILAST, INFO 
      LOGICAL             ILQ, ILZ, ILSCHR
      

*     ..
*     .. Array Arguments .. 
*     ..
      DOUBLE PRECISION    ALPHAI( * ), ALPHAR( * ), BETA( * ), H( * ), 
     $   T( * ), Q( * ), Z( * ), DWORK( * )
      INTEGER             DESCH( * ), DESCT( * ), DESCQ( * ), 
     $   DESCZ( * ), SELECT( * )

*  Purpose
*  =======
*  PDHGEQZ7 works on H(IFIRST:ILAST, IFIRST:ILAST) 
*  and T(IFIRST:MYILAST, IFIRST:MYILAST)
*  The aim is to push 0 elements in the diagonal of T to the 
*  bottom or top f T, while keeping H, T, Q and Z updated.
*  First it inspects the bottom and top of the T matrix, and deflate
*  occurences of infinte eigenvalues. Afterward, it calls PDHGEQZ8 
*  to push zeros on diagonal of T downwards, and PDHGEQ9 to push
*  zeros upwards. 
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
*
*  ILSCHR  (global input) LOGICAL
*          = TRUE: Apply all transformations to A and B. 
*          = FALSE: A and B are not referenced.
*      
*  ILQ     (input) LOGICAL
*          = TRUE: Compute Q.
*          = FALSE: Q is not referenced.
*      
*  ILZ     (global input) LOGICAL
*          = TRUE: Compute Z.
*          = FALSE: Z is not referenced.
*
*  IFIRST  (global input) INTEGRE
*          Row/column index for where to start the process.
*           
*  ILAST   (global input) INTEGER 
*          Row/column index for where to stop the process. 
*
*  ILO     (global input) INTEGER
*          Updates wont be applied below this row/column index.
*
*  IHI     (global input) INTEGER
*          Updates wont be applied to indexes greater than this.
*
*  H       (local input/output) DOUBLE PRECISION array, dimension (LLD_H, LOCc(N)).
*          On entry, H contains the upper Hessenberg matrix H.
*          On exit, if ILSCHR=TRUE, H is quasi-upper triangular in rows and columns 
*          (ILO:IHI), with 1 x 1 and 2 x 2 blocks on the diagonal where 
*          the 2 x 2 blocks corresponds to complex conjugated pairs of 
*          eigenvalues. If ILSCHR=FALSE, H is unspecified on exit.      
*
*  DESCH   (global and local input) INTEGER array of dimension DLEN_.
*          The array descriptor for the distributed matrix H. 
*      
*  T       (local input/output) DOUBLE PRECISION array, dimension (LLD_T, LOCc(N)).
*          On entry, the N-by-N upper triangular matrix T. 
*          On exit, if ILSHUR the upper triangular matrix the updated tringular 
*          matrix T.
*      
*  DESCT   (global and local input) INTEGER array of dimension DLEN_.
*          The array descriptor for the distributed matrix T. 
*      
*  ALPHAR  (global output) DOUBLE PRECISION array, dimension (N)
*          ALPHAR(1:N) will be set to real parts of the diagonal 
*          elements of H that would result from reducing H and T to
*          Schur form and then further reducing them both to triangular
*          form using unitary transformations s.t. the diagonal of T
*          was non-negative real.  Thus, if H(j,j) is in a 1-by-1 block
*          (i.e., H(j+1,j)=H(j,j+1)=0), then ALPHAR(j)=H(j,j). 
*          Note that the (real or complex) values 
*          (ALPHAR(j) + i*ALPHAI(j))/BETA(j), j=1,...,N, are the 
*          generalized eigenvalues of the matrix pencil H - wT. 
*      
*  ALPHAI  (global output) DOUBLE PRECISION array, dimension (N)
*          ALPHAI(1:N) will be set to imaginary parts of the diagonal
*          elements of H that would result from reducing H and T to
*          Schur form and then further reducing them both to triangular
*          form using unitary transformations s.t. the diagonal of T
*          was non-negative real.  Thus, if H(j,j) is in a 1-by-1 block
*          (i.e., H(j+1,j)=J(j,j+1)=0), then ALPHAI(j)=0. 
*          Note that the (real or complex) values 
*          (ALPHAR(j) + i*ALPHAI(j))/BETA(j), j=1,...,N, are the 
*          generalized eigenvalues of the matrix pencil H - wT. 
*      
*  BETA    (global output) DOUBLE PRECISION array, dimension (N)
*          BETA(1:N) will be set to the (real) diagonal elements of T
*          that would result from reducing H and T to Schur form and
*          then further reducing them both to triangular form using
*          unitary transformations s.t. the diagonal of T was 
*          non-negative real.  Thus, if H(j,j) is in a 1-by-1 block
*          (i.e., H(j+1,j)=H(j,j+1)=0), then BETA(j)=T(j,j). 
*          Note that the (real or complex) values 
*          (ALPHAR(j) + i*ALPHAI(j))/BETA(j), j=1,...,N, are the 
*          generalized eigenvalues of the matrix pencil H - wT. 
*          (Note that BETA(1:N) will always be non-negative, and no
*          BETAI is necessary.)
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
*  GMAXWIN (global input) INTEGER
*          Total number of windows to setup and use - have affect on avialble workspace
*          since these are allocated using ALLOCATE
*
*  ND      (global output) INTEGER
*          The number of converged eigenvalues discovered by this
*          subroutine.
*     
*  DWORK   (local workspace/global output) DOUBLE PRECISION array,
*          dimension (LDWORK)
*     
*  LDWORK  (global input) INTEGER
*          The dimension of the array DWORK.
*          If LDWORK = - 1 then a workspace query is assumed. Optimal workspace
*          is then returned in DWORK(1)
*
*  SELECT  (global input) INTEGER, array dimension N
*          Placeholder for index of discovered infite eigenvalues
*
*  INFO    (global output) INTEGER
*          = 0 : Succesful call
*
*  =========================================

      
      
 

      
*     ..
*     .. Parameters ..
*     ..
      DOUBLE PRECISION    ZERO, ONE
      PARAMETER           ( ZERO = 0.0D+0, ONE = 1.0D+0 )      
      INTEGER            BLOCK_CYCLIC_2D, CSRC_, CTXT_, DLEN_, DT_, 
     $     LLD_, MB_, M_, NB_, N_, RSRC_
      PARAMETER          ( BLOCK_CYCLIC_2D = 1, DLEN_ = 9, DT_ = 1,
     $     CTXT_ = 2, M_ = 3, N_ = 4, MB_ = 5, NB_ = 6,
     $     RSRC_ = 7, CSRC_ = 8, LLD_ = 9 )
      
*     ..
*     .. Precision paramters ..
*     ..
      DOUBLE PRECISION    SMLNUM, ULP, TFNORM, T1NORM, HFNORM, HTOL,
     $                    TTOL
      COMMON /PREC/       SMLNUM, ULP, TFNORM, T1NORM, HFNORM, HTOL, 
     $                    TTOL 

*     ..
*     .. Local Scalars ..
*     ..
      INTEGER             J, CTX, NB, MYILAST, MYIFIRST, N, NPROW, 
     $   NPCOL, MYROW, MYCOL, IAM, NPROCS, NDTOP, NDBOT, BSIZE,
     $   LWKOPT, NDUP, NDDW, STARTI, ENDI, RSRC1, CSRC1, NUM0ZEL, 
     $   I, LIR, LIC, NH, LEFT, RIGHT, UP, DOWN, WINSZ, MID

      LOGICAL             LQUERY, TEST, INFFOUND
      DOUBLE PRECISION    S, C, TEMP, LVAL
      

*     ..  
*     ..  Local Arrays
*     ..

      DOUBLE PRECISION TMPH( 2, 2 )       

*     ..
*     .. Externals
*     ..
      EXTERNAL            INDXG2P, ICEIL, INDXG2L, NUMROC
      INTEGER             INDXG2P, ICEIL, INDXG2L, NUMROC
      
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
*     Set current window sizes. 
*      
      WINSZ = ILAST - IFIRST + 1
      N = IHI - ILO + 1
*     
*     Calulate maximum needed work, for dgemm updates only      
*     
      LWKOPT = NB*NB
      IF ( LQUERY ) THEN
         DWORK(1) = LWKOPT
         RETURN
      END IF
      
      

      
*
*     Get grid column indices for neighbours.
*     
      LEFT = MODULO( MYCOL - 1, NPCOL )
      RIGHT = MODULO( MYCOL + 1, NPCOL )      
      UP = MODULO( MYROW - 1, NPROW )
      DOWN = MODULO( MYROW + 1, NPROW )
      
      
      
      
*
*     Set the number of found zeros in the diagonal of B to none.
*     
      SELECT( 1 : IHI ) = 0
      
*              
*     Check (locally) for zeros in the diagonal of B.       
*                 
      NUM0ZEL = 0
      I = ILAST
      DO WHILE ( I .GE. IFIRST )
          BSIZE = 1 + MOD( I - 1, NB )
*
*         If last block, recalc the block size.
*          
          IF (I - BSIZE + 1 .LE. IFIRST ) THEN
              BSIZE = I - IFIRST + 1
          END IF
          
          RSRC1 = INDXG2P( I, NB, 0, DESCT( RSRC_ ), NPROW )
          CSRC1 = INDXG2P( I, NB, 0, DESCT( CSRC_ ), NPCOL )

*
*         Work on locally on block T( I - BSIZE + 1 : I, I - BSIZE + 1 : I )
*
          IF ( MYROW .EQ. RSRC1 .AND. MYCOL .EQ. CSRC1 ) THEN
                            
              LIR = INDXG2L( I, NB, 0, 0, NPROW )
              LIC = INDXG2L( I, NB, 0, 0, NPCOL )
                      
              DO J = 0, BSIZE - 1
                  LVAL = T( LIR - J + ( LIC - J - 1 ) * DESCT( LLD_ ) ) 

                  TEST = ABS( LVAL ).LE.MAX( TTOL , SMLNUM ) 
                  IF ( TEST ) THEN
                     T( LIR - J + ( LIC - J - 1 ) * DESCT( LLD_ ) ) = 
     $                    ZERO
                     SELECT( I - J ) = 1
                     NUM0ZEL = NUM0ZEL + 1                      
                  END IF                      
              END DO
          END IF
*         
*         Goto next block
*      
          I = I - BSIZE
      END DO
      
      
*     
*     All check if infite eigenvalues exists
*              
      CALL IGSUM2D( CTX, 'A', ' ', 1, 1, 
     $    NUM0ZEL, 1, -1, -1 )  
            
*
*     Return if nothing to do
*              
      IF ( NUM0ZEL .EQ. 0 ) THEN
         RETURN
      END IF
      MYILAST = ILAST
      MYIFIRST = IFIRST
      
      
*
*     Let everyone have the same startinformation 
*     about where zeros are. This is used to deflate
*     at top and bottom.
*          
*     
      
      CALL IGAMX2D( CTX, 'A', ' ', IHI, 1, 
     $    SELECT, 1, -1, -1, -1, -1, -1 ) 

      NDBOT = 0
*     
*     Precheck 1
*     Check B( MYILAST, MYILAST ) 
*     
      I = MYILAST
      DO WHILE ( I .GT. MYIFIRST ) 
         IF ( SELECT( I ) .NE. 1) EXIT
*        
*        Construct a column rotation that annihilates H( I, I - 1 )     
*        and updates H( I, I )
         CALL PDLACP3( 2, I - 1, H, DESCH, TMPH, 2, -1, -1, 0 )
         TEMP = TMPH( 2, 2 )                  
         CALL DLARTG( TEMP, TMPH( 2, 1 ), C, S, TMPH( 2, 2 ) )
         TMPH( 2, 1 ) = ZERO                  
*        
*        Apply the column rotation to H, T, and Z. This won't 
*        cause any fill-in in T (since T(MYILAST, MYILAST) is zero).
*        
         CALL PDROT( I - ILO, H, ILO, I, DESCH, 1,
     $      H, ILO, I - 1, DESCH, 1, C, S, DWORK, LWORK, INFO )
         
         CALL PDROT( I - ILO, T, ILO, I, DESCT, 1,
     $      T, ILO, I - 1, DESCT, 1, C, S, DWORK, LWORK, INFO )
         
         IF (ILZ) THEN
            CALL PDROT( IHIZ - ILOZ + 1, Z, ILOZ, I, DESCZ, 1, 
     $         Z, ILOZ, I - 1, DESCZ, 1, C, S, DWORK, LWORK, INFO )
         END IF
*        
*        Update the 2 the global elements affected by the local DLARTG
*        
         CALL PDELSET( H, I, I - 1, DESCH, ZERO )
         CALL PDELSET( H, I, I, DESCH, TMPH( 2, 2 ) )
         ALPHAR( I ) = TMPH( 2, 2 )
         ALPHAI( I ) = ZERO
*
*        It's enough to set the BETA to zero to indicate an
*        infinite eigenvalue
*
         BETA( I ) = ZERO
         CALL PDELGET('A', ' ', TEMP, T, I - 1, I - 1, DESCT)
         IF ( ABS(TEMP) .GT. MAX( TTOL , SMLNUM )) THEN
            SELECT( I - 1 ) = 0
         ELSE
            SELECT( I - 1 ) = 1
            CALL PDELSET( T, I - 1, I - 1, DESCT, ZERO )
         END IF

         MYILAST = MYILAST - 1
         I = I - 1       
         NDBOT = NDBOT + 1
      END DO

*     
*     Precheck 2
*     Check B(MYIFIRST, MYIFIRST) 
*     
      NDTOP = 0
      I = MYIFIRST
      DO WHILE ( I .LT. MYILAST ) 
         IF ( SELECT( I ) .NE. 1) EXIT
*        
*        Construct a row rotation that annihilates H( I + 1, I )  
*        and updates H( I, I )
         CALL PDLACP3( 2, I, H, DESCH, TMPH, 2, -1, -1, 0 )
         TEMP = TMPH( 1, 1 )
         CALL DLARTG( TEMP, TMPH( 2, 1 ), C, S,
     $      TMPH( 1, 1 ) )
         TMPH( 2, 1 ) = ZERO
         
*        
*        Update the 2 the global elements affected by the local DLARTG
*        
         CALL PDELSET( H, I + 1, I, DESCH, ZERO )
         CALL PDELSET( H, I, I, DESCH, TMPH( 1, 1 ) )
*
*        Set the delfated element in ALPHAR

         ALPHAR( I ) = TMPH( 1, 1)
         ALPHAI( I ) = ZERO
*        
*        Apply the row rotation to H, T, and Q. This won't 
*        cause any fill-in in T (since T(I, I) is zero).                  
        
         CALL PDROT( IHI - I, H, I, I + 1, DESCH, DESCH( N_ ),
     $      H, I + 1, I + 1, DESCH, DESCH( N_ ), C, S, DWORK, LWORK, 
     $      INFO )
         CALL PDROT( IHI - I, T, I, I + 1, DESCT, DESCT( N_ ),
     $      T, I + 1, I + 1, DESCT, DESCT( N_ ), C, S, DWORK, LWORK, 
     $      INFO )
         IF ( ILQ ) THEN
            CALL PDROT( IHIQ - ILOQ + 1, Q, ILOQ, I, DESCQ, 1, 
     $         Q, ILOQ, I + 1, DESCQ, 1, C, S, DWORK, LWORK, INFO )     
         END IF
*         
*        Does H(I, I - 1) exists ? If so apply C to it!

*
*        Set BETA to zero to indicate an
*        infinite eigenvalue
*         
         BETA( I ) = ZERO
         CALL PDELGET( 'A', ' ', TEMP, T, I + 1, I + 1, DESCT )
*
*        Check if this pertubation had impact on next row in T
         IF ( ABS( TEMP ) .GT. MAX( TTOL , SMLNUM )) THEN
            SELECT( I + 1 ) = 0
         ELSE
            SELECT( I + 1 ) = 1
            CALL PDELSET( T, I + 1, I + 1, DESCT, ZERO )
         END IF

         MYIFIRST = MYIFIRST + 1
         I = I + 1
         NDTOP = NDTOP + 1
         
      END DO
 
      ND = NDBOT + NDTOP
      NDDW = 0
      NDUP = 0
      NH = MYILAST - MYIFIRST + 1
*
*     Exit if only small problem remains
*
      IF (NH.LE.1) GOTO 200
 
      IF ( NH .LE. NB ) THEN
         MID = MYIFIRST + 1
         STARTI = MYIFIRST      
         ENDI = MYILAST
      ELSE
         MID = MYIFIRST + ( NH / 2 )
         STARTI = MID - 1
         ENDI = MID + 1
      END IF


      INFFOUND = .FALSE.
      DO WHILE ( .NOT. INFFOUND .AND. STARTI .LE. MYILAST ) 
         INFFOUND = (SELECT( STARTI ) .EQ. 1 )
         IF ( .NOT. INFFOUND ) THEN
            STARTI = STARTI + 1
         END IF
      END DO
      IF ( INFFOUND ) THEN         
*        
*        Chase to bottom, MYILAST might be updated
*        
         STARTI = MAX( STARTI - 1, MID - 1 )
         CALL PDHGEQZ8( ILSCHR, ILQ, ILZ, 
     $      STARTI, MYILAST, ILO, IHI, 
     $      H, DESCH, T, DESCT,
     $      ALPHAR, ALPHAI, BETA, 
     $      Q, DESCQ, Z, DESCZ, 
     $      ILOQ, IHIQ, ILOZ, IHIZ,
     $      GMAXWIN, NDDW, SELECT, DWORK, LWORK, INFO )
         
      END IF
      
*      
*     Make sure all have same (min)value of MYILAST
*
      CALL IGAMN2D( CTX, 'A', ' ', 1, 1, 
     $    MYILAST, 1, -1, -1, -1, -1, -1 ) 
      IF ( NH.LE.NB ) GOTO 200
      

      
*
*     Set the number of found zeros in the diagonal of B to none.
*     
      SELECT( 1 : IHI ) = 0
      
*              
*     Check (locally) for zeros in the diagonal of B, PDHGEQZ8 might have
*     perturbed some elements to or from zero.       
*                 
      NUM0ZEL = 0
      I = ENDI
      DO WHILE ( I .GE. MYIFIRST )
         BSIZE = 1 + MOD( I - 1, NB )
*        
*        If last block, recalc the block size.
*        
         IF (I - BSIZE + 1 .LE. MYIFIRST ) THEN
            BSIZE = I - MYIFIRST + 1
         END IF
         
         RSRC1 = INDXG2P( I, NB, 0, DESCT( RSRC_ ), NPROW )
         CSRC1 = INDXG2P( I, NB, 0, DESCT( CSRC_ ), NPCOL )
         
*
*        Work on locally on block T( I - BSIZE + 1 : I, I - BSIZE + 1 : I )
*        
         IF ( MYROW .EQ. RSRC1 .AND. MYCOL .EQ. CSRC1 ) THEN
            
            LIR = INDXG2L( I, NB, 0, 0, NPROW )
            LIC = INDXG2L( I, NB, 0, 0, NPCOL )
            
            DO J = 0, BSIZE - 1
               LVAL = T( LIR - J + ( LIC - J - 1 ) * DESCT( LLD_ ) ) 
               
               TEST = ABS( LVAL ).LE.MAX( TTOL , SMLNUM ) 
               IF ( TEST ) THEN
                  T( LIR - J + ( LIC - J - 1 ) * DESCT( LLD_ ) ) = 
     $               ZERO
                  SELECT( I - J ) = 1
                  NUM0ZEL = NUM0ZEL + 1                      
               END IF                      
            END DO
         END IF
*        
*        Goto next block
*        
         I = I - BSIZE
      END DO
      
      
*     
*     All check if infite eigenvalues exists
*     
      CALL IGSUM2D( CTX, 'A', ' ', 1, 1, NUM0ZEL, 1, -1, -1)  
      
*
*     Return if nothing to do
*              
      IF ( NUM0ZEL .EQ. 0 ) THEN
         ILAST = MYILAST
         IFIRST = MYIFIRST
         ND = ND + NDDW
         RETURN
      END IF      
*     
*     Let everyone have the same information 
*     about where the zeros are.
*           
      CALL IGAMX2D( CTX, 'A', ' ', IHI, 1, SELECT, 1, -1, -1, -1, 
     $   -1, -1 ) 
      

      INFFOUND = .FALSE.
      DO WHILE ( .NOT. INFFOUND .AND. ENDI .GE. MYIFIRST ) 
         INFFOUND = ( SELECT( ENDI ) .EQ. 1 )          
         IF ( .NOT. INFFOUND ) THEN
            ENDI = ENDI - 1
         END IF
      END DO  
      IF ( INFFOUND ) THEN   
         
*
*        Chase to top, MYIFIRST might be updated
*        
         ENDI = MIN( ENDI + 1, MID + 1 )
         
         CALL PDHGEQZ9( ILSCHR, ILQ, ILZ, 
     $      MYIFIRST, ENDI, ILO, IHI, 
     $      H, DESCH, T, DESCT,
     $      ALPHAR, ALPHAI, BETA, 
     $      Q, DESCQ, Z, DESCZ,
     $      ILOQ, IHIQ, ILOZ, IHIZ,
     $      GMAXWIN, NDUP, SELECT, DWORK, LWORK, INFO)
      END IF
      
*     
*     Let everyone have the same (max) value of MYIFIRST
*     
      CALL IGAMX2D( CTX, 'A', ' ', 1, 1, MYIFIRST, 1, -1, -1, -1, 
     $   -1, -1 )

 200  CONTINUE
      ND = ND + NDDW + NDUP
      ILAST = MYILAST
      IFIRST = MYIFIRST
      
      RETURN
*     
*     End of PDHGEQZ7
*     
      END


