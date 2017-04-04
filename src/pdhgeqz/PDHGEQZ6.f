***********************************************************************
*                                                                     *
*     PDHGEQZ6.f:                                                     *
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
      SUBROUTINE PDHGEQZ6(DIR, ILSCHR, ILQ, ILZ, A, DESCA, B, DESCB, Q,
     $   DESCQ, Z, DESCZ, U, V, LDUV, LTOP, ILO, IHI, N, NW, ILOQ, IHIQ,
     $   ILOZ, IHIZ, LDWORK, DWORK, INFO)



      IMPLICIT NONE
*     ..
*     .. Scalar Arguments ..
*     ..
      INTEGER             LTOP, ILO, IHI, N, NW, LDWORK, ILOZ, IHIZ,
     $   ILOQ, IHIQ, LDUV
      LOGICAL             ILQ, ILZ, ILSCHR
*     ..
*     .. Array Arguments ..
*     ..
      DOUBLE PRECISION    A( * ), B( * ), Q( * ), Z( * ), DWORK( * ) ,U(
     $   LDUV, * ), V( LDUV, * )
      INTEGER             DESCA( 9 ), DESCB( 9 ), DESCQ( 9 ), DESCZ( 9 )
      CHARACTER           DIR
      
*     Purpose
*     =======
*     
*     
*     PDHGEQZ6 is an auxiliary routine used to update off diagonal blocks
*     in A, B, Q, and Z with the orthognal matrices U, from left and V, from right.
*     
*     
*     All the inputs are assumed to be valid without 
*     checking.
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
*     DIR     (global input) CHARACTER*1
*     ='A': Update both from left and right
*     ='L': Update only from left
*     ='R': Update only from right
*     
*     ILSCHR  (global input) LOGICAL
*     = TRUE: Apply all transformations to A and B. 
*     = FALSE: A and B are not referenced.
*     
*     ILQ     (input) LOGICAL
*     = TRUE: Compute Q.
*     = FALSE: Q is not referenced.
*     
*     ILZ     (global input) LOGICAL
*     = TRUE: Compute Z.
*     = FALSE: Z is not referenced.
*     
*     A       (local input/output) DOUBLE PRECISION array, dimension (LLD_A, LOCc(N)).
*     On entry, the matrix A whose off-diagonal blocks are to be updated.
*     On exit, if ILSCHR is .TRUE., A is updated w.r.t to U, from left and 
*     V from right. 
*     If ILSCHR is .FALSE. A is not referenced.
*     
*     DESCA   (global and local input) INTEGER array of dimension DLEN_.
*     The array descriptor for the distributed matrix A. 
*     
*     B       (local input/output) DOUBLE PRECISION array, dimension (LLD_B, LOCc(N)).
*     On entry, the matrix B whose off-diagonal blocks are to be updated.
*     On exit, if ILSCHR is .TRUE., B is updated w.r.t to U, from left and 
*     V from right. 
*     If ILSCHR is .FALSE. B is not referenced.
*     
*     DESCB   (global and local input) INTEGER array of dimension DLEN_.
*     The array descriptor for the distributed matrix B. 
*     
*     Q       (local input/output) DOUBLE PRECISION array, dimension (LLD_Q, LOCc(N)).
*     If ILQ = FALSE:  Q is not referenced. 
*     If ILQ = TRUE:  on entry, Q must contain an orthogonal matrix.
*     On exit, Q is multiplied on the right with 
*     the product of transformations applied on H,T from the left.
*     
*     DESCQ   (global and local input) INTEGER array of dimension DLEN_.
*     The array descriptor for the distributed matrix Q. 
*     
*     Z       (local input/output) DOUBLE PRECISION array, dimension (LLD_Z, LOCc(N)).
*     If ILZ = FALSE:  Z is not referenced. 
*     If ILZ = TRUE:  on entry, Z must contain an orthogonal matrix.
*     On exit, Z is multiplied on the right with 
*     the product of transformations applied on H,T from the right. 
*     
*     DESCZ   (global and local input) INTEGER array of dimension DLEN_.
*     The array descriptor for the distributed matrix Z. 
*     
*     U       (local input) DOUBLE PRECSION array, dimension(LDUV, NW)
*     U holds the transformations to be applied from the left on A,B and on 
*     the right on Q.
*     
*     V       (local input) DOUBLE PRECSION array, dimension(LDUV, NW)
*     V holds the transformations to be applied from the right on A,B, and Z. 
*     
*     LDUV    (local input) INTEGER
*     Leading dimension for the matrices U and V.
*     
*     LTOP    (global input) INTEGER
*     Defines the row index where the diagonal block begins. 
*     Row updates on A and B are restricted to columns > LTOP+NW-1
*     Column updates on A and B are restricted < LTOP
*     
*     ILO     (global input) INTERGER
*     Updates are not applied to rows(from left)/cols(from right) on A and B below this index. 
*     
*     IHI     (global input) INTERGER
*     Updates are not applied to rows(from left)/cols(from right) on A and B greater than this index. 
*     
*     N       (global input) INTEGER 
*     The order of the matrices A,B,Q, and Z.  N >= 0. 
*     
*     NW      (global input) INTEGER 
*     Size of the diagonal block and transformation 
*     matrices U and V (NW * NW)
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
*     
*     DWORK   (local workspace/global output) DOUBLE PRECISION array,
*     dimension (LWORK) 
*     
*     LDWORK  (global input) INTEGER 
*     The dimension of the array WORK. 
*     Vital that all participating processors have the same value for LDWORK
*     
*     INFO    (global output) INTEGER
*     = 0:  successful exit.
*     ===========================================================

      
*     ..
*     .. Parameters ..
*     ..
      DOUBLE PRECISION   HALF, ZERO, ONE, SAFETY
      PARAMETER          ( HALF = 0.5D+0, ZERO = 0.0D+0, ONE = 1.0D+0,
     $   SAFETY = 1.0D+2 )
      INTEGER            BLOCK_CYCLIC_2D, CSRC_, CTXT_, DLEN_, DT_, 
     $   LLD_, MB_, M_, NB_, N_, RSRC_
      PARAMETER          ( BLOCK_CYCLIC_2D = 1, DLEN_ = 9, DT_ = 1,
     $   CTXT_ = 2, M_ = 3, N_ = 4, MB_ = 5, NB_ = 6,
     $   RSRC_ = 7, CSRC_ = 8, LLD_ = 9 )
*     ..
*     .. Local Scalars ..
*     ..
      INTEGER             LDA, LDB, LDQ, LDZ, NB, ICTXT, MYROW, MYCOL,
     $   NPROW, NPCOL, LLDTMP, NODE, LEFT, RIGHT, UP, DOWN, HSTEP, VSTEP
     $   , II, JJ, KLN, ITMP1, ITMP2, ICOL, ICOL1, KKCOL, IROW, IROW1,
     $   KKROW, D1, D2, INFO, IAM, NPROCS, RSRC, CSRC
      LOGICAL             DOLEFT, DORIGHT

*     ..
*     .. Local Arrays
*     ..
      INTEGER             DESCUV(9), DESCWH(9), DESCWV(9)

*     ..
*     .. Externals ..
*     ..
      EXTERNAL            NUMROC, INDXG2L, INDXG2P, LSAME
      INTEGER             NUMROC, INDXG2L, INDXG2P
      LOGICAL             LSAME
      
*     ..
*     .. Intrinsic Functions ..
*     ..
      INTRINSIC           ABS, DBLE, MAX, MIN, SQRT, MOD


      INFO = 0

*     Extract current communication context and get grid parameters
      ICTXT = DESCA( CTXT_ )
      CALL BLACS_GRIDINFO( ICTXT, NPROW, NPCOL, MYROW, MYCOL )  
      CALL BLACS_PINFO( IAM, NPROCS )

      
      LDA = DESCA( LLD_ )
      LDB = DESCB( LLD_ )
      LDQ = DESCQ( LLD_ )
      LDZ = DESCZ( LLD_ )

*     NB_ .EQ. MB_ is assumed. Get current blocking factor and store in NB.
      NB = DESCA( NB_ )

      NODE = MYROW * NPCOL + MYCOL
      LEFT = MOD( MYCOL + NPCOL- 1, NPCOL )
      RIGHT = MOD( MYCOL + 1, NPCOL )
      UP = MOD( MYROW + NPROW - 1, NPROW )
      DOWN = MOD( MYROW + 1, NPROW )
      DOLEFT = LSAME( DIR, 'L' ) .OR. LSAME( DIR, 'A' )
      DORIGHT = LSAME( DIR, 'R' ) .OR. LSAME( DIR, 'A' )
      
      RSRC = DESCA( RSRC_ )
      CSRC = DESCA( CSRC_ )
      
      
      
      IF( MOD( ILO - 1, NB ) + NW .LE. NB .OR. NPROW * NPCOL .EQ. 1 )
     $   THEN
         
         HSTEP = LDWORK / NW
         VSTEP = HSTEP          

*        Simplest case: the deflation window is located on one
*        processor.
*        Call DGEMM directly to perform the update.
         IF ( ILSCHR ) THEN
            IF ( DOLEFT ) THEN

*              Update horizontal slab in A and B with U
               CALL INFOG2L( ILO, IHI + 1, DESCA, NPROW, NPCOL, MYROW,
     $            MYCOL, IROW, ICOL, II, JJ )
               IF( MYROW .EQ. II) THEN
                  ICOL1 = NUMROC( N, NB, MYCOL, CSRC, NPCOL )
                  DO 10 KKCOL = ICOL, ICOL1, HSTEP
                     KLN = MIN( HSTEP, ICOL1 - KKCOL + 1 )   
*                    A
                     CALL DGEMM( 'T', 'N', NW, KLN, NW, ONE, U,
     $                  LDUV, A( IROW + ( KKCOL - 1) * LDA ), LDA, ZERO,
     $                  DWORK, NW )    
                     CALL DLACPY( 'A', NW, KLN, DWORK, NW,
     $                  A( IROW + ( KKCOL - 1 ) * LDA ), LDA )
*                    B 
                     CALL DGEMM( 'T', 'N', NW, KLN, NW, ONE, U,
     $                  LDUV, B( IROW + ( KKCOL - 1) * LDB ), LDB, ZERO,
     $                  DWORK, NW )    
                     CALL DLACPY( 'A', NW, KLN, DWORK, NW,
     $                  B( IROW + ( KKCOL - 1) * LDB ), LDB )
 10               CONTINUE
               END IF          
            END IF
         END IF
 

         IF ( ILSCHR ) THEN
            IF ( DORIGHT ) THEN
*              Update vertical slab in A and B with V
               CALL INFOG2L( LTOP, ILO, DESCA, NPROW, NPCOL, MYROW,
     $            MYCOL, IROW, ICOL, II, JJ )
               IF( MYCOL .EQ. JJ ) THEN
                  CALL INFOG2L( ILO-1, ILO, DESCA, NPROW, NPCOL,
     $               MYROW, MYCOL, IROW1, ICOL1, ITMP1, ITMP2 )
                  IF( MYROW .NE. ITMP1 ) IROW1 = IROW1 - 1
                  
                  DO 20 KKROW = IROW, IROW1, VSTEP
                     KLN = MIN( VSTEP, IROW1 - KKROW + 1 )
*                    A
                     CALL DGEMM( 'N', 'N', KLN, NW, NW, ONE,
     $                  A( KKROW + ( ICOL - 1 ) * LDA ), LDA, V, 
     $                  LDUV,  ZERO, DWORK, KLN )
                     CALL DLACPY( 'A', KLN, NW, DWORK, KLN,
     $                  A( KKROW + ( ICOL - 1 ) * LDA ), LDA )
*                    B
                     CALL DGEMM( 'N', 'N', KLN, NW, NW, ONE,
     $                  B( KKROW +( ICOL - 1) * LDB ), LDB, V, 
     $                  LDUV,  ZERO, DWORK, KLN )
                     CALL DLACPY( 'A', KLN, NW, DWORK, KLN,
     $                  B( KKROW + ( ICOL - 1) * LDB ), LDB )
 20               CONTINUE
               END IF
            END IF
         END IF      

*        Update vertical slab in Q with U
         IF( ILQ .AND. DORIGHT ) THEN
            CALL INFOG2L( ILOQ, ILO, DESCQ, NPROW, NPCOL, MYROW,
     $         MYCOL, IROW, ICOL, II, JJ )
            IF( MYCOL .EQ. JJ ) THEN
               CALL INFOG2L( IHIQ, ILO, DESCQ, NPROW, NPCOL,
     $            MYROW, MYCOL, IROW1, ICOL1, ITMP1, ITMP2 )
               IF( MYROW .NE. ITMP1 ) IROW1 = IROW1-1
               DO 30 KKROW = IROW, IROW1, VSTEP
                  KLN = MIN( VSTEP, IROW1-KKROW+1 )
                  CALL DGEMM( 'N', 'N', KLN, NW, NW, ONE,
     $               Q( KKROW+(ICOL-1)*LDQ ), LDQ, U, LDUV, 
     $               ZERO, DWORK, KLN )
                  CALL DLACPY( 'A', KLN, NW, DWORK, KLN,
     $               Q( KKROW+(ICOL-1)*LDQ ), LDQ )
                  
 30            CONTINUE
            END IF
         END IF 

*        Update vertical slab in Z with V
         IF( ILZ .AND. DORIGHT ) THEN
            CALL INFOG2L( ILOZ, ILO, DESCZ, NPROW, NPCOL, MYROW,
     $         MYCOL, IROW, ICOL, II, JJ )
            IF( MYCOL .EQ. JJ ) THEN
               CALL INFOG2L( IHIZ, ILO, DESCZ, NPROW, NPCOL,
     $            MYROW, MYCOL, IROW1, ICOL1, ITMP1, ITMP2 )
               IF( MYROW .NE. ITMP1 ) IROW1 = IROW1 - 1
               DO 35 KKROW = IROW, IROW1, VSTEP
                  KLN = MIN( VSTEP, IROW1-KKROW+1 )
                  CALL DGEMM( 'N', 'N', KLN, NW, NW, ONE,
     $               Z( KKROW + ( ICOL - 1 ) * LDZ ), LDZ, V, LDUV, 
     $               ZERO, DWORK, KLN )
                  CALL DLACPY( 'A', KLN, NW, DWORK, KLN,
     $               Z( KKROW + ( ICOL - 1 )  *LDZ ), LDZ )
 35            CONTINUE
            END IF
            
         END IF
         
      ELSE IF( MOD( ILO - 1, NB ) + NW .LE. 2 * NB ) THEN

*        More complicated case: the deflation window lay on a 2x2
*        processor mesh
*        Call DGEMM locally and communicate by pair

         D1 = NB - MOD( ILO-1, NB )
         D2 = NW - D1
         HSTEP = LDWORK / NW
         VSTEP = HSTEP     
         

         IF( ILSCHR ) THEN
            
*           Update horizontal slab in A and B with U
            IF ( DOLEFT ) THEN
               CALL INFOG2L( ILO,  IHI + 1, DESCA, NPROW, NPCOL, MYROW,
     $            MYCOL, IROW, ICOL, II, JJ )
               IF( MYROW .EQ. UP ) THEN
                  IF( MYROW .EQ. II ) THEN
                     ICOL1 = NUMROC( N, NB, MYCOL, CSRC, NPCOL )
                     DO 40 KKCOL = ICOL, ICOL1, HSTEP
                        KLN = MIN( HSTEP, ICOL1 - KKCOL + 1 )

*                       A
                        CALL DGEMM( 'T', 'N', NW, KLN, NW, ONE, U, LDUV,
     $                     A( IROW + ( KKCOL - 1 ) *LDA ), LDA, ZERO,
     $                     DWORK, NW )
                        CALL DLACPY( 'A', NW, KLN, DWORK, NW,
     $                     A( IROW + ( KKCOL - 1) * LDA ), LDA )

*                       B
                        CALL DGEMM( 'T', 'N', NW, KLN, NW, ONE, U, LDUV,
     $                     B( IROW + ( KKCOL - 1) * LDB ), LDB, ZERO,
     $                     DWORK, NW )
                        CALL DLACPY( 'A', NW, KLN, DWORK, NW,
     $                     B( IROW + ( KKCOL - 1 ) * LDB ), LDB )
 40                  CONTINUE
                  END IF
               ELSE
                  IF( MYROW .EQ. II ) THEN
                     ICOL1 = NUMROC( N, NB, MYCOL, CSRC, NPCOL )
                     DO 50 KKCOL = ICOL, ICOL1, HSTEP
                        KLN = MIN( HSTEP, ICOL1-KKCOL+1 )

*                       A 
                        CALL DGEMM( 'T', 'N', D2, KLN, D1, ONE,
     $                     U( 1, D1+1 ), LDUV, 
     $                     A( IROW+(KKCOL-1)*LDA ), LDA, ZERO, 
     $                     DWORK( D1+1 ), NW )
                        CALL DGESD2D( ICTXT, D2, KLN, DWORK( D1+1 ),
     $                     NW, DOWN, MYCOL )
                        CALL DGERV2D( ICTXT, D1, KLN, DWORK, NW, 
     $                     DOWN, MYCOL )
                        CALL DGEMM( 'T', 'N', D1, KLN, D1, ONE,
     $                     U, LDUV, A( IROW+(KKCOL-1)*LDA ), LDA, ONE,
     $                     DWORK, NW )
                        CALL DLACPY( 'A', D1, KLN, DWORK, NW,
     $                     A( IROW+(KKCOL-1)*LDA ), LDA )

*                       B 
                        CALL DGEMM( 'T', 'N', D2, KLN, D1, ONE,
     $                     U( 1, D1+1 ), LDUV, 
     $                     B( IROW+(KKCOL-1)*LDB ), LDB, ZERO, 
     $                     DWORK( D1+1 ), NW )
                        CALL DGESD2D( ICTXT, D2, KLN, DWORK( D1+1 ),
     $                     NW, DOWN, MYCOL )
                        CALL DGERV2D( ICTXT, D1, KLN, DWORK, NW, 
     $                     DOWN, MYCOL )
                        CALL DGEMM( 'T', 'N', D1, KLN, D1, ONE,
     $                     U, LDUV, B( IROW+(KKCOL-1)*LDB ), LDB, ONE,
     $                     DWORK, NW )
                        CALL DLACPY( 'A', D1, KLN, DWORK, NW,
     $                     B( IROW+(KKCOL-1)*LDB ), LDB )
 50                  CONTINUE
                  ELSE IF( UP .EQ. II ) THEN
                     ICOL1 = NUMROC( N, NB, MYCOL, CSRC, NPCOL )
                     DO 60 KKCOL = ICOL, ICOL1, HSTEP
                        KLN = MIN( HSTEP, ICOL1-KKCOL+1 )

*                       A
                        CALL DGEMM( 'T', 'N', D1, KLN, D2, ONE,
     $                     U( D1+1, 1 ), LDUV, 
     $                     A( IROW+(KKCOL-1)*LDA ), LDA, ZERO, 
     $                     DWORK, NW )
                        CALL DGESD2D( ICTXT, D1, KLN, DWORK, NW, UP,
     $                     MYCOL )
                        CALL DGERV2D( ICTXT, D2, KLN, DWORK( D1+1 ),
     $                     NW, UP, MYCOL )
                        CALL DGEMM( 'T', 'N', D2, KLN, D2, ONE,
     $                     U( D1+1, D1+1 ), LDUV,
     $                     A( IROW+(KKCOL-1)*LDA ), LDA, ONE,
     $                     DWORK( D1+1 ), NW )
                        CALL DLACPY( 'A', D2, KLN, DWORK( D1+1 ), NW,
     $                     A( IROW+(KKCOL-1)*LDA ), LDA )

*                       B
                        CALL DGEMM( 'T', 'N', D1, KLN, D2, ONE,
     $                     U( D1+1, 1 ), LDUV, 
     $                     B( IROW+(KKCOL-1)*LDB ), LDB, ZERO, 
     $                     DWORK, NW )
                        CALL DGESD2D( ICTXT, D1, KLN, DWORK, NW, UP,
     $                     MYCOL )
                        CALL DGERV2D( ICTXT, D2, KLN, DWORK( D1+1 ),
     $                     NW, UP, MYCOL )
                        CALL DGEMM( 'T', 'N', D2, KLN, D2, ONE,
     $                     U( D1+1, D1+1 ), LDUV,
     $                     B( IROW+(KKCOL-1)*LDB ), LDB, ONE,
     $                     DWORK( D1+1 ), NW )
                        CALL DLACPY( 'A', D2, KLN, DWORK( D1+1 ), NW,
     $                     B( IROW+(KKCOL-1)*LDB ), LDB )
 60                  CONTINUE
                  END IF
               END IF
            END IF

*           .. Update vertical slab in A and B with V
            IF ( DORIGHT ) THEN
               CALL INFOG2L( LTOP, ILO, DESCA, NPROW, NPCOL, MYROW,
     $            MYCOL, IROW, ICOL, II, JJ )
               IF( MYCOL .EQ. LEFT ) THEN
                  IF( MYCOL .EQ. JJ ) THEN
                     CALL INFOG2L( ILO-1, ILO, DESCA, NPROW, NPCOL,
     $                  MYROW, MYCOL, IROW1, ICOL1, ITMP1, ITMP2 )
                     IF( MYROW .NE. ITMP1 ) IROW1 = IROW1-1
                     DO 70 KKROW = IROW, IROW1, VSTEP
                        KLN = MIN( VSTEP, IROW1-KKROW+1 )

*                       A
                        CALL DGEMM( 'N', 'N', KLN, NW, NW, ONE,
     $                     A( KKROW+(ICOL-1)*LDA ), LDA, 
     $                     V, LDUV, ZERO,
     $                     DWORK, KLN )
                        CALL DLACPY( 'A', KLN, NW, DWORK, KLN,
     $                     A( KKROW+(ICOL-1)*LDA ), LDA )

*                       B                      
                        CALL DGEMM( 'N', 'N', KLN, NW, NW, ONE,
     $                     B( KKROW+(ICOL-1)*LDB ), LDB, V, LDUV, 
     $                     ZERO, DWORK, KLN )
                        CALL DLACPY( 'A', KLN, NW, DWORK, KLN,
     $                     B( KKROW+(ICOL-1)*LDB ), LDB )
                        
 70                  CONTINUE
                  END IF
               ELSE
                  IF( MYCOL .EQ. JJ ) THEN
                     CALL INFOG2L( ILO-1, ILO, DESCA, NPROW, NPCOL,
     $                  MYROW, MYCOL, IROW1, ICOL1, ITMP1, ITMP2 )
                     IF( MYROW .NE. ITMP1 ) IROW1 = IROW1-1
                     DO 80 KKROW = IROW, IROW1, VSTEP
                        KLN = MIN( VSTEP, IROW1-KKROW+1 )

*                       A   
                        CALL DGEMM( 'N', 'N', KLN, D2, D1, ONE,
     $                     A( KKROW+(ICOL-1)*LDA ), LDA,
     $                     V( 1, D1+1 ), LDUV, ZERO, 
     $                     DWORK( 1+D1*KLN ), KLN )
                        CALL DGESD2D( ICTXT, KLN, D2, 
     $                     DWORK( 1+D1*KLN ), KLN, MYROW, RIGHT )
                        CALL DGERV2D( ICTXT, KLN, D1, DWORK, KLN, 
     $                     MYROW, RIGHT )
                        CALL DGEMM( 'N', 'N', KLN, D1, D1, ONE,
     $                     A( KKROW+(ICOL-1)*LDA ), LDA, V, LDUV, 
     $                     ONE, DWORK, KLN )
                        CALL DLACPY( 'A', KLN, D1, DWORK, KLN,
     $                     A( KKROW+(ICOL-1)*LDA ), LDA )

*                       B  
                        CALL DGEMM( 'N', 'N', KLN, D2, D1, ONE,
     $                     B( KKROW+(ICOL-1)*LDB ), LDB,
     $                     V( 1, D1+1 ), LDUV, ZERO, 
     $                     DWORK( 1+D1*KLN ), KLN )
                        CALL DGESD2D( ICTXT, KLN, D2, 
     $                     DWORK( 1+D1*KLN ), KLN, MYROW, RIGHT )
                        CALL DGERV2D( ICTXT, KLN, D1, DWORK, KLN, 
     $                     MYROW, RIGHT )
                        CALL DGEMM( 'N', 'N', KLN, D1, D1, ONE,
     $                     B( KKROW+(ICOL-1)*LDB ), LDB, V, LDUV, 
     $                     ONE, DWORK, KLN )
                        CALL DLACPY( 'A', KLN, D1, DWORK, KLN,
     $                     B( KKROW+(ICOL-1)*LDB ), LDB )
 80                  CONTINUE
                  ELSE IF ( LEFT .EQ. JJ ) THEN
                     CALL INFOG2L(ILO-1, ILO, DESCA, NPROW, NPCOL,
     $                  MYROW, MYCOL, IROW1, ICOL1, ITMP1, ITMP2 )
                     IF( MYROW .NE. ITMP1 ) IROW1 = IROW1-1
                     DO 90 KKROW = IROW, IROW1, VSTEP
                        KLN = MIN( VSTEP, IROW1-KKROW+1 )

*                       A   
                        CALL DGEMM( 'N', 'N', KLN, D1, D2, ONE,
     $                     A( KKROW+(ICOL-1)*LDA ), LDA, 
     $                     V( D1+1, 1 ), LDUV, ZERO, DWORK, KLN )
                        CALL DGESD2D( ICTXT, KLN, D1, DWORK, KLN, 
     $                     MYROW, LEFT )
                        CALL DGERV2D( ICTXT, KLN, D2, 
     $                     DWORK( 1+D1*KLN ), KLN, MYROW, LEFT )
                        CALL DGEMM( 'N', 'N', KLN, D2, D2, ONE,
     $                     A( KKROW+(ICOL-1)*LDA ), LDA, 
     $                     V( D1+1, D1+1 ), LDUV, ONE, 
     $                     DWORK( 1+D1*KLN ), KLN )
                        CALL DLACPY( 'A', KLN, D2, DWORK( 1+D1*KLN ), 
     $                     KLN, A( KKROW+(ICOL-1)*LDA ), LDA )

*                       B    
                        CALL DGEMM( 'N', 'N', KLN, D1, D2, ONE,
     $                     B( KKROW+(ICOL-1)*LDB ), LDB, 
     $                     V( D1+1, 1 ), LDUV, ZERO, DWORK, KLN )
                        CALL DGESD2D( ICTXT, KLN, D1, DWORK, KLN,
     $                     MYROW, LEFT )
                        CALL DGERV2D( ICTXT, KLN, D2, 
     $                     DWORK( 1+D1*KLN ), KLN, MYROW, LEFT )
                        CALL DGEMM( 'N', 'N', KLN, D2, D2, ONE,
     $                     B( KKROW+(ICOL-1)*LDB ), LDB, 
     $                     V( D1+1, D1+1 ), LDUV, ONE, 
     $                     DWORK( 1+D1*KLN ), KLN )
                        CALL DLACPY( 'A', KLN, D2, DWORK( 1+D1*KLN ), 
     $                     KLN, B( KKROW+(ICOL-1)*LDB ), LDB )
 90                  CONTINUE
                  END IF
               END IF
            END IF
         END IF

*        Update vertical slab in Q with U
         IF( ILQ .AND. DORIGHT ) THEN
            CALL INFOG2L( ILOQ, ILO, DESCQ, NPROW, NPCOL, MYROW,
     $         MYCOL, IROW, ICOL, II, JJ )
            IF( MYCOL .EQ. LEFT ) THEN
               IF( MYCOL .EQ. JJ ) THEN
                  CALL INFOG2L( IHIZ, ILO, DESCQ, NPROW, NPCOL,
     $               MYROW, MYCOL, IROW1, ICOL1, ITMP1, ITMP2 )
                  IF( MYROW .NE. ITMP1 ) IROW1 = IROW1-1
                  DO 100 KKROW = IROW, IROW1, VSTEP
                     KLN = MIN( VSTEP, IROW1-KKROW+1 )
                     CALL DGEMM( 'N', 'N', KLN, NW, NW, ONE,
     $                  Q( KKROW+(ICOL-1)*LDQ ), LDQ, U, LDUV, 
     $                  ZERO, DWORK, KLN )
                     CALL DLACPY( 'A', KLN, NW, DWORK, KLN,
     $                  Q( KKROW+(ICOL-1)*LDQ ), LDQ )
 100              CONTINUE
               END IF
            ELSE
               IF( MYCOL .EQ. JJ ) THEN
                  CALL INFOG2L( IHIQ, ILO, DESCQ, NPROW, NPCOL,
     $               MYROW, MYCOL, IROW1, ICOL1, ITMP1, ITMP2 )
                  IF( MYROW .NE. ITMP1 ) IROW1 = IROW1-1
                  DO 110 KKROW = IROW, IROW1, VSTEP
                     KLN = MIN( VSTEP, IROW1-KKROW+1 )
                     CALL DGEMM( 'N', 'N', KLN, D2, D1, ONE,
     $                  Q( KKROW+(ICOL-1)*LDZ ), LDQ,
     $                  U( 1, D1+1 ), LDUV, ZERO, DWORK( 1+D1*KLN ),
     $                  KLN )
                     CALL DGESD2D( ICTXT, KLN, D2, 
     $                  DWORK( 1+D1*KLN ),KLN, MYROW, RIGHT )
                     CALL DGERV2D( ICTXT, KLN, D1, DWORK, KLN, 
     $                  MYROW, RIGHT )
                     CALL DGEMM( 'N', 'N', KLN, D1, D1, ONE,
     $                  Q( KKROW+(ICOL-1)*LDQ ), LDQ, U, LDUV, ONE,
     $                  DWORK, KLN )
                     CALL DLACPY( 'A', KLN, D1, DWORK, KLN,
     $                  Q( KKROW+(ICOL-1)*LDQ ), LDQ )
 110              CONTINUE
               ELSE IF( LEFT .EQ. JJ ) THEN
                  CALL INFOG2L( IHIQ, ILO, DESCQ, NPROW, NPCOL,
     $               MYROW, MYCOL, IROW1, ICOL1, ITMP1, ITMP2 )
                  IF( MYROW .NE. ITMP1 ) IROW1 = IROW1-1
                  DO 120 KKROW = IROW, IROW1, VSTEP
                     KLN = MIN( VSTEP, IROW1-KKROW+1 )
                     CALL DGEMM( 'N', 'N', KLN, D1, D2, ONE,
     $                  Q( KKROW+(ICOL-1)*LDQ ), LDQ,
     $                  U( D1+1, 1 ), LDUV, ZERO, DWORK, KLN )
                     CALL DGESD2D( ICTXT, KLN, D1, DWORK, KLN, 
     $                  MYROW, LEFT )
                     CALL DGERV2D( ICTXT, KLN, D2,
     $                  DWORK( 1+D1*KLN ), KLN, MYROW, LEFT )
                     CALL DGEMM( 'N', 'N', KLN, D2, D2, ONE,
     $                  Q( KKROW+(ICOL-1)*LDQ ), LDQ,
     $                  U( D1+1, D1+1 ), LDUV, ONE,
     $                  DWORK( 1+D1*KLN ), KLN )
                     CALL DLACPY( 'A', KLN, D2, DWORK( 1+D1*KLN ),
     $                  KLN, Q( KKROW+(ICOL-1)*LDQ ), LDQ )
 120              CONTINUE
               END IF
            END IF
         END IF

*        Update vertical slab in Z with V
         IF( ILZ .AND. DORIGHT) THEN
            CALL INFOG2L( ILOZ, ILO, DESCZ, NPROW, NPCOL, MYROW,
     $         MYCOL, IROW, ICOL, II, JJ )
            IF( MYCOL .EQ. LEFT ) THEN
               IF( MYCOL .EQ. JJ ) THEN
                  CALL INFOG2L( IHIZ,ILO, DESCA, NPROW, NPCOL,
     $               MYROW, MYCOL, IROW1, ICOL1, ITMP1, ITMP2 )
                  IF( MYROW .NE. ITMP1 ) IROW1 = IROW1-1
                  DO 200 KKROW = IROW, IROW1, VSTEP
                     KLN = MIN( VSTEP, IROW1-KKROW+1 )
                     CALL DGEMM( 'N', 'N', KLN, NW, NW, ONE,
     $                  Z( KKROW+(ICOL-1)*LDZ ), LDZ, V, LDUV, 
     $                  ZERO, DWORK, KLN )
                     CALL DLACPY( 'A', KLN, NW, DWORK, KLN,
     $                  Z( KKROW+(ICOL-1)*LDZ ), LDZ )
 200              CONTINUE
               END IF
            ELSE
               IF( MYCOL .EQ. JJ ) THEN
                  CALL INFOG2L( IHIZ, ILO, DESCZ, NPROW, NPCOL,
     $               MYROW, MYCOL, IROW1, ICOL1, ITMP1, ITMP2 )
                  IF( MYROW .NE. ITMP1 ) IROW1 = IROW1-1
                  DO 210 KKROW = IROW, IROW1, VSTEP
                     KLN = MIN( VSTEP, IROW1-KKROW+1 )
                     CALL DGEMM( 'N', 'N', KLN, D2, D1, ONE,
     $                  Z( KKROW+(ICOL-1)*LDZ ), LDZ,
     $                  V( 1, D1+1 ), LDUV, ZERO, 
     $                  DWORK( 1+D1*KLN ), KLN )
                     CALL DGESD2D( ICTXT, KLN, D2, 
     $                  DWORK( 1+D1*KLN ),KLN, MYROW, RIGHT )
                     CALL DGERV2D( ICTXT, KLN, D1, DWORK, KLN, 
     $                  MYROW, RIGHT )
                     CALL DGEMM( 'N', 'N', KLN, D1, D1, ONE,
     $                  Z( KKROW+(ICOL-1)*LDZ ), LDZ, 
     $                  V, LDUV, ONE,
     $                  DWORK, KLN )
                     CALL DLACPY( 'A', KLN, D1, DWORK, KLN,
     $                  Z( KKROW+(ICOL-1)*LDZ ), LDZ )
 210              CONTINUE
               ELSE IF( LEFT .EQ. JJ ) THEN
                  CALL INFOG2L( IHIZ, ILO, DESCA, NPROW, NPCOL,
     $               MYROW, MYCOL, IROW1, ICOL1, ITMP1, ITMP2 )
                  IF( MYROW .NE. ITMP1 ) IROW1 = IROW1-1
                  DO 220 KKROW = IROW, IROW1, VSTEP
                     KLN = MIN( VSTEP, IROW1-KKROW+1 )
                     CALL DGEMM( 'N', 'N', KLN, D1, D2, ONE,
     $                  Z( KKROW+(ICOL-1)*LDZ ), LDZ,
     $                  V( D1+1, 1 ), LDUV, ZERO, DWORK, KLN )
                     CALL DGESD2D( ICTXT, KLN, D1, DWORK, KLN, 
     $                  MYROW, LEFT )
                     CALL DGERV2D( ICTXT, KLN, D2,
     $                  DWORK( 1+D1*KLN ), KLN, MYROW, LEFT )
                     CALL DGEMM( 'N', 'N', KLN, D2, D2, ONE,
     $                  Z( KKROW+(ICOL-1)*LDZ ), LDZ,
     $                  V( D1+1, D1+1 ), LDUV, ONE,
     $                  DWORK( 1+D1*KLN ), KLN )
                     CALL DLACPY( 'A', KLN, D2, DWORK( 1+D1*KLN ),
     $                  KLN, Z( KKROW+(ICOL-1)*LDZ ), LDZ )
 220              CONTINUE
               END IF
            END IF
         END IF        
      ELSE
         WRITE(*,*)'Shouldnt be here...'
         STOP

*        Most complicated case: the deflation window lay across the ..
*        border of the processor mesh ..
*        Treat U and V as a distributed matrix and call PDGEMM ..
         HSTEP = LDWORK / NW * NPCOL
         VSTEP = LDWORK / NW * NPROW
         LLDTMP = NUMROC( NW, NW, MYROW, 0, NPROW )
         LLDTMP = MAX( 1, LLDTMP )
         
         
         CALL DESCINIT( DESCUV, NW, NW, NW, NW, 0, 0, ICTXT,
     $      LLDTMP, INFO )
         CALL DESCINIT( DESCWH, NW, HSTEP, NW, LDWORK / NW, 0,
     $      0, ICTXT, LLDTMP, INFO )
         
         IF( ILSCHR ) THEN
            IF (DOLEFT) THEN
*              Update horizontal slab in A and B with U
               DO 230 KKCOL = IHI+1, N, HSTEP
                  KLN = MIN( HSTEP, N-KKCOL+1 )
*                 A                      
                  CALL PDGEMM( 'T', 'N', NW, KLN, NW, ONE, U, 1, 1,
     $               DESCUV, A, ILO, KKCOL, DESCA, ZERO, DWORK, 1,
     $               1, DESCWH )
                  CALL PDGEMR2D( NW, KLN, DWORK, 1, 1, DESCWH, A,
     $               ILO, KKCOL, DESCA, ICTXT )
*                 B      
                  CALL PDGEMM( 'T', 'N', NW, KLN, NW, ONE, U, 1, 1,
     $               DESCUV, B, ILO, KKCOL, DESCB, ZERO, DWORK, 1,
     $               1, DESCWH )
                  CALL PDGEMR2D( NW, KLN, DWORK, 1, 1, DESCWH, B,
     $               ILO, KKCOL, DESCB, ICTXT )
                  
 230           CONTINUE
            END IF
*           Update vertical slab in A and B with V
            IF (DORIGHT) THEN
               DO 240 KKROW = LTOP, ILO-1, VSTEP
                  KLN = MIN( VSTEP, ILO-KKROW  )
                  LLDTMP = NUMROC( KLN, LDWORK / NW, MYROW, 0, NPROW )
                  LLDTMP = MAX( 1, LLDTMP )              
                  CALL DESCINIT( DESCWV, KLN, NW, LDWORK / NW, NW, 0,
     $               0, ICTXT, LLDTMP, INFO )
*                 A                      
                  CALL PDGEMM( 'N', 'N', KLN, NW, NW, ONE, A, KKROW,
     $               ILO, DESCA, V, 1, 1, DESCUV, ZERO, DWORK, 1, 1,
     $               DESCWV )
                  CALL PDGEMR2D( KLN, NW, DWORK, 1, 1, DESCWV, A, KKROW,
     $               ILO, DESCA, ICTXT )              
*                 B                      
                  CALL PDGEMM( 'N', 'N', KLN, NW, NW, ONE, B, KKROW,
     $               ILO, DESCB, V, 1, 1, DESCUV, ZERO, DWORK, 1, 1,
     $               DESCWV )
                  CALL PDGEMR2D( KLN, NW, DWORK, 1, 1, DESCWV, B, KKROW,
     $               ILO, DESCB, ICTXT )              
 240           CONTINUE
            END IF
         END IF

*        Update vertical slab in Q with U
         IF( ILQ .AND. DORIGHT) THEN
            DO 250 KKROW = ILOQ, IHIQ, VSTEP
               KLN = MIN( VSTEP, IHIQ-KKROW+1 )
               LLDTMP = NUMROC( KLN, LDWORK / NW, MYROW, 0, NPROW )
               LLDTMP = MAX( 1, LLDTMP )
               CALL DESCINIT( DESCWV, KLN, NW, LDWORK / NW, NW,
     $            0, 0, ICTXT, LLDTMP, INFO )
               CALL PDGEMM( 'N', 'N', KLN, NW, NW, ONE, Q, KKROW,
     $            ILO, DESCQ, U, 1, 1, DESCUV, ZERO, DWORK, 1,
     $            1, DESCWV )
               CALL PDGEMR2D( KLN, NW, DWORK, 1, 1, DESCWV, Q,
     $            KKROW, ILO, DESCQ, ICTXT )
 250        CONTINUE
         END IF  
         
*        Update vertical slab in Z with V
         IF( ILZ .AND. DORIGHT ) THEN
            DO 260 KKROW = ILOZ, IHIZ, VSTEP
               KLN = MIN( VSTEP, IHIZ-KKROW+1 )
               LLDTMP = NUMROC( KLN, LDWORK / NW, MYROW, 0, NPROW )
               LLDTMP = MAX( 1, LLDTMP )
               CALL DESCINIT( DESCWV, KLN, NW, LDWORK / NW, NW,
     $            0, 0, ICTXT, LLDTMP, INFO )
               CALL PDGEMM( 'N', 'N', KLN, NW, NW, ONE, Z, KKROW,
     $            ILO, DESCZ, V, 1, 1, DESCUV, ZERO, DWORK, 1,
     $            1, DESCWV )
               CALL PDGEMR2D( KLN, NW, DWORK, 1, 1, DESCWV, Z,
     $            KKROW, ILO, DESCZ, ICTXT )
 260        CONTINUE
         END IF       
      END IF       

*     End of PDHGEQZ6
      END
