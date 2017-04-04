***********************************************************************
*                                                                     *
*     PDLACP4.f:                                                      *
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
      SUBROUTINE PDLACP4( J, ILAST, WSIZE, H, T, DESC_HT, TMPH, TMPT, 
     $   LDTMP_HT, TYPE, DIR )
      IMPLICIT NONE
*     ..
*     .. Scalar Arguments ..
*     ..
      INTEGER             J, WSIZE, ILAST, LDTMP_HT, DIR, TYPE
*     ..
*     .. Array Arguments ..
*     ..
      DOUBLE PRECISION    T( * ) , H( * ), TMPH( LDTMP_HT, * ), 
     $   TMPT( LDTMP_HT, * )
      INTEGER             DESC_HT( * )
*     
*     Purpose
*     =======
*     PDLACP4 either copies or restores the submatrix (J:J+WSIZE-1, J:J+WSIZE-1) of 
*     H and T using TMPH and TMPT as destination or source. DIR = 0 copies to
*     TMPH and TMPT. If DIR = 1 TMPH and TMPT are stored in H and T. 
*     
*     Only the diagonal processes (at most 2) will receive data in TMPH and TMPT.   
*     
*     The submatrices of H and T are shared, by at most, by 4 processes.
*     H and T are assumed to have the same distribution.      
*     
*     If bulges, 3x3 at most, are present on the diagonal of H and T,
*     this is indicated by setting TYPE=1, 
*     else H is assumed to be upper Hessenberg, and T upper triangular.  
*     
*     Auxiliary routine - all inputs are assumed valid without checking 
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
*     J       (global input) INTEGER
*     Holds the row/column index of H and T where to copy from.
*     
*     ILAST   (global input) INTEGER
*     Holds last row/column of H and T to work with.
*     
*     WSIZE   (global input) INTEGER
*     Holds the size of the window to copy/restore.
*     
*     H       (local input/output) DOUBLE PRECISION array, dimension (LLD_HT, LOCc(N)). 
*     On entry, H contains the upper Hessenberg matrix H to copy/restore data
*     from.
*     
*     T       (local input/output) DOUBLE PRECISION array, dimension (LLD_HT, LOCc(N)). 
*     On entry, T contains the upper triangular matrix T to copy/restore data
*     from.
*     
*     DESC_HT (global and local input) INTEGER array of dimension DLEN_.
*     The array descriptor for the distributed matrices H and T.
*     
*     TMPH   (local input/output) DOUBLE PRECISION array, dimension (LDTMP_HT,*)
*     Local storage to either store or restore part of matrix H.
*     
*     TMPT   (local input/output) DOUBLE PRECISION array, dimension (LDTMP_HT,*)
*     Local storage to either store or restore part of matrix T.
*     
*     LDTMP_HT (local input) INTEGER
*     Leading dimension for the matrices TMPH and TMPT
*     
*     Type   (global input) INTEGER
*     =1 : Hessenberg/triangular matrices with "bumps" are assumed, 
*     otherwise regular Hessenberg/triangular matrices are assumed.
*     
*     Dir    (global input) INTEGER
*     = 0 : Store parts of distributed H,T into TMPH, TMPT.
*     = 1 : Restore parts of distributed H,T from TMPH, TMPT.
*     
*     =========================================
*     
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

*     ..
*     .. Local Arrays ..
*     ..      


*     ..
*     .. Local Scalars ..
*     ..
      INTEGER             ICTXT, MYROW, MYCOL, NPROW, NPCOL, LDH, LDT, 
     $   NB, RSRC1, CSRC1, RSRC2, CSRC2, RSRC3, CSRC3, RSRC4, CSRC4, 
     $   ILOC, JLOC, DIM1, DIM4, RWS3, CLS3
      LOGICAL             PARTOWNER, HASBULGES
*     ..
*     .. Intrinsic Functions ..
*     ..
      INTRINSIC           ABS, DBLE, MAX, MIN, SQRT, MOD
*     ..
*     .. Externals
*     ..
      INTEGER             INDXG2L, INDXG2P
      EXTERNAL            INDXG2L, INDXG2P



*     
*     Extract current communication context and get grid parameters
*     
      ICTXT = DESC_HT( CTXT_ )
      CALL BLACS_GRIDINFO( ICTXT, NPROW, NPCOL, MYROW, MYCOL )


*     
*     NB_ .EQ. MB_ is assumed, we use NB as blocking factor.
*     
      NB = DESC_HT(NB_)
      
      
      
      NB = DESC_HT(NB_)
      LDH = DESC_HT(LLD_)
      LDT = DESC_HT(LLD_)

      HASBULGES = TYPE .EQ . 1
      
*     Extract processors for current window,
*     as below:
*     
*     1 | 2
*     -----
*     3 | 4
*     
      RSRC1 = INDXG2P( J, NB, -1, DESC_HT(RSRC_), NPROW)
      CSRC1 = INDXG2P( J, NB, -1, DESC_HT(CSRC_), NPCOL)
      RSRC2 = RSRC1
      CSRC2 = INDXG2P( J + WSIZE - 1, NB, -1, DESC_HT(CSRC_), 
     $   NPCOL)
      RSRC3 = INDXG2P( J + WSIZE - 1, NB, -1, DESC_HT(RSRC_), 
     $   NPROW)
      CSRC3 = CSRC1
      RSRC4 = RSRC3
      CSRC4 = CSRC2
      
      
      PARTOWNER = ( MYROW .EQ. RSRC1. AND. MYCOL. EQ. CSRC1 ) .OR.
     $   ( MYROW .EQ. RSRC2 .AND. MYCOL .EQ. CSRC2) .OR.
     $   ( MYROW .EQ. RSRC3 .AND. MYCOL .EQ. CSRC3) .OR.
     $   ( MYROW .EQ. RSRC4 .AND. MYCOL .EQ. CSRC4)
      
      DIM1 = NB - MOD( J - 1, NB )
      IF ( J + DIM1 -1 .GE. ILAST ) THEN
         DIM1 = WSIZE
      END IF
      
      DIM4 = WSIZE - DIM1
      IF ( NPROW * NPCOL .EQ. 1 ) THEN
         DIM1 = WSIZE
         DIM4 = 0
      END IF

      RWS3 = 1
      CLS3 = 1
      IF ( HASBULGES ) THEN
         RWS3 = MIN( 3, DIM4 )
         CLS3 = MIN( 3, DIM1 )
      END IF
      
      IF ( .NOT. PARTOWNER ) GOTO 200
*     DIR.EQ.0 means store H,T in TMPH, TMPT      
*     DIR.EQ.1 means restore  
*     If restore, jump to that part of code      
*     
      IF ( DIR .EQ. 1 ) GOTO 150

      
*     
*     Store distributed H and T in TMPH and TMPT
*     
      CALL DLASET( 'A', WSIZE, WSIZE, ZERO, ZERO, TMPH,
     $   LDTMP_HT )
      CALL DLASET( 'A', WSIZE, WSIZE, ZERO, ZERO, TMPT,
     $   LDTMP_HT )    

*     Let process 1 and 4 exchange data
      IF ( MYROW .EQ. RSRC1 .AND. MYCOL .EQ. CSRC1 ) THEN
         ILOC = INDXG2L( J, NB, MYROW, DESC_HT( RSRC_ ), NPROW )
         JLOC = INDXG2L( J, NB, MYCOL, DESC_HT( CSRC_ ), NPCOL )
         CALL DLACPY( 'All', DIM1, DIM1,
     $      H( ( JLOC - 1 ) * LDH + ILOC ), LDH, TMPH,
     $      LDTMP_HT )
         CALL DLACPY( 'All', DIM1, DIM1,
     $      T( ( JLOC - 1 ) * LDT + ILOC ), LDT, TMPT,
     $      LDTMP_HT )            
         IF ( DIM4 .EQ. 0 ) GOTO 100
         IF( RSRC1 .NE. RSRC4 .OR. CSRC1 .NE. CSRC4 ) THEN
*           Exch #1 <==> #4
            CALL DGESD2D( ICTXT, DIM1, DIM1,
     $         TMPH, LDTMP_HT, RSRC4, CSRC4 )
            CALL DGESD2D( ICTXT, DIM1, DIM1,
     $         TMPT, LDTMP_HT, RSRC4, CSRC4 )

            CALL DGERV2D( ICTXT, DIM4, DIM4,
     $         TMPH( DIM1 + 1, DIM1 + 1 ),
     $         LDTMP_HT, RSRC4, CSRC4 )
            CALL DGERV2D( ICTXT, DIM4, DIM4,
     $         TMPT( DIM1 + 1, DIM1 + 1 ),
     $         LDTMP_HT, RSRC4, CSRC4 )      
         END IF
      END IF
      
      IF ( DIM4 .EQ. 0 ) GOTO 100
      
      IF ( MYROW.EQ.RSRC4 .AND. MYCOL.EQ.CSRC4 ) THEN
         ILOC = INDXG2L( J + DIM1, NB, MYROW,
     $      DESC_HT( RSRC_ ), NPROW )
         JLOC = INDXG2L( J + DIM1, NB, MYCOL,
     $      DESC_HT( CSRC_ ), NPCOL )    
         CALL DLACPY( 'All', DIM4, DIM4,
     $      H( ( JLOC - 1 ) * LDH + ILOC ), LDH,
     $      TMPH( DIM1 + 1, DIM1 + 1 ),
     $      LDTMP_HT )
         CALL DLACPY( 'All', DIM4, DIM4,
     $      T( ( JLOC - 1 ) * LDT + ILOC ), LDT,
     $      TMPT( DIM1 + 1, DIM1 + 1 ),
     $      LDTMP_HT )
         
         IF( RSRC4.NE.RSRC1 .OR. CSRC4.NE.CSRC1 ) THEN
*           Exch #4 <==> #1
            CALL DGESD2D( ICTXT, DIM4, DIM4,
     $         TMPH( DIM1 + 1, DIM1 + 1 ),
     $         LDTMP_HT, RSRC1, CSRC1 )
            CALL DGESD2D( ICTXT, DIM4, DIM4,
     $         TMPT( DIM1 + 1, DIM1 + 1 ),
     $         LDTMP_HT, RSRC1, CSRC1 )

            CALL DGERV2D( ICTXT, DIM1, DIM1,
     $         TMPH, LDTMP_HT, RSRC1, CSRC1 )
            
            CALL DGERV2D( ICTXT, DIM1, DIM1,
     $         TMPT, LDTMP_HT, RSRC1, CSRC1 )   
         END IF
         
      END IF

*     Retrive data from process 2 and 3, send to 1 and 4
      IF ( MYROW.EQ.RSRC2 .AND. MYCOL.EQ.CSRC2 ) THEN          
         ILOC = INDXG2L( J, NB, MYROW,
     $      DESC_HT( RSRC_ ), NPROW )
         JLOC = INDXG2L( J+DIM1, NB, MYCOL,
     $      DESC_HT( CSRC_ ), NPCOL )
         CALL DLACPY( 'All', DIM1, DIM4,
     $      H( ( JLOC - 1 ) * LDH + ILOC ), LDH,
     $      TMPH( 1, DIM1 + 1 ), LDTMP_HT )
         CALL DLACPY( 'All', DIM1, DIM4,
     $      T( ( JLOC - 1 ) * LDT + ILOC ), LDT,
     $      TMPT( 1, DIM1 + 1 ), LDTMP_HT )   
         
         IF( RSRC2.NE.RSRC1 .OR. CSRC2.NE.CSRC1 ) THEN
*           Send #2 ==> #1
            CALL DGESD2D( ICTXT, DIM1, DIM4,
     $         TMPH( 1, DIM1 + 1 ),
     $         LDTMP_HT, RSRC1, CSRC1 )   
            CALL DGESD2D( ICTXT, DIM1, DIM4,
     $         TMPT( 1, DIM1 + 1 ),
     $         LDTMP_HT, RSRC1, CSRC1 )                   
         END IF
         IF( RSRC2.NE.RSRC4 .OR. CSRC2.NE.CSRC4 ) THEN
*           Send #2 ==> #4 
            CALL DGESD2D( ICTXT, DIM1, DIM4,
     $         TMPH( 1, DIM1+ 1 ),
     $         LDTMP_HT, RSRC4, CSRC4 )
            CALL DGESD2D( ICTXT, DIM1, DIM4,
     $         TMPT( 1, DIM1 + 1 ),
     $         LDTMP_HT, RSRC4, CSRC4 )                  
         END IF
      END IF

      IF ( MYROW .EQ. RSRC3 .AND. MYCOL .EQ. CSRC3 ) THEN
         ILOC = INDXG2L( J + DIM1, NB, MYROW,
     $      DESC_HT (RSRC_ ), NPROW )
         JLOC = INDXG2L( J + DIM1 - CLS3, NB, MYCOL,
     $      DESC_HT( CSRC_ ), NPCOL )
         CALL DLACPY( 'All', RWS3, CLS3,
     $      H( ( JLOC - 1) * LDH + ILOC ), LDH,
     $      TMPH( DIM1 + 1, DIM1-CLS3+1),
     $      LDTMP_HT )    
         CALL DLACPY( 'All', RWS3, CLS3,
     $      T( ( JLOC - 1) * LDT + ILOC ), LDT,
     $      TMPT( DIM1 + 1, DIM1 - CLS3 + 1 ),
     $      LDTMP_HT )                 
         IF( RSRC3.NE.RSRC1 .OR. CSRC3.NE.CSRC1 ) THEN
*           Send #3 ==> #1
            CALL DGESD2D( ICTXT, RWS3, CLS3,
     $         TMPH ( DIM1 + 1, DIM1 - CLS3 + 1 ),
     $         LDTMP_HT, RSRC1, CSRC1 )
            IF (HASBULGES) THEN
               CALL DGESD2D( ICTXT, RWS3, CLS3,
     $            TMPT( DIM1 + 1, DIM1 - CLS3 + 1 ),
     $            LDTMP_HT, RSRC1, CSRC1 )                  
            END IF
         END IF
         
         IF( RSRC3 .NE. RSRC4 .OR. CSRC3 .NE. CSRC4 ) THEN
*           Send #3 ==> #4
            CALL DGESD2D( ICTXT, RWS3, CLS3,
     $         TMPH( DIM1 + 1, DIM1-CLS3+1),
     $         LDTMP_HT, RSRC4, CSRC4 )
            IF ( HASBULGES ) THEN
               CALL DGESD2D( ICTXT, RWS3, CLS3,
     $            TMPT( DIM1 + 1, DIM1 - CLS3 + 1 ),
     $            LDTMP_HT, RSRC4, CSRC4 ) 
            END IF
         END IF
         
      END IF

      IF( MYROW .EQ. RSRC1 .AND. MYCOL .EQ. CSRC1 ) THEN
         IF( RSRC1 .NE. RSRC2 .OR. CSRC1 .NE. CSRC2 ) THEN
*           Recv #1 <== #2
            CALL DGERV2D( ICTXT, DIM1, DIM4,
     $         TMPH( 1,  DIM1 + 1 ),
     $         LDTMP_HT, RSRC2, CSRC2 )
            CALL DGERV2D( ICTXT, DIM1, DIM4,
     $         TMPT( 1 ,DIM1 + 1 ),
     $         LDTMP_HT, RSRC2, CSRC2 ) 
         END IF
         IF( RSRC1 .NE. RSRC3 .OR. CSRC1 .NE. CSRC3 ) THEN
*           Recv #1 <== #3
            CALL DGERV2D( ICTXT, RWS3, CLS3,
     $         TMPH( DIM1 + 1, DIM1 - CLS3 + 1 ),
     $         LDTMP_HT, RSRC3, CSRC3 )
            IF ( HASBULGES ) THEN
               CALL DGERV2D( ICTXT, RWS3, CLS3,
     $            TMPT( DIM1 + 1, DIM1 - CLS3 + 1 ),
     $            LDTMP_HT, RSRC3, CSRC3 )
            END IF

         END IF

      END IF
      
      IF( MYROW .EQ. RSRC4 .AND. MYCOL .EQ. CSRC4 ) THEN
         IF( RSRC4 .NE. RSRC2 .OR. CSRC4 .NE. CSRC2 ) THEN
*           Recv #4 <== #2  
            CALL DGERV2D( ICTXT, DIM1, DIM4,
     $         TMPH( 1, DIM1 + 1 ),
     $         LDTMP_HT, RSRC2, CSRC2 )
            CALL DGERV2D( ICTXT, DIM1, DIM4,
     $         TMPT( 1, DIM1 + 1),
     $         LDTMP_HT, RSRC2, CSRC2 )    
         END IF
         IF( RSRC4 .NE. RSRC3 .OR. CSRC4 .NE. CSRC3 ) THEN
*           Recv #4 <== #3
            CALL DGERV2D( ICTXT, RWS3, CLS3,
     $         TMPH( DIM1 + 1, DIM1 - CLS3 + 1 ),
     $         LDTMP_HT, RSRC3, CSRC3 )     
            IF ( HASBULGES ) THEN
               CALL DGERV2D( ICTXT, RWS3, CLS3,
     $            TMPT( DIM1 + 1, DIM1 - CLS3 + 1 ),
     $            LDTMP_HT, RSRC3, CSRC3 )     
            END IF
         END IF
      END IF
      
 100  CONTINUE  
      
      GOTO 200

*     Restore H and T using TMPH and TMPT as input
 150  CONTINUE
      
*     Restore parts of H and T using process 1 and 4
      IF( MYROW .EQ. RSRC1 .AND. MYCOL .EQ. CSRC1 ) THEN
         ILOC = INDXG2L( J, NB, MYROW,
     $      DESC_HT( RSRC_ ), NPROW )
         JLOC = INDXG2L( J, NB, MYCOL,
     $      DESC_HT( CSRC_ ), NPCOL )
         
         CALL DLACPY( 'All', DIM1, DIM1, TMPH,
     $      LDTMP_HT, H( ( JLOC - 1 ) * LDH + ILOC ),
     $      LDH )
         CALL DLACPY( 'All', DIM1, DIM1, TMPT,
     $      LDTMP_HT, T( ( JLOC - 1 ) * LDT + ILOC ),
     $      LDT )              
      END IF
      
      IF ( DIM4 .EQ. 0) GOTO 200
      
      IF ( MYROW .EQ. RSRC4 .AND. MYCOL .EQ. CSRC4 ) THEN
         ILOC = INDXG2L( J + DIM1, NB, MYROW,
     $      DESC_HT( RSRC_ ), NPROW )
         JLOC = INDXG2L( J + DIM1, NB, MYCOL,
     $      DESC_HT( CSRC_ ), NPCOL )
         CALL DLACPY( 'All', DIM4, DIM4,
     $      TMPH( DIM1 + 1, DIM1 + 1 ),
     $      LDTMP_HT, H( ( JLOC - 1 ) * LDH + ILOC ),
     $      LDH )
         CALL DLACPY( 'All', DIM4, DIM4,
     $      TMPT( DIM1 + 1, DIM1 + 1 ),
     $      LDTMP_HT, T( ( JLOC - 1 ) * LDT + ILOC ),
     $      LDT )              
      END IF
      
      

*     Return data to process 2 and 3 and update H and T
      RWS3 = MIN( 3, DIM4 )
      CLS3 = MIN( 3, DIM1 )
      IF ( MYROW .EQ. RSRC1 .AND. MYCOL .EQ. CSRC1 ) THEN
         IF ( RSRC1 .NE. RSRC3 .OR. CSRC1 .NE. CSRC3 ) THEN
*           Send #1 ==> #3                  
            CALL DGESD2D( ICTXT, RWS3, CLS3,
     $         TMPH( DIM1 + 1, DIM1 - CLS3 + 1 ),
     $         LDTMP_HT, RSRC3, CSRC3 )
            IF ( HASBULGES ) THEN
               CALL DGESD2D( ICTXT, RWS3, CLS3,
     $            TMPT( DIM1 + 1, DIM1 - CLS3 + 1 ),
     $            LDTMP_HT, RSRC3, CSRC3 )                  
            END IF
         END IF
      END IF

      IF ( MYROW .EQ. RSRC4 .AND. MYCOL .EQ. CSRC4 ) THEN
         IF( RSRC4 .NE. RSRC2 .OR. CSRC4 .NE. CSRC2 ) THEN
*           Send #4 ==> #2
            CALL DGESD2D( ICTXT, DIM1, DIM4,
     $         TMPH(1, DIM1 + 1 ),
     $         LDTMP_HT, RSRC2, CSRC2 )
            CALL DGESD2D( ICTXT, DIM1, DIM4,
     $         TMPT( 1, DIM1 + 1 ),
     $         LDTMP_HT, RSRC2, CSRC2 )
         END IF
      END IF

      IF ( MYROW .EQ. RSRC2 .AND. MYCOL .EQ. CSRC2 ) THEN
         ILOC = INDXG2L( J, NB, MYROW,
     $      DESC_HT( RSRC_ ), NPROW )
         JLOC = INDXG2L( J + DIM1, NB, MYCOL,
     $      DESC_HT( CSRC_ ), NPCOL )
         IF( RSRC2 .NE. RSRC4 .OR. CSRC2 .NE. CSRC4 ) THEN
*           Recv #2 <== #4
            CALL DGERV2D( ICTXT, DIM1, DIM4,
     $         TMPH( 1, DIM1 + 1 ),
     $         LDTMP_HT, RSRC4, CSRC4 )
            CALL DGERV2D( ICTXT, DIM1, DIM4,
     $         TMPT(1, DIM1 + 1 ),
     $         LDTMP_HT, RSRC4, CSRC4 )   
         END IF
         CALL DLACPY( 'All', DIM1, DIM4,
     $      TMPH( 1, DIM1 + 1 ), LDTMP_HT,
     $      H( ( JLOC - 1 ) * LDH + ILOC ), LDH )
         CALL DLACPY( 'All', DIM1, DIM4,
     $      TMPT( 1, DIM1 + 1 ), LDTMP_HT,
     $      T( ( JLOC -1 )* LDT + ILOC ), LDT )              
      END IF                  

      IF( MYROW .EQ. RSRC3 .AND. MYCOL .EQ. CSRC3 ) THEN
         ILOC = INDXG2L( J + DIM1, NB, MYROW,
     $      DESC_HT( RSRC_ ), NPROW )
         JLOC = INDXG2L( J + DIM1 - CLS3, NB, MYCOL,
     $      DESC_HT( CSRC_ ), NPCOL )
         IF( RSRC3 .NE. RSRC1 .OR. CSRC3 .NE. CSRC1 ) THEN
*           Recv #3 <== #1
            CALL DGERV2D( ICTXT, RWS3, CLS3,
     $         TMPH( DIM1 + 1, DIM1 - CLS3 + 1 ),
     $         LDTMP_HT, RSRC1, CSRC1 )
            IF ( HASBULGES ) THEN
               CALL DGERV2D( ICTXT, RWS3, CLS3,
     $            TMPT( DIM1 + 1, DIM1 - CLS3 + 1 ),
     $            LDTMP_HT, RSRC1, CSRC1 )                  
            END IF
         END IF
         CALL DLACPY( 'A', RWS3, CLS3,
     $      TMPH( DIM1 + 1,DIM1 - CLS3 + 1 ),
     $      LDTMP_HT, H( ( JLOC - 1) * LDH + ILOC ),
     $      LDH )
         IF (HASBULGES) THEN
            CALL DLACPY( 'A', RWS3, CLS3,
     $         TMPT( DIM1 + 1,DIM1 - CLS3 + 1 ),
     $         LDTMP_HT, T( ( JLOC - 1 ) * LDT + ILOC ),
     $         LDT )     
         END IF

      END IF

 200  CONTINUE

      
      
      RETURN
*     End of PDLACP4

      END

