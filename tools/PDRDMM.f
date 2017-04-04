********************************************************************
*                                                                  *
* PDRDMM - Routine for parsing an open Matrix Market file          *
*                                                                  *
********************************************************************
      SUBROUTINE PDRDMM( NIN, REP, FIELD, SYMM, NROWS, NCOLS, 
     $     TOTNNZ, A, DESCA, LDWORK, DWORK, INFO ) 

      IMPLICIT NONE
*     Scalar arguments
      INTEGER           NIN, NROWS, NCOLS, LDWORK, INFO
      INTEGER*8         TOTNNZ
      CHARACTER         REP*10
      CHARACTER         FIELD*7
      CHARACTER         SYMM*19

*     Arrays arguments
      DOUBLE PRECISION  A( * ), DWORK( * )
      INTEGER           DESCA(*)
      

*     .. Parameters
      DOUBLE PRECISION  ZERO, ONE
      PARAMETER         ( ZERO = 0.0D0, ONE = 1.0D0 )

      INTEGER            BLOCK_CYCLIC_2D, CSRC_, CTXT_, DLEN_, DTYPE_,
     $                   LLD_, MB_, M_, NB_, N_, RSRC_
      PARAMETER          ( BLOCK_CYCLIC_2D = 1, DLEN_ = 9, DTYPE_ = 1,
     $                     CTXT_ = 2, M_ = 3, N_ = 4, MB_ = 5, NB_ = 6,
     $                     RSRC_ = 7, CSRC_ = 8, LLD_ = 9 )


*     Local Scalars

      INTEGER           N, I
      INTEGER*8         READNNZ
      INTEGER           BLOCKNNZ
      INTEGER           ICTXT, NPROW, NPCOL, MYROW, MYCOL
      INTEGER           IAM, NPROCS
      INTEGER           NB


*     Local arrays
      COMPLEX           CVAL( 1 )
      INTEGER           IVAL( 1 )
      INTEGER, ALLOCATABLE :: INDX( : ), JNDX ( : )

      INFO = 0
*     Extract context and gridinformation
      ICTXT = DESCA( CTXT_ )            
      CALL BLACS_GRIDINFO( ICTXT, NPROW, NPCOL, MYROW, MYCOL )
      CALL BLACS_PINFO(IAM, NPROCS)

      NB = DESCA( NB_ )
      N = NROWS

      BLOCKNNZ = NB * NB
*     Allocate storage for index pointers
      ALLOCATE( INDX( BLOCKNNZ ) )
      ALLOCATE( JNDX( BLOCKNNZ ) )

      CALL PDLASET('A', N, N, ZERO, ZERO, A, 1, 1, DESCA )

*     process the A file in chunks of NB*NB
      READNNZ = 0
      DO WHILE ( READNNZ .LT. TOTNNZ )

         BLOCKNNZ = min(BLOCKNNZ, TOTNNZ-READNNZ)
         
         CALL mmread2( NIN, REP, FIELD, SYMM, NROWS, NCOLS,
     $        BLOCKNNZ, BLOCKNNZ,
     $        INDX, JNDX, IVAL, DWORK, CVAL)
         DO I = 1, BLOCKNNZ
            CALL PDELSET( A, INDX(I), JNDX(I), DESCA, DWORK(I))
         END DO

         IF (SYMM .EQ. 'symmetric') THEN
            DO I = 1, BLOCKNNZ
               CALL PDELSET( A, JNDX(I), INDX(I), DESCA, DWORK(I))
            END DO
         END IF

         READNNZ = READNNZ + BLOCKNNZ

      END DO
      
      DEALLOCATE( INDX, JNDX )


C
      END
