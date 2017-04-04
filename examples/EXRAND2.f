***********************************************************************
*                                                                     *
*     Example driver for ScaLAPACK-style routine PDHGEQZ              *
*     (random matrix pair)                                            *
*                                                                     *
*     Authors: Bjorn Adlerborn                                        *
*              Bo Kagstrom                                            *
*              Daniel Kressner                                        *
*                                                                     *
*     Department of Computing Science and HPC2N, Umea University      *
*     MATHICSE ANCHP, EPF Lausanne                                    *
*                                                                     *
***********************************************************************
*
      PROGRAM EXRAND2
*     ..
*     .. Declarations
*     ..
      IMPLICIT NONE
*
*     .. Parameters
      INTEGER            BLOCK_CYCLIC_2D, CSRC_, CTXT_, DLEN_, DTYPE_,
     $     LLD_, MB_, M_, NB_, N_, RSRC_
      PARAMETER          ( BLOCK_CYCLIC_2D = 1, DLEN_ = 9, DTYPE_ = 1,
     $     CTXT_ = 2, M_ = 3, N_ = 4, MB_ = 5, NB_ = 6,
     $     RSRC_ = 7, CSRC_ = 8, LLD_ = 9 )

      INTEGER           N
      PARAMETER         (
*     Problem size.
     $                    N = 1000 )

      DOUBLE PRECISION    ZERO, ONE
      PARAMETER           ( ZERO = 0.D0, ONE = 1.0D+0 )   
*     
*     .. Local Scalars
      DOUBLE PRECISION    ANORM, BNORM, QNORM, ZNORM
      INTEGER*8           LDWORK, LIWORK
      INTEGER             ICTXT, NPROW, NPCOL, MYROW, MYCOL
      INTEGER             IAM, SYSPROCS, NPROCS
      INTEGER             ILO, IHI, NP, NQ, NB
      INTEGER             I, J, INFO 
      DOUBLE PRECISION    EPS
      DOUBLE PRECISION    TOTTIME, QZTIME, RESTIME
      DOUBLE PRECISION    AB2HTTIME, TIMEDORMQR, TIMEQRF
      DOUBLE PRECISION    TRANSPTIME

*     
*     .. Local Arrays 
      INTEGER             DESCA( DLEN_ ), DESCB( DLEN_ )
      INTEGER             DESCQ( DLEN_ ), DESCZ( DLEN_ )
           
      DOUBLE PRECISION    DDUM(1)
      
*     Placeholders for the computed eigenvalues
      DOUBLE PRECISION, ALLOCATABLE :: AR(:)
      DOUBLE PRECISION, ALLOCATABLE :: AI(:) 
      DOUBLE PRECISION, ALLOCATABLE :: BETA(:)
*     Need two copies of A and B to be able to perform residual check 
*     after computation.     
      DOUBLE PRECISION, ALLOCATABLE :: A1(:), A2(:)
      DOUBLE PRECISION, ALLOCATABLE :: B1(:), B2(:)
      DOUBLE PRECISION, ALLOCATABLE :: Q1(:), Q2(:)
      DOUBLE PRECISION, ALLOCATABLE :: Z1(:)
*     Workspace
      INTEGER, ALLOCATABLE :: IWORK(:)
      DOUBLE PRECISION, ALLOCATABLE :: DWORK(:)  
*     
*     .. Externals 
      EXTERNAL          CHISQ, MPI_WTIME, PDLAMCH
      DOUBLE PRECISION  CHISQ, MPI_WTIME, PDLAMCH
      EXTERNAL          NUMROC
      INTEGER           NUMROC

*     
*     .. Intrinsic Functions 
      INTRINSIC           DBLE, MAX

       
*     
*     .. Executable Statements


      INFO = 0


*     Get number of participating proceccsors spawned by mpirun
*     and form the grid

      CALL BLACS_PINFO( IAM, SYSPROCS )
      NPROW = INT( SQRT( DBLE(SYSPROCS) ) )
      NPCOL = SYSPROCS / NPROW
      CALL BLACS_GET( 0, 0, ICTXT )
      CALL BLACS_GRIDINIT( ICTXT, '2D', NPROW, NPCOL )
      CALL BLACS_GRIDINFO( ICTXT, NPROW, NPCOL, MYROW, MYCOL )
      NPROCS = NPROW * NPCOL
*     Configure parameters for the underlying KKQZ 
      EPS = PDLAMCH( ICTXT, 'PRECISION' )
      CALL KKQZCONF( EPS )

      TOTTIME = MPI_WTIME()

*     Hard coded workspace for integers, needed by infinite eigenvalue computation
*     as well as reordering.
*     N more than 2 * N should be needed at any point
      LIWORK = 4 * N 
      ALLOCATE( IWORK( LIWORK ), STAT = INFO )
      IF( INFO.NE.0 ) THEN
          WRITE(*,*) 'Could not allocate IWORK. INFO = ', INFO
          GOTO 777 
      END IF  

*     Initialize the random seed generator
      CALL RANDOM_SEED ()
 
*     Define the blocking factor, select depending upon N
      IF ( N .LE. 2000 ) THEN
          NB = 60
      ELSE
          NB = 130
      END IF


*     Print welcoming message.

      IF( IAM .EQ. 0 ) THEN
         WRITE(*,*)
         WRITE(*,*) '%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%'
         WRITE(*,*) '%%      EXAMPLE PROGRAM 2 FOR PDHGEQZ       %%'
         WRITE(*,*) '%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%'
         WRITE(*,*)
      END IF

*     We dont use any balancing routine so set ILO to 1 and IHI to N
      ILO = 1 
      IHI = N

      NP = NUMROC( N, NB, MYROW, 0, NPROW )
      NQ = NUMROC( N, NB, MYCOL, 0, NPCOL )
*     Initialize the array descriptor for the matrix A, B
      CALL DESCINIT( DESCA, N, N, NB, NB, 0, 0, ICTXT, 
     $     MAX( 1, NP ), INFO )
      CALL DESCINIT( DESCB, N, N, NB, NB, 0, 0, ICTXT, 
     $     MAX( 1, NP ), INFO )          

*     Initialize the array descriptor for the matrix Q, Z
      CALL DESCINIT( DESCQ, N, N, NB, NB, 0, 0, ICTXT, 
     $     MAX( 1, NP ), INFO )     
      CALL DESCINIT( DESCZ, N, N, NB, NB, 0, 0, ICTXT, 
     $     MAX( 1, NP ), INFO ) 
         
      IF ( IAM .EQ. 0 ) THEN
         WRITE(*,*) '%% Allocating memory for matrices and eigenvalues'
      END IF
*     Allocate needed space for the matrices A*2,B*2,Q and Z.
      ALLOCATE ( A1( DESCA( LLD_ )*NQ), STAT = INFO )
      IF( INFO.NE.0 ) THEN
           WRITE(*,*) '% Could not allocate A1. INFO = ', INFO
           GOTO 777 
      END IF           
           
      ALLOCATE ( B1( DESCA( LLD_ )*NQ), STAT = INFO )
      IF( INFO.NE.0 ) THEN
          WRITE(*,*) '% Could not allocate B1. INFO = ', INFO
          GOTO 777 
      END IF           
          
      ALLOCATE ( Q1( DESCA( LLD_ )*NQ), STAT = INFO )
      IF( INFO.NE.0 ) THEN
           WRITE(*,*) '% Could not allocate Q1. INFO = ', INFO
           GOTO 777 
      END IF    

      ALLOCATE ( Q2( DESCA( LLD_ )*NQ), STAT = INFO )
      IF( INFO.NE.0 ) THEN
           WRITE(*,*) '% Could not allocate Q2. INFO = ', INFO
           GOTO 777
      END IF

          
      ALLOCATE ( Z1( DESCA( LLD_ )*NQ), STAT = INFO )
      IF( INFO.NE.0 ) THEN
           WRITE(*,*) '% Could not allocate Z1. INFO = ', INFO
           GOTO 777 
      END IF    

      ALLOCATE ( A2( DESCA( LLD_ )*NQ), STAT = INFO )
      IF( INFO.NE.0 ) THEN
           WRITE(*,*) '% Could not allocate A2. INFO = ', INFO
           GOTO 777 
      END IF           
          
      ALLOCATE ( B2( DESCA( LLD_ )*NQ), STAT = INFO )
      IF( INFO.NE.0 ) THEN
           WRITE(*,*) '% Could not allocate B2. INFO = ', INFO
           GOTO 777 
      END IF               
              
         
*     Allocate space for eigenvalues
      ALLOCATE ( AR( N ), STAT = INFO )
      IF( INFO.NE.0 ) THEN
            WRITE(*,*) '% Could not allocate AR. INFO = ', INFO
           GOTO 777 
      END IF
      ALLOCATE ( AI( N ), STAT = INFO )
      IF( INFO.NE.0 ) THEN
           WRITE(*,*) '% Could not allocate AI. INFO = ', INFO
           GOTO 777 
      END IF
      ALLOCATE ( BETA( N ), STAT = INFO )
      IF( INFO.NE.0 ) THEN
           WRITE(*,*) '% Could not allocate BETA. INFO = ', INFO
           GOTO 777 
      ENDIF
*     BETA should never be negativ, so set a negtive value for 
*     being 'non set'.
      BETA ( 1 : N ) = -1 
*
      LDWORK = 0
*     Perform a workspace query for PDGEQRF
      CALL PDGEQRF(N, N, B1, 1, 1, DESCB, DDUM,
     $     DDUM, -1, INFO)
      LDWORK = MAX(LDWORK, INT(DDUM(1)))

      CALL PDORMQR( 'L', 'T', N, N, N, B1, 1, 1, DESCB,
     $     DDUM, A1, 1, 1, DESCA, DDUM, -1,
     $     INFO ) 
      LDWORK = MAX(LDWORK, INT(DDUM(1)+N))

*     Perform a workspace query to PDGGHRD
      CALL PDGGHRD( 'I', 'I', N, ILO, IHI,
     $     A1, DESCA, B1, DESCB,
     $     Q1, DESCQ, Z1, DESCZ,
     $     DDUM, -1, INFO )
      LDWORK = MAX(LDWORK, INT(DDUM(1)+N))

          
*     Perform a workspace query to PDHGEQZ      
      CALL PDHGEQZ( 'S', 'V', 'V', N, ILO, IHI, 
     $     A1, DESCA, B1, DESCB,
     $     AR, AI, BETA,
     $     Q1, DESCQ, Z1, DESCZ, 
     $     DDUM, -1,
     $     IWORK, LIWORK, INFO )
      LDWORK = MAX(LDWORK, INT(DDUM(1)))

         
*     Make sure all processors allocate the same amount of work          
      CALL IGAMX2D( ICTXT, 'All', ' ', 1, 1, LDWORK, 1, -1, -1, -1,
     $     -1, -1 )   
           
      IF ( IAM .EQ. 0 ) WRITE(*,*)'%% Allocating workspace :' , 
     $     LDWORK
*     Allocate workspace          
      ALLOCATE ( DWORK( LDWORK ), STAT = INFO )
      IF( INFO.NE.0 ) THEN
          WRITE(*,*) 'Could not allocate Workspace. INFO = ', INFO
          GOTO 777 
      END IF 

      IF ( IAM .EQ. 0) WRITE(*,*)'%% Done Allocating.'
      CALL FLUSH(6)
    
*     Generate Data
                               
      IF ( IAM .EQ. 0 ) WRITE(*,*)'%% Generating data'


      CALL PDMATGEN2( ICTXT, 'Random', 'NoDiagDominant',
     $     DESCA( N_ ), DESCA( N_ ),
     $     DESCA( NB_ ), DESCA( NB_), A1,
     $     DESCA( LLD_ ), 0, 0, 6,
     $     0, NP, 0, NQ, MYROW, MYCOL, NPROW, NPCOL )
           
      CALL PDMATGEN2( ICTXT, 'Random', 'NoDiagDominant',
     $     DESCA( N_ ), DESCA( N_ ),
     $     DESCA( NB_ ), DESCA( NB_), B1,
     $     DESCA( LLD_ ), 0, 0, 7,
     $     0, NP, 0, NQ, MYROW, MYCOL, NPROW, NPCOL )
      
*     Initiliaze Q and Z to the identity      
      CALL PDLASET( 'Full', N, N, ZERO, ONE, Q1, 1, 1, DESCQ )
      CALL PDLASET( 'Full', N, N, ZERO, ONE, Z1, 1, 1, DESCZ )

*     Take a copy of A1 and B1, store in A2 and B2
      CALL PDLACPY( 'All', N, N, A1, 1, 1, DESCA, A2,
     $           1, 1, DESCA )
      CALL PDLACPY( 'All', N, N, B1, 1, 1, DESCA, B2,
     $           1, 1, DESCA )


*     QR factorize B
      IF ( IAM .EQ. 0 ) WRITE(*,*)'%% QR factor B'
      TIMEQRF = MPI_WTIME()
      CALL PDGEQRF(N, N, B1, 1, 1, DESCB, DWORK,
     $     DWORK( N + 1 ), LDWORK - N, INFO)
      TIMEQRF = MPI_WTIME() - TIMEQRF
      IF ( IAM .EQ. 0 ) WRITE(*,*)'%% PDGEQRF took ', TIMEQRF, 
     $   ' seconds to complete.'

      IF ( INFO .NE. 0 ) THEN
         WRITE(*,*) '% PDGEQRF: INFO =', INFO
      END IF
      LDWORK = MAX(LDWORK, INT(DDUM(1)))
      INFO = 0

*     Apply the Q factor on A
      TIMEDORMQR = MPI_WTIME()
      IF ( IAM .EQ. 0 ) WRITE(*,*)'%% Apply QR factor on A'
      CALL PDORMQR( 'L', 'T', N, N, N, B1, 1, 1, DESCB,
     $     DWORK, A1, 1, 1, DESCA, DWORK( N + 1 ), LDWORK - N,
     $     INFO )
*     Apply the Q factor on Q 
      IF ( IAM .EQ. 0 ) WRITE(*,*)'%% Apply QR factor on Q'
      CALL PDORMQR( 'L', 'T', N, N, N, B1, 1, 1, DESCB,
     $     DWORK, Q1, 1, 1, DESCQ, DWORK( N + 1 ), LDWORK - N,
     $     INFO )
      TIMEDORMQR = MPI_WTIME() - TIMEDORMQR

      IF ( IAM .EQ. 0 ) WRITE(*,*)'%% PDORMQR took ', TIMEDORMQR,
     $   ' seconds to complete.'
 
*     Set lower part of B to zero      
      DO J = 1, N
          DO I = J + 1, N
              CALL PDELSET( B1, I, J, DESCB, ZERO)
          END DO
      END DO
      INFO = 0
*     Reduce the matrix pair to Hessenberg-triangular form
      IF ( IAM .EQ. 0 ) WRITE(*,*)'%% Calling PDGGHRD'
      AB2HTTIME = MPI_WTIME()
      CALL PDGGHRD( 'V', 'V', N, ILO, IHI,
     $     A1, DESCA, B1, DESCB, 
     $     Q1, DESCQ, Z1, DESCZ,
     $     DWORK, LDWORK, INFO )
      AB2HTTIME = MPI_WTIME() - AB2HTTIME
      IF ( IAM .EQ. 0 ) WRITE(*,*)'%% PDGGHRD took ', AB2HTTIME,
     $   ' seconds to complete.'

      IF ( INFO .NE. 0 ) THEN
         WRITE(*,*) '% PDHGGHRD: INFO =', INFO
      END IF
*
*     Transpose Q1 into Q2
*
      TRANSPTIME = MPI_WTIME()
      IF ( IAM .EQ. 0 ) THEN
         WRITE(*,*)'%% Transposing Q1 into Q2'
      END IF
      CALL PDGEADD('T', N, N, ONE, Q1, 1, 1, DESCQ,
     $    ZERO, Q2, 1, 1, DESCQ)
*     Copy Q2 to Q1
      CALL PDGEADD('N', N, N, ONE, Q2, 1, 1, DESCQ,
     $    ZERO, Q1, 1, 1, DESCQ)
      TRANSPTIME = MPI_WTIME() - TRANSPTIME
      IF ( IAM .EQ. 0 ) THEN
         WRITE(*,*)'%% Transposing Q1 into Q2 took', TRANSPTIME,
     $       'seconds'
         CALL FLUSH(6)
      END IF

*     QZ : Reduction to (S,T) from (H,T)
      INFO = 0
      IF ( IAM .EQ. 0 ) WRITE(*,*)'%% Calling PDHGEQZ'

      QZTIME = MPI_WTIME()
      CALL PDHGEQZ( 'S', 'V', 'V', N, ILO, IHI, 
     $     A1, DESCA, B1, DESCB,
     $     AR, AI, BETA,
     $     Q1, DESCQ, Z1, DESCZ, 
     $     DWORK, LDWORK,
     $     IWORK, LIWORK, INFO )              
      QZTIME = MPI_WTIME() - QZTIME
        
      IF ( IAM .EQ. 0 ) WRITE(*,*)'%% PDHGEQZ took ', QZTIME,
     $   ' seconds to complete.'
      IF ( INFO .NE. 0 ) THEN
         WRITE(*,*) '% PDHGEQZ: INFO =', INFO
      END IF                     
      
      RESTIME = MPI_WTIME()
      IF ( IAM .EQ. 0 ) WRITE(*,*)'%% Computing residuals.'
              
      CALL PCOMPRES( A2, B2, Q1, Z1, A1, B1,DESCA,DESCB,
     $     DESCZ, DESCQ, DESCA, DESCB,
     $     N, ANORM, BNORM, QNORM, ZNORM, .TRUE.)
      RESTIME = MPI_WTIME() - RESTIME 
         
      IF ( IAM .EQ. 0 ) WRITE(*,*)'%% Residual comp. took ', RESTIME,
     $   ' seconds to complete.' 

      TOTTIME = MPI_WTIME() - TOTTIME

      IF( IAM.EQ.0 ) WRITE(*,*)
     $   '%% Total execution time in seconds:', TOTTIME
      IF ( IAM .EQ. 0 ) THEN 
          WRITE(*,*)'% ||(S-Q^(A)Z)|| / ||A||=', ANORM 
          WRITE(*,*)'% ||(T-Q^(B)Z)|| / ||B||=', BNORM 
          WRITE(*,*)'% ||(QQ^ - I)|| / (N*eps)=', QNORM 
          WRITE(*,*)'% ||(ZZ^ - I)|| / (N*eps)=', ZNORM
      END IF               
      IF (ANORM.LT.1E-14.AND.BNORM.LT.1E-14.AND.QNORM.LT.10.AND.
     $      ZNORM.LT.10) THEN
          IF (IAM.EQ.0) WRITE(*,*)'PDHGEQZ END TEST ALL OK'
      END IF
      DEALLOCATE ( DWORK )         
      DEALLOCATE ( AR, AI, BETA )
                     
      DEALLOCATE( A1, B1, Q1, Z1, Q2 )
      DEALLOCATE( A2, B2 )

      DEALLOCATE( IWORK )

      CALL BLACS_GRIDEXIT( ICTXT )

      
 777  CONTINUE  
      CALL BLACS_EXIT( 0 )
*     
*     End of EXRAND2
*     
      END
 
