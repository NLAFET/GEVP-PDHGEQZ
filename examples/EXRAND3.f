***********************************************************************
*                                                                     *
*     Example driver for ScaLAPACK-style routine PDHGEQZ              *
*     (reads a matrix pair from matrix market files)                  *
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
      PROGRAM EXRAND3
*     ..
*     .. Declarations
*     ..
      IMPLICIT NONE
*     
*     .. Parameters
      INTEGER            BLOCK_CYCLIC_2D, CSRC_, CTXT_, DLEN_, DTYPE_,
     $   LLD_, MB_, M_, NB_, N_, RSRC_
      PARAMETER          ( BLOCK_CYCLIC_2D = 1, DLEN_ = 9, DTYPE_ = 1,
     $   CTXT_ = 2, M_ = 3, N_ = 4, MB_ = 5, NB_ = 6,
     $   RSRC_ = 7, CSRC_ = 8, LLD_ = 9 )

*     Matrix Market files
      CHARACTER(40)       AMATRIX, BMATRIX
c     PARAMETER           (AMATRIX = 'mhd4800a.mtx', 
c     $     BMATRIX='mhd4800b.mtx')
      PARAMETER         (AMATRIX = 'sparseA.mtx', BMATRIX='sparseB.mtx')

      DOUBLE PRECISION    ZERO, ONE
      PARAMETER           ( ZERO = 0.D0, ONE = 1.0D+0 )  

      INTEGER             NIN1, NIN2
      PARAMETER           ( NIN1 = 11, NIN2 = 12 ) 
*     
*     .. Local Scalars
      INTEGER             N
      INTEGER*8           LDWORK, LIWORK
      INTEGER             ICTXT, NPROW, NPCOL, MYROW, MYCOL
      INTEGER             IAM, SYSPROCS, NPROCS
      INTEGER             ILO, IHI, NP, NQ, NB
      INTEGER             I, J, INFO 
      INTEGER*8           TOTNNZA, TOTNNZB
      INTEGER             NROWSA, NCOLSA, NROWSB, NCOLSB
      INTEGER             FCNT, M, IDUM
      DOUBLE PRECISION    ANORM, BNORM, QNORM, ZNORM
      DOUBLE PRECISION    EPS
      DOUBLE PRECISION    TOTTIME, QZTIME, RESTIME
      DOUBLE PRECISION    AB2HTTIME, TIMEDORMQR, TIMEQRF, REORDERTIME

      CHARACTER           REPA*10, REPB*10
      CHARACTER           FIELDA*7, FIELDB*7
      CHARACTER           SYMMA*19, SYMMB*19
      LOGICAL             CHKR, CHKI
*     
*     .. Local Arrays 
      INTEGER             DESCA( DLEN_ ), DESCB( DLEN_ )
      INTEGER             DESCQ( DLEN_ ), DESCZ( DLEN_ )
      INTEGER             PARA(6)
      DOUBLE PRECISION    DDUM(1)
      
*     Placeholders for the computed eigenvalues
      DOUBLE PRECISION, ALLOCATABLE :: AR(:)
      DOUBLE PRECISION, ALLOCATABLE :: AI(:) 
      DOUBLE PRECISION, ALLOCATABLE :: BETA(:)
*     Need two copies of A and B to be able to perform residual check 
*     after computation.     
      DOUBLE PRECISION, ALLOCATABLE :: A1(:)
      DOUBLE PRECISION, ALLOCATABLE :: A2(:)
      DOUBLE PRECISION, ALLOCATABLE :: B1(:)
      DOUBLE PRECISION, ALLOCATABLE :: B2(:)
      DOUBLE PRECISION, ALLOCATABLE :: Q1(:)
      DOUBLE PRECISION, ALLOCATABLE :: Q2(:)
      DOUBLE PRECISION, ALLOCATABLE :: Z1(:)
*     Workspace
      INTEGER, ALLOCATABLE :: IWORK(:)
      DOUBLE PRECISION, ALLOCATABLE :: DWORK(:)  
*     
*     .. Externals 
      EXTERNAL          CHISQ, MPI_WTIME, PDLAMCH
      DOUBLE PRECISION  CHISQ, MPI_WTIME, PDLAMCH
      EXTERNAL          NUMROC, PILAENVX
      INTEGER           NUMROC, PILAENVX

*     
*     .. Intrinsic Functions 
      INTRINSIC         DBLE, MAX

      
*     
*     .. Executable Statements

      INFO = 0
      TOTTIME = MPI_WTIME()


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


*     Read information from the matrix market file
*     First open the files
      OPEN( UNIT = NIN1, FILE = trim(AMATRIX), STATUS = 'OLD')
      OPEN( UNIT = NIN2, FILE = trim(BMATRIX), STATUS = 'OLD')

*     Use the version that does not rewind the file
      CALL mminfo2( NIN1, REPA, FIELDA, SYMMA, NROWSA, NCOLSA, 
     $   TOTNNZA )
      CALL mminfo2( NIN2, REPB, FIELDB, SYMMB, NROWSB, NCOLSB, 
     $   TOTNNZB )

*     Assume NROWSA = NCOLSA = NROWSB = NCOLSB
      N = NROWSA

*     Hard coded workspace for integers, needed by infinite eigenvalue computation
*     as well as reordering.
*     N more than 4 * N should be needed at any point
      LIWORK = 4 * N 
      ALLOCATE( IWORK( LIWORK ), STAT = INFO )
      IF( INFO.NE.0 ) THEN
         WRITE(*,*) 'Could not allocate IWORK. INFO = ', INFO
         GOTO 777 
      END IF  

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
         WRITE(*,*) '%%      EXAMPLE PROGRAM 3 FOR PDHGEQZ       %%'
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
     $   MAX( 1, NP ), INFO )
      CALL DESCINIT( DESCB, N, N, NB, NB, 0, 0, ICTXT, 
     $   MAX( 1, NP ), INFO )          

*     Initialize the array descriptor for the matrix Q, Z
      CALL DESCINIT( DESCQ, N, N, NB, NB, 0, 0, ICTXT, 
     $   MAX( 1, NP ), INFO )     
      CALL DESCINIT( DESCZ, N, N, NB, NB, 0, 0, ICTXT, 
     $   MAX( 1, NP ), INFO ) 
      
      IF ( IAM .EQ. 0 ) THEN
         WRITE(*,*) '%% Allocating memory for matrices and eigenvalues'
      END IF
*     Allocate needed space for the matrices A*2,B*2,Q and Z.
      ALLOCATE ( A1( DESCA( LLD_ )*NQ), STAT = INFO )
      IF( INFO.NE.0 ) THEN
         WRITE(*,*) '% Could not allocate A1. INFO = ', INFO
         GOTO 777 
      END IF           
      
      ALLOCATE ( B1( DESCB( LLD_ )*NQ))
      IF( INFO.NE.0 ) THEN
         WRITE(*,*) '% Could not allocate B1. INFO = ', INFO
         GOTO 777 
      END IF           
      
      ALLOCATE ( Q1( DESCQ( LLD_ )*NQ))
      IF( INFO.NE.0 ) THEN
         WRITE(*,*) '% Could not allocate Q1. INFO = ', INFO
         GOTO 777 
      END IF    
      
      ALLOCATE ( Q2( DESCQ( LLD_ )*NQ))
      IF( INFO.NE.0 ) THEN
         WRITE(*,*) '% Could not allocate Q2. INFO = ', INFO
         GOTO 777 
      END IF  

      ALLOCATE ( Z1( DESCZ( LLD_ )*NQ))
      IF( INFO.NE.0 ) THEN
         WRITE(*,*) '% Could not allocate Z1. INFO = ', INFO
         GOTO 777 
      END IF    

      ALLOCATE ( A2( DESCA( LLD_ )*NQ))
      IF( INFO.NE.0 ) THEN
         WRITE(*,*) '% Could not allocate A2. INFO = ', INFO
         GOTO 777 
      END IF           
      
      ALLOCATE ( B2( DESCB( LLD_ )*NQ))
      IF( INFO.NE.0 ) THEN
         WRITE(*,*) '% Could not allocate B2. INFO = ', INFO
         GOTO 777 
      END IF               
            
*     Allocate space for eigenvalues
      ALLOCATE ( AR( N ))
      IF( INFO.NE.0 ) THEN
         WRITE(*,*) '% Could not allocate AR. INFO = ', INFO
         GOTO 777 
      END IF
      ALLOCATE ( AI( N ))
      IF( INFO.NE.0 ) THEN
         WRITE(*,*) '% Could not allocate AI. INFO = ', INFO
         GOTO 777 
      END IF
      ALLOCATE ( BETA( N ))
      IF( INFO.NE.0 ) THEN
         WRITE(*,*) '% Could not allocate BETA. INFO = ', INFO
         GOTO 777 
      ENDIF
*     BETA should never be negativ, so set a negtive value for 
*     being "non set".
      Q2(N*N) = ZERO
      A2(N*N) = ZERO
      B2(N*N) = ZERO
      BETA ( 1 : N ) = -1 
*     
      LDWORK = 0
*     Perform a workspace query for PDGEQRF
      CALL PDGEQRF(N, N, B1, 1, 1, DESCB, DDUM,
     $   DDUM, -1, INFO)
      LDWORK = MAX(LDWORK, INT(DDUM(1)))
*     Perform a workspace query for PDORMQR
      CALL PDORMQR( 'L', 'T', N, N, N, B1, 1, 1, DESCB,
     $   DDUM, A1, 1, 1, DESCA, DDUM, -1,
     $   INFO ) 
      LDWORK = MAX(LDWORK, INT(DDUM(1)+N))


*     Perform a workspace query to PDGGHRD
      CALL PDGGHRD( 'I', 'I', N, ILO, IHI,
     $   A1, DESCA, B1, DESCB,
     $   Q1, DESCQ, Z1, DESCZ,
     $   DDUM, -1, INFO )
      LDWORK = MAX(LDWORK, INT(DDUM(1)+N))

      
*     Perform a workspace query to PDHGEQZ      
      CALL PDHGEQZ( 'S', 'V', 'V', N, ILO, IHI, 
     $   A1, DESCA, B1, DESCB,
     $   AR, AI, BETA,
     $   Q1, DESCQ, Z1, DESCZ, 
     $   DDUM, -1,
     $   IWORK, LIWORK, INFO )
      LDWORK = MAX(LDWORK, INT(DDUM(1)))

      
*     Make sure all processors allocate the same amount of work          
      CALL IGAMX2D( ICTXT, 'All', ' ', 1, 1, LDWORK, 1, -1, -1, -1,
     $   -1, -1 )   

      LDWORK = LDWORK + 2000000

      IF ( IAM .EQ. 0 ) WRITE(*,*)'%% Allocating workspace :' , 
     $   LDWORK
*     Allocate workspace          
      ALLOCATE ( DWORK( LDWORK ))
      IF( INFO.NE.0 ) THEN
         WRITE(*,*) 'Could not allocate Workspace. INFO = ', INFO
         GOTO 777 
      END IF 

      IF ( IAM .EQ. 0) WRITE(*,*)'%% Done Allocating.'
      CALL FLUSH( 6 )
      
*     Read the data
      
      IF ( IAM .EQ. 0 ) WRITE(*,*)'%% Reading data. N=', N
      CALL FLUSH( 6 )
      CALL PDRDMM( NIN1, REPA, FIELDA, SYMMA, NROWSA, NCOLSA, 
     $   TOTNNZA, A1, DESCA, LDWORK, DWORK, INFO )

      CALL PDRDMM( NIN2, REPB, FIELDB, SYMMB, NROWSB, NCOLSB, 
     $   TOTNNZB, B1, DESCB, LDWORK, DWORK, INFO )

      IF ( IAM .EQ. 0) WRITE(*,*)'%% Done Reading data.'
      

*     Close the files since they are no longer needed
      CLOSE( NIN1 )
      CLOSE( NIN2 )

*     Initiliaze Q and Z to the identity      
      CALL PDLASET( 'Full', N, N, ZERO, ONE, Q1, 1, 1, DESCQ )
      CALL PDLASET( 'Full', N, N, ZERO, ONE, Z1, 1, 1, DESCZ )

*     Take a copy of A1 and B1, store in A2 and B2      
      CALL PDLACPY( 'All', N, N, A1, 1, 1, DESCA, A2,
     $   1, 1, DESCA )
      CALL PDLACPY( 'All', N, N, B1, 1, 1, DESCB, B2,
     $   1, 1, DESCB )

*     QR factorize B
      IF ( IAM .EQ. 0 ) WRITE(*,*)'%% QR factor B'
      CALL FLUSH( 6 )
      TIMEQRF = MPI_WTIME()
      CALL PDGEQRF(N, N, B1, 1, 1, DESCB, DWORK,
     $   DWORK( N + 1 ), LDWORK - N, INFO)
      TIMEQRF = MPI_WTIME() - TIMEQRF
      IF ( IAM .EQ. 0 ) WRITE(*,*)'%% PDGEQRF took ', TIMEQRF, 
     $   ' seconds to complete.'

      IF ( INFO .NE. 0 ) THEN
         WRITE(*,*) '% PDGEQRF: INFO =', INFO
      END IF
      CALL FLUSH( 6 ) 
      LDWORK = MAX(LDWORK, INT(DDUM(1)))
      INFO = 0

*     Apply the Q factor on A
      TIMEDORMQR = MPI_WTIME()
      IF ( IAM .EQ. 0 ) WRITE(*,*)'%% Apply QR factor on A'
      CALL PDORMQR( 'L', 'T', N, N, N, B1, 1, 1, DESCB,
     $   DWORK, A1, 1, 1, DESCA, DWORK( N + 1 ), LDWORK - N,
     $   INFO )
*     Apply the Q factor on Q 
      IF ( IAM .EQ. 0 ) WRITE(*,*)'%% Apply QR factor on Q'
      CALL PDORMQR( 'L', 'T', N, N, N, B1, 1, 1, DESCB,
     $   DWORK, Q1, 1, 1, DESCQ, DWORK( N + 1 ), LDWORK - N,
     $   INFO )
      TIMEDORMQR = MPI_WTIME() - TIMEDORMQR

      IF ( IAM .EQ. 0 ) WRITE(*,*)'%% PDORMQR took ', TIMEDORMQR,
     $   ' seconds to complete.'
      CALL FLUSH( 6 ) 
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
     $   A1, DESCA, B1, DESCB, 
     $   Q1, DESCQ, Z1, DESCZ,
     $   DWORK, LDWORK, INFO )

      IF ( IAM .EQ. 0 ) THEN
         WRITE(*,*)'% Transposing Q1 into Q2'
      END IF
      CALL PDGEADD('T', N, N, ONE, Q1, 1, 1, DESCQ, 
     $   ZERO, Q2, 1, 1, DESCQ)
*     Copy Q2 to Q1
      CALL PDGEADD('N', N, N, ONE, Q2, 1, 1, DESCQ, 
     $   ZERO, Q1, 1, 1, DESCQ)

      AB2HTTIME = MPI_WTIME() - AB2HTTIME
      IF ( IAM .EQ. 0 ) WRITE(*,*)'%% PDGGHRD took ', AB2HTTIME,
     $   ' seconds to complete.'

      IF ( INFO .NE. 0 ) THEN
         WRITE(*,*) '% PDHGGHRD: INFO =', INFO
      END IF
      CALL FLUSH( 6 ) 
      GOTO 600
*     QZ : Reduction to (S,T) from (H,T)
      INFO = 0
      IF ( IAM .EQ. 0 ) WRITE(*,*)'%% Calling PDHGEQZ'

      QZTIME = MPI_WTIME()
      CALL PDHGEQZ( 'S', 'V', 'V', N, ILO, IHI, 
     $   A1, DESCA, B1, DESCB,
     $   AR, AI, BETA,
     $   Q1, DESCQ, Z1, DESCZ, 
     $   DWORK, LDWORK,
     $   IWORK, LIWORK, INFO )              
      QZTIME = MPI_WTIME() - QZTIME
      
      IF ( IAM .EQ. 0 ) WRITE(*,*)'%% PDHGEQZ took ', QZTIME,
     $   ' seconds to complete.'
      IF ( INFO .NE. 0 ) THEN
         WRITE(*,*) '% PDHGEQZ: INFO =', INFO
      END IF                     
      CALL FLUSH( 6 )
      GOTO 400
*     
*     Now reorder the matrix pair such that eigenvalues belonging to 
*     the unit circle is found at the top
*     
      
*     Use first N rows of IWORK as SELECT vector, just clear all first
      IWORK (1 : LIWORK ) = 0

*     Scan the eigenvalues and mark in the SELECT vector to move
*     up if they match our condition abs(AR(I)/BETA(I))<1 and
*     abs(AI(I)/BETA(I))< 1.
*     FCNT will reflect how many eigenvalues that met the condition
      


      FCNT = 1
      DO I = 1, N
         IF ( BETA( I ) .NE. 0 ) THEN
            CHKR = ABS(AR ( I ) / BETA ( I ) ) .LT. ONE
            CHKI = ABS(AI ( I ) / BETA ( I ) ) .LT. ONE
            CHKR = CHKR .AND. ( ( AR( I ) / BETA ( I ) ) .GE. 0 )
            CHKI = CHKI .AND. ( ( AI( I ) / BETA ( I ) ) .GE. 0 )
            IF ( CHKR .AND. CHKI ) THEN
               IWORK( I ) = 1
               FCNT = FCNT + 1
            END IF
         END IF
      END DO
      
      IF ( IAM .EQ. 0 ) WRITE(*,*)'%% Found ', FCNT, ' eigenvalues
     $   that lies within the unit circle.'
      CALL FLUSH( 6 ) 

      
*     Call the parallel reordering routine to move the selected
*     eigenvalues in place
*     
      
*     Perpare the call by setting parameters
*     Maximum number of independent computational windows
      PARA( 1 ) = PILAENVX(ICTXT, 80, 'PDHGEQZ', '', N, NB,
     $   IDUM, IDUM)
*     Number of eigenvalues in each window
      PARA( 2 ) = PILAENVX(ICTXT, 81, 'PDHGEQZ', '', IDUM, NB,
     $   IDUM, IDUM)
*     Computational window size
      PARA( 3 ) = PILAENVX(ICTXT, 82, 'PDHGEQZ', '', IDUM, NB,
     $   IDUM, IDUM)
*     Minimal percentage of flops required for
*     performing matrix-matrix multiplications instead of
*     pipelined orthogonal transformations
      PARA( 4 ) = PILAENVX(ICTXT, 83, 'PDHGEQZ', '', IDUM, IDUM,
     $   IDUM, IDUM)
*     Width of block column slabs for row-wise
*     application of pipelined orthogonal transformations in
*     their factorized form
      PARA( 5 ) = PILAENVX(ICTXT, 84, 'PDHGEQZ', '', IDUM, NB,
     $   IDUM, IDUM)
*     Maximum number of eigenvalues to bring over
*     the block border 
      PARA( 6 ) = PILAENVX(ICTXT, 85, 'PDHGEQZ', '', IDUM, NB,
     $   IDUM, IDUM)

      INFO = 0

      IF ( FCNT .GT. 0 ) THEN
         REORDERTIME = MPI_WTIME()
         CALL PDTGORD( .TRUE. , .TRUE., IWORK, PARA, N,
     $      A1, 1, 1, DESCA, B1, 1, 1, DESCB, 
     $      Q1, 1, 1, DESCQ, Z1, 1, 1, DESCZ,
     $      AR, AI, BETA,
     $      M, DWORK, LDWORK, IWORK( N + 1 ), LIWORK - N, INFO )
         REORDERTIME = MPI_WTIME() - REORDERTIME

         IF ( IAM .EQ. 0 ) WRITE(*,*)'%% Reordering took ', REORDERTIME,
     $      ' seconds to complete. Defl. subspace M=', M
         IF ( INFO .NE. 0 ) THEN
            WRITE(*,*) '%% PDTGORD: INFO =', INFO
         END IF   
         CALL FLUSH( 6 ) 
      END IF

*     Print the FCNT eigenvalues, in Matlab style 
*     if no failure occured
      IF ( INFO. EQ. 0 .AND. FCNT .GT. 0 ) THEN
         IF ( IAM .EQ. 0 ) THEN
            WRITE(*,*)'Eigs=['
            DO I = 1, FCNT               
               WRITE(*,*) AR( I ), AI( I ), BETA( I )
            END DO
            WRITE(*,*)'];'
         END IF
      END IF
      GOTO 600
 400  CONTINUE
      IF (IAM.EQ.0) THEN
         WRITE(*,*)'Eigs=['
         DO I = 1, N
            WRITE(*,*) AR( I ), AI( I ), BETA(I)
            CALL FLUSH(6)
         END DO
         WRITE(*,*)'];'    
         CALL FLUSH( 6 ) 
      END IF
 600  CONTINUE

      
      RESTIME = MPI_WTIME()
      IF ( IAM .EQ. 0 ) WRITE(*,*)'%% Computing residuals.'
      
      CALL PCOMPRES( A2, B2, Q1, Z1, A1, B1, DESCA, DESCB,
     $   DESCZ, DESCQ, DESCA, DESCB, N, ANORM, BNORM, QNORM, ZNORM, 
     $   .TRUE.)
      RESTIME = MPI_WTIME() - RESTIME 
      
      IF ( IAM .EQ. 0 ) WRITE(*,*)'%% Residual comp. took ', RESTIME,
     $   ' seconds to complete.' 
         
      CALL FLUSH( 6 ) 
      
      TOTTIME = MPI_WTIME() - TOTTIME
     
      IF ( IAM .EQ. 0 ) THEN 
         WRITE(*,*)'% ||(S-Q^(A)Z)|| / ||A||=', ANORM 
         WRITE(*,*)'% ||(T-Q^(B)Z)|| / ||B||=', BNORM 
         WRITE(*,*)'% ||(QQ^ - I)|| / (N*eps)=', QNORM 
         WRITE(*,*)'% ||(ZZ^ - I)|| / (N*eps)=', ZNORM
      END IF  
 
      IF( IAM.EQ.0 ) WRITE(*,*)
     $   '%% Total execution time in seconds:', TOTTIME

      DEALLOCATE ( DWORK )         
      DEALLOCATE ( AR, AI, BETA )
      
      DEALLOCATE( A1, B1, Q1, Q2, Z1 )
      DEALLOCATE( A2, B2 )
      DEALLOCATE( IWORK )

      CALL BLACS_GRIDEXIT( ICTXT )
      
      GOTO 999

 777  CONTINUE  
      CLOSE( NIN1 )
      CLOSE( NIN2 )
      
 999  CONTINUE

      CALL BLACS_EXIT( 0 )
*     
*     End of EXRAND3
*     
      END
 
