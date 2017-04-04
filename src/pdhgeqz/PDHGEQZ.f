***********************************************************************
*                                                                     *
*     PDHGEQZ.f:                                                      *
*         Driver routine in the package PDHGEQZ.                      *
*                                                                     *
*     Contributors: Bjorn Adlerborn                                   *
*                   Bo Kagstrom                                       *
*                   Daniel Kressner                                   *
*                                                                     *
*     Department of Computing Science and HPC2N, Umea University      *
*     MATHICSE ANCHP, EPF Lausanne                                    *
*	                                                              * 
***********************************************************************

      SUBROUTINE PDHGEQZ( JOB, COMPQ, COMPZ, N, ILO, IHI, H, DESCH, T, 
     $   DESCT, ALPHAR, ALPHAI, BETA, Q, DESCQ, Z, DESCZ, DWORK, LDWORK, 
     $   IWORK, LIWORK, INFO )

      IMPLICIT NONE
*     ..
*     .. Scalar Arguments ..
*     ..
      CHARACTER          COMPQ, COMPZ, JOB
      INTEGER            IHI, ILO, INFO, LDWORK, N, LIWORK
*     ..
*     .. Array Arguments .. 
*     ..
      DOUBLE PRECISION   H( * ), T( * ), Q( * ), Z( * ), ALPHAI( * ), 
     $   ALPHAR( * ), BETA( * ), DWORK( * )
      INTEGER            IWORK( * ), DESCH( 9 ), DESCT( 9 ), DESCQ( 9 ), 
     $   DESCZ( 9 )
      
      

*     Purpose
*     =======
*     
*     PDHGEQZ : implements a parallel single-/double-shift version of the QZ method for
*     finding the generalized eigenvalues
*     
*     w(j)=(ALPHAR(j) + i*ALPHAI(j))/BETAR(j)   of the equation
*     
*     det( H - w(i) T ) = 0
*     
*     In addition, the pair H,T may be reduced to generalized Schur form:
*     T is upper triangular, and H is block upper triangular, where the
*     diagonal blocks are either 1-by-1 or 2-by-2, the 2-by-2 blocks having
*     complex generalized eigenvalues (see the description of the argument
*     JOB.)
*     
*     If JOB='S', then the pair (H,T) is simultaneously reduced to Schur
*     form by applying one orthogonal tranformation (usually called Q) on
*     the left and another (usually called Z) on the right.  The 2-by-2
*     upper-triangular diagonal blocks of T corresponding to 2-by-2 blocks
*     of H will be reduced to positive diagonal matrices.  (I.e.,
*     if H(j+1,j) is non-zero, then T(j+1,j)=T(j,j+1)=0 and T(j,j) and
*     T(j+1,j+1) will be positive.)
*     
*     If JOB='E', then at each iteration, the same transformations
*     are computed, but they are only applied to those parts of H and T
*     which are needed to compute ALPHAR, ALPHAI, and BETAR.
*     
*     If JOB='S' and COMPQ and COMPZ are 'V' or 'I', then the orthogonal
*     transformations used to reduce (H,T) are accumulated into the arrays
*     Q and Z s.t.:
*     
*     Q(in) H(in) Z(in)* = Q(out) H(out) Z(out)*
*     Q(in) T(in) Z(in)* = Q(out) T(out) Z(out)*
*     
*     Ref: C.B. Moler & G.W. Stewart, "An Algorithm for Generalized Matrix
*     Eigenvalue Problems", SIAM J. Numer. Anal., 10(1973),
*     pp. 241--256.
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
*     Arguments
*     =========
*     
*     JOB     (input) CHARACTER*1
*     = 'E': compute only ALPHAR, ALPHAI, and BETA.  A and B will
*     not necessarily be put into generalized Schur form.
*     = 'S': put A and B into generalized Schur form, as well
*     as computing ALPHAR, ALPHAI, and BETA.
*     
*     COMPQ   (global input) CHARACTER*1
*     = 'N': do not compute Q;
*     = 'I': Q is initialized to the unit matrix, and the
*     orthogonal matrix Q is returned;
*     = 'V': Q must contain an orthogonal matrix Q1 on entry,
*     and the product Q1*Q is returned.
*     
*     COMPZ   (global input) CHARACTER*1
*     = 'N': do not compute Z;
*     = 'I': Z is initialized to the unit matrix, and the
*     orthogonal matrix Z is returned;
*     = 'V': Z must contain an orthogonal matrix Z1 on entry,
*     and the product Z1*Z is returned.
*     
*     N       (global input) INTEGER
*     The order of the matrices H,T,Q, and Z.  N >= 0.
*     
*     ILO     (global input) INTEGER
*     IHI     (global input) INTEGER
*     It is assumed that A is already upper triangular in rows and
*     columns 1:ILO-1 and IHI+1:N.  ILO and IHI are normally set
*     by a previous call to DGGBAL; otherwise they should be set
*     to 1 and N respectively.
*     1 <= ILO <= IHI <= N, if N > 0; ILO=1 and IHI=0, if N=0.
*     
*     H       (local input/output) DOUBLE PRECISION array, dimension (LLD_H, LOCc(N)). 
*     On entry, H contains the upper Hessenberg matrix H.
*     On exit, if JOB = 'S', H is quasi-upper triangular in rows and columns 
*     (ILO:IHI), with 1 x 1 and 2 x 2 blocks on the diagonal where 
*     the 2 x 2 blocks corresponds to complex conjugated pairs of 
*     eigenvalues. If JOB='E', H is unspecified on exit.      
*     
*     DESCH   (global and local input) INTEGER array of dimension DLEN_.
*     The array descriptor for the distributed matrix H.
*     
*     T       (local input/output) DOUBLE PRECISION array, dimension (LLD_T, LOCc(N)).
*     On entry, the N-by-N upper triangular matrix B.
*     On exit, the upper triangular matrix T = Q' B Z.  The
*     elements below the diagonal are set to zero.
*     
*     DESCT   (global and local input) INTEGER array of dimension DLEN_.
*     The array descriptor for the distributed matrix T.
*     
*     ALPHAR  (global output) DOUBLE PRECISION array, dimension (N)
*     ALPHAR(1:N) will be set to real parts of the diagonal
*     elements of H that would result from reducing H and T to
*     Schur form and then further reducing them both to triangular
*     form using unitary transformations s.t. the diagonal of T
*     was non-negative real.  Thus, if H(j,j) is in a 1-by-1 block
*     (i.e., H(j+1,j)=H(j,j+1)=0), then ALPHAR(j)=H(j,j).
*     Note that the (real or complex) values
*     (ALPHAR(j) + i*ALPHAI(j))/BETA(j), j=1,...,N, are the
*     generalized eigenvalues of the matrix pencil H - wT.
*     
*     ALPHAI  (global output) DOUBLE PRECISION array, dimension (N)
*     ALPHAI(1:N) will be set to imaginary parts of the diagonal
*     elements of A that would result from reducing H and T to
*     Schur form and then further reducing them both to triangular
*     form using unitary transformations s.t. the diagonal of T
*     was non-negative real.  Thus, if H(j,j) is in a 1-by-1 block
*     (i.e., H(j+1,j)=H(j,j+1)=0), then ALPHAI(j)=0.
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
*     Q       (local input/output) DOUBLE PRECISION array, dimension (LLD_Q, LOCc(N)).
*     If COMPQ='N':  Q is not referenced.
*     If COMPQ='I':  on entry, Q need not be set, and on exit it
*     contains the orthogonal matrix Q, where Q'
*     is the product of the Givens transformations
*     which are applied to A and B on the left.
*     If COMPQ='V':  on entry, Q must contain an orthogonal matrix
*     Q1, and on exit this is overwritten by Q1*Q.
*     
*     DESCQ   (global and local input) INTEGER array of dimension DLEN_.
*     The array descriptor for the distributed matrix Q.
*     
*     Z       (local input/output) DOUBLE PRECISION array, dimension (LLD_Z, LOCc(N)).
*     If COMPZ='N':  Z is not referenced.
*     If COMPZ='I':  on entry, Z need not be set, and on exit it
*     contains the orthogonal matrix Z, which is
*     the product of the Givens transformations
*     which are applied to A and B on the right.
*     If COMPZ='V':  on entry, Z must contain an orthogonal matrix
*     Z1, and on exit this is overwritten by Z1*Z.
*     
*     DESCZ   (global and local input) INTEGER array of dimension DLEN_.
*     The array descriptor for the distributed matrix Z.
*     
*     DWORK   (local workspace/global output) DOUBLE PRECISION array,
*     dimension (LDWORK) 
*     
*     LDWORK  (global input) INTEGER 
*     The dimension of the array DWORK. 
*     If LWORK = -1, then a workspace query is assmed, and optimal 
*     workspace is returned in DWORK(1)
*     
*     IWORK   (local workspace/global output) INTEGER array, 
*     dimension (LIWORK)  
*     
*     LIWORK  (global input) INTEGER  
*     The dimension of the array IWORK. 
*     If LIWORK = -1 then workspace query is assumed, and optimal  
*     workspace is returned in IWORK(1)
*     
*     INFO    (global output) INTEGER
*     = 0:  successful exit.
*     < 0:  if INFO = -i, the i-th argument had an illegal value.
*     
*     
*     =====================================================================



*     ..
*     .. Parameters ..
*     ..
      DOUBLE PRECISION   ZERO, ONE, TWO
      PARAMETER          ( ZERO = 0.0D+0, ONE = 1.0D+0, TWO = 2.0D+0 )

      INTEGER            BLOCK_CYCLIC_2D, CSRC_, CTXT_, DLEN_, DT_, 
     $   LLD_, MB_, M_, NB_, N_, RSRC_
      PARAMETER          ( BLOCK_CYCLIC_2D = 1, DLEN_ = 9, DT_ = 1,
     $   CTXT_ = 2, M_ = 3, N_ = 4, MB_ = 5, NB_ = 6,
     $   RSRC_ = 7, CSRC_ = 8, LLD_ = 9 )

      
*     ..
*     .. Local Scalars ..
*     ..      
      LOGICAL             DEBUG, ILQ, ILSCHR, ILZ, LQUERY
      INTEGER             ICOMPQ, ICOMPZ, ISCHUR, ICTXT, NPROW, NPCOL, 
     $   MYROW, MYCOL, IERR, NMIN1, IAM, NRPROCS, IDUM, LDWKOPT, 
     $   LIWKOPT, NH
      DOUBLE PRECISION    UNFL, OVFL


*     ..
*     .. Precision paramters ..
*     ..
      DOUBLE PRECISION    SMLNUM, ULP, TFNORM, T1NORM, HFNORM, HTOL,
     $   TTOL
      COMMON /PREC/       SMLNUM, ULP, TFNORM, T1NORM, HFNORM, HTOL, 
     $   TTOL
      
*     ..
*     .. Local Arrays ..
*     ..
      
*     ..
*     .. Externals ..
*     ..      
      EXTERNAL            LSAME, PDLAMCH, PDLANHS, PDLANTR,
     $   PDLASET, PXERBLA, ICEIL, NUMROC, INDXG2L, INDXG2P, 
     $   PDLANGE, PILAENVX
      DOUBLE PRECISION    PDLAMCH, PDLANHS, PDLANTR, PDLANGE
      INTEGER             INDXG2L, INDXG2P, ICEIL, NUMROC, PILAENVX
      LOGICAL             LSAME

*     ..
*     .. Intrinsic Functions ..
*     ..
      INTRINSIC           ABS, DBLE, MAX, MIN, SQRT, MOD, LOG, NINT
      
*     ..
*     .. Executable Statements ..
*     ..

*     Quick return if possible
      IF( N .LE. 0 .OR. IHI .LT. ILO ) THEN
         RETURN
      END IF

*     Extract current communication context and get grid parameters            
      ICTXT = DESCH(CTXT_)
      CALL BLACS_PINFO(IAM, NRPROCS)
      CALL BLACS_GRIDINFO( ICTXT, NPROW, NPCOL, MYROW, MYCOL )

      DEBUG = .FALSE.
      DEBUG = DEBUG .AND. IAM .EQ. 0
*     Check for workspace query     
      LQUERY = LDWORK .EQ. -1 .OR. LIWORK .EQ. -1


*     Define current problem size
      NH = IHI - ILO + 1


*     NMIN1 defines the minimum problem size to be solved by PDHGEQZ0. 
*     If problem size is less than this value, PDHGEQZ1 is called instead
      NMIN1 = PILAENVX(ICTXT, 50, 'PDHGEQZ', '', IDUM, IDUM, IDUM, IDUM)
      IF ( DEBUG ) WRITE(*,*)'% PDHGEQZ : NMIN1 = ', NMIN1
*     Decode JOB
      IF( LSAME( JOB, 'E' ) ) THEN
         ILSCHR = .FALSE.
         ISCHUR = 1
      ELSE IF( LSAME( JOB, 'S' ) ) THEN
         ILSCHR = .TRUE.
         ISCHUR = 2
      ELSE
         ISCHUR = 0
      END IF      

*     Decode COMPQ
      IF( LSAME( COMPQ, 'N' ) ) THEN
         ILQ = .FALSE.
         ICOMPQ = 1
      ELSE IF( LSAME( COMPQ, 'V' ) ) THEN
         ILQ = .TRUE.
         ICOMPQ = 2
      ELSE IF( LSAME( COMPQ, 'I' ) ) THEN
         ILQ = .TRUE.
         ICOMPQ = 3
      ELSE
         ICOMPQ = 0
      END IF

*     Decode COMPZ
      IF( LSAME( COMPZ, 'N' ) ) THEN
         ILZ = .FALSE.
         ICOMPZ = 1
      ELSE IF( LSAME( COMPZ, 'V' ) ) THEN
         ILZ = .TRUE.
         ICOMPZ = 2
      ELSE IF( LSAME( COMPZ, 'I' ) ) THEN
         ILZ = .TRUE.
         ICOMPZ = 3
      ELSE
         ICOMPZ = 0
      END IF

*     Check Argument Values
      INFO = 0
      IF( ISCHUR.EQ.0 ) THEN
         INFO = -1
      ELSE IF( ICOMPQ.EQ.0 ) THEN
         INFO = -2
      ELSE IF( ICOMPZ.EQ.0 ) THEN
         INFO = -3
      ELSE IF( N.LT.0 ) THEN
         INFO = -4
      ELSE IF( ILO.LT.1 ) THEN
         INFO = -5
      ELSE IF( IHI.GT.N .OR. IHI.LT.ILO-1 ) THEN
         INFO = -6
      ELSE IF( DESCH(NB_) .LT. 5 ) THEN
         INFO = -7
      ELSE IF ( DESCH(NB_) .NE. DESCH(MB_) ) THEN
         INFO = -8
      END IF

      IF( INFO.NE.0 ) THEN
         CALL PXERBLA(ICTXT, 'PDHGEQZ', -INFO )
         RETURN
      END IF



*     Workspace query call to PDHGEQZ0/PDHGEQZ1
      IF (NH .LE. NMIN1) THEN
         CALL PDHGEQZ1(ILSCHR, ILQ, ILZ, N, ILO, IHI, 
     $      H, DESCH, T, DESCT,
     $      ALPHAR, ALPHAI, BETA, Q, DESCQ, Z, DESCZ,
     $      DWORK, -1, IWORK, -1, IERR )
      ELSE
         CALL PDHGEQZ0(ILSCHR, ILQ, ILZ, N, ILO, IHI, 
     $      H, DESCH, T, DESCT,
     $      ALPHAR, ALPHAI, BETA, Q, DESCQ, Z, DESCZ,
     $      DWORK, -1, IWORK, -1, IERR, 0 )
      END IF

      LDWKOPT = INT( DWORK( 1 ) )
      LIWKOPT = IWORK( 1 )

*     Quick return if workspace query
      IF ( LQUERY ) THEN
         DWORK( 1 ) = LDWKOPT
         IWORK( 1 ) = LIWKOPT
         RETURN
      END IF
      
      
*     Check Argument Values (workspace)
      INFO = 0
      IF( LDWORK .LT. MAX( 1, LDWKOPT ) ) THEN
         INFO = -9
      ELSE IF( LIWORK .LT. MAX( 1, LIWKOPT ) ) THEN
         INFO = -10
      END IF

      IF( INFO.NE.0 ) THEN
         CALL PXERBLA(ICTXT, 'PDHGEQZ', -INFO )
         RETURN
      END IF

*     Machine constants
*     SMLNUM, ULP, TFNORM, T1NORM, HFNORM, HTOL, AND TTOL
*     are all common variables
      UNFL = PDLAMCH( ICTXT, 'SAFE MINIMUM' )
      OVFL = ONE / UNFL
      CALL PDLABAD( ICTXT, UNFL, OVFL )
      ULP = PDLAMCH( ICTXT, 'PRECISION' )
      SMLNUM = UNFL*( NH / ULP )
      HFNORM = PDLANGE( 'F', NH, NH, H, 1, 1, DESCH, DWORK )
      HTOL = MAX(UNFL, ULP*HFNORM)
      TFNORM = PDLANGE( 'F', NH, NH, T, 1, 1, DESCT, DWORK)
      T1NORM = TFNORM
      TTOL = MAX(UNFL, ULP*TFNORM)

      IF ( DEBUG ) THEN
         WRITE(*,*)'% Machine const: '
         WRITE(*,*)'% EPS =', ULP
         WRITE(*,*)'% HTOL =', HTOL
         WRITE(*,*)'% TTOL =', TTOL 
         WRITE(*,*)'% SMLNUM = ', SMLNUM
         WRITE(*,*)'% End Const.'         
      END IF
*     Initilize computed eigenvalues as not set
      ALPHAR( ILO : IHI ) = -ONE
      ALPHAI( ILO : IHI ) = -ONE
      BETA( ILO : IHI ) = -ONE     
      
*     Initialize Q and Z if desired
      IF( ICOMPQ.EQ.3 ) THEN
         CALL PDLASET( 'A', N, N, ZERO, ONE, Q, 1, 1, DESCQ )
      END IF
      IF( ICOMPZ.EQ.3 ) THEN
         CALL PDLASET( 'A', N, N, ZERO, ONE, Z, 1, 1, DESCZ )
      END IF

*     Call the appropriate routine
      IF (NH .LE. NMIN1) THEN
         IF ( DEBUG ) WRITE(*,*)'% PDHGEQZ : Calling PDHGEQZ1'
         CALL PDHGEQZ1( ILSCHR, ILQ, ILZ, N, ILO, IHI, 
     $      H, DESCH, T, DESCT,
     $      ALPHAR, ALPHAI, BETA, Q, DESCQ, Z, DESCZ,
     $      DWORK, LDWORK, IWORK, LIWORK, IERR ) 
      ELSE
         IF ( DEBUG ) WRITE(*,*)'% PDHGEQZ : Calling PDHGEQZ0'
         CALL PDHGEQZ0( ILSCHR, ILQ, ILZ, N, ILO, IHI, 
     $      H, DESCH, T, DESCT,
     $      ALPHAR, ALPHAI, BETA, Q, DESCQ, Z, DESCZ,
     $      DWORK, LDWORK, IWORK, LIWORK, IERR, 0 )
      END IF

      IF ( IERR .NE. 0 ) THEN
*        Abonormal termination
         INFO = IERR
      ELSE     
*        Normal Termination
         INFO = 0
      END IF
      RETURN
*     
*     End of PDHGEQZ
*     
      END SUBROUTINE





