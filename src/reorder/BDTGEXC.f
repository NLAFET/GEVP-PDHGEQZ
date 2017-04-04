      SUBROUTINE BDTGEXC( N, A, LDA, B, LDB, IFST, ILST, NITRAF, ITRAF,
     $                    NDTRAF, QTRAF, ZTRAF, WORK, LWORK, INFO )
	IMPLICIT NONE
*
*  -- LAPACK routine (version 3.0) --
*     Univ. of Tennessee, Univ. of California Berkeley, NAG Ltd.,
*     Courant Institute, Argonne National Lab, and Rice University
*     June 30, 1999
*
*     .. Scalar Arguments ..
      INTEGER            IFST, ILST, INFO, LDA, LDB, LWORK, N, NDTRAF,
     $                   NITRAF
*     ..
*     .. Array Arguments ..
      INTEGER            ITRAF( * )
      DOUBLE PRECISION   A( LDA, * ), B( LDB, * ), QTRAF( * ),
     $                   WORK( * ), ZTRAF( * )
      LOGICAL            DBG
*     ..
*
*  Purpose
*  =======
*
*  BDTGEXC reorders the generalized real Schur decomposition of a real
*  matrix pair (A,B) using an orthogonal equivalence transformation
*
*                 (A, B) = Q * (A, B) * Z'',
*
*  so that the diagonal block of (A, B) with row index IFST is moved
*  to row ILST.
*
*  (A, B) must be in generalized real Schur canonical form (as returned
*  by DGGES), i.e. A is block upper triangular with 1-by-1 and 2-by-2
*  diagonal blocks. B is upper triangular.
*
*
*  Arguments
*  =========
*
*  N       (input) INTEGER
*          The order of the matrices A and B. N >= 0.
*
*  A       (input/output) DOUBLE PRECISION array, dimension (LDA,N)
*          On entry, the matrix A in generalized real Schur canonical
*          form.
*          On exit, the updated matrix A, again in generalized
*          real Schur canonical form.
*
*  LDA     (input)  INTEGER
*          The leading dimension of the array A. LDA >= max(1,N).
*
*  B       (input/output) DOUBLE PRECISION array, dimension (LDB,N)
*          On entry, the matrix B in generalized real Schur canonical
*          form (A,B).
*          On exit, the updated matrix B, again in generalized
*          real Schur canonical form (A,B).
*
*  LDB     (input)  INTEGER
*          The leading dimension of the array B. LDB >= max(1,N).
*
*  IFST    (input/output) INTEGER
*  ILST    (input/output) INTEGER
*          Specify the reordering of the diagonal blocks of (A, B).
*          The block with row index IFST is moved to row ILST, by a
*          sequence of swapping between adjacent blocks.
*          On exit, if IFST pointed on entry to the second row of
*          a 2-by-2 block, it is changed to point to the first row;
*          ILST always points to the first row of the block in its
*          final position (which may differ from its input value by
*          +1 or -1). 1 <= IFST, ILST <= N.
*
*  NITRAF  (input/output) INTEGER
*          On entry, length of the array ITRAF.
*          As a minimum requirement, NITRAF >= max(1,|ILST-IFST|).
*          If there are 2-by-2 blocks in (A,B) then NITRAF should be
*          larger; a safe choice is NITRAF >= max(1,2*|ILST-IFST|).
*          On exit, actual length of the array ITRAF.
*
*  ITRAF   (output) INTEGER array, length NITRAF
*          List of parameters for representing the transformation
*          matrices Q and Z, see further details.
*
*  NDTRAF  (input/output) INTEGER
*          On entry, length of the arrays QTRAF and ZTRAF.
*          As a minimum requirement, NDTRAF >= max(1,2*|ILST-IFST|).
*          If there are 2-by-2 blocks in (A,B) then NDTRAF must be
*          larger; a safe choice is NDTRAF >= max(1,18*|ILST-IFST|).
*          On exit, actual length of the arrays QTRAF and ZTRAF.
*
*  QTRAF   (output) DOUBLE PRECISION array, length NDTRAF
*          List of parameters for representing the transformation
*          matrix Q, see further details.
*
*  ZTRAF   (output) DOUBLE PRECISION array, length NDTRAF
*          List of parameters for representing the transformation
*          matrix Z, see further details.
*
*  WORK    (workspace/output) DOUBLE PRECISION array, dimension (LWORK)
*          On exit, if INFO = 0, WORK(1) returns the optimal LWORK.
*
*  LWORK   (input) INTEGER
*          The dimension of the array WORK. LWORK >= 4*N + 16.
*
*          If LWORK = -1, then a workspace query is assumed; the routine
*          only calculates the optimal size of the WORK array, returns
*          this value as the first entry of the WORK array, and no error
*          message related to LWORK is issued by XERBLA.
*
*  INFO    (output) INTEGER
*           =0:  successful exit.
*           <0:  if INFO = -i, the i-th argument had an illegal value.
*           =1:  The transformed matrix pair (A, B) would be too far
*                from generalized Schur form; the problem is ill-
*                conditioned. (A, B) may have been partially reordered,
*                and ILST points to the first row of the current
*                position of the block being moved.
*
*  Further Details
*  ===============
*
*  The orthogonal transformation matrices Q and Z are products of
*  NITRAF orthogonal transformations; either Givens rotations or
*  3-by-3 / 4-by-4 orthogonal matrices. The parameters defining these
*  transformations are stored in the arrays ITRAF and DTRAF as follows:
*
*  Consider the i-th transformation acting on rows/columns POS,
*  POS+1, ... If this transformation is
*
*     (1) a Givens rotation with cosine C and sine S, then
*
*           ITRAF(i) = POS,
*           DTRAF(j) = C,    DTRAF(j+1) = S;
*
*     (2) a 3-by-3 orthogonal matrix L, then
*
*           ITRAF(i)    = N + POS,
*           DTRAF(j)    = L(1,1),     DTRAF(j+1) = L(2,1), ...,
*           DTRAF(j+8)  = L(3,3).
*
*     (3) a 4-by-4 orthogonal matrix L, then
*
*           ITRAF(i)    = 2*N + POS,
*           DTRAF(j)    = L(1,1),     DTRAF(j+1) = L(2,1), ...,
*           DTRAF(j+15) = L(4,4).
*
*  Note that the parameters in QTRAF and ZTRAF are stored
*  consecutively.
*
*
*  Based on contributions by
*     Bo Kagstrom and Peter Poromaa, Department of Computing Science,
*     Umea University, S-901 87 Umea, Sweden.
*
*  [1] B. Kagstrom; A Direct Method for Reordering Eigenvalues in the
*      Generalized Real Schur Form of a Regular Matrix Pair (A, B), in
*      M.S. Moonen et al (eds), Linear Algebra for Large Scale and
*      Real-Time Applications, Kluwer Academic Publ. 1993, pp 195-218.
*
*  =====================================================================
*
*     .. Parameters ..
      DOUBLE PRECISION   ZERO
      PARAMETER          ( ZERO = 0.0D+0 )
      INTEGER            DLNGTH(3)
*     ..
*     .. Local Scalars ..
      LOGICAL            LQUERY
      INTEGER            CDTRAF, CITRAF, LDTRAF, HERE, I, LWMIN,
     $                   NBF, NBL, NBNEXT
*     ..
*     .. External Subroutines ..
      EXTERNAL           BDTGEX2, XERBLA
*     ..
*     .. Intrinsic Functions ..
      INTRINSIC          MAX
*     .. Data Statements ..
      DATA ( DLNGTH(I), I = 1, 3 ) / 2, 9, 16 / 
*     ..
*     .. Executable Statements ..
*
*     Decode and test input arguments.
*
      DBG = .FALSE.
      INFO = 0
      LWMIN = MAX( 1, 4*N+16 )
      LQUERY = ( LWORK.EQ.-1 )
      IF( N.LT.0 ) THEN
         INFO = -1
      ELSE IF( LDA.LT.MAX( 1, N ) ) THEN
         INFO = -3
      ELSE IF( LDB.LT.MAX( 1, N ) ) THEN
         INFO = -5
      ELSE IF( IFST.LT.1 .OR. IFST.GT.N ) THEN
         INFO = -6
      ELSE IF( ILST.LT.1 .OR. ILST.GT.N ) THEN
         INFO = -7
      ELSE IF ( NITRAF.LT.MAX( 1, ABS( ILST-IFST ) ) ) THEN
         INFO = -8
      ELSE IF ( NDTRAF.LT.MAX( 1, 2*ABS( ILST-IFST ) ) ) THEN
         INFO = -10
      ELSE IF( LWORK.LT.LWMIN .AND. .NOT.LQUERY ) THEN
         INFO = -14
      END IF
*
      IF( INFO.EQ.0 ) THEN
         WORK( 1 ) = LWMIN
      END IF
*
      IF( INFO.NE.0 ) THEN
         CALL XERBLA( 'BDTGEXC', -INFO )
         RETURN
      ELSE IF( LQUERY ) THEN
         RETURN
      END IF
*
*     Quick return if possible
*
      IF( N.LE.1 )
     $   RETURN
      CITRAF = 1
      CDTRAF = 1
*
*     Determine the first row of the specified block and find out
*     if it is 1-by-1 or 2-by-2.
*
      IF( IFST.GT.1 ) THEN
         IF( A( IFST, IFST-1 ).NE.ZERO )
     $      IFST = IFST - 1
      END IF
      NBF = 1
      IF( IFST.LT.N ) THEN
         IF( A( IFST+1, IFST ).NE.ZERO )
     $      NBF = 2
      END IF
*
*     Determine the first row of the final block
*     and find out if it is 1-by-1 or 2-by-2.
*
      IF( ILST.GT.1 ) THEN
         IF( A( ILST, ILST-1 ).NE.ZERO )
     $      ILST = ILST - 1
      END IF
      NBL = 1
      IF( ILST.LT.N ) THEN
         IF( A( ILST+1, ILST ).NE.ZERO )
     $      NBL = 2
      END IF
      IF( IFST.EQ.ILST )
     $   RETURN
*

      IF (DBG) WRITE(*,*)IFST, ILST, NBF, NBL
      IF( IFST.LT.ILST ) THEN
*
*        Update ILST.
*
         IF( NBF.EQ.2 .AND. NBL.EQ.1 )
     $      ILST = ILST - 1
         IF( NBF.EQ.1 .AND. NBL.EQ.2 )
     $      ILST = ILST + 1
*
         HERE = IFST
*
   10    CONTINUE
         IF (DBG) WRITE(*,*)'10', HERE, ILST, NBF, NBL 
*
*        Swap with next one below.
*
         IF( NBF.EQ.1 .OR. NBF.EQ.2 ) THEN
*
*           Current block either 1-by-1 or 2-by-2.
*
            NBNEXT = 1
            IF( HERE+NBF+1.LE.N ) THEN
               IF( A( HERE+NBF+1, HERE+NBF ).NE.ZERO )
     $            NBNEXT = 2
            END IF
*
            LDTRAF = DLNGTH(NBF+NBNEXT-1)
            IF( CITRAF.GT.NITRAF ) THEN
               INFO = -8
               CALL XERBLA( 'BDTGEXC', -INFO )
               RETURN
            END IF
            IF( CDTRAF+LDTRAF-1.GT.NDTRAF ) THEN
               INFO = -10
               CALL XERBLA( 'BDTGEXC', -INFO )
               RETURN
            END IF
            CALL BDTGEX2( N, A, LDA, B, LDB, HERE, NBF, NBNEXT,
     $                    ITRAF(CITRAF), QTRAF(CDTRAF), ZTRAF(CDTRAF),
     $                    WORK, LWORK, INFO )
            IF( INFO.NE.0 ) THEN
               ILST = HERE
               NITRAF = CITRAF - 1
               NDTRAF = CDTRAF - 1
               RETURN
            END IF
            CITRAF = CITRAF + 1
            CDTRAF = CDTRAF + LDTRAF
            HERE = HERE + NBNEXT
*
*           Test if 2-by-2 block breaks into two 1-by-1 blocks.
*
            IF( NBF.EQ.2 ) THEN
               IF( A( HERE+1, HERE ).EQ.ZERO )
     $            NBF = 3
            END IF
*
         ELSE
*
*           Current block consists of two 1-by-1 blocks, each of which
*           must be swapped individually.
*
            NBNEXT = 1
            IF( HERE+3.LE.N ) THEN
               IF( A( HERE+3, HERE+2 ).NE.ZERO )
     $            NBNEXT = 2
            END IF
            LDTRAF = DLNGTH(NBNEXT)
            IF( CITRAF.GT.NITRAF ) THEN
               INFO = -8
               CALL XERBLA( 'BDTGEXC', -INFO )
               RETURN
            END IF
            IF( CDTRAF+LDTRAF-1.GT.NDTRAF ) THEN
               INFO = -10
               CALL XERBLA( 'BDTGEXC', -INFO )
               RETURN
            END IF
            CALL BDTGEX2( N, A, LDA, B, LDB, HERE+1, 1, NBNEXT,
     $                    ITRAF(CITRAF), QTRAF(CDTRAF), ZTRAF(CDTRAF),
     $                    WORK, LWORK, INFO )
            IF( INFO.NE.0 ) THEN
               NITRAF = CITRAF - 1
               NDTRAF = CDTRAF - 1
               ILST = HERE
               RETURN
            END IF
            CITRAF = CITRAF + 1
            CDTRAF = CDTRAF + LDTRAF
*
            IF( NBNEXT.EQ.1 ) THEN
*
*              Swap two 1-by-1 blocks.
*
               LDTRAF = DLNGTH(1)
               IF( CITRAF.GT.NITRAF ) THEN
                  INFO = -8
                  CALL XERBLA( 'BDTGEXC', -INFO )
                  RETURN
               END IF
               IF( CDTRAF+LDTRAF-1.GT.NDTRAF ) THEN
                  INFO = -10
                  CALL XERBLA( 'BDTGEXC', -INFO )
                  RETURN
               END IF
               CALL BDTGEX2( N, A, LDA, B, LDB, HERE, 1, 1,
     $                       ITRAF(CITRAF), QTRAF(CDTRAF),
     $                       ZTRAF(CDTRAF), WORK, LWORK, INFO )
               IF( INFO.NE.0 ) THEN
                  NITRAF = CITRAF - 1
                  NDTRAF = CDTRAF - 1
	            INFO = 2
                  ILST = HERE
                  RETURN
               END IF
               CITRAF = CITRAF + 1
               CDTRAF = CDTRAF + LDTRAF
               HERE = HERE + 1
*
            ELSE
*
*              Recompute NBNEXT in case of 2-by-2 split.
*
               IF( A( HERE+2, HERE+1 ).EQ.ZERO )
     $            NBNEXT = 1
               IF( NBNEXT.EQ.2 ) THEN
*
*                 2-by-2 block did not split.
*
                  LDTRAF = DLNGTH(2)
                  IF( CITRAF.GT.NITRAF ) THEN
                     INFO = -8
                     CALL XERBLA( 'BDTGEXC', -INFO )
                     RETURN
                  END IF
                  IF( CDTRAF+LDTRAF-1.GT.NDTRAF ) THEN
                     INFO = -10
                     CALL XERBLA( 'BDTGEXC', -INFO )
                     RETURN
                  END IF
                  CALL BDTGEX2( N, A, LDA, B, LDB, HERE, 1, NBNEXT,
     $                          ITRAF(CITRAF), QTRAF(CDTRAF),
     $                          ZTRAF(CDTRAF), WORK, LWORK, INFO )
                  IF( INFO.NE.0 ) THEN
                     NITRAF = CITRAF - 1
                     NDTRAF = CDTRAF - 1
	               INFO = 2
                     ILST = HERE
                     RETURN
                  END IF
                  CITRAF = CITRAF + 1
                  CDTRAF = CDTRAF + LDTRAF
                  HERE = HERE + 2
               ELSE
*
*                 2-by-2 block did split.
*
                  LDTRAF = DLNGTH(1)
                  IF( CITRAF.GT.NITRAF ) THEN
                     INFO = -8
                     CALL XERBLA( 'BDTGEXC', -INFO )
                     RETURN
                  END IF
                  IF( CDTRAF+LDTRAF-1.GT.NDTRAF ) THEN
                     INFO = -10
                     CALL XERBLA( 'BDTGEXC', -INFO )
                     RETURN
                  END IF
                  CALL BDTGEX2( N, A, LDA, B, LDB, HERE, 1, 1,
     $                          ITRAF(CITRAF), QTRAF(CDTRAF),
     $                          ZTRAF(CDTRAF), WORK, LWORK, INFO )
                  IF( INFO.NE.0 ) THEN
                     NITRAF = CITRAF - 1
                     NDTRAF = CDTRAF - 1
	               INFO = 2
                     ILST = HERE
                     RETURN
                  END IF
                  CITRAF = CITRAF + 1
                  CDTRAF = CDTRAF + LDTRAF
                  HERE = HERE + 1
*
                  LDTRAF = DLNGTH(1)
                  IF( CITRAF.GT.NITRAF ) THEN
                     INFO = -8
                     CALL XERBLA( 'BDTGEXC', -INFO )
                     RETURN
                  END IF
                  IF( CDTRAF+LDTRAF-1.GT.NDTRAF ) THEN
                     INFO = -10
                     CALL XERBLA( 'BDTGEXC', -INFO )
                     RETURN
                  END IF
                  CALL BDTGEX2( N, A, LDA, B, LDB, HERE, 1, 1,
     $                          ITRAF(CITRAF), QTRAF(CDTRAF),
     $                          ZTRAF(CDTRAF), WORK, LWORK, INFO )
                  IF( INFO.NE.0 ) THEN
                     NITRAF = CITRAF - 1
                     NDTRAF = CDTRAF - 1
	               INFO = 3
                     ILST = HERE
                     RETURN
                  END IF
                  CITRAF = CITRAF + 1
                  CDTRAF = CDTRAF + LDTRAF
                  HERE = HERE + 1
               END IF
*
            END IF
         END IF
         IF( HERE.LT.ILST )
     $      GO TO 10
      ELSE
         HERE = IFST
*
   20    CONTINUE
         IF (DBG) WRITE(*,*)'20', HERE, ILST, NBF, NBL
*
*        Swap with next one below.
*
         IF( NBF.EQ.1 .OR. NBF.EQ.2 ) THEN
*
*           Current block either 1-by-1 or 2-by-2.
*
            NBNEXT = 1
            IF( HERE.GE.3 ) THEN
               IF( A( HERE-1, HERE-2 ).NE.ZERO )
     $            NBNEXT = 2
            END IF
*           
            IF (DBG) WRITE(*,*)'NBF=',NBF, NBNEXT
            LDTRAF = DLNGTH(NBF+NBNEXT-1)
            IF( CITRAF.GT.NITRAF ) THEN
               INFO = -8
               CALL XERBLA( 'BDTGEXC', -INFO )
               RETURN
            END IF
            IF( CDTRAF+LDTRAF-1.GT.NDTRAF ) THEN
               INFO = -10
               CALL XERBLA( 'BDTGEXC', -INFO )
               RETURN
            END IF
            CALL BDTGEX2( N, A, LDA, B, LDB, HERE-NBNEXT, NBNEXT, NBF,
     $                    ITRAF(CITRAF), QTRAF(CDTRAF), ZTRAF(CDTRAF),
     $                    WORK, LWORK, INFO )
            IF( INFO.NE.0 ) THEN
               ILST = HERE
               NITRAF = CITRAF - 1
               NDTRAF = CDTRAF - 1
               RETURN
            END IF
            CITRAF = CITRAF + 1
            CDTRAF = CDTRAF + LDTRAF

            HERE = HERE - NBNEXT
*
*           Test if 2-by-2 block breaks into two 1-by-1 blocks.
*
            IF( NBF.EQ.2 ) THEN
               IF( A( HERE+1, HERE ).EQ.ZERO )
     $            NBF = 3
            END IF
*
         ELSE
*
*           Current block consists of two 1-by-1 blocks, each of which
*           must be swapped individually.
*
            IF (DBG) WRITE(*,*)'NBF=',NBF, NBNEXT
            NBNEXT = 1
            IF( HERE.GE.3 ) THEN
               IF( A( HERE-1, HERE-2 ).NE.ZERO )
     $            NBNEXT = 2
            END IF
*
            LDTRAF = DLNGTH(NBNEXT)
            IF( CITRAF.GT.NITRAF ) THEN
               INFO = -8
               CALL XERBLA( 'BDTGEXC', -INFO )
               RETURN
            END IF
            IF( CDTRAF+LDTRAF-1.GT.NDTRAF ) THEN
               INFO = -10
               CALL XERBLA( 'BDTGEXC', -INFO )
               RETURN
            END IF
            CALL BDTGEX2( N, A, LDA, B, LDB, HERE-NBNEXT, NBNEXT, 1,
     $                    ITRAF(CITRAF), QTRAF(CDTRAF), ZTRAF(CDTRAF),
     $                    WORK, LWORK, INFO )
            IF( INFO.NE.0 ) THEN
               ILST = HERE
               NITRAF = CITRAF - 1
               NDTRAF = CDTRAF - 1
               RETURN
            END IF
            CITRAF = CITRAF + 1
            CDTRAF = CDTRAF + LDTRAF
*
            IF( NBNEXT.EQ.1 ) THEN
*
*              Swap two 1-by-1 blocks.
*
               IF (DBG) WRITE(*,*)'Swap 1by1 blocks'  
               LDTRAF = DLNGTH(1)
               IF( CITRAF.GT.NITRAF ) THEN
                  INFO = -8
                  CALL XERBLA( 'BDTGEXC', -INFO )
                  RETURN
               END IF
               IF( CDTRAF+LDTRAF-1.GT.NDTRAF ) THEN
                  INFO = -10
                  CALL XERBLA( 'BDTGEXC', -INFO )
                  RETURN
               END IF
               CALL BDTGEX2( N, A, LDA, B, LDB, HERE, NBNEXT, 1,
     $                       ITRAF(CITRAF), QTRAF(CDTRAF),
     $                       ZTRAF(CDTRAF), WORK, LWORK, INFO )
               IF( INFO.NE.0 ) THEN
                  NITRAF = CITRAF - 1
                  NDTRAF = CDTRAF - 1
                  ILST = HERE
                  RETURN
               END IF
               CITRAF = CITRAF + 1
               CDTRAF = CDTRAF + LDTRAF
               HERE = HERE - 1
            ELSE
*
*             Recompute NBNEXT in case of 2-by-2 split.
*
               IF (DBG) WRITE(*,*)'Recomp NBNEXT'  
               IF( A( HERE, HERE-1 ).EQ.ZERO )
     $            NBNEXT = 1
               IF (DBG) WRITE(*,*)'Recomp NBNEXT',NBNEXT 
               IF( NBNEXT.EQ.2 ) THEN
*
*                 2-by-2 block did not split.
*
                  IF (DBG) WRITE(*,*)'No split, Here=',HERE  
                  LDTRAF = DLNGTH(2)
                  IF( CITRAF.GT.NITRAF ) THEN
                     INFO = -8
                     CALL XERBLA( 'BDTGEXC', -INFO )
                     RETURN
                  END IF
                  IF( CDTRAF+LDTRAF-1.GT.NDTRAF ) THEN
                     INFO = -10
                     CALL XERBLA( 'BDTGEXC', -INFO )
                     RETURN
                  END IF
                  CALL BDTGEX2( N, A, LDA, B, LDB, HERE-1, 2, 1,
     $                          ITRAF(CITRAF), QTRAF(CDTRAF),
     $                          ZTRAF(CDTRAF), WORK, LWORK, INFO )
                  IF( INFO.NE.0 ) THEN
	               INFO = 2
                     ILST = HERE
                     NITRAF = CITRAF - 1
                     NDTRAF = CDTRAF - 1
                     RETURN
                  END IF
                  CITRAF = CITRAF + 1
                  CDTRAF = CDTRAF + LDTRAF
                  HERE = HERE - 2
               ELSE
*
*                 2-by-2 block did split.
*
                  IF (DBG) WRITE(*,*)'Split, Here=',HERE  
                  LDTRAF = DLNGTH(1)
                  IF( CITRAF.GT.NITRAF ) THEN
                     INFO = -8
                     CALL XERBLA( 'BDTGEXC', -INFO )
                     RETURN
                  END IF
                  IF( CDTRAF+LDTRAF-1.GT.NDTRAF ) THEN
                     INFO = -10
                     CALL XERBLA( 'BDTGEXC', -INFO )
                     RETURN
                  END IF
                  CALL BDTGEX2( N, A, LDA, B, LDB, HERE, 1, 1,
     $                          ITRAF(CITRAF), QTRAF(CDTRAF),
     $                          ZTRAF(CDTRAF), WORK, LWORK, INFO )
                  IF( INFO.NE.0 ) THEN
	               INFO = 2
                     NITRAF = CITRAF - 1
                     NDTRAF = CDTRAF - 1
                     ILST = HERE
                     RETURN
                  END IF
                  CITRAF = CITRAF + 1
                  CDTRAF = CDTRAF + LDTRAF
                  HERE = HERE - 1
*
                  IF( CITRAF.GT.NITRAF ) THEN
                     INFO = -8
                     CALL XERBLA( 'BDTGEXC', -INFO )
                     RETURN
                  END IF
                  IF( CDTRAF+LDTRAF-1.GT.NDTRAF ) THEN
                     INFO = -10
                     CALL XERBLA( 'BDTGEXC', -INFO )
                     RETURN
                  END IF
                  CALL BDTGEX2( N, A, LDA, B, LDB, HERE, 1, 1,
     $                          ITRAF(CITRAF), QTRAF(CDTRAF),
     $                          ZTRAF(CDTRAF), WORK, LWORK, INFO )
                  IF( INFO.NE.0 ) THEN
	               INFO = 3
                     NITRAF = CITRAF - 1
                     NDTRAF = CDTRAF - 1
                     ILST = HERE
                     RETURN
                  END IF
                  CITRAF = CITRAF + 1
                  CDTRAF = CDTRAF + LDTRAF
                  HERE = HERE - 1
               END IF
            END IF
         END IF
         IF( HERE.GT.ILST )
     $      GO TO 20
      END IF
      ILST = HERE
      NITRAF = CITRAF - 1
      NDTRAF = CDTRAF - 1
      WORK( 1 ) = LWMIN
      IF (DBG) STOP
      RETURN
*
*     End of BDTGEXC
*
      END
