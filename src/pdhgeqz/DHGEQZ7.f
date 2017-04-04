***********************************************************************
*                                                                     *
*     DHGEQZ7.f:                                                      *
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
      SUBROUTINE DHGEQZ7( DIR, M, N, H, LDH, T, LDT, Q, Z, LDQZ, SELECT, 
     $   FIRSTLASTBLK, MOVECNT )
      IMPLICIT NONE
*     ..
*     .. Scalar Arguments ..
*     ..
      CHARACTER           DIR
      INTEGER             LDH, LDT, N, M, LDQZ, MOVECNT
      LOGICAL             FIRSTLASTBLK
*     ..
*     .. Array Arguments ..
*     ..
      DOUBLE PRECISION    H( LDH, * ), T( LDT, * ), Q( LDQZ, * ), 
     $   Z( LDQZ, * )
      INTEGER             SELECT( M )

*     Purpose
*     =======
*     
*     DHGEQZ7 chases diagonal zeros elements in T, indicated by SELECT(i) = 1, 
*     down or up to the end/begining of T. H, T, Q, and Z are of size M*N.
*     On exit, MOVECNT is set to the number of moved zeros.a
*     
*     The SELECT vector is updated to reflect moved elements.
*     
*     Auxiliary routine - all inputs are assumed valid without checking
*     
*     Arguments 
*     ========= 
*     
*     DIR     (input) CHARACTER*1
*     = 'U', zeros on the diagonal of T are chased up.
*     <>'U', zeros on the diaognal of T are chased down.
*     
*     M       (input) INTEGER
*     Number of rows in H and T.
*     Number of columns in Q and Z.
*     N       (input) INTEGER
*     Number of columns in H and T.
*     Number of rows in Q and Z.
*     H       (input) DOUBLE PRECISION array of dimension ( LDH, * )
*     Hessenberg matrix which is updated when a zero in T is 
*     moved up or down.
*     
*     LDH     (input) INTEGER
*     Leading dimension of the matrix H.
*     
*     T       (input/output) DOUBLE PRECISION array of dimension (LDT,*)
*     Triangular matrix with potential zeros on the diagonal 
*     to be chased either up or down.
*     
*     LDT     (input) INTEGER
*     Leading dimension of the matrix T.
*     
*     Q       (input/output) DOUBLE PRECISION array of dimension ( LDQZ, * )
*     Q is multiplied, from right, with all transformations
*     applied from left on H and T.          
*     
*     Z       (input/output) DOUBLE PRECISION array of dimension ( LDQZ, * )
*     Z is multiplied, from right, with all transformations
*     applied from right on H and T.   
*     
*     LDQZ    (input) INTEGER 
*     Leading dimension of the matrices Q and Z.
*     
*     SELECT  (input/output) INTEGER array of dimension (M)
*     if SELECT (I) = 1, then a zero is found at T(I,I)
*     SELECT is maintained to reflect any movements.
*     
*     FIRSTLASTBLK (input) LOGICAL
*     If .TRUE. deflation is performed, 
*     otherwise zeros are moved as close to matrix boundaries 
*     as possible.
*     
*     MOVECNT (output) INTEGER
*     On exit, MOVECNT holds the number of moved zeros. A zero
*     value can mean no zeros have been found, or they cant 
*     be moved due to they are already packed as much as possible.
*     
*     =============================================================== 

*     ..
*     .. Parameters ..
*     ..
      DOUBLE PRECISION    ZERO, ONE
      PARAMETER           ( ZERO = 0.0D+0, ONE = 1.0D+0 )
*     ..
*     .. Local Scalars ..
*     ..
      INTEGER             I, J, JLAST, JFIRST, LASTSAFEENTRY, 
     $   FIRSTSAFEENTRY
      DOUBLE PRECISION    C, S, TEMP

*     ..
*     .. Local Arrays ..
*     ..

*     ..
*     .. Externals 
*     ..      
      EXTERNAL            LSAME
      LOGICAL             LSAME
*     ..
*     .. Precision paramters ..
*     ..
      DOUBLE PRECISION    SMLNUM, ULP, TFNORM, T1NORM, HFNORM, HTOL,
     $   TTOL
      COMMON /PREC/       SMLNUM, ULP, TFNORM, T1NORM, HFNORM, HTOL,
     $   TTOL

      CALL DLASET( 'A', M, N, ZERO, ONE, Q, LDQZ )
      CALL DLASET( 'A', M, N, ZERO, ONE, Z, LDQZ )

      JFIRST = 1   
      JLAST = M
      MOVECNT = 0
      
*     
*     DIR, direction to push zeros, is assumed to be either U=Up or D=Down
*     
      IF ( LSAME( DIR, 'U' ) ) GOTO 40
      
      LASTSAFEENTRY = JLAST - 2
      IF ( FIRSTLASTBLK ) THEN
         LASTSAFEENTRY = JLAST - 1
      END IF 
      
      DO 20 I = M, 2, -1
*        
*        Find the last zero element in T above last element
*        
         IF ( SELECT( I ) .EQ. 1 .AND. I .LE. LASTSAFEENTRY) THEN
*           
*           Keep track of how many elements we have moved
*           
            MOVECNT = MOVECNT + 1
            
            

*           
*           Move T(I, I) == 0 to T(JLAST, JLAST)
*           
            DO 10 J = I, LASTSAFEENTRY              
*              
*              Create a row rotation that annihilates diagonal entry T( J + 1, J + 1 ) 
*              
               TEMP = T( J, J + 1 )
               CALL DLARTG( TEMP, T( J + 1, J + 1 ), C, S,
     $            T( J, J + 1 ) )
               T( J + 1, J + 1 ) = ZERO
*              
*              Update the select vector to indicate the newly created 0
*              and mark previous entry as to be set as non zero
*              IS THIS ALWAYS TRUE ?
*              
               SELECT ( J ) = 0
               SELECT ( J + 1 ) = 1
*              
*              Apply the row rotation to T, H and Q.
*              
               IF( J .LT. N - 1 ) THEN
                  CALL DROT( N - J - 1, T( J, J + 2 ), LDT,
     $               T( J + 1, J + 2 ), LDT, C, S )
               END IF
               
               CALL DROT( N - J + 2, H( J, J - 1 ), LDH,
     $            H( J + 1, J - 1 ), LDH, C, S )

               CALL DROT( M, Q( 1, J ), 1, Q( 1, J + 1 ), 1,
     $            C, S )
               
*              
*              Create a column rotation that annihilates the fill-in 
*              in H created by the row rotation above.
*              Annihilated element is H (J + 1, J - 1)
*              
               TEMP = H( J + 1, J )
               CALL DLARTG( TEMP, H( J + 1, J - 1 ), C, S,
     $            H( J + 1, J ) )                  
               H( J + 1, J - 1 ) = ZERO
*              
*              Apply the column rotation to H, T and Z. This won't cause 
*              any new fill-in outside diagonal of T.                  
*              
               CALL DROT( J , H( 1, J ), 1,
     $            H( 1, J - 1 ), 1, C, S )
               
               CALL DROT( J - 1, T( 1, J ), 1,
     $            T( 1, J - 1 ), 1, C, S )
               IF ( ABS(T( J - 1, J - 1 )) .GT. TTOL ) SELECT(J - 1) = 0
               
               CALL DROT( M, Z( 1, J ), 1, Z( 1, J - 1 ), 1,
     $            C, S )
 10         CONTINUE
            
*           
*           If this is the last block in the distributed (H,T), 
*           annihilate the last subdiagonal in H and apply the column-rotation 
*           on remaining local (H, T). This leads to a converged infinite 
*           eigenvalue.              
*           
            IF ( FIRSTLASTBLK ) THEN       
*              
*              Construct a column rotation that annihilates H( M, M -  1)                  
*              
               TEMP = H( JLAST, JLAST )                  
               CALL DLARTG( TEMP, H( JLAST, JLAST - 1 ), C, S,
     $            H( JLAST, JLAST ) )                  
               H( JLAST, JLAST-1 ) = ZERO
*              
*              Apply the column rotation to H, T, and Z. This won't 
*              cause any fill-in in T (since T(ILAST, ILAST) is zero). 
*              
               CALL DROT( JLAST - 1, H( 1, JLAST ), 1,
     $            H( 1, JLAST - 1 ), 1, C, S )
               
               CALL DROT( JLAST - 1, T( 1, JLAST ), 1,
     $            T( 1, JLAST - 1 ), 1, C, S )
               
               CALL DROT( M, Z( 1, JLAST ), 1, 
     $            Z( 1, JLAST - 1 ), 1, C, S )
               JLAST = JLAST - 1   
               LASTSAFEENTRY = LASTSAFEENTRY - 1
            ELSE 
               LASTSAFEENTRY = LASTSAFEENTRY - 2
            END IF
         ELSE IF ( SELECT( I ) .EQ. 1 ) THEN
            LASTSAFEENTRY = LASTSAFEENTRY - 2
         END IF
 20   CONTINUE
      GOTO 100

      
 40   CONTINUE

      FIRSTSAFEENTRY = JFIRST + 2
      IF ( FIRSTLASTBLK ) THEN
         FIRSTSAFEENTRY = JFIRST + 1
      END IF

      DO 50 I = 2, M - 1
*        
*        Find the first zero element in T above first element
*        
         IF ( SELECT( I ) .EQ. 1 .AND. I .GE. FIRSTSAFEENTRY) THEN
*           
*           Keep track of how many elements we have moved
*           
            MOVECNT = MOVECNT + 1
*           
*           Move T(I, I) == 0 to T(JFIRST, JFIRST)
*           
            DO 60 J = I, FIRSTSAFEENTRY, -1
*              
*              Create a column rotation that annihilates 
*              diagonal entry T( J - 1, J - 1 ) 
*              
               TEMP = T( J - 1, J )
               CALL DLARTG( TEMP, T( J - 1, J - 1 ), C, S,
     $            T( J - 1, J ) )
               T( J - 1, J - 1 ) = ZERO
*              
*              Update the select vector to indicate the newly created 0
*              and mark previous entry as to be set as non zero
*              IS THIS ALWAYS TRUE ?
*              
               SELECT ( J ) = 0
               SELECT ( J - 1 ) = 1
*              
*              Apply the column rotation to T, H and Z.
*              
               IF( J .GT. 2 ) THEN
                  CALL DROT( J - 2, T( 1, J ), 1,
     $               T( 1, J - 1 ), 1, C, S )
               END IF
               CALL DROT( MIN( J + 1, M ), H( 1, J ), 1,
     $            H( 1, J - 1 ), 1, C, S )

               CALL DROT( M, Z( 1, J ), 1, Z( 1, J - 1 ), 1,
     $            C, S )
*              
*              Create a row rotation that annihilates the fill-in 
*              in H created by the column rotation above.
*              Annihilated element is H (J + 1, J - 1)
*              
               TEMP = H( J, J - 1)
               CALL DLARTG( TEMP, H( J + 1, J - 1 ), C, S,
     $            H( J, J - 1 ) )
               H( J + 1, J - 1 ) = ZERO
*              
*              Apply the row rotation to H, T and Q. This won't cause 
*              any new fill-in outside diagonal of T.                
*              
               CALL DROT( N - J + 1, H( J, J ), LDH,
     $            H( J + 1, J ), LDH, C, S )

               CALL DROT( N - J, T( J, J + 1 ), LDT,
     $            T( J + 1, J + 1 ), LDT, C, S )
               IF ( ABS(T( J + 1, J + 1)) .GT. TTOL) SELECT (J + 1) = 0
               CALL DROT( M, Q( 1, J ), 1, Q( 1, J + 1 ), 1,
     $            C, S )
 60         CONTINUE
*           
*           If this is the first block in the distributed (H,T), 
*           annihilate the first subdiagonal in H and apply the row-rotation 
*           on remaining local (H, T). This leads to a converged infinite 
*           eigenvalue.              
*           
            IF ( FIRSTLASTBLK ) THEN
*              
*              Construct a row rotation that annihilates H( JFIRST + 1, JFIRST)                  
*              
               TEMP = H( JFIRST, JFIRST )
               CALL DLARTG( TEMP, H( JFIRST + 1, JFIRST ), C, S,
     $            H( JFIRST, JFIRST ) )
               H( JFIRST + 1, JFIRST ) = ZERO
*              
*              Apply the row rotation to H, T, and Q. This won't 
*              cause any fill-in in T (since T(JFIRST, JFIRST) is zero).                  
*              
               CALL DROT( N - JFIRST, H( JFIRST, JFIRST + 1), LDH,
     $            H( JFIRST + 1, JFIRST + 1 ), LDH, C, S )

               CALL DROT( N - JFIRST, T( JFIRST, JFIRST + 1 ), LDT,
     $            T( JFIRST + 1, JFIRST + 1 ), LDT, C, S )

               CALL DROT( M, Q( 1, JFIRST ), 1,
     $            Q( 1, JFIRST + 1 ), 1, C, S )
               JFIRST = JFIRST + 1
               FIRSTSAFEENTRY = FIRSTSAFEENTRY + 1
            ELSE
               FIRSTSAFEENTRY = FIRSTSAFEENTRY + 2
            END IF
         ELSE IF ( SELECT( I ) .EQ. 1 ) THEN
            FIRSTSAFEENTRY = FIRSTSAFEENTRY + 2
         END IF
 50   CONTINUE


 100  CONTINUE
      
      RETURN
      END
     


