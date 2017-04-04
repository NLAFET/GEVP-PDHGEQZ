      SUBROUTINE QZCONF( FNAME )
      IMPLICIT NONE
C     
C     Reads the parameters for the QZ algorithm from a text file called
C     "FNAME" and stores them in COMMON blocks called "PARA_IT",
C     "PARA_DEF"´ and "PARA_INF".
C     
      CHARACTER         FNAME*255
      INTEGER           IOS
C     
C     Parameters
C     
      LOGICAL           PARA_ITTRI
      INTEGER           PARA_ITASH, PARA_ITBULG, PARA_ITDACK, 
     $   PARA_ITDANB, PARA_ITEXC, PARA_ITEXT, PARA_ITKK, PARA_ITMAX, 
     $   PARA_ITNS, PARA_ITOPP
      COMMON /PARA_IT/  PARA_ITASH, PARA_ITBULG, PARA_ITDACK, 
     $   PARA_ITDANB, PARA_ITEXC, PARA_ITEXT, PARA_ITKK, PARA_ITMAX, 
     $   PARA_ITNS, PARA_ITOPP, PARA_ITTRI
      
      INTEGER           PARA_AGGCRIT, PARA_AGGSHF, PARA_AGGMIN,
     $   PARA_AGGWIN, PARA_HDEF
      
      LOGICAL           PARA_AGGRQUP
      DOUBLE PRECISION  PARA_AGGREP   
      
      COMMON /PARA_DEF/ PARA_AGGCRIT,
     $   PARA_AGGSHF, PARA_AGGMIN, PARA_AGGWIN, PARA_HDEF, PARA_AGGRQUP,
     $   PARA_AGGREP    

      LOGICAL           PARA_INFALL, PARA_INFBLOCK
      INTEGER           PARA_INFDEF, PARA_INFDIR, PARA_INFNZ,
     $   PARA_INFWDW
      DOUBLE PRECISION  PARA_INFTOL
      COMMON /PARA_INF/ PARA_INFALL, PARA_INFBLOCK, PARA_INFDEF, 
     $   PARA_INFDIR, PARA_INFNZ, PARA_INFWDW, PARA_INFTOL

      INTEGER           PARA_LDADD
      COMMON /PARA_MIS/ PARA_LDADD
C     
      OPEN( 1, IOSTAT = IOS, STATUS = 'OLD', FILE = FNAME )
      IF ( IOS.NE.0 ) THEN
         PRINT*, 'Error in QZCONF: Parameter file not found'
      END IF
C     
C     Skip first line.
C     
      READ ( 1, FMT = '()', IOSTAT = IOS )
C     
C     Settings for the multishift QZ algorithm.
C     
C     PARA_ITKK       = minimal size addressed by chains of bulges
      READ ( 1, FMT = '()', IOSTAT = IOS )
      READ ( 1, FMT = *,  IOSTAT = IOS )   PARA_ITKK
      IF ( PARA_ITKK.LE.2 ) THEN
         PRINT*, 'Error in QZCONF: PARA_ITKK must be an integer > 2.'
      END IF
C     PARA_ITDACK     = minimal size addressed by Dackland/Kagstrom QZ
      READ ( 1, FMT = '()', IOSTAT = IOS )
      READ ( 1, FMT = *,  IOSTAT = IOS )   PARA_ITDACK
      IF ( PARA_ITDACK.LE.2 ) THEN
         PRINT*, 'Error in QZCONF: PARA_ITDACK must be an integer > 2.'
      END IF
      IF ( PARA_ITDACK.GT.PARA_ITKK ) THEN
         PRINT*, 'Error in QZCONF: PARA_ITDACK should not be larger
     $      than PARA_ITKK.'
      END IF
C     PARA_ITASH      = usage of shifts from early deflation
      READ ( 1, FMT = '()', IOSTAT = IOS )
      READ ( 1, FMT = *,  IOSTAT = IOS )   PARA_ITASH
      IF ( PARA_ITASH.LT.0 .OR. PARA_ITASH.GT.2 ) THEN
         PRINT*, 'Error in QZCONF: PARA_ITASH must be in [0,2].'
      END IF
C     PARA_ITBULG     = number of shifts in each bulge
      READ ( 1, FMT = '()', IOSTAT = IOS )
      READ ( 1, FMT = *,  IOSTAT = IOS )   PARA_ITBULG
      IF ( PARA_ITBULG.LT.2 .OR. MOD( PARA_ITBULG, 2 ).NE.0 ) THEN
         PRINT*, 'Error in QZCONF: PARA_ITBULG must be an even integer >
     $      =2.'
      END IF
C     PARA_ITNS       = number of shifts in block QZ iteration
      READ ( 1, FMT = '()', IOSTAT = IOS )
      READ ( 1, FMT = *,  IOSTAT = IOS )   PARA_ITNS
      IF ( PARA_ITNS.LT.2 .OR. MOD( PARA_ITNS, PARA_ITBULG ).NE.0 ) THEN
         PRINT*, 'Error in QZCONF: PARA_ITNS must be an integer multiple
     $      of PARA_ITBULG.'
      END IF
      IF ( PARA_ITNS.LT.PARA_ITBULG*PARA_ITBULG ) THEN
         PRINT*, 'Error in QZCONF: PARA_ITNS must not be smaller than
     $      PARA_ITBULG*PARA_ITBULG.'
      END IF
C     PARA_ITEXT      = size of extension window in block QZ iteration
      READ ( 1, FMT = '()', IOSTAT = IOS )
      READ ( 1, FMT = *,  IOSTAT = IOS )   PARA_ITEXT
      IF ( PARA_ITEXT.LT.1 ) THEN
         PRINT*, 'Error in QZCONF: PARA_ITEXT must be a positive
     $      integer' 
      END IF
C     PARA_ITTRI      = exploit block triangular structure
      READ ( 1, FMT = '()', IOSTAT = IOS )
      READ ( 1, FMT = *,  IOSTAT = IOS )   PARA_ITTRI
C     PARA_ITDANB     = block parameter in Dackland/Kagstrom QZ
      READ ( 1, FMT = '()', IOSTAT = IOS )
      READ ( 1, FMT = *,  IOSTAT = IOS )   PARA_ITDANB
      IF ( PARA_ITDANB.LE.4 ) THEN
         PRINT*, 'Error in QZCONF: PARA_ITDANB must be larger than 4'
      END IF
C     PARA_ITEXC      - number of QZ iterations before exceptional shift
      READ ( 1, FMT = '()', IOSTAT = IOS )
      READ ( 1, FMT = *,  IOSTAT = IOS )   PARA_ITEXC
      IF ( PARA_ITEXC.LT.1 ) THEN
         PRINT*, 'Error in QZCONF: PARA_ITEXC must be a positive
     $      integer'
      END IF
C     PARA_ITMAX      = maximum number of shifts / N overally applied
      READ ( 1, FMT = '()', IOSTAT = IOS )
      READ ( 1, FMT = *,  IOSTAT = IOS )   PARA_ITMAX
      IF ( PARA_ITMAX.LT.1 ) THEN
         PRINT*, 'Error in QZCONF: PARA_ITMAX must be a positive
     $      integer'
      END IF
C     PARA_ITOPP      = how to compute opposite Householder
      READ ( 1, FMT = '()', IOSTAT = IOS )
      READ ( 1, FMT = *,  IOSTAT = IOS )   PARA_ITOPP
      IF ( PARA_ITOPP.LT.1 .OR. PARA_ITOPP.GT.2 ) THEN
         PRINT*, 'Error in QZCONF: PARA_ITOPP must be in [1,2].'
      END IF
C     
C     Settings for (aggressive early) deflation.
C     
C     PARA_HDEF       = deflation criterion for subdiagonal elements
      READ ( 1, FMT = '()', IOSTAT = IOS )
      READ ( 1, FMT = *,  IOSTAT = IOS )   PARA_HDEF
      IF ( PARA_HDEF.LT.1 .OR. PARA_HDEF.GT.2 ) THEN
         PRINT*, 'Error in QZCONF: PARA_HDEF must be in [1,2].'
      END IF
C     PARA_AGGCRIT    = aggressive early deflation criterion
      READ ( 1, FMT = '()', IOSTAT = IOS )
      READ ( 1, FMT = *,  IOSTAT = IOS )   PARA_AGGCRIT
      IF ( PARA_AGGCRIT.LT.1 .OR. PARA_AGGCRIT.GT.5 ) THEN
         PRINT*, 'Error in QZCONF: PARA_AGGCRIT must be in [1,5].'
      END IF
C     PARA_AGGSHF   = number of shifts after each early deflation
      READ ( 1, FMT = '()', IOSTAT = IOS )
      READ ( 1, FMT = *,  IOSTAT = IOS )   PARA_AGGSHF
      IF ( PARA_AGGSHF.LT.1 ) THEN
         PRINT*, 'Error in QZCONF: PARA_AGGSHF must be a positive
     $      integer.'
      END IF
C     PARA_AGGREP     = ratio for repeating aggressive early deflation
      READ ( 1, FMT = '()', IOSTAT = IOS )
      READ ( 1, FMT = *,  IOSTAT = IOS )   PARA_AGGREP
      IF ( PARA_AGGREP.LT.0.0D0 .OR. PARA_AGGREP.GT.1.0D0 ) THEN
         PRINT*, 'Error in QZCONF: PARA_AGGREP must be in [0.0,1.0].'
      END IF
C     PARA_AGGMIN     = minimal size of active submatrix pair
      READ ( 1, FMT = '()', IOSTAT = IOS )
      READ ( 1, FMT = *,  IOSTAT = IOS )   PARA_AGGMIN
      IF ( PARA_AGGMIN.LT.0 ) THEN
         PRINT*, 'Error in QZCONF: PARA_AGGMIN must be an integer.'
      END IF
C     PARA_AGGWIN     = window size of aggressive early deflation
      READ ( 1, FMT = '()', IOSTAT = IOS )
      READ ( 1, FMT = *,  IOSTAT = IOS )   PARA_AGGWIN
      IF ( PARA_AGGWIN.GT.PARA_AGGMIN ) THEN
         PRINT*, 'Error in QZCONF: PARA_AGGWIN must be an integer not
     $      greater than PARA_AGGMIN.'
      END IF
C     PARA_AGGRQUP    = use RQ update after deflation
      READ ( 1, FMT = '()', IOSTAT = IOS )
      READ ( 1, FMT = *,  IOSTAT = IOS )   PARA_AGGRQUP
C     
C     Settings for deflating infinite eigenvalues.
C     
C     PARA_INFDEF     = deflation criterion for infinite eigenvalues
      READ ( 1, FMT = '()', IOSTAT = IOS )
      READ ( 1, FMT = *,  IOSTAT = IOS )   PARA_INFDEF
      IF ( PARA_INFDEF.LT.0 .OR. PARA_INFDEF.GT.2 ) THEN
         PRINT*, 'Error in QZCONF: PARA_INFDEF is not in [0,2].'
      END IF
C     PARA_INFTOL     = deflation tolerance for infinite eigenvalues
      READ ( 1, FMT = '()', IOSTAT = IOS )
      READ ( 1, FMT = *,  IOSTAT = IOS )   PARA_INFTOL
C     PARA_INFDIR    = direction to chase infinite eigenvalues
      READ ( 1, FMT = '()', IOSTAT = IOS )
      READ ( 1, FMT = *,  IOSTAT = IOS )   PARA_INFDIR
      IF ( PARA_INFDEF.LT.1 .OR. PARA_INFDEF.GT.3 ) THEN
         PRINT*, 'Error in QZCONF: PARA_INFDIR is not in [1,3].'
      END IF
C     PARA_INFALL     = location of infinite eigenvalues to be deflated
      READ ( 1, FMT = '()', IOSTAT = IOS )
      READ ( 1, FMT = *,  IOSTAT = IOS )   PARA_INFALL
C     PARA_INFBLOCK   = use of block algorithms for deflating inf. ev's
      READ ( 1, FMT = '()', IOSTAT = IOS )
      READ ( 1, FMT = *,  IOSTAT = IOS )   PARA_INFBLOCK
C     PARA_INFNZ      = maxium number of zeros to be chased in one block chase
      READ ( 1, FMT = '()', IOSTAT = IOS )
      READ ( 1, FMT = *,  IOSTAT = IOS )   PARA_INFNZ
      IF ( PARA_INFNZ.LT.1 ) THEN
         PRINT*, 'Error in QZCONF: PARA_INFNZ is not a positive',
     $      ' integer.'
      END IF
C     PARA_INFWDW     = window size for one block sweep ( > 2*PARA_INFNZ )
      READ ( 1, FMT = '()', IOSTAT = IOS )
      READ ( 1, FMT = *,  IOSTAT = IOS )   PARA_INFWDW
      IF ( PARA_INFWDW.LE.2*PARA_INFNZ ) THEN
         PRINT*, 'Error in QZCONF: PARA_INFWDW is not greater than',
     $      ' 2*PARA_INFNZ.'
      END IF
C     
C     Other settings.
C     
C     PARA_LDADD      = extra leading dimension in workspace arrays
      READ ( 1, FMT = '()', IOSTAT = IOS )
      READ ( 1, FMT = *,  IOSTAT = IOS )   PARA_LDADD
      IF ( PARA_LDADD.LT.0 ) THEN
         PRINT*, 'Error in QZCONF: PARA_LDADD must be an integer.'
      END IF
C     
      RETURN
C     *** Last line of QZCONF ***
      END
