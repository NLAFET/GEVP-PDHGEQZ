      SUBROUTINE KKQZCONF(eps)
      IMPLICIT NONE
      DOUBLE PRECISION eps
      
C
C     Sets parameters for the KKQZ algorithm instead of 
C     reading them from a file (see QZCONF.f)
C     and stores them in COMMON blocks called "PARA_IT",
C     "PARA_DEF"´ and "PARA_INF".
C
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
C     Settings for the multishift QZ algorithm.
C
C     PARA_ITKK       = minimal size addressed by chains of bulges
      PARA_ITKK = 60
C     PARA_ITDACK     = minimal size addressed by Dackland/Kagstrom QZ
      PARA_ITDACK = 4
      IF ( PARA_ITDACK.GT.PARA_ITKK ) THEN
         PRINT*, 'Error in KKQZCONF: PARA_ITDACK should not be larger
     $ than PARA_ITKK.'
      END IF
C     PARA_ITASH      = usage of shifts from early deflation
      PARA_ITASH = 1
      IF ( PARA_ITASH.LT.0 .OR. PARA_ITASH.GT.2 ) THEN
         PRINT*, 'Error in QZCONF: PARA_ITASH must be in [0,2].'
      END IF
C     PARA_ITBULG     = number of shifts in each bulge
      PARA_ITBULG = 2
      IF ( PARA_ITBULG.LT.2 .OR. MOD( PARA_ITBULG, 2 ).NE.0 ) THEN
         PRINT*, 'Error in QZCONF: PARA_ITBULG must be an even integer >
     $=2.'
      END IF
C     PARA_ITNS       = number of shifts in block QZ iteration 
      PARA_ITNS = 40
      IF ( PARA_ITNS.LT.2 .OR. MOD( PARA_ITNS, PARA_ITBULG ).NE.0 ) THEN
         PRINT*, 'Error in QZCONF: PARA_ITNS must be an integer multiple
     $ of PARA_ITBULG.'
      END IF
      IF ( PARA_ITNS.LT.PARA_ITBULG*PARA_ITBULG ) THEN
           PRINT*, 'Error in QZCONF: PARA_ITNS must not be smaller than
     $ PARA_ITBULG*PARA_ITBULG.'
      END IF
C     PARA_ITEXT      = size of extension window in block QZ iteration
      PARA_ITEXT = 50
      IF ( PARA_ITEXT.LT.1 ) THEN
         PRINT*, 'Error in QZCONF: PARA_ITEXT must be a positive
     $ integer' 
      END IF
C     PARA_ITTRI      = exploit block triangular structure
      PARA_ITTRI = .TRUE.
C     PARA_ITDANB     = block parameter in Dackland/Kagstrom QZ 
      PARA_ITDANB = 32
      IF ( PARA_ITDANB.LE.4 ) THEN
         PRINT*, 'Error in QZCONF: PARA_ITDANB must be larger than 4'
      END IF
C     PARA_ITEXC      - number of QZ iterations before exceptional shift
      PARA_ITEXC = 10
      IF ( PARA_ITEXC.LT.1 ) THEN
         PRINT*, 'Error in QZCONF: PARA_ITEXC must be a positive
     $ integer'
      END IF
C     PARA_ITMAX      = maximum number of shifts / N overally applied
      PARA_ITMAX = 30
      IF ( PARA_ITMAX.LT.1 ) THEN
         PRINT*, 'Error in QZCONF: PARA_ITMAX must be a positive
     $ integer'
      END IF
C     PARA_ITOPP      = how to compute opposite Householder
      PARA_ITOPP = 2
      IF ( PARA_ITOPP.LT.1 .OR. PARA_ITOPP.GT.2 ) THEN
         PRINT*, 'Error in QZCONF: PARA_ITOPP must be in [1,2].'
      END IF
C
C     Settings for (aggressive early) deflation.
C
C     PARA_HDEF       = deflation criterion for subdiagonal elements
      PARA_HDEF = 2
      IF ( PARA_HDEF.LT.1 .OR. PARA_HDEF.GT.2 ) THEN
         PRINT*, 'Error in QZCONF: PARA_HDEF must be in [1,2].'
      END IF
C     PARA_AGGCRIT    = aggressive early deflation criterion
      PARA_AGGCRIT = 5
      IF ( PARA_AGGCRIT.LT.1 .OR. PARA_AGGCRIT.GT.5 ) THEN
         PRINT*, 'Error in QZCONF: PARA_AGGCRIT must be in [1,5].'
      END IF
C     PARA_AGGSHF   = number of shifts after each early deflation
      PARA_AGGSHF = 16
      IF ( PARA_AGGSHF.LT.1 ) THEN
         PRINT*, 'Error in QZCONF: PARA_AGGSHF must be a positive
     $ integer.'
      END IF
C     PARA_AGGREP     = ratio for repeating aggressive early deflation
      PARA_AGGREP = 0.2D0
      IF ( PARA_AGGREP.LT.0.0D0 .OR. PARA_AGGREP.GT.1.0D0 ) THEN
         PRINT*, 'Error in QZCONF: PARA_AGGREP must be in [0.0,1.0].'
      END IF
C     PARA_AGGMIN     = minimal size of active submatrix pair
      PARA_AGGMIN = 60
      IF ( PARA_AGGMIN.LT.0 ) THEN
         PRINT*, 'Error in QZCONF: PARA_AGGMIN must be an integer.'
      END IF
C     PARA_AGGWIN     = window size of aggressive early deflation
      PARA_AGGWIN = 60
      IF ( PARA_AGGWIN.GT.PARA_AGGMIN ) THEN
         PRINT*, 'Error in QZCONF: PARA_AGGWIN must be an integer not
     $ greater than PARA_AGGMIN.'
      END IF
C     PARA_AGGRQUP    = use RQ update after deflation
      PARA_AGGRQUP = .TRUE.
C
C     Settings for deflating infinite eigenvalues.
C
C     PARA_INFDEF     = deflation criterion for infinite eigenvalues
      PARA_INFDEF = 1
      IF ( PARA_INFDEF.LT.0 .OR. PARA_INFDEF.GT.2 ) THEN
         PRINT*, 'Error in QZCONF: PARA_INFDEF is not in [0,2].'
      END IF
C     PARA_INFTOL     = deflation tolerance for infinite eigenvalues
      PARA_INFTOL = eps
C     PARA_INFDIR    = direction to chase infinite eigenvalues
      PARA_INFDIR = 3
      IF ( PARA_INFDEF.LT.1 .OR. PARA_INFDEF.GT.3 ) THEN
         PRINT*, 'Error in QZCONF: PARA_INFDIR is not in [1,3].'
      END IF
C     PARA_INFALL     = location of infinite eigenvalues to be deflated
      PARA_INFALL = .TRUE.
C     PARA_INFBLOCK   = use of block algorithms for deflating inf. ev''s
      PARA_INFBLOCK = .TRUE.
C     PARA_INFNZ      = maxium number of zeros to be chased in one block chasea
      PARA_INFNZ = 3
      IF ( PARA_INFNZ.LT.1 ) THEN
         PRINT*, 'Error in QZCONF: PARA_INFNZ is not a positive',
     $           ' integer.'
      END IF
C     PARA_INFWDW     = window size for one block sweep ( > 2*PARA_INFNZ )
      PARA_INFWDW= 13
      IF ( PARA_INFWDW.LE.2*PARA_INFNZ ) THEN
         PRINT*, 'Error in QZCONF: PARA_INFWDW is not greater than',
     $           ' 2*PARA_INFNZ.'
      END IF
C
C     Other settings.
C
C     PARA_LDADD      = extra leading dimension in workspace arrays
      PARA_LDADD = 0
      IF ( PARA_LDADD.LT.0 ) THEN
         PRINT*, 'Error in QZCONF: PARA_LDADD must be an integer.'
      END IF
C
      RETURN
C *** Last line of KKQZCONF ***
      END
