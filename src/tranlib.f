C%Z%===================================================================
C%Z%
C%Z%
C%Z%                   FILE =  %M%
C%Z%
C%Z%  ---------------  VERSION = %I%
C%Z%  |  SCCS  FILE |
C%Z%  |   SUMMARY   |  CURRENT CHECKOUT DATE = %H%
C%Z%  ---------------                           at %T%
C%Z%                   DATE OF NEWEST DELTA = %G%
C%Z%                                            at %U%
C%Z%  SCCS file name = %F%
C%Z%===================================================================
C23456789012345678901234567890123456789012345678901234567890123456789012
C
C     Version 3.5
C
         SUBROUTINE MCABS
C
C        CHANGES FROM VERSION 1.0:
C        1. Changed REAL*8 to DOUBLE PRECISION
C        CHANGES FROM VERSION 1.1:
C        1. Eliminated many GOTO's
C        CHANGES FOR VERSION 1.3
C        1. SUBROUTINE MCLEN
C        CHANGES FOR VERSION 1.4
C        1. Additional record to binary file indicates 
C           version, machine precision, and error status
C        2. Additional record to binary file has required lengths for 
C           integer, real work arrays
C        3. New Subroutines MCPNT, MCSAVE read, write binary
C           file information, work arrays
C        CHANGES FOR VERSION 1.6
C        1. Versions 1.8 and 1.9 added (TRANLIB V.1.5
C           and TRANFIT V.1.8 were intermediate versions which may
C           not be legitimate; TRANFIT V.1.9 is actually a
C           correction to V.1.7, and TRANLIB 1.6 is an update of
C           V.1.4)
C
C        CHANGES FOR VERSION 1.7 (10/1/92 F. Rupley per M. Coltrin)
C        1. Created MCABS to hold version and change information
C        2. COMMON /MCCONS/ VERS, PREC, KERR, LENI, LENR eliminates
C           the need for argument LINKMC in the MCSAVE call list
C     CHANGES FOR VERSION 3.0 (3/15/94 F. Rupley)
C     1.  DOS/PC compatibility effort includes adding file names to
C         OPEN statements, removing unused variables in CALL lists,
C         unusued but possibly initialized variables.
C     CHANGES FOR VERSION 3.1 (4/5/94 M. Coltrin)
C     1.  Fix error in MCLMDT around line 1674 ("ZERO" replaced by
C         "ONE").
C     2.  Indices in formula for BINDIF(I,I) in loop 150 should have
C         been (K,K).
C     3.  A factor of PFAC was omitted from the statement before
C         statement number 1600 in MCLMDT.
C     4.  The factor PIFAC and the statement in which it was used in
C         loop 1600 in MCLMDT were both incorrect.
C     CHANGES FOR VERSION 3.2 (6/8/94 H.K. Moffat
C     1.  Made MCLMDT roughly 30 times faster on a sample problem.
C         Opcount in the function now scales like KK**2, instead
C         of the previous KK**3.
C     2.  No longer need a delta function in the library -
C         CHANGES to MCLMDT and MCORDF
C     CHANGES FOR VERSION 3.3 (6/9/94) H.K. Moffat
C     1.  Changed the Gas constant to agree with 1986 CODATA value.
C         (two extra significant digits of accuracy)
C     2.  Changed parameter values for pi and its powers to
C         so that the values are accurate to machine precision.
C     CHANGES FOR VERSION 3.4 (6/12/94) H.K. Moffat
C     1.  Added change blocks for LAPACK linear algebra
C     CHANGES FOR VERSION 3.5 (8/10/94 H.K. Moffat
C     1.  Accepts to 3.2 version number
C
C
      RETURN
      END
C---------------------------------------------------------------------
C
      SUBROUTINE MCINIT (LINKMC, LOUT, LENIMC, LENRMC, IMCWRK, RMCWRK)
C
C*****precision > double
      IMPLICIT DOUBLE PRECISION (A-H, O-Z), INTEGER (I-N)
C*****END precision > double
C
C*****precision > single
C      IMPLICIT REAL (A-H, O-Z), INTEGER (I-N)
C*****END precision > single
C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C
C  SUBROUTINE MCINIT (LINKMC, LOUT, LENIMC, LENRMC, IMCWRK, RMCWRK)
C
C    THIS SUBROUTINE SERVES TO READ THE LINKING FILE FROM THE FITTING
C    CODE AND TO CREATE THE INTERNAL STORAGE AND WORK ARRAYS, IMCWRK(*)
C    AND RMCWRK(*).  MCINIT MUST BE CALLED BEFORE ANY OTHER TRANSPORT
C    SUBROUTINE IS CALLED.  IT MUST BE CALLED AFTER THE CHEMKIN PACKAGE
C    IS INITIALIZED.
C
C  INPUT-
C    LINKMC  - LOGICAL UNIT NUMBER OF THE LINKING FILE.
C                  FITTING CODE WRITE TO DEFAULT UNIT 35
C    LOUT    - LOGICAL UNIT NUMBER FOR PRINTED OUTPUT.
C    LENIMC  - ACTUAL DIMENSION OF THE INTEGER STORAGE AND WORKING
C              SPACE, ARRAY IMCWRK(*).  LENIMC MUST BE AT LEAST:
C                LENIMC = 4*KK + NLITE
C                 WHERE, KK    = NUMBER OF SPECIES.
C                        NLITE = NUMBER OF SPECIES WITH MOLECULAR WEIGHT
C                                LESS THAN 5.
C    LENRMC  - ACTUAL DIMENSION OF THE FLOATING POINT STORAGE AND
C              WORKING SPACE, ARRAY RMCWRK(*).  LENRMC MUST BE AT LEAST:
C                LENRMC = KK*(19 + 2*NO + NO*NLITE) + (NO+15)*KK**2
C                 WHERE, KK    = NUMBER OF SPECIES.
C                        NO    = ORDER OF THE POLYNOMIAL FITS,
C                                DEFAULT, NO=4.
C                        NLITE = NUMBER OF SPECIES WITH MOLECULAR WEIGHT
C                                LESS THAN 5.
C
C  WORK-
C    IMCWRK  - ARRAY OF INTEGER STORAGE AND WORK SPACE.  THE STARTING
C              ADDRESSES FOR THE IMCWRK SPACE ARE STORED IN
C              COMMON /MCMCMC/.
C                  DIMENSION IMCWRK(*) AT LEAST LENIMC.
C    RMCWRK  - ARRAY OF FLOATING POINT STORAGE AND WORK SPACE.  THE
C              STARTING ADDRESSES FOR THE RMCWRK SPACE ARE STORED IN
C              COMMON /MCMCMC/.
C                  DIMENSION RMCWRK(*) AT LEAST LENRMC.
C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C
      DIMENSION IMCWRK(*), RMCWRK(*)
      CHARACTER*16 VERS, PREC
      LOGICAL IOK, ROK, KERR
      COMMON /MCCONS/ VERS, PREC, KERR, LENI, LENR
C
      COMMON /MCMCMC/ RU, PATMOS, SMALL, NKK, NO, NLITE, INLIN, IKTDIF,
     1                IPVT, NWT, NEPS, NSIG, NDIP, NPOL, NZROT, NLAM,
     2                NETA, NDIF, NTDIF, NXX, NVIS, NXI, NCP,
     3                NCROT, NCINT, NBIND, NEOK, NSGM,
     4                NAST, NBST, NCST, NXL, NR, NWRK, K3
C
C
C         THE FOLLOWING NUMBER SMALL IS USED IN THE MIXTURE DIFFUSION
C        COEFFICIENT CALCULATION.  ITS USE ALLOWS A SMOOTH AND WELL
C        DEFINED DIFFUSION COEFFICIENT AS THE MIXTURE APPROACHES A
C        PURE SPECIES, EVEN THOUGH STRICTLY SPEAKING THERE DOES NOT
C        EXIST A DIFFUSION COEFFICIENT IN THIS CASE.  THE VALUE OF
C        "SMALL" SHOULD BE SMALL RELATIVE TO ANY SPECIES MOLE FRACTION
C        OF IMPORTANCE, BUT LARGE ENOUGH TO BE REPRESENTED ON THE
C        COMPUTER.
C
      SMALL = 1.0E-20
C
C          Gas constant as reported in 1993 CRC, (J. Research of 
C        National Bureau of Standards, 92, 85, 1987).
C         ( 8.314510(70)E+07 Joules mol-1 K-1)
C
      RU    = 8.314510E+07
C
C          Standard atmosphere (defined as an exact quantity)
C
      PATMOS= 1.01325E+06
C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C
C          WRITE VERSION NUMBER
C
!       WRITE (LOUT, 15)
!    15 FORMAT(
!      1/' TRANLIB:  Multicomponent transport library,',
!      2/'           CHEMKIN-II Version 3.5, August 1994',
! C*****precision > double
!      3/'           DOUBLE PRECISION')
C*****END precision > double
C*****precision > single
C     3/'           SINGLE PRECISION')
C*****END precision > single
C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C
C        READ THE PROBLEM SIZE
C
      CALL MCLEN (LINKMC, LOUT, LI, LR)
      IOK = (LENIMC .GE. LI)
      ROK = (LENRMC .GE. LR)
C
      IF (.NOT.IOK .OR. .NOT.ROK) THEN
         IF (.NOT. IOK) WRITE (LOUT, 300) LI
         IF (.NOT. ROK) WRITE (LOUT, 350) LR
         STOP
      ENDIF
C
      REWIND LINKMC
      READ (LINKMC, ERR=999) VERS, PREC, KERR
      READ (LINKMC, ERR=999) LI, LR, NO, NKK, NLITE
C
      NK  = NO*NKK
      NK2 = NO*NKK*NKK
      K2  = NKK*NKK
      K3  = 3*NKK
      K32 = K3*K3
      NKT = NO*NKK*NLITE
C
C        APPORTION THE REAL WORK SPACE
C           THE POINTERS HAVE THE FOLLOWING MEANINGS:
C             NWT   - THE SPECIES MOLECULAR WEIGHTS.
C             NEPS  - THE EPSILON/K WELL DEPTH FOR THE SPECIES.
C             NSIG  - THE COLLISION DIAMETER FOR THE SPECIES.
C             NDIP  - THE DIPOLE MOMENTS FOR THE SPECIES.
C             NPOL  - THE POLARIZABILITIES FOR THE SPECIES.
C             NZROT - THE ROTATIONAL RELAXATION COLLISION NUMBERS.
C             NLAM  - THE COEFFICIENTS FOR THE CONDUCTIVITY FITS.
C             NETA  - THE COEFFICIENTS FOR THE VISCOSITY FITS.
C             NTDIF - THE COEFFICIENTS FOR THE THERMAL DIFFUSION
C                     RATIO FITS.
C             NXX   - THE MOLE FRACTIONS.
C             NVIS  - THE SPECIES VISCOSITIES.
C             NXI   - THE ROTATIONAL RELAXATION COLLISION NUMBERS BEFORE
C                     THE PARKER COFFECTION.
C             NCP   - THE SPECIES SPECIFIC HEATS.
C             NCROT - THE ROTATIONAL PARTS OF THE SPECIFIC HEATS.
C             NCINT - THE INTERNAL PARTS OF THE SPECIFIC HEATS.
C             NBIND - THE BINARY DIFFUSION COEFFICIENTS.
C             NEOK  - THE MATRIX OF REDUCED WELL DEPTHS.
C             NSGM  - THE MATRIX OF REDUCED COLLISION DIAMETERS.
C             NAST  - THE MATRIX OF A* COLLISION INTEGRALS FOR EACH
C                     SPECIES PAIR.
C             NBST  - THE MATRIX OF B* COLLISION INTEGRALS FOR EACH
C                     SPECIES PAIR.
C             NCST  - THE MATRIX OF C* COLLISION INTEGRALS FOR EACH
C                     SPECIES PAIR.
C             NXL   - THE "L" MATRIX.
C             NR    - THE RIGHT HAND SIDES OF THE LINEAR SYSTEM
C                     INVOLVING THE "L" MATRIX.
C             NWRK  - THE WORK SPACE NEEDED BY LINPACK TO SOLVE THE
C                     "L" MATRIX LINEAR SYSTEM.
C
      NWT  = 1
      NEPS = NWT + NKK
      NSIG = NEPS + NKK
      NDIP = NSIG + NKK
      NPOL = NDIP + NKK
      NZROT= NPOL + NKK
C
      NLAM = NZROT + NKK
      NETA = NLAM + NK
      NDIF = NETA + NK
      NTDIF= NDIF + NK2
C
      NXX  = NTDIF + NO*NKK*NLITE
      NVIS = NXX + NKK
      NXI  = NVIS + NKK
      NCP  = NXI + NKK
      NCROT= NCP + NKK
      NCINT= NCROT + NKK
C
      NBIND= NCINT + NKK
      NEOK = NBIND + K2
      NSGM = NEOK + K2
      NAST = NSGM + K2
      NBST = NAST + K2
      NCST = NBST + K2
C
      NXL  = NCST + K2
C
      NR   = NXL + K32
      NWRK = NR + K3
C      NTOT = NWRK + K3 - 1
C
C           APPORTION THE INTEGER WORK SPACE
C              THE POINTERS HAVE THE FOLLOWING MEANING:
C
C                INLIN - THE INDICATORS FOR THE MOLECULE LINEARITY.
C                IKTDIF- THE SPECIES INDICIES FOR THE "LIGHT" SPECIES.
C                IPVT  - THE PIVOT INDICIES FOR LINPACK CALLS.
C
      INLIN = 1
      IKTDIF= INLIN + NKK
      IPVT  = IKTDIF + NLITE
C      ITOT  = IPVT + K3 - 1
C
C        READ THE DATA FROM THE LINK FILE
C

      READ (LINKMC, ERR=999) PATMOS, (RMCWRK(NWT+N-1),
     1   RMCWRK(NEPS+N-1), RMCWRK(NSIG+N-1),
     2   RMCWRK(NDIP+N-1), RMCWRK(NPOL+N-1), RMCWRK(NZROT+N-1),
     3   IMCWRK(INLIN+N-1), N=1,NKK),
     4   (RMCWRK(NLAM+N-1), N=1,NK), (RMCWRK(NETA+N-1), N=1,NK),
     5   (RMCWRK(NDIF+N-1), N=1,NK2),
     6   (IMCWRK(IKTDIF+N-1), N=1,NLITE), (RMCWRK(NTDIF+N-1), N=1,NKT)
C
C        SET EPS/K AND SIG FOR ALL I,J PAIRS
C
      CALL MCEPSG (NKK, RMCWRK(NEPS), RMCWRK(NSIG), RMCWRK(NDIP),
     1            RMCWRK(NPOL), RMCWRK(NEOK), RMCWRK(NSGM) )
C
  300 FORMAT (10X,'IMCWRK MUST BE DIMENSIONED AT LEAST ', I5)
  350 FORMAT (10X,'RMCWRK MUST BE DIMENSIONED AT LEAST ', I5)
      RETURN
  999 CONTINUE
      WRITE (LOUT, *) ' Error reading Transport binary file...'
      STOP
      END
C
C---------------------------------------------------------------------
C
      SUBROUTINE MCLEN (LINKMC, LOUT, LI, LR)
C
C*****precision > double
      IMPLICIT DOUBLE PRECISION (A-H, O-Z), INTEGER (I-N)
C*****END precision > double
C
C*****precision > single
C      IMPLICIT REAL (A-H, O-Z), INTEGER (I-N)
C*****END precision > single
C
      PARAMETER (NLIST = 3)
      LOGICAL KERR, VOK, POK
      CHARACTER*16 LIST(NLIST), VERS, PREC
      COMMON /MCCONS/ VERS, PREC, KERR, LENI, LENR
      DATA LIST/'3.0', '3.1', '3.2'/
C
      VERS = ' '
      PREC = ' '
      LENI = 0
      LENR = 0
      LI   = LENI
      LR   = LENR
C
      REWIND (LINKMC)
      READ (LINKMC, ERR=999) VERS, PREC, KERR
C
      VOK = .FALSE.
      DO 5 N = 1, NLIST
         IF (VERS .EQ. LIST(N)) VOK = .TRUE.
    5 CONTINUE
C
      POK = .FALSE.
C*****precision > double
      IF (INDEX(PREC, 'DOUB') .GT. 0) POK = .TRUE.
C*****END precision > double
C*****precision > single
C      IF (INDEX(PREC, 'SING') .GT. 0) POK = .TRUE.
C*****END precision > single
C
      IF (KERR .OR. (.NOT.POK) .OR. (.NOT.VOK)) THEN
         IF (KERR) THEN
            WRITE (LOUT,'(/A,/A)')
     1      ' There is an error in the transport binary file...',
     2      ' Check TRANFIT output for error conditions.'
         ENDIF
         IF (.NOT. VOK) THEN
            WRITE (LOUT,'(/A,A)')
     1      ' Transport binary file is incompatible with Transport',
     2      ' Library Version 1.7'
         ENDIF
         IF (.NOT. POK) THEN
            WRITE (LOUT, '(/A,A)')
     1      ' Precision of Transport binary file does not agree with',
     2      ' precision of Transport Library'
         ENDIF
         STOP
      ENDIF
C
      READ (LINKMC, ERR=999) LENIMC, LENRMC, NO, NKK, NLITE
      REWIND (LINKMC)
      LENI = LENIMC
      LENR = LENRMC
      LI   = LENI
      LR   = LENR
      RETURN
C
  999 CONTINUE
      WRITE (LOUT, 50)
   50 FORMAT
     1 (' Error reading Multi-component Transport binary file.')
      STOP
      END
C
C----------------------------------------------------------------------C
C
      SUBROUTINE MCPNT (LSAVE, LOUT, NPOINT, V, P, LI, LR, IERR)
C
C  START PROLOGUE
C
C  SUBROUTINE MCPNT (LSAVE, LOUT, NPOINT, VERS, PREC, LENI, LENR,
C                    KERR)
C     Reads from a binary file information about a Transport
C     binary file, pointers for the Transport Library, and
C     returns lengths of work arrays.
C
C  INPUT
C     LSAVE  - Integer input unit for binary data file.
C                   Data type - integer scalar
C     LOUT   - Integer output unit for printed diagnostics.
C                   Data type - integer scalar
C
C  OUTPUT
C     NPOINT - Total number of pointers.
C                   Data type - integer scalar
C     VERS   - Version number of the Transport binary file.
C                   Data type - real scalar
C     PREC   - Machine precision of the Transport binary file.
C                   Data type - character string
C     LENI   - Minimum length required for the integer work array.
C                   Data type - integer scalar
C     LENR   - Minimum length required for the real work array.
C                   Data type - integer scalar
C     KERR   - Logical error flag.
C
C  END PROLOGUE
C
C*****precision > double
        IMPLICIT DOUBLE PRECISION (A-H, O-Z), INTEGER (I-N)
C*****END precision > double
C*****precision > single
C        IMPLICIT REAL (A-H, O-Z), INTEGER (I-N)
C*****END precision > single
C
      COMMON /MCMCMC/ RU, PATMOS, SMALL, NKK, NO, NLITE, INLIN, IKTDIF,
     1                IPVT, NWT,  NEPS,  NSIG, NDIP, NPOL, NZROT, NLAM,
     2                NETA, NDIF,  NTDIF, NXX, NVIS, NXI,  NCP,  NCROT,
     3                NCINT, NBIND, NEOK, NSGM, NAST, NBST, NCST,
     4                NXL, NR, NWRK, K3
      LOGICAL KERR, IERR
      CHARACTER PREC*16, VERS*16, P*16, V*16
      COMMON /MCCONS/ VERS, PREC, KERR, LENI, LENR
C
      KERR = .FALSE.
      IERR = .FALSE.
C
      READ (LSAVE, ERR=100) NPOINT, VERS, PREC, LENI, LENR,
     *                RU, PATMOS, SMALL, NKK, NO, NLITE, INLIN, IKTDIF,
     1                IPVT, NWT, NEPS, NSIG, NDIP, NPOL, NZROT, NLAM,
     2                NETA, NDIF, NTDIF, NXX, NVIS, NXI, NCP,
     3                NCROT, NCINT, NBIND, NEOK, NSGM,
     4                NAST, NBST, NCST, NXL, NR, NWRK, K3
      LI = LENI
      LR = LENR
      IERR = KERR
      V = VERS
      P = PREC
      RETURN
C
  100 CONTINUE
      WRITE (LOUT, *) ' Error reading Transport binary file data...'
      KERR   = .TRUE.
      IERR   = KERR
      NPOINT = 0
      VERS   = ' '
      V      = VERS
      PREC   = ' '
      P      = PREC
      RETURN
      END
C
      SUBROUTINE MCPRAM (IMCWRK, RMCWRK, EPS, SIG, DIP, POL, ZROT,
     1                   NLIN)
C
C*****precision > double
      IMPLICIT DOUBLE PRECISION (A-H, O-Z), INTEGER (I-N)
C*****END precision > double
C
C*****precision > single
C      IMPLICIT REAL (A-H, O-Z), INTEGER (I-N)
C*****END precision > single
C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C
C  SUBROUTINE MCPRAM (IMCWRK, RMCWRK, EPS, SIG, DIP, POL, ZROT, NLIN)
C
C    THIS SUBROUTINE IS CALLED TO RETURN THE ARRAYS OF MOLECULAR
C    PARAMETERS AS READ FROM THE TRANSPORT DATA BASE.
C
C  WORK-
C    IMCWRK  - ARRAY OF INTEGER STORAGE AND WORK SPACE.  THE STARTING
C              ADDRESSES FOR THE IMCWRK SPACE ARE STORED IN
C              COMMON /MCMCMC/.
C                  DIMENSION IMCWRK(*) AT LEAST LENIMC.
C    RMCWRK   - ARRAY OF FLOATING POINT STORAGE AND WORK SPACE.  THE
C              STARTING ADDRESSES FOR THE RMCWRK SPACE ARE STORED IN
C              COMMON /MCMCMC/.
C                  DIMENSION RMCWRK(*) AT LEAST LENRMC.
C
C  OUTPUT-
C    EPS     - ARRAY OF LENNARD-JONES POTENTIAL WELL DEPTHS.
C                  CGS UNITS - K.
C                  DIMENSION EPS(*) AT LEAST KK
C    SIG     - ARRAY OF LENNARD-JONES COLLISION DIAMETERS.
C                  UNITS - ANGSTROMS.
C                  DIMENSION SIG(*) AT LEAST KK
C    DIP     - ARRAY OF DIPOLE MOMENTS
C                  UNITS - DEBYE
C                  DIMENSION DIP(*) AT LEAST KK
C    POL     - ARRAY OF POLARIZABILITIES.
C                  UNITS - ANGSTROMS**3.
C                  DIMENSION POL(*) AT LEAST KK
C    ZROT    - ARRAY OF ROTATIONAL COLLISION NUMBERS EVALUATED
C              AT 298K.
C                  UNITS - NONE
C                  DIMENSION ZROT(*) AT LEAST KK
C    NLIN    - ARRAY OF FLAGS INDICATING WHETHER THE MOLECULE
C              LINEAR OR NOT.
C              NLIN=0, SINGLE ATOM.
C              NLIN=1, LINEAR MOLECULE.
C              NLIN=2, NONLINEAR MOLECULE.
C                  UNITS - NONE.
C                  DIMENSION NLIN(*) AT LEAST KK
C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C
      DIMENSION IMCWRK(*), RMCWRK(*), EPS(*), SIG(*), DIP(*), POL(*),
     1          ZROT(*), NLIN(*)
C
      COMMON /MCMCMC/ RU, PATMOS, SMALL, NKK, NO, NLITE, INLIN, IKTDIF,
     1                IPVT, NWT, NEPS, NSIG, NDIP, NPOL, NZROT, NLAM,
     2                NETA, NDIF, NTDIF, NXX, NVIS, NXI, NCP, NCROT,
     3                NCINT, NBIND, NEOK, NSGM, NAST, NBST,
     4                NCST, NXL, NR, NWRK, K3
C
      DO 100 K = 1, NKK
         EPS(K) = RMCWRK(NEPS+K-1)
         SIG(K) = RMCWRK(NSIG+K-1)
         DIP(K) = RMCWRK(NDIP+K-1)
         POL(K) = RMCWRK(NPOL+K-1)
         ZROT(K)= RMCWRK(NZROT+K-1)
         NLIN(K)= IMCWRK(INLIN+K-1)
  100 CONTINUE
C
      RETURN
      END
C
C----------------------------------------------------------------------C
C
      SUBROUTINE MCSAVE (LOUT, LSAVE, IMCWRK, RMCWRK)
C
C  START PROLOGUE
C
C  SUBROUTINE MCSAVE (LOUT, LSAVE, IMCWRK, RMCWRK)
C     Writes to a binary file information about a Trasport
C     binary file, pointers for the Transport Library, and
C     Transport work arrays.
C
C  INPUT
C     LOUT   - Output file for printed diagnostics.
C                   Data type - integer scalar
C     LSAVE  - Integer output unit.
C                   Data type - integer scalar
C     IMCWRK - Array of integer workspace containing integer data.
C                   Data type - integer array
C     RMCWRK - Array of real workspace containing real data.
C                   Data type - real array
C
C  END PROLOGUE
C
C*****precision > double
        IMPLICIT DOUBLE PRECISION (A-H, O-Z), INTEGER (I-N)
C*****END precision > double
C*****precision > single
C        IMPLICIT REAL (A-H, O-Z), INTEGER (I-N)
C*****END precision > single
C
      DIMENSION IMCWRK(*), RMCWRK(*)
      COMMON /MCMCMC/ RU, PATMOS, SMALL, NKK, NO, NLITE, INLIN, IKTDIF,
     1                IPVT, NWT, NEPS, NSIG, NDIP, NPOL, NZROT, NLAM,
     2                NETA, NDIF, NTDIF, NXX, NVIS, NXI, NCP,
     3                NCROT, NCINT, NBIND, NEOK, NSGM,
     4                NAST, NBST, NCST, NXL, NR, NWRK, K3
      CHARACTER VERS*16, PREC*16
      LOGICAL KERR
      COMMON /MCCONS/ VERS, PREC, KERR, LENI, LENR
C
      NPOINT = 41
      WRITE (LSAVE,ERR=999)   NPOINT, VERS,   PREC,   LENI,   LENR,
     *                RU, PATMOS, SMALL, NKK, NO, NLITE, INLIN, IKTDIF,
     1                IPVT, NWT, NEPS, NSIG, NDIP, NPOL, NZROT, NLAM,
     2                NETA, NDIF, NTDIF, NXX, NVIS, NXI, NCP,
     3                NCROT, NCINT, NBIND, NEOK, NSGM,
     4                NAST, NBST, NCST, NXL, NR, NWRK, K3
      WRITE (LSAVE,ERR=999) (IMCWRK(L), L = 1, LENI)
      WRITE (LSAVE,ERR=999) (RMCWRK(L), L = 1, LENR)
C
      RETURN
  999 CONTINUE
      WRITE (LOUT, *)
     1    ' Error writing Transport binary file information...'
      RETURN
      END
C
C
C
C
      SUBROUTINE MCSVIS (T, RMCWRK, VIS)
C
C*****precision > double
      IMPLICIT DOUBLE PRECISION (A-H, O-Z), INTEGER (I-N)
C*****END precision > double
C
C*****precision > single
C      IMPLICIT REAL (A-H, O-Z), INTEGER (I-N)
C*****END precision > single
C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C
C  SUBROUTINE MCSVIS (T, RMCWRK, VIS)
C
C    THIS SUBROUTINE COMPUTES THE ARRAY OF PURE SPECIES VISCOSITIES,
C    GIVEN THE TEMPERATURE.
C
C  INPUT-
C    T       - TEMPERATURE
C                  CGS UNITS - K.
C
C  WORK-
C    RMCWRK   - ARRAY OF FLOATING POINT STORAGE AND WORK SPACE.  THE
C              STARTING ADDRESSES FOR THE RMCWRK SPACE ARE STORED IN
C              COMMON /MCMCMC/.
C                  DIMENSION RMCWRK(*) AT LEAST LENRMC.
C
C  OUTPUT-
C    VIS     - ARRAY OF SPECIES VISCOSITIES.
C                  CGS UNITS - GM/CM*S.
C                  DIMENSION VIS(*) AT LEAST KK.
C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C
      DIMENSION VIS(*), RMCWRK(*)
C
      COMMON /MCMCMC/ RU, PATMOS, SMALL, NKK, NO, NLITE, INLIN, IKTDIF,
     1                IPVT, NWT, NEPS, NSIG, NDIP, NPOL, NZROT, NLAM,
     2                NETA, NDIF, NTDIF, NXX, NVIS, NXI, NCP, NCROT,
     3                NCINT, NBIND, NEOK, NSGM, NAST, NBST,
     4                NCST, NXL, NR, NWRK, K3
C
      ALOGT = LOG(T)
      CALL MCEVAL (ALOGT, NKK, NO, RMCWRK(NETA), VIS)
      DO 25 K = 1, NKK
         VIS(K) = EXP(VIS(K))
   25 CONTINUE
C
      RETURN
      END
C
C---------------------------------------------------------------------
C
      SUBROUTINE MCAVIS (T, X, RMCWRK, VISMIX)
C
C*****precision > double
      IMPLICIT DOUBLE PRECISION (A-H, O-Z), INTEGER (I-N)
C*****END precision > double
C
C*****precision > single
C      IMPLICIT REAL (A-H, O-Z), INTEGER (I-N)
C*****END precision > single
C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C
C  SUBROUTINE MCAVIS (T, X, RMCWRK, VISMIX)
C
C    THIS SUBROUTINE COMPUTES THE MIXTURE VISCOSITY, GIVEN
C    THE TEMPERATURE AND THE SPECIES MOLE FRACTIONS.  IT USES MODIFICATI
C    OF THE WILKE SEMI-EMPIRICAL FORMULAS.
C
C  INPUT-
C    T       - TEMPERATURE
C                  CGS UNITS - K.
C    X       - ARRAY OF MOLE FRACTIONS OF THE MIXTURE.
C                  DIMENSION Y(*) AT LEAST KK.
C
C  WORK-
C    RMCWRK  - ARRAY OF FLOATING POINT STORAGE AND WORK SPACE.  THE
C              STARTING ADDRESSES FOR THE RMCWRK SPACE ARE STORED IN
C              COMMON /MCMCMC/.
C                  DIMENSION RMCWRK(*) AT LEAST LENRMC.
C
C  OUTPUT-
C    VISMIX  - MIXTURE VISCOSITY
C                  CGS UNITS - GM/CM*S.
C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C
      DIMENSION X(*), RMCWRK(*)
C
      COMMON /MCMCMC/ RU, PATMOS, SMALL, NKK, NO, NLITE, INLIN, IKTDIF,
     1                IPVT, NWT, NEPS, NSIG, NDIP, NPOL, NZROT, NLAM,
     2                NETA, NDIF, NTDIF, NXX, NVIS, NXI, NCP, NCROT,
     3                NCINT, NBIND, NEOK, NSGM, NAST, NBST,
     4                NCST, NXL, NR, NWRK, K3
C
      SQRT8R = 0.3535534
      ZERO   = 0.0
      ONE    = 1.0
C
C       IN THE FOLLOWING CALL:
C         THE SPECIES MOLECULAR WEIGHTS ARE STORED IN RMCWRK(NWT)
C         THE PURE SPECIES VISCOSITIES ARE IN RMCWRK(NVIS)
C
      ALOGT = LOG(T)
      CALL MCEVAL (ALOGT, NKK, NO, RMCWRK(NETA), RMCWRK(NVIS))
      DO 25 K = 1, NKK
         RMCWRK(NVIS+K-1) = EXP(RMCWRK(NVIS+K-1))
   25 CONTINUE
C
      SUMO = ZERO
      DO 200 K = 1, NKK
C
         SUMI = ZERO
         DO 100 J = 1, NKK
            TOP = (ONE + SQRT(RMCWRK(NVIS+K-1)/RMCWRK(NVIS+J-1)) *
     1               (RMCWRK(NWT+J-1)/RMCWRK(NWT+K-1))**0.25E0 )**2
            BOT = SQRT(ONE + RMCWRK(NWT+K-1) / RMCWRK(NWT+J-1))
            PHIKJ = SQRT8R*TOP/BOT
            SUMI = SUMI + X(J)*PHIKJ
  100    CONTINUE
C
         SUMO = SUMO + X(K)*RMCWRK(NVIS+K-1) / SUMI
  200 CONTINUE
C
      VISMIX = SUMO
C
      RETURN
      END
C
C---------------------------------------------------------------------
C
      SUBROUTINE MCSCON (T, RMCWRK, CON)
C
C*****precision > double
      IMPLICIT DOUBLE PRECISION (A-H, O-Z), INTEGER (I-N)
C*****END precision > double
C
C*****precision > single
C      IMPLICIT REAL (A-H, O-Z), INTEGER (I-N)
C*****END precision > single
C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C
C  SUBROUTINE MCSCON (T, RMCWRK, CON)
C
C    THIS SUBROUTINE COMPUTES THE ARRAY OF PURE SPECIES CONDUCTIVITIES
C    GIVEN THE TEMPERATURE.
C
C  INPUT-
C    T       - TEMPERATURE
C                  CGS UNITS - K.
C
C  WORK-
C    RMCWRK  - ARRAY OF FLOATING POINT STORAGE AND WORK SPACE.  THE
C              STARTING ADDRESSES FOR THE RMCWRK SPACE ARE STORED IN
C              COMMON /MCMCMC/.
C                  DIMENSION RMCWRK(*) AT LEAST LENRMC.
C
C  OUTPUT-
C    CON     - ARRAY OF SPECIES THERMAL CONDUCTIVITIES.
C                  CGS UNITS - ERG/CM*K*S.
C                  DIMENSION CON(*) AT LEAST KK.
C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C
      DIMENSION CON(*), RMCWRK(*)
C
      COMMON /MCMCMC/ RU, PATMOS, SMALL, NKK, NO, NLITE, INLIN, IKTDIF,
     1                IPVT, NWT, NEPS, NSIG, NDIP, NPOL, NZROT, NLAM,
     2                NETA, NDIF, NTDIF, NXX, NVIS, NXI, NCP, NCROT,
     3                NCINT, NBIND, NEOK, NSGM, NAST, NBST,
     4                NCST, NXL, NR, NWRK, K3
C
      ALOGT = LOG(T)
      CALL MCEVAL (ALOGT, NKK, NO, RMCWRK(NLAM), CON)
C
      DO 25 K = 1, NKK
         CON(K) = EXP(CON(K))
   25 CONTINUE
C
      RETURN
      END
C
C------------------------------------------------------------------
C
      SUBROUTINE MCACON (T, X, RMCWRK, CONMIX)
C
C*****precision > double
      IMPLICIT DOUBLE PRECISION (A-H, O-Z), INTEGER (I-N)
C*****END precision > double
C
C*****precision > single
C      IMPLICIT REAL (A-H, O-Z), INTEGER (I-N)
C*****END precision > single
C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C
C  SUBROUTINE MCACON (T, X, RMCWRK, CONMIX)
C
C    THIS SUBROUTINE COMPUTES THE MIXTURE THERMAL CONDUCTIVITY, GIVEN
C    THE TEMPERATURE AND THE SPECIES MOLE FRACTIONS.
C
C  INPUT-
C    T       - TEMPERATURE
C                  CGS UNITS - K.
C    X       - ARRAY OF MOLE FRACTIONS OF THE MIXTURE.
C                  DIMENSION X(*) AT LEAST KK.
C
C  WORK-
C    RMCWRK  - ARRAY OF FLOATING POINT STORAGE AND WORK SPACE.  THE
C              STARTING ADDRESSES FOR THE RMCWRK SPACE ARE STORED IN
C              COMMON /MCMCMC/.
C                  DIMENSION RMCWRK(*) AT LEAST LENRMC.
C
C  OUTPUT-
C    CONMIX  - MIXTURE THERMAL CONDUCTIVITY
C                  CGS UNITS - ERG/CM*K*S.
C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C
      DIMENSION X(*), RMCWRK(*)
C
      COMMON /MCMCMC/ RU, PATMOS, SMALL, NKK, NO, NLITE, INLIN, IKTDIF,
     1                IPVT, NWT, NEPS, NSIG, NDIP, NPOL, NZROT, NLAM,
     2                NETA, NDIF, NTDIF, NXX, NVIS, NXI, NCP, NCROT,
     3                NCINT, NBIND, NEOK, NSGM, NAST, NBST,
     4                NCST, NXL, NR, NWRK, K3
C
      ONE  = 1.0d0
      ZERO = 0.0d0
C
C       IN THE FOLLOWING CALL:
C         THE PURE SPECIES CONDUCTIVITIES ARE IN RMCWRK(NXI)
C
      ALOGT = dLOG(T)
      CALL MCEVAL (ALOGT, NKK, NO, RMCWRK(NLAM), RMCWRK(NXI))
C
      dSUM = ZERO
      SUMR = ZERO
      DO 100 K = 1, NKK
         RMCWRK(NXI+K-1) = dEXP(RMCWRK(NXI+K-1))
         dSUM =  dSUM  + X(K)*RMCWRK(NXI+K-1)
         SUMR = SUMR + X(K)/RMCWRK(NXI+K-1)
  100 CONTINUE
C
      CONMIX = 0.5d0 * (dSUM + ONE/SUMR)
C
      RETURN
      END
C
C---------------------------------------------------------------------
C
      SUBROUTINE MCSDIF (P, T, KDIM, RMCWRK, DJK)
C
C*****precision > double
      IMPLICIT DOUBLE PRECISION (A-H, O-Z), INTEGER (I-N)
C*****END precision > double
C
C*****precision > single
C      IMPLICIT REAL (A-H, O-Z), INTEGER (I-N)
C*****END precision > single
C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C
C  SUBROUTINE MCSDIF (P, T, KDIM, RMCWRK, DJK)
C
C    THIS SUBROUTINE COMPUTES THE BINARY DIFFUSION COEFFICIENTS, GIVEN
C    THE PRESSURE AND TEMPERATURE.
C
C  INPUT-
C    P       - PRESSURE
C                  CGS UNITS - DYNES/CM**2.
C    T       - TEMPERATURE
C                  CGS UNITS - K.
C    KDIM    - ACTUAL FIRST DIMENSION OF DJK(KDIM,KK)
C
C  WORK-
C    RMCWRK   - ARRAY OF FLOATING POINT STORAGE AND WORK SPACE.  THE
C              STARTING ADDRESSES FOR THE RMCWRK SPACE ARE STORED IN
C              COMMON /MCMCMC/.
C                  DIMENSION RMCWRK(*) AT LEAST LENRMC.
C
C  OUTPUT-
C    DJK     - MATRIX OF BINARY DIFFUSION COEFFICIENTS.  DJK(J,K) IS
C              DIFFUSION COEFFICIENT OF SPECIES J IN SPECIES K.
C                  CGS UNITS - CM**2/S
C                  DIMENSION DJK(KDIM,*) EXACTLY KDIM FOR THE FIRST
C                   DIMENSION AND AT LEAST KK FOR THE SECOND.
C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C
      DIMENSION RMCWRK(*), DJK(KDIM,*)
C
      COMMON /MCMCMC/ RU, PATMOS, SMALL, NKK, NO, NLITE, INLIN, IKTDIF,
     1                IPVT, NWT, NEPS, NSIG, NDIP, NPOL, NZROT, NLAM,
     2                NETA, NDIF, NTDIF, NXX, NVIS, NXI, NCP, NCROT,
     3                NCINT, NBIND, NEOK, NSGM, NAST, NBST,
     4                NCST, NXL, NR, NWRK, K3
C
      PFAC = PATMOS/P
      ALOGT = LOG(T)
C
      DO 100 K = 1, NKK
         ISTART = NDIF + (K-1)*NO*NKK
         CALL MCEVAL (ALOGT, NKK, NO, RMCWRK(ISTART), DJK(1,K))
         DO 100 J = 1, NKK
            DJK(J,K) = EXP(DJK(J,K)) * PFAC
  100 CONTINUE
C
      RETURN
      END
C
C----------------------------------------------------------------------
C
      SUBROUTINE MCADIF (P, T, X, RMCWRK, D)
C
C*****precision > double
      IMPLICIT DOUBLE PRECISION (A-H, O-Z), INTEGER (I-N)
C*****END precision > double
C
C*****precision > single
C      IMPLICIT REAL (A-H, O-Z), INTEGER (I-N)
C*****END precision > single
C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C
C  SUBROUTINE MCADIF (P, T, X, RMCWRK, D)
C
C    THIS SUBROUTINE COMPUTES MIXTURE-AVERAGED DIFFUSION COEFFICIENTS
C    GIVEN THE PRESSURE, TEMPERATURE, AND SPECIES MASS FRACTIONS.
C
C  INPUT-
C    P       - PRESSURE
C                  CGS UNITS - DYNES/CM**2.
C    T       - TEMPERATURE
C                  CGS UNITS - K.
C    X       - ARRAY OF MOLE FRACTIONS OF THE MIXTURE.
C                  DIMENSION X(*) AT LEAST KK.
C
C  WORK-
C    RMCWRK  - ARRAY OF FLOATING POINT STORAGE AND WORK SPACE.  THE
C              STARTING ADDRESSES FOR THE RMCWRK SPACE ARE STORED IN
C              COMMON /MCMCMC/.
C                  DIMENSION RMCWRK(*) AT LEAST LENRMC.
C
C  OUTPUT-
C    D       - ARRAY OF MIXTURE DIFFUSION COEFFICIENTS
C                  CGS UNITS - CM**2/S.
C                  DIMENSION D(*) AT LEAST KK.
C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C
      DIMENSION X(*), D(*), RMCWRK(*)
C
      COMMON /MCMCMC/ RU, PATMOS, SMALL, NKK, NO, NLITE, INLIN, IKTDIF,
     1                IPVT, NWT, NEPS, NSIG, NDIP, NPOL, NZROT, NLAM,
     2                NETA, NDIF, NTDIF, NXX, NVIS, NXI, NCP, NCROT,
     3                NCINT, NBIND, NEOK, NSGM, NAST, NBST,
     4                NCST, NXL, NR, NWRK, K3
C
      CALL MCEDIF (T, NO, NKK, X, RMCWRK(NDIF), RMCWRK(NWT), SMALL,
     1             RMCWRK(NXX), RMCWRK(NBIND), D)
C
      DO 100 K = 1, NKK
         D(K) = D(K) * PATMOS/P
  100 CONTINUE
C
      RETURN
      END
C
C----------------------------------------------------------------------
C
      SUBROUTINE MCEDIF(T, NO, KK, X, COFD, WT, SMALL, XX, DJK, D)
C
C*****precision > double
      IMPLICIT DOUBLE PRECISION (A-H, O-Z), INTEGER (I-N)
C*****END precision > double
C
C*****precision > single
C      IMPLICIT REAL (A-H, O-Z), INTEGER (I-N)
C*****END precision > single
C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C
C  SUBROUTINE MCEDIF(T, NO, KK, X, COFD, WT, SMALL, XX, DJK, D)
C
C    THIS SUBROUTINE IS USED INTERNALLY TO COMPUTE THE MIXTURE
C    DIFFUSION COEFFICIENTS.  NORMALLY NOT CALLED BY THE PACKAGE USER.
C
C  INPUT-
C    T       - TEMPERATURE
C                  CGS UNITS - K.
C    NO      - ORDER OF FIT.
C    KK      - NUMBER OF SPECIES.
C    X       - ARRAY OF MOLE FRACTIONS OF THE MIXTURE.
C                  DIMENSION X(*) AT LEAST KK.
C    COFD    - COEFFICIENTS OF THE FITS FOR THE BINARY DIFFUSION
C              COEFFICIENTS.
C                  DIMENSION COFD(NO,KK,*) EXACTLY NO FOR
C                   FIRST DIMENSION, KK FOR THE SECOND, AND AT
C                   LEAST KK FOR THE THIRD.
C    WT      - ARRAY OF SPECIES MOLECULAR WEIGHTS.
C                   DIMENSION WT(*) AT LEAST KK.
C    SMALL   - A SMALL NUMBER ADDED TO ALL MOLE FRACTIONS BEFORE
C                COMPUTING THE MIXTURE DIFFUSION COEFFICIENTS.
C                THIS PROCESS AVOIDS AN UNDEFINED SITUATION WHEN
C                A PURE SPECIES CONDITION IS APPROACHED.
C    XX      - ARRAY OF MOLE FRACTIONS PLUS "SMALL," TO AVOID THE
C              PROBLEM OF A PURE SPECIES.
C                  DIMENSION XX(*) AT LEAST KK.
C
C  WORK-
C    RMCWRK  - ARRAY OF FLOATING POINT STORAGE AND WORK SPACE.  THE
C              STARTING ADDRESSES FOR THE RMCWRK SPACE ARE STORED IN
C              COMMON /MCMCMC/.
C                  DIMENSION RMCWRK(*) AT LEAST LENRMC.
C
C  OUTPUT-
C    D       - ARRAY OF MIXTURE DIFFUSION COEFFICIENTS
C                  CGS UNITS - CM**2/S.
C                  DIMENSION D(*) AT LEAST KK.
C    DJK     - MATRIX OF BINARY DIFFUSION COEFFICIENTS.  DJK(J,K) IS
C              DIFFUSION COEFFICIENT OF SPECIES J IN SPECIES K.
C                  CGS UNITS - CM**2/S
C                  DIMENSION DJK(KDIM,*) EXACTLY KDIM FOR THE FIRST
C                   DIMENSION AND AT LEAST KK FOR THE SECOND.
C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C
      DIMENSION X(*), COFD(NO,KK,*), WT(*), XX(*), DJK(KK,*), D(*)
C
      ZERO = 0.0
      ALOGT = LOG(T)
      DO 100 K = 1, KK
         CALL MCEVAL (ALOGT, KK, NO, COFD(1,1,K), DJK(1,K) )
  100 CONTINUE
C
      DO 150 K = 1, KK
         DO 150 J = 1, KK
            DJK(J,K) = EXP(DJK(J,K))
  150 CONTINUE
C
      WTM = ZERO
      DO 175 K = 1, KK
         WTM = WTM + WT(K)*X(K)
         XX(K) = X(K) + SMALL
  175 CONTINUE
C
      DO 300 K = 1, KK
C
         SUMXW = ZERO
         SUMXOD = ZERO
C
         DO 200 J = 1, KK
            IF (J .NE. K) THEN
               SUMXW  = SUMXW  + XX(J)*WT(J)
               SUMXOD = SUMXOD + XX(J)/DJK(J,K)
            ENDIF
  200    CONTINUE
C
         D(K) = SUMXW/(WTM*SUMXOD)
C
  300 CONTINUE
C
      RETURN
      END
C
C---------------------------------------------------------------------
C
      SUBROUTINE MCMDIF (P, T, X, KDIM, IMCWRK, RMCWRK, D)
C
C*****precision > double
      IMPLICIT DOUBLE PRECISION (A-H, O-Z), INTEGER (I-N)
C*****END precision > double
C
C*****precision > single
C      IMPLICIT REAL (A-H, O-Z), INTEGER (I-N)
C*****END precision > single
C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C
C  SUBROUTINE MCMDIF (P, T, X, KDIM, IMCWRK, RMCWRK, D)
C
C    THIS SUBROUTINE COMPUTES THE ORDINARY MULTICOMPONENT DIFFUSION
C    COEFFICIENTS, GIVEN THE PRESSURE, TEMPERATURE, AND MOLE FRACTIONS.
C
C  INPUT-
C    P       - PRESSURE
C                  CGS UNITS - DYNES/CM**2.
C    T       - TEMPERATURE
C                  CGS UNITS - K.
C    X       - ARRAY OF MOLE FRACTIONS OF THE MIXTURE.
C                  DIMENSION X(*) AT LEAST KK.
C    KDIM    - ACTUAL FIRST DIMENSION OF D(KDIM,KK).  KDIM MUST BE AT
C                LEAST THE NUMBER OF SPECIES, KK.
C
C  WORK-
C    IMCWRK  - ARRAY OF INTEGER STORAGE AND WORK SPACE.  THE
C              STARTING ADDRESSES FOR THE IMCWRK SPACE ARE STORED IN
C              COMMON /MCMCMC/.
C                  DIMENSION IMCWRK(*) AT LEAST LENIMC.
C    RMCWRK  - ARRAY OF FLOATING POINT STORAGE AND WORK SPACE.  THE
C              STARTING ADDRESSES FOR THE RMCWRK SPACE ARE STORED IN
C              COMMON /MCMCMC/.
C                  DIMENSION RMCWRK(*) AT LEAST LENRMC.
C
C  OUTPUT-
C    D       - MATRIX OF ORDINARY MULTICOMPONENT DIFFUSION COEFFICIENTS.
C                  CGS UNITS - CM**2/S
C                  DIMENSION DJK(KDIM,*) EXACTLY KDIM FOR THE FIRST
C                   DIMENSION AND AT LEAST KK FOR THE SECOND.
C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C
      DIMENSION X(*), IMCWRK(*), RMCWRK(*), D(KDIM,*)
C
      COMMON /MCMCMC/ RU, PATMOS, SMALL, NKK, NO, NLITE, INLIN, IKTDIF,
     1                IPVT, NWT, NEPS, NSIG, NDIP, NPOL, NZROT, NLAM,
     2                NETA, NDIF, NTDIF, NXX, NVIS, NXI, NCP, NCROT,
     3                NCINT, NBIND, NEOK, NSGM, NAST, NBST,
     4                NCST, NXL, NR, NWRK, K3
C
      CALL MCORDF (P, T, X, NKK, KDIM, SMALL, RMCWRK(NWT), RMCWRK,
     1             RMCWRK(NXX), RMCWRK(NBIND), RMCWRK(NXL),
     2             RMCWRK(NWRK), IMCWRK(IPVT), D)
C
      RETURN
      END
C
C---------------------------------------------------------------------
C
      SUBROUTINE MCMCDT (P, T, X, IMCWRK, RMCWRK, ICKWRK, CKWRK,
     1                   DT, COND)
C
C*****precision > double
      IMPLICIT DOUBLE PRECISION (A-H, O-Z), INTEGER (I-N)
C*****END precision > double
C
C*****precision > single
C      IMPLICIT REAL (A-H, O-Z), INTEGER (I-N)
C*****END precision > single
C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C
C   SUBROUTINE MCMCDT (P, T, X, KDIM, IMCWRK, RMCWRK, ICKWRK, CKWRK,
C  1                   DT, COND)
C
C    THIS SUBROUTINE COMPUTES THE THERMAL DIFFUSION COEFFICIENTS, AND
C    MIXTURE THERMAL CONDUCTIVITIES, GIVEN THE PRESSURE, TEMPERATURE,
C    AND MOLE FRACTIONS.
C
C  INPUT-
C    P       - PRESSURE
C                  CGS UNITS - DYNES/CM**2.
C    T       - TEMPERATURE
C                  CGS UNITS - K.
C    X       - ARRAY OF MOLE FRACTIONS OF THE MIXTURE.
C                  DIMENSION X(*) AT LEAST KK.
C
C  WORK-
C    IMCWRK  - ARRAY OF INTEGER STORAGE AND WORK SPACE.  THE
C              STARTING ADDRESSES FOR THE IMCWRK SPACE ARE STORED IN
C              COMMON /MCMCMC/.
C                  DIMENSION IMCWRK(*) AT LEAST LENIMC.
C    RMCWRK  - ARRAY OF FLOATING POINT STORAGE AND WORK SPACE.  THE
C              STARTING ADDRESSES FOR THE RMCWRK SPACE ARE STORED IN
C              COMMON /MCMCMC/.
C                  DIMENSION RMCWRK(*) AT LEAST LENRMC.
C    ICKWRK  - ARRAY OF INTEGER CHEMKIN STORAGE AND WORK SPACE.  SEE
C                  CHEMKIN DOCUMENTATION.
C                  DIMENSION ICKWRK(*) AT LEAST LENIWRK.
C    CKWRK   - ARRAY OF FLOATING POINT CHEMKIN STORAGE AND WORK SPACE.
C              SEE CHEMKIN DOCUMENTATION.
C                  DIMENSION CKWRK(*) AT LEAST LENWRK.
C
C  OUTPUT-
C    DT      - VECTOR OF THERMAL MULTICOMPONENT DIFFUSION COEFFICIENTS.
C                  CGS UNITS - GM/(CM*SEC)
C    COND    - MIXTURE THERMAL CONDUCTIVITY
C                  CGS UNITS - ERG/(CM*K*S).
C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C
      DIMENSION X(*), IMCWRK(*), RMCWRK(*), ICKWRK(*), CKWRK(*), DT(*)
C
      COMMON /MCMCMC/ RU, PATMOS, SMALL, NKK, NO, NLITE, INLIN, IKTDIF,
     1                IPVT, NWT, NEPS, NSIG, NDIP, NPOL, NZROT, NLAM,
     2                NETA, NDIF, NTDIF, NXX, NVIS, NXI, NCP, NCROT,
     3                NCINT, NBIND, NEOK, NSGM, NAST, NBST,
     4                NCST, NXL, NR, NWRK, K3
C
      CALL MCLMDT (P, T, X, NKK, K3, SMALL, RMCWRK(NWT), RMCWRK(NEOK),
     1             RMCWRK(NZROT), IMCWRK(INLIN), RMCWRK(NEPS),
     2             ICKWRK, CKWRK, RMCWRK, RMCWRK(NXX), RMCWRK(NVIS),
     3             RMCWRK(NAST), RMCWRK(NBST), RMCWRK(NCST),
     4             RMCWRK(NXI),  RMCWRK(NCP), RMCWRK(NCROT),
     5             RMCWRK(NCINT), RMCWRK(NXL), RMCWRK(NR), 
     6             RMCWRK(NBIND), IMCWRK(IPVT), DT, COND)
C
      RETURN
      END
C
C---------------------------------------------------------------------
C
      SUBROUTINE MCEVAL (TF, KK, NO, COF, VAL)
C
C*****precision > double
      IMPLICIT DOUBLE PRECISION (A-H, O-Z), INTEGER (I-N)
C*****END precision > double
C
C*****precision > single
C      IMPLICIT REAL (A-H, O-Z), INTEGER (I-N)
C*****END precision > single
C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C
C  SUBROUTINE MCEVAL (TF, KK, NO, COF, VAL)
C
C    THIS SUBROUTINE USES HORNERS ALGORITHM TO EVALUATE A POLYNOMIAL
C    FIT.  THIS ROUTINE IS NOT NORMALLY CALLED BY THE PACKAGE USER.
C
C  INPUT-
C    TF      - INDEPENDENT VARIABLE OF FIT. EITHER TEMPERATURE
C              OR LOG TEMPERATURE.
C    KK      - NUMBER OF SPECIES.
C    NO      - ORDER OF FIT.
C    COF     - MATRIX OF FIT COEFFICIENTS.  COF(N,K) IS THE NTH
C              COEFFICIENT OF A FIT FOR KTH SPECIES PROPERTY.
C                 DIMENSION COF(NO,*) EXACTLY NO FOR THE FIRST
C                 DIMENSION AND AT LEAST KK FOR THE SECOND.
C
C  OUTPUT-
C    VAL     - ARRAY OF VALUES, EVALUATED FROM THE FIT AT TF.
C                 DIMENSION VAL(*) AT KEAST KK.
C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C
      DIMENSION COF(NO,*), VAL(*)
C
      NOM1 = NO-1
C
      DO 200 K = 1, KK
         B = COF(NO,K)
         DO 100 I = 1, NOM1
            B = COF(NO-I,K) + B*TF
  100    CONTINUE
         VAL(K) = B
  200 CONTINUE
C
      RETURN
      END
C
C---------------------------------------------------------------------
C
      SUBROUTINE MCEPSG (KK, EPS, SIG, DIP, POL, EOK, SGM)
C
C*****precision > double
      IMPLICIT DOUBLE PRECISION (A-H, O-Z), INTEGER (I-N)
C*****END precision > double
C
C*****precision > single
C      IMPLICIT REAL (A-H, O-Z), INTEGER (I-N)
C*****END precision > single
C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C
C  SUBROUTINE MCEPSG (KK, EPS, SIG, DIP, POL, EOK, SGM)
C
C    THIS SUBROUTINE COMPUTES THE REDUCED WELL DEPTH EOK(I,J) AND
C    COLLISION DIAMETER SGM(I,J) FOR EACH I,J SPECIES PAIR.  THE
C    ROUTINE IS CALLED ONLY ONCE BY THE INITIALIZATION SUBROUTINE MCINIT
C    THIS ROUTINE IS NORMALLY NOT CALLED BY THE USER.
C
C  INPUT-
C    KK      - NUMBER OF SPECIES
C    EPS     - ARRAY OF LENNARD-JONES POTENTIAL WELL DEPTHS.
C                  CGS UNITS - K.
C                  DIMENSION EPS(*) AT LEAST KK
C    SIG     - ARRAY OF LENNARD-JONES COLLISION DIAMETERS.
C                  UNITS - ANGSTROMS.
C                  DIMENSION SIG(*) AT LEAST KK
C    DIP     - ARRAY OF DIPOLE MOMENTS
C                  UNITS - DEBYE
C                  DIMENSION DIP(*) AT LEAST KK
C    POL     - ARRAY OF POLARIZABILITIES.
C                  UNITS - ANGSTROMS**3.
C                  DIMENSION POL(*) AT LEAST KK
C
C  OUTPUT-
C    EOK     - MATRIX OF REDUCED WELL DEPTHS FOR EACH SPECIES PAIR.
C                  UNITS - K
C                  DIMENSION EOK(KDIM,*) EXACTLY KDIM FOR THE FIRST
C                    DIMENSION, AND AT LEAST KK FOR THE SECOND.
C    SGM     - MATRIX OF REDUCED COLLISION DIAMETERS FOR EACH SPECIES
C              PAIR.
C                  UNITS - ANGSTROMS.
C                  DIMENSION SGM(KDIM,*) EXACTLY KDIM FOR THE FIRST
C                    DIMENSION, AND AT LEAST KK FOR THE SECOND.
C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C
      DIMENSION EPS(*), SIG(*), DIP(*), POL(*), EOK(KK,*), SGM(KK,*)
C
      DATA FDTCGS/1.0E-18/, FATCM/1.0E8/,
     1     DIPMIN/1.0E-20/, BOLTZ/1.38056E-16/
C
      ONE = 1.0
C         COMPUTE AND STORE EPS/K AND SIGMA FOR ALL PAIRS
C
      DO 1000 J = 1, KK
C
         DO 1000 K = 1, J
C
           IF((DIP(J).LT.DIPMIN .AND. DIP(K).GT.DIPMIN)) THEN
C
C                K IS POLAR, J IS NONPOLAR
C
              XI = ONE + 0.25*(POL(J)/SIG(J)**3) *
     1                     (FDTCGS**2*FATCM**3/BOLTZ) *
     2                     (DIP(K)**2/(EPS(K)*SIG(K)**3)) *
     3                      SQRT(EPS(K)/EPS(J))
              SGM(K,J) = 0.5 * (SIG(J)+SIG(K)) * XI**(-ONE/6.0)
              SGM(J,K) = SGM(K,J)
              EOK(K,J) = SQRT(EPS(J)*EPS(K)) * XI**2
              EOK(J,K) = EOK(K,J)
C
          ELSE IF((DIP(J).GT.DIPMIN .AND. DIP(K).LT.DIPMIN)) THEN
C
C             J IS POLAR, K IS NONPOLAR
C
              XI = ONE + 0.25*(POL(K)/SIG(K)**3) *
     1                     (FDTCGS**2*FATCM**3/BOLTZ) *
     2                     (DIP(J)**2/(EPS(J)*SIG(J)**3)) *
     3                      SQRT(EPS(J)/EPS(K))
              SGM(K,J) = 0.5 * (SIG(J)+SIG(K)) * XI**(-ONE / 6.0)
              SGM(J,K) = SGM(K,J)
              EOK(K,J) = SQRT(EPS(J)*EPS(K)) * XI**2
              EOK(J,K) = EOK(K,J)
C
          ELSE
C
C              NORMAL CASE, EITHER BOTH POLAR OR BOTH NONPOLAR
C
              SGM(K,J) = 0.5 * (SIG(J) + SIG(K))
              SGM(J,K) = SGM(K,J)
              EOK(K,J) = SQRT(EPS(J)*EPS(K))
              EOK(J,K) = EOK(K,J)
C
          ENDIF
1000  CONTINUE
C
      RETURN
      END
C
C---------------------------------------------------------------------
C
      SUBROUTINE MCORDF (P, T, X, KK, KDIM, SMALL, WT, RMCWRK, XX,
     1                   BINDIF, XL0000, WORK, IPVT, D)
C
C*****precision > double
      IMPLICIT DOUBLE PRECISION (A-H, O-Z), INTEGER (I-N)
      DOUBLE PRECISION ONE,         ZERO
      PARAMETER       (ONE = 1.0D0, ZERO = 0.0D0)
C*****END precision > double
C
C*****precision > single
C      IMPLICIT REAL (A-H, O-Z), INTEGER (I-N)
C      REAL       ONE,       ZERO
C      PARAMETER (ONE = 1.0, ZERO = 0.0)
C*****END precision > single
C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C
C  SUBROUTINE MCORDF (P, T, X, KK, KDIM, SMALL, WT, RMCWRK, XX,
C 1                   BINDIF, XL0000, WORK, IPVT, D)
C
C    THIS SUBROUTINE COMPUTES ORDINARY MULTICOMPONENT DIFFUSION COEFFICI
C    COEFFICIENT MATRIX.  IT DOES SO BY COMPUTING THE INVERSE OF THE
C    L00,00 MATRIX.  THIS ROUTINE IS NOT NORMALLY CALLED DIRECTLY BY THE
C    USER; THE USER CALLS MCMDIF, WHICH IN TURN CALLS MCORDF.
C
C  INPUT-
C    P       - PRESSURE
C                  CGS UNITS - DYNES/CM**2.
C    T       - TEMPERATURE
C                  CGS UNITS - K.
C    X       - ARRAY OF MOLE FRACTIONS OF THE MIXTURE.
C                  DIMENSION X(*) AT LEAST KK.
C    KK      - NUMBER OF SPECIES.
C    KDIM    - ACTUAL FIRST DIMENSION OF D(KDIM,KK).  KDIM MUST BE AT
C                LEAST THE NUMBER OF SPECIES, KK.
C    SMALL   - THE MOLE FRACTIONS USED IN THE TRANSPORT COMPUTATION
C              ARE GIVEN BY XX(K) = X(K) + SMALL.
C    WT      - THE ARRAY OF SPECIES MOLECULAR WEIGHTS.
C                   DIMENSION WT(*) AT LEAST KK.
C
C  WORK AND SCRATCH SPACE
C    RMCWRK  - ARRAY OF FLOATING POINT STORAGE AND WORK SPACE.  THE
C              STARTING ADDRESSES FOR THE RMCWRK SPACE ARE STORED IN
C              COMMON /MCMCMC/.
C                  DIMENSION RMCWRK(*) AT LEAST LENRMC.
C    XX      - THE SPECIES MOLE FRACTION ARRAY THAT IS USED IN THE
C              TRANSPORT COMPUTATION.  XX(K) = X(K) + SMALL.
C    BINDIF  - MATRIX OF BINARY DIFFUSION COEFFICIENTS.
C                  CGS UNITS - CM**2/S
C                  DIMENSION BINDIF(KDIM,*) EXACTLY KDIM FOR THE FIRST
C                   DIMENSION AND AT LEAST KK FOR THE SECOND.
C    XL0000  - THE L00,00 MATRIX.
C                  DIMENSION L0000(KDIM,*) EXACTLY KDIM FOR THE FIRST
C                   DIMENSION AND AT LEAST KK FOR THE SECOND.
C    WORK    - ARRAY OF WORK SPACE FOR THE INVERSION OF THE L00,00
C              MATRIX BY THE LINPACK ROUTINES SGEFA AND SGEDI.
C                  DIMENSION WORK(*) AT LEAST KK.
C    IPVT    - ARRAY OF PIVOT INDICES FOR THE INVERSION OF THE L00,00
C              MATRIX BY THE LINPACK ROUTINES SGEFA AND SGEDI.
C                  DIMENSION IPVT(*) AT LEAST KK.
C
C  OUTPUT-
C    D       - MATRIX OF ORDINARY MULTICOMPONENT DIFFUSION COEFFICIENTS.
C                  CGS UNITS - CM**2/S
C                  DIMENSION DJK(KDIM,*) EXACTLY KDIM FOR THE FIRST
C                   DIMENSION AND AT LEAST KK FOR THE SECOND.
C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C
      DIMENSION X(*), WT(*), BINDIF(KK,*), XL0000(KK,*), IPVT(*),
     1          DET(2), WORK(*), XX(*), RMCWRK(*), D(KDIM,*)

C local variables
C
      INTEGER   I, J, K
C
      SAVE      JOB
      DATA      JOB /1/
C
C
C         SET MINIMUM MOLE FRACTION TO SMALL
C
      DO 50 I = 1, KK
         XX(I) = MAX( X(I) , SMALL)
   50 CONTINUE
C
C        EVALUATE THE BINARY DIFFUSION COEFFICIENTS
C
      CALL MCSDIF (P, T, KK, RMCWRK, BINDIF)
C
C        ASSEMBLE L00,00
C
      PFAC = 16.0 * T / (25.0 * P)
      DO 200 I = 1, KK
        SUM = -XX(I) / BINDIF(I,I)
        DO 90 J = 1, KK
          SUM = SUM + XX(J) / BINDIF(I,J)
  90    CONTINUE 
        SUM = SUM / WT(I)
        DO 100 J = 1, KK
          XL0000(I,J) = PFAC * XX(J) *
     $                  (WT(J) * SUM + XX(I) / BINDIF(I,J))
 100    CONTINUE
        XL0000(I,I) = ZERO
 200  CONTINUE

C
C        INVERT L00,00 USING LAPACK or LINPACK
C
*
C*****precision > double - lapack
      CALL DGETRF (KK, KK, XL0000, KK, IPVT, INFO)
      IF (INFO .NE. 0) THEN
          WRITE (6,*) ' ERROR IN DGETRF, INFO = ', INFO
          STOP
      ENDIF
      CALL DGETRI(KK, XL0000, KK, IPVT, WORK, KK, INFO)
      IF (INFO .NE. 0) THEN
          WRITE (6,*) ' ERROR IN DGETRI, INFO = ', INFO
          STOP
      ENDIF
C*****END precision > double - lapack
*
C*****precision > double - linpack
C      CALL DGEFA (XL0000, KK, KK, IPVT, INFO)
C      IF (INFO .NE. 0) THEN
C        WRITE (6, *) ' ERROR IN DGEFA, INFO = ', INFO
C        STOP
C      ENDIF
C      CALL DGEDI (XL0000, KK, KK, IPVT, DET, WORK, JOB)
C*****END precision > double - linpack
*
C*****precision > single - lapack
C      CALL SGETRF (KK, KK, XL0000, KK, IPVT, INFO)
C      IF (INFO .NE. 0) THEN
C          WRITE (6,*) ' ERROR IN SGETRF, INFO = ', INFO
C          STOP
C      ENDIF
C      CALL SGETRI(KK, XL0000, KK, IPVT, WORK, KK, INFO)
C      IF (INFO .NE. 0) THEN
C          WRITE (6,*) ' ERROR IN SGETRI, INFO = ', INFO
C          STOP
C      ENDIF
C*****END precision > single - lapack
*
C*****precision > single - linpack
C      CALL SGEFA (XL0000, KK, KK, IPVT, INFO)
C      IF (INFO .NE. 0) THEN
C        WRITE (6, *) ' ERROR IN SGEFA, INFO = ', INFO
C        STOP
C      ENDIF
C      CALL SGEDI (XL0000, KK, KK, IPVT, DET, WORK, JOB)
C*****END precision > single - linpack
*
C
C        COMPUTE THE ORDINARY MULTICOMPONENT DIFFUSION COEFFICIENTS
C
      SUM = ZERO
      DO 400 I = 1, KK
        SUM = SUM + WT(I) * X(I)
  400 CONTINUE
      DO 500 J = 1, KK
         PFAC_J = PFAC * SUM / WT(J)
         DO 500 I = 1, KK
            D(I,J) = PFAC_J * XX(I) * (XL0000(I,J)-XL0000(I,I))
  500 CONTINUE
C
      RETURN
      END
C
C---------------------------------------------------------------------
C
      SUBROUTINE MCLMDT (P, T, X, KK, KK3, SMALL, WT, EOK, ZROT, LIN,
     1                   EPS, ICKWRK, CKWRK, RMCWRK, XX, VIS, ASTAR,
     2                   BSTAR, CSTAR, XI, CPOR, CROTOR, CINTOR, XL,
     3                   R, BINDIF, IPVT, DT, COND)
C
C*****precision > double
      IMPLICIT DOUBLE PRECISION (A-H, O-Z), INTEGER (I-N)
      DOUBLE PRECISION ONE,         ZERO
      PARAMETER       (ONE = 1.0D0, ZERO = 0.0D0)
      DOUBLE PRECISION RU, PI, PI32O2, P2O4P2,  PI32
      PARAMETER (RU=8.314510D+07, PI= 3.1415926535897932D+00,
     $                        PI32O2= 2.7841639984158539D+00,
     $                        P2O4P2= 4.4674011002723397D+00,
     $                          PI32= 5.5683279968317078D+00)
C*****END precision > double
C
C*****precision > single
C      IMPLICIT REAL (A-H, O-Z), INTEGER (I-N)
C      REAL       ONE,       ZERO
C      PARAMETER (ONE = 1.0E0, ZERO = 0.0E0)
C      REAL       RU, PI, PI32O2, P2O4P2,  PI32
C      PARAMETER (RU=8.314510E+07, PI= 3.1415926535897932E+00,
C     $                        PI32O2= 2.7841639984158539E+00,
C     $                        P2O4P2= 4.4674011002723397E+00,
C     $                          PI32= 5.5683279968317078E+00)
C*****END precision > single
C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C
C  SUBROUTINE MCLMDT (P, T, X, KK, KK3, SMALL, WT, EOK, ZROT,
C 1                   LIN, EPS, ICKWRK, CKWRK, RMCWRK, XX, VIS,
C 2                   ASTAR, BSTAR, CSTAR, XI, CPOR, CROTOR,
C 3                   CINTOR, XL, R, BINDIF, IPVT,
C 4                   DT, COND)
C
C
C    THIS SUBROUTINE COMPUTES THE THERMAL CONDUCTIVITY, AND THE THERMAL
C    DIFFUSION COEFFICIENT ARRAY.  IT DOES SO BY FIRST FORMING THE L
C    MATRIX, AND THEN SOLVING EQ. 24A.  THIS ROUTINE IS NOT NORMALLY CAL
C    DIRECTLY BY THE USER; THE USER CALLS MCMCDT, WHICH IN TURN CALLS
C
C  INPUT-
C    P       - PRESSURE
C                  CGS UNITS - DYNES/CM**2.
C    T       - TEMPERATURE
C                  CGS UNITS - K.
C    X       - ARRAY OF MOLE FRACTIONS OF THE MIXTURE.
C                  DIMENSION X(*) AT LEAST KK.
C    KK      - NUMBER OF SPECIES.
C    KK3     - THREE TIMES THE NUMBER OF SPECIES.  THE SIZE OF THE L
C              MATRIX IS KK3 * KK3.
C    SMALL   - THE MOLE FRACTIONS USED IN THE TRANSPORT COMPUTATION
C              ARE GIVEN BY XX(K) = X(K) + SMALL.
C    WT      - THE ARRAY OF SPECIES MOLECULAR WEIGHTS.
C                   DIMENSION WT(*) AT LEAST KK.
C    EOK     - MATRIX OF REDUCED WELL DEPTHS FOR EACH SPECIES PAIR.
C                  UNITS - K
C                  DIMENSION EOK(KK,*) EXACTLY KDIM FOR THE FIRST
C                    DIMENSION, AND AT LEAST KK FOR THE SECOND.
C    ZROT    - ARRAY OF ROTATIONAL COLLISION NUMBERS EVALUATED
C              AT 298K.
C                  UNITS - NONE
C                  DIMENSION ZROT(*) AT LEAST KK
C    LIN     - ARRAY OF FLAGS INDICATING WHETHER THE MOLECULE
C              LINEAR OR NOT.
C              NLIN=0, SINGLE ATOM.
C              NLIN=1, LINEAR MOLECULE.
C              NLIN=2, NONLINEAR MOLECULE.
C                  UNITS - NONE.
C                  DIMENSION NLIN(*) AT LEAST KK
C    EPS     - ARRAY OF LENNARD-JONES POTENTIAL WELL DEPTHS.
C                  CGS UNITS - K.
C                  DIMENSION EPS(*) AT LEAST KK
C
C  WORK AND SCRATCH SPACE
C    ICKWRK  - ARRAY OF CHEMKIN INTEGER WORK SPACE.
C    CKWRK   - ARRAY OF CHEMKIN REAL WORK SPACE.
C    RMCWRK  - ARRAY OF FLOATING POINT STORAGE AND WORK SPACE.  THE
C              STARTING ADDRESSES FOR THE RMCWRK SPACE ARE STORED IN
C              COMMON /MCMCMC/.
C                  DIMENSION RMCWRK(*) AT LEAST LENRMC.
C    XX      - THE SPECIES MOLE FRACTION ARRAY THAT IS USED IN THE
C              TRANSPORT COMPUTATION.  XX(K) = X(K) + SMALL.
C    VIS     - ARRAY OF SPECIES VISCOSITIES.  EVALUATED FROM MCSVIS.
C                  CGS UNITS - GM/CM-S
C                  DIMENSION VIS(*) AT LEAST KK.
C    ASTAR   - MATRIX OF COLLISION INTEGRALS A*, FOR EACH SPECIES PAIR.
C                  DIMENSION ASTAR(KDIM,*) EXACTLY KDIM FOR THE FIRST
C                   DIMENSION AND AT LEAST KK FOR THE SECOND.
C    BSTAR   - MATRIX OF COLLISION INTEGRALS B*, FOR EACH SPECIES PAIR.
C                  DIMENSION BSTAR(KDIM,*) EXACTLY KDIM FOR THE FIRST
C                   DIMENSION AND AT LEAST KK FOR THE SECOND.
C    CSTAR   - MATRIX OF COLLISION INTEGRALS C*, FOR EACH SPECIES PAIR.
C                  DIMENSION CSTAR(KDIM,*) EXACTLY KDIM FOR THE FIRST
C                   DIMENSION AND AT LEAST KK FOR THE SECOND.
C    XI      - ARRAY COLLISION NUMBERS FOR THE TRANSFER OF ROTATIONAL
C              ENERGY OF SPECIES I INTO TRANSLATIONAL ENERGY OF
C              SPECIES J (EQ. 42).  WE ASSUME THAT ALL XI(I,J) = XI(I,I)
C              SEE P. 132 FOR DISCUSSION.
C    CPOR    - ARRAY OF DIMENSIONLESS SPECIFIC HEATS, CP/R.  EVALUATED
C              FROM CKCPOR.
C                   DIMENSION CPOR(*) AT LEAST KK.
C    CROT    - ARRAY OF DIMENSIONLESS ROTATIONAL CONTRIBUTIONS TO THE
C              SPECIES SPECIFIC HEATS.
C                   DIMENSION CROT(*) AT LEAST KK.
C    CINT    - ARRAY OF DIMENSIONLESS INTERNAL CONTRIBUTIONS TO THE
C              SPECIES SPECIFIC HEATS.
C                   DIMENSION CINT(*) AT LEAST KK.
C    XL      - THE L MATRIX, EQ. 43 AND 49.
C                  DIMENSION XL(3*KK,*) EXACTLY KDIM FOR THE FIRST
C                   DIMENSION AND AT LEAST 3*KK FOR THE SECOND.
C    R       - ARRAY OF RIGHT HAND SIDES OF EQ. 24A.
C                   DIMENSION R(*) AT LEAST 3*KK.
C    BINDIF  - MATRIX OF BINARY DIFFUSION COEFFICIENTS.
C                  CGS UNITS - CM**2/S
C                  DIMENSION BINDIF(KK,*) EXACTLY KDIM FOR THE FIRST
C                   DIMENSION AND AT LEAST KK FOR THE SECOND.
C    IPVT    - ARRAY OF PIVOT INDICES FOR THE INVERSION OF THE XL
C              MATRIX BY THE LINPACK ROUTINES SGEFA AND SGEDI.
C                  DIMENSION IPVT(*) AT LEAST 3*KK.
C
C  OUTPUT-
C    DT      - VECTOR OF THERMAL MULTICOMPONENT DIFFUSION COEFFICIENTS.
C                  CGS UNITS - GM/(CM*SEC)
C                  DIMENSION DT(*) AT LEAST KK.
C    COND    - MIXTURE THERMAL CONDUCTIVITY
C                  CGS UNITS - ERG/(CM*K*S).
C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C

      DIMENSION X(*), ICKWRK(*), CKWRK(*), RMCWRK(*), WT(*), XX(*),
     1          VIS(*), EOK(KK,*), ZROT(*), LIN(*), EPS(*),
     2          ASTAR(KK,*), BSTAR(KK,*), CSTAR(KK,*), XI(*), CPOR(*),
     3          CROTOR(*), CINTOR(*), XL(KK3,*), R(*), BINDIF(KK,*),
     4          IPVT(*), DT(*), FITAST(7), FITBST(7), FITCST(7)
      SAVE      FITAST, FITBST, FITCST
C
C            FITS OF A*, B*, AND C* AS FUNCTIONS OF LN(T*)
C
      DATA FITAST / .1106910525E+01, -.7065517161E-02,
     1             -.1671975393E-01,  .1188708609E-01,
     2              .7569367323E-03, -.1313998345E-02,
     3              .1720853282E-03/
C
      DATA FITBST / .1199673577E+01, -.1140928763E+00,
     1             -.2147636665E-02,  .2512965407E-01,
     2             -.3030372973E-02, -.1445009039E-02,
     3              .2492954809E-03/
C
      DATA FITCST / .8386993788E+00,  .4748325276E-01,
     1              .3250097527E-01, -.1625859588E-01,
     2             -.2260153363E-02,  .1844922811E-02,
     3             -.2115417788E-03/
C
C
C         SET MINIMUM MOLE FRACTION TO SMALL, 
C            (Note, the possibility of negative mole fractions 
C             necessitates the use of the MAX function ).
C
      DO 50 K = 1, KK
        XX(K) = MAX( X(K) ,  SMALL)
   50 CONTINUE
C
C          DETERMINE A*, B*, AND C* FOR EACH SPECIES
C              Note, these are symmetric matrices because EOK(I,J)
C                    is symmetric
C
      DO 100 J = 1, KK
         DO 100 I = 1, J
C
         TSLOG = LOG ( T/EOK(I,J) )

         T1 = TSLOG
         T2 = TSLOG*T1
         T3 = TSLOG*T2
         T4 = TSLOG*T3
         T5 = TSLOG*T4
         T6 = TSLOG*T5
         ASTAR(I,J) = FITAST(1)    + FITAST(2)*T1 + FITAST(3)*T2 +
     1                FITAST(4)*T3 + FITAST(5)*T4 + FITAST(6)*T5 +
     2                FITAST(7)*T6
         ASTAR(J,I) = ASTAR(I,J)
         BSTAR(I,J) = FITBST(1)    + FITBST(2)*T1 + FITBST(3)*T2 +
     1                FITBST(4)*T3 + FITBST(5)*T4 + FITBST(6)*T5 +
     2                FITBST(7)*T6
         BSTAR(J,I) = BSTAR(I,J)
         CSTAR(I,J) = FITCST(1)    + FITCST(2)*T1 + FITCST(3)*T2 +
     1                FITCST(4)*T3 + FITCST(5)*T4 + FITCST(6)*T5 +
     2                FITCST(7)*T6
         CSTAR(J,I) = CSTAR(I,J)
  100 CONTINUE
C
C        EVALUATE THE BINARY DIFFUSION COEFFICIENTS AND VISCOSITY
C
      CALL MCSDIF (P, T, KK, RMCWRK, BINDIF)
      CALL MCSVIS (T, RMCWRK, VIS)
C
      PFAC = 1.2 * RU * T / P
      DO 150 K = 1, KK
C
C        EVALUATE BINARY SELF-DIFFUSION COEFFICIENTS FROM VISCOSITY
C
         BINDIF(K,K) = PFAC * ASTAR(K,K) * VIS(K) / WT(K)
C
C         COMPUTE PARKER CORRECTION FOR ZROT
C
         DD = EPS(K) / T
         DR = EPS(K) / 298.0
         SQRTDD = SQRT(DD)
         SQRTDR = SQRT(DR)
         DD32 = SQRTDD*DD
         DR32 = SQRTDR*DR
         XI(K) = ( (ONE + PI32O2*SQRTDR + P2O4P2*DR + PI32*DR32) /
     1             (ONE + PI32O2*SQRTDD + P2O4P2*DD + PI32*DD32)  )
     2            * MAX(ONE, ZROT(K))
  150 CONTINUE
C
C         ROTATIONAL AND INTERNAL PARTS OF SPECIFIC HEAT
C
      CALL CKCPOR (T, ICKWRK, CKWRK, CPOR)
      DO 400 K = 1, KK
         IF (LIN(K) .EQ. 0) THEN
            CROTOR(K) = ZERO
            CINTOR(K) = ZERO
         ELSEIF (LIN(K) .EQ. 1) THEN
            CROTOR(K) = ONE
            CINTOR(K) = CPOR(K) - 2.5
         ELSEIF (LIN(K) .EQ. 2) THEN
            CROTOR(K) = 1.5
            CINTOR(K) = CPOR(K) - 2.5
         ENDIF
  400 CONTINUE
C
C        ASSEMBLE L00,00
C
      PFAC = 16.0 * T / (25.0 * P)
      DO 600 I = 1, KK
         SUM = - XX(I) / BINDIF(I,I)
         DO 450 K = 1, KK
           SUM = SUM + XX(K) / BINDIF(I,K)
  450    CONTINUE
         SUM = SUM / WT(I)
         DO 500 J = 1, KK
           XL(I,J) =   PFAC * XX(J) * 
     $                        (WT(J) * SUM + XX(I) / BINDIF(J,I))
500      CONTINUE
         XL(I,I) = ZERO
  600 CONTINUE
C
C         ASSEMBLE L00,10 and L10,00
C
      PFAC = 8.0 * T / (5.0 * P)
      DO 1200 J = 1, KK
         WTJ_TMP = WT(J)
         XJ_TMP  = X(J)
         SUM     = ZERO 
         DO 1150 I = 1, KK
            XL(I, J+KK) = -  PFAC * XX(I) * XJ_TMP * WT(I) 
     1                       * (1.2*CSTAR(J,I) - ONE) /
     2                         ((WTJ_TMP + WT(I)) * BINDIF(J,I))
            XL(J+KK, I) = XL(I, J+KK)
            SUM = SUM   + XL(I, J+KK)
 1150    CONTINUE
         XL(J, J+KK) = XL(J, J+KK) - SUM 
         XL(J+KK, J) = XL(J, J+KK)
 1200 CONTINUE
C
C
C         ASSEMBLE  L01,00 AND L00,01
C
      DO 1400 J = 1, KK
         DO 1400 I = 1, KK
            XL(2*KK+I, J) = ZERO
            XL(I, 2*KK+J) = ZERO
 1400 CONTINUE
C
C         ASSEMBLE diagonal and off-diagonal elements of L10,10
C
      PFAC = 16.0D0 * T / (25.0 * P)
      PIFAC = 5.0 / (3.0*PI)
      DO 1600 J = 1, KK
        WTJ_TMP = WT(J)
        CROT_J  = CROTOR(J) / XI(J)
        PFAC_J  = PFAC * XX(J) * WTJ_TMP
        SUM     = ZERO
        DO 1550 I = 1, KK

          FAC_1 = XX(I) / ((WT(I) + WTJ_TMP)**2 * BINDIF(I,J))

          FAC_2 = 4.0*ASTAR(I,J)*
     $              (ONE + PIFAC*(CROTOR(I)/XI(I) + CROT_J))

          XL(I+KK, J+KK) = PFAC_J * WT(I) * FAC_1
     $                    * ( 13.75 - 3.0*BSTAR(I,J) - FAC_2 )

          SUM = SUM + FAC_1 
     $              * (   7.5*WTJ_TMP**2
     $                  + WT(I)*(  WT(I)*(6.25 - 3.0*BSTAR(J,I)) 
     $                           + WTJ_TMP * FAC_2 )
     $                )
 1550   CONTINUE
        XL(J+KK, J+KK) = XL(J+KK, J+KK) - PFAC*XX(J)*SUM
 1600 CONTINUE

C
C         ASSEMBLE L10,01 AND L01,10, both the off-diagonal entries
C         and the on-diagonal entries.
C
      NN = 2*KK
      PFAC = 32.0 * T / (5.0 * PI * P)
      DO 1850 J = 1, KK
         IF (LIN(J) .NE. 0) THEN
            NN = NN + 1
            SUM = ZERO
            WTJ_TMP = WT(J)
            PFAC_J =   ( PFAC * WTJ_TMP * XX(J) * CROTOR(J) )
     $               / ( CINTOR(J) * XI(J) )
            DO 1800 I = 1, KK 
C                             The L10,01 term:     
              XL(I+KK, NN) = ( PFAC_J * ASTAR(J,I) * XX(I)     )
     $                     / ( (WTJ_TMP + WT(I)) * BINDIF(J,I) )
C                             The L01,10 term:
              XL(NN, I+KK) =  XL(I+KK, NN)
C                             The extra term that get's stuck
C                             on the diagonal:
              SUM    = SUM +  XL(I+KK, NN)
 1800       CONTINUE
C
C           Extra diagonal entries:
C               (These use the viscosity, eq. 49, in their formulation,
C                because the self-diffusion coefficient has been
C                reevaluated to be consistent with the viscosity.)      
C
            XL(J+KK, NN) = XL(J+KK, NN) + SUM
            XL(NN, J+KK) = XL(NN, J+KK) + SUM

         ENDIF
 1850 CONTINUE
C
C        ASSEMBLE L01,01, USING VISCOSITY EQ. (49)
C
      DO 2000 J = 1, KK
         DO 2000 I = 1, KK
            XL(2*KK+I, 2*KK+J) = ZERO
 2000 CONTINUE
C
      NN = 2*KK
      PFAC  = 4.0 * T / P
      PIFAC = 12.0 / (5.0 * PI)
      PIRU  = - 8.0 / (PI * RU)
      DO 2200 I = 1, KK
         IF (LIN(I) .NE. 0) THEN
            NN = NN + 1
            SUM = ZERO
            FAC_1 = ( PIFAC * WT(I) * CROTOR(I) )
     $            / ( CINTOR(I) * XI(I)         )
            DO 2100 K = 1, KK
               FAC_2 = XX(K) / BINDIF(I,K)
               SUM = SUM + FAC_2
               IF (I .NE. K) THEN
                  SUM = SUM + (FAC_1 * FAC_2 * ASTAR(I,K)) / WT(K)
               ENDIF
 2100       CONTINUE
            FAC_2 = XX(I) / CINTOR(I)
            XL(NN, NN) =     ( PIRU * WT(I) * FAC_2 * CROTOR(I) )
     1                     / ( VIS(I) * XI(I)                   )
     2                   - PFAC * SUM
            XL(NN, NN) = FAC_2 * XL(NN, NN)
         ENDIF
 2200 CONTINUE
C
C          ASSEMBLE THE RIGHT HAND SIDE FOR SOLVING EQ. (24)
C
      NN = 2*KK
      DO 3300 I = 1, KK
         R(I)    = ZERO
         R(I+KK) = XX(I)
         IF (LIN(I) .NE. 0) THEN
            NN  = NN + 1
            R(NN) = XX(I)
         ENDIF
 3300 CONTINUE
C
C          FACTOR AND SOLVE EQ. (24)
C
*
C*****precision > double - lapack
      CALL DGETRF (NN, NN, XL, KK3, IPVT, INFO)
      IF (INFO .NE. 0) THEN
          WRITE (6,*) ' ERROR IN DGETRF, INFO = ', INFO
          STOP
      ENDIF
      CALL DGETRS('N', NN, 1, XL, KK3, IPVT, R, NN, INFO)
      IF (INFO .NE. 0) THEN
          WRITE (6,*) ' ERROR IN DGETRS, INFO = ', INFO
          STOP
      ENDIF
C*****END precision > double - lapack
*
C*****precision > double - linpack
C      CALL DGEFA (XL, KK3, NN, IPVT, INFO)
C      IF (INFO .NE. 0) THEN
C          WRITE (6,*) ' ERROR IN DGEFA, INFO = ', INFO
C          STOP
C      ENDIF
C      CALL DGESL (XL, KK3, NN, IPVT, R, 0)
C*****END precision > double - linpack
*
C*****precision > single - lapack
C      CALL SGETRF (NN, NN, XL, KK3, IPVT, INFO)
C      IF (INFO .NE. 0) THEN
C          WRITE (6,*) ' ERROR IN SGETRF, INFO = ', INFO
C          STOP
C      ENDIF
C      CALL SGETRS('N', NN, 1, XL, KK3, IPVT, R, NN, INFO)
C      IF (INFO .NE. 0) THEN
C          WRITE (6,*) ' ERROR IN SGETRS, INFO = ', INFO
C          STOP
C      ENDIF
C*****END precision > single - lapack
*
C*****precision > single - linpack
C      CALL SGEFA (XL, KK3, NN, IPVT, INFO)
C      IF (INFO .NE. 0) THEN
C          WRITE (6,*) ' ERROR IN SGEFA, INFO = ', INFO
C          STOP
C      ENDIF
C      CALL SGESL (XL, KK3, NN, IPVT, R, 0)
C*****END precision > single - linpack
C
C
C          FORM THERMAL DIFFUSION COEFFICIENTS
C
      PFAC = 1.6 / RU
      DO 4000 K = 1, KK
         DT(K)  = PFAC * WT(K) * XX(K) * R(K)
 4000 CONTINUE
C
C          FORM THE THERMAL CONDUCTIVITY
C
      CONDTR = ZERO
      DO 4100 K = 1, KK
         CONDTR = CONDTR + X(K) * R(KK+K)
 4100 CONTINUE
C
      NN = 2*KK
      CONDIN = ZERO
      DO 4200 K = 1, KK
         IF (LIN(K) .NE. 0) THEN
            NN = NN + 1
            CONDIN = CONDIN + X(K) * R(NN)
         ENDIF
 4200 CONTINUE
C
      COND = -4.0 * (CONDTR + CONDIN)
C
      RETURN
      END
C
C----------------------------------------------------------------------
C
      SUBROUTINE MCATDR (T, X, IMCWRK, RMCWRK, TDR)
C
C*****precision > double
      IMPLICIT DOUBLE PRECISION (A-H, O-Z), INTEGER (I-N)
C*****END precision > double
C
C*****precision > single
C      IMPLICIT REAL (A-H, O-Z), INTEGER (I-N)
C*****END precision > single
C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C
C  SUBROUTINE MCATDR (T, X, IMCWRK, RMCWRK, TDR)
C
C    THIS SUBROUTINE COMPUTES THE THERMAL DIFFUSION RATIOS FOR THE
C    LIGHT SPECIES INTO THE MIXTURE.
C
C  INPUT-
C    T       - TEMPERATURE
C                  CGS UNITS - K.
C    X       - ARRAY OF MOLE FRACTIONS OF THE MIXTURE.
C                  DIMENSION X(*) AT LEAST KK.
C
C  WORK-
C    IMCWRK  - ARRAY OF INTEGER STORAGE AND WORK SPACE.  THE STARTING
C              ADDRESSES FOR THE IMCWRK SPACE ARE STORED IN
C              COMMON /MCMCMC/.
C                  DIMENSION IMCWRK(*) AT LEAST LENIMC.
C    RMCWRK  - ARRAY OF FLOATING POINT STORAGE AND WORK SPACE.  THE
C              STARTING ADDRESSES FOR THE RMCWRK SPACE ARE STORED IN
C              COMMON /MCMCMC/.
C                  DIMENSION RMCWRK(*) AT LEAST LENRMC.
C
C  OUTPUT-
C    TDR     - ARRAY OF THERMAL DIFFUSION RATIOS FOR THE KK SPECIES.
C              TDR(K) = 0 FOR ANY SPECIES WITH MOLECULAR WEIGHT LESS
C              THAN 5.
C                  DIMENSION TDR(*) AT LEAST KK.
C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C
      DIMENSION X(*), IMCWRK(*), RMCWRK(*), TDR(*)
C
      COMMON /MCMCMC/ RU, PATMOS, SMALL, NKK, NO, NLITE, INLIN, IKTDIF,
     1                IPVT, NWT, NEPS, NSIG, NDIP, NPOL, NZROT, NLAM,
     2                NETA, NDIF, NTDIF, NXX, NVIS, NXI, NCP, NCROT,
     3                NCINT, NBIND, NEOK, NSGM, NAST, NBST,
     4                NCST, NXL, NR, NWRK, K3
C
C       IN THIS SUBROUTINE, TEMPORARY STORAGE IS ASSIGNED AS:
C         A VECTOR OF THE "FITTED" PARTS OF THE THERMAL DIFUSSION
C           RATIOS ARE STORED IN RMCWRK(NXI).  SPECIFICALLY, THE
C           VECTOR REPRESENTS THE J COMPONENTS OF TDR(J,K), WHERE
C           K IS THE LIGHT SPECIES.
C
      ZERO = 0.0
C
      DO 100 K = 1, NKK
         TDR(K) = ZERO
  100 CONTINUE
C
      DO 500 L = 1, NLITE
         K = IMCWRK(IKTDIF+L-1)
         IF (K.eq.0) THEN
            WRITE(*,*) 'IKTDIF', IKTDIF
            WRITE(*,*) 'L', L
            WRITE(*,*) 'IKTDIF+L-1', IKTDIF+L-1
            WRITE(*,*) 'K', K
            call exit(1)
         END IF
C
         ISTRT = NTDIF + (L-1)*NO*NKK
         CALL MCEVAL (T, NKK, NO, RMCWRK(ISTRT), RMCWRK(NXI))
C
         DO 500 J = 1, NKK
            TDR(K) = TDR(K) + RMCWRK(NXI+J-1)*X(K)*X(J)
  500 CONTINUE
C
      RETURN
      END


