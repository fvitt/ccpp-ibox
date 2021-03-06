C-----------------------------------------------------------------------
      SUBROUTINE SPTGPSD(IROMB,MAXWV,KMAX,NPS,
     &                   KWSKIP,KGSKIP,NISKIP,NJSKIP,
     &                   TRUE,XMESH,ORIENT,WAVE,XN,YN,XS,YS)
C$$$  SUBPROGRAM DOCUMENTATION BLOCK
C
C SUBPROGRAM:  SPTGPSD    TRANSFORM SPECTRAL TO POLAR STEREO. GRADIENTS
C   PRGMMR: IREDELL       ORG: W/NMC23       DATE: 96-02-29
C
C ABSTRACT: THIS SUBPROGRAM PERFORMS A SPHERICAL TRANSFORM
C           FROM SPECTRAL COEFFICIENTS OF SCALAR FIELDS
C           TO GRADIENT FIELDS ON A PAIR OF POLAR STEREOGRAPHIC GRIDS.
C           THE WAVE-SPACE CAN BE EITHER TRIANGULAR OR RHOMBOIDAL.
C           THE WAVE AND GRID FIELDS MAY HAVE GENERAL INDEXING,
C           BUT EACH WAVE FIELD IS IN SEQUENTIAL 'IBM ORDER',
C           I.E. WITH ZONAL WAVENUMBER AS THE SLOWER INDEX.
C           THE TWO SQUARE POLAR STEREOGRAPHIC GRIDS ARE CENTERED
C           ON THE RESPECTIVE POLES, WITH THE ORIENTATION LONGITUDE
C           OF THE SOUTHERN HEMISPHERE GRID 180 DEGREES OPPOSITE
C           THAT OF THE NORTHERN HEMISPHERE GRID.
C           THE VECTORS ARE AUTOMATICALLY ROTATED TO BE RESOLVED
C           RELATIVE TO THE RESPECTIVE POLAR STEREOGRAPHIC GRIDS.
C
C           THE TRANSFORM IS MADE EFFICIENT             \ 4 | 5 /
C           BY COMBINING POINTS IN EIGHT SECTORS         \  |  /
C           OF EACH POLAR STEREOGRAPHIC GRID,           3 \ | / 6
C           NUMBERED AS IN THE DIAGRAM AT RIGHT.           \|/
C           THE POLE AND THE SECTOR BOUNDARIES          ----+----
C           ARE TREATED SPECIALLY IN THE CODE.             /|\
C           UNFORTUNATELY, THIS APPROACH INDUCES        2 / | \ 7
C           SOME HAIRY INDEXING AND CODE LOQUACITY,      /  |  \
C           FOR WHICH THE DEVELOPER APOLOGIZES.         / 1 | 8 \
C
C           THE TRANSFORMS ARE ALL MULTIPROCESSED OVER SECTOR POINTS.
C           TRANSFORM SEVERAL FIELDS AT A TIME TO IMPROVE VECTORIZATION.
C           SUBPROGRAM CAN BE CALLED FROM A MULTIPROCESSING ENVIRONMENT.
C
C PROGRAM HISTORY LOG:
C   96-02-29  IREDELL
C 1998-12-15  IREDELL  OPENMP DIRECTIVES INSERTED
C
C USAGE:    CALL SPTGPSD(IROMB,MAXWV,KMAX,NPS,
C    &                   KWSKIP,KGSKIP,NISKIP,NJSKIP,
C    &                   TRUE,XMESH,ORIENT,WAVE,XP,YP)
C   INPUT ARGUMENTS:
C     IROMB    - INTEGER SPECTRAL DOMAIN SHAPE
C                (0 FOR TRIANGULAR, 1 FOR RHOMBOIDAL)
C     MAXWV    - INTEGER SPECTRAL TRUNCATION
C     KMAX     - INTEGER NUMBER OF FIELDS TO TRANSFORM.
C     NPS      - INTEGER ODD ORDER OF THE POLAR STEREOGRAPHIC GRIDS
C     KWSKIP   - INTEGER SKIP NUMBER BETWEEN WAVE FIELDS
C                (DEFAULTS TO (MAXWV+1)*((IROMB+1)*MAXWV+2) IF KWSKIP=0)
C     KGSKIP   - INTEGER SKIP NUMBER BETWEEN GRID FIELDS
C                (DEFAULTS TO NPS*NPS IF KGSKIP=0)
C     NISKIP   - INTEGER SKIP NUMBER BETWEEN GRID I-POINTS
C                (DEFAULTS TO 1 IF NISKIP=0)
C     NJSKIP   - INTEGER SKIP NUMBER BETWEEN GRID J-POINTS
C                (DEFAULTS TO NPS IF NJSKIP=0)
C     TRUE     - REAL LATITUDE AT WHICH PS GRID IS TRUE (USUALLY 60.)
C     XMESH    - REAL GRID LENGTH AT TRUE LATITUDE (M)
C     ORIENT   - REAL LONGITUDE AT BOTTOM OF NORTHERN PS GRID
C                (SOUTHERN PS GRID WILL HAVE OPPOSITE ORIENTATION.)
C     WAVE     - REAL (*) WAVE FIELDS
C   OUTPUT ARGUMENTS:
C     XN       - REAL (*) NORTHERN POLAR STEREOGRAPHIC X-GRADIENTS
C     YN       - REAL (*) NORTHERN POLAR STEREOGRAPHIC Y-GRADIENTS
C     XS       - REAL (*) SOUTHERN POLAR STEREOGRAPHIC X-GRADIENTS
C     YS       - REAL (*) SOUTHERN POLAR STEREOGRAPHIC Y-GRADIENTS
C
C SUBPROGRAMS CALLED:
C   SPWGET       GET WAVE-SPACE CONSTANTS
C   SPLAPLAC     COMPUTE LAPLACIAN IN SPECTRAL SPACE
C   SPTGPSV      TRANSFORM SPECTRAL VECTOR TO POLAR STEREO.
C
C ATTRIBUTES:
C   LANGUAGE: FORTRAN 77
C
C$$$
      REAL WAVE(*),XN(*),YN(*),XS(*),YS(*)
      REAL EPS((MAXWV+1)*((IROMB+1)*MAXWV+2)/2),EPSTOP(MAXWV+1)
      REAL ENN1((MAXWV+1)*((IROMB+1)*MAXWV+2)/2)
      REAL ELONN1((MAXWV+1)*((IROMB+1)*MAXWV+2)/2)
      REAL EON((MAXWV+1)*((IROMB+1)*MAXWV+2)/2),EONTOP(MAXWV+1)
      REAL WD((MAXWV+1)*((IROMB+1)*MAXWV+2)/2*2+1,KMAX)
      REAL WZ((MAXWV+1)*((IROMB+1)*MAXWV+2)/2*2+1,KMAX)
C - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
C  CALCULATE PRELIMINARY CONSTANTS
      CALL SPWGET(IROMB,MAXWV,EPS,EPSTOP,ENN1,ELONN1,EON,EONTOP)
      MX=(MAXWV+1)*((IROMB+1)*MAXWV+2)/2
      MDIM=2*MX+1
      KW=KWSKIP
      IF(KW.EQ.0) KW=2*MX
C - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
C  CALCULATE GRADIENTS
C$OMP PARALLEL DO PRIVATE(KWS)
      DO K=1,KMAX
        KWS=(K-1)*KW
        CALL SPLAPLAC(IROMB,MAXWV,ENN1,WAVE(KWS+1),WD(1,K),1)
        WZ(1:2*MX,K)=0.
      ENDDO
      CALL SPTGPSV(IROMB,MAXWV,KMAX,NPS,MDIM,KGSKIP,NISKIP,NJSKIP,
     &             TRUE,XMESH,ORIENT,WD,WZ,XN,YN,XS,YS)
C - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
      END
