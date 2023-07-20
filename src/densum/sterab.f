      SUBROUTINE STERAB( N ,JMAX ,DELE ,SUM, DENS )
C    CALCULATES COMPLEXIONS, DENSITIES, AND SUMS OF STATES ACCORDING TO
C     S.E. STEIN AND B.S. RABINOVITCH, J.CHEM.PHYS. 58, 2438 (1973)
C
C     RETURNS SUMS AND DENSITIES OF RO/VIB STATES IN
C           ARRAYS SUM AND DENS
C
C     DELE = GRAIN SIZE (CM-1)
C
C     ALL ENERGIES AND FREQUENCIES IN WAVENUMBERS
C
C	N	= degrees of freedom
C	JMAX	= index of energy bin corresponding to Emax
C	DELE	= energy grain (cm-1)
C
C
C	OTHER VARIABLES:
C	IDOF	= type of degree of freedom
C		= VIB (Morse or harmonic vibration)
C		= ROT (classical free rotation)
C		= QRO (quantized free rotation)
C                   = BOX (quantized particle-in-a-box)
C                   = HRx (hindered rotation; x=a,b,c; see below)
C                   = TRN (relative translation)
C	WE	= vibrational frequency (cm-1)
C         AMOM      = moment of inertia (amu Ang^2)
C	ANH	= vector of anharmonicities, or rotor symmetry number
C	NG	= vib. degeneracy, or rotor dimension
C	SUM	= sums of states
C	DENS	= densities of states
C         zzpe      = zero point energy (vibrations and hindered rotors)
C
C     Free Rotations:
C         AMOM(I) is moment of inertia (amu Ang^2)
C         ANH(I) is rotor foldedness (symmetry number)
C         NG(I) is rotor dimension (renamed NDIM)
C
C     1-D Hindered Rotations:
C         B(I) is rotational constant (cm-1)
C         W(I) small vibration frequency (cm-1)
C         NG(I) is internal rotor potential energy periodicity (foldedness)
C         NSIG(I) is internal rotor symmetry number
C         VV(I) is hindrance barrier height (cm-1)
C
C     Vibrations:
C         IF ANH=0.0,  HARMONIC OSCILLATORS ASSUMED
C         NG(I) IS OSCILLATOR DEGENERACY (OR ROTOR DIMENSION);
C         BOX:  NG(I) is oscillator degeneracy
C
C      JMAX IS NUMBER OF GRAINS AND CANNOT EXCEED 20001
C
C     THIS ROUTINE USES SUBROUTINES MORLEV(....)
C                                   BOXLEV(...)
C                                   GAM(X)
C                                   ROTLEV
C                                   HRLEVC
!
      USE SEPSUBS                   ! Declarations and subroutines for separable d.of.f.
      USE HRMOD                     ! separable 1D Hindered rotor subroutines
      IMPLICIT NONE

      INTEGER, INTENT(IN) :: N                           ! degrees of freedom
      INTEGER, INTENT(IN) :: JMAX                        ! index of energy bin corresponding to Emax
      REAL(8), INTENT(IN) :: DELE                        ! energy grain (cm-1)
      REAL(8), INTENT(OUT), DIMENSION(JMAX) :: SUM, DENS

      REAL(8), DIMENSION(JMAX) :: T                      ! Beyer-Swinehart array
      INTEGER, PARAMETER :: IMAX=501                     ! hindered rotor grid points and maximum number of eigenvalues
      INTEGER :: NTRN, NTRNFLAG, JK
      REAL(8) :: FACTR, atr, BB, BG, E, FAC
      INTEGER :: I, K, L
      INTEGER :: NDIM, NGG, NR, NSYMM
      REAL(8) :: R2, RG, VVV, zap
C
C          INITIALIZE T
C
      T(1) = 1.0D+00
      DO I = 2 , JMAX
         T(I) = 0.0D+00
      END DO
C
C          TROE METHOD FOR CLASSICAL ROTATIONS, EXTENDED TO TRANSLATIONS
C          [Astholz, Troe, and Wieters, J.Chem.Phys., 70, 5107, 1979]
C          [Equation numbers from Robinson & Holbrook, Unimolecular Reactions (Wiley, 1973)]
c
      NR = 0
      NTRN = 0
      NTRNFLAG = 0
      FAC = 1.D+00
      NSYMM = 1
      DO I = 1 , N                                        ! TEST FOR CLASS. ROTORS AND TRANSLATION
         IF ( IDOF(I).EQ.ROT ) THEN
            IF ( NG(I).LE.0 ) NG(I) = 1   ! free rotor dimension
            NDIM = NG(I)                  ! free rotor dimension
            BB = 16.85763D+00/AMOM(I)     ! rotational constant (cm-1)
            NR = NR + NDIM                ! number of classical free rotor degrees of freedom
            RG = DBLE(NDIM)/2.D+00        ! in Eq. 5.16
            BG = SQRT( BB**NDIM )         ! in Eq. 5.16
            FAC = FAC*GAM(RG)/(BG*ANH(I)) ! in Eq. 5.16; ANH(I) is free rotor foldedness (symmetry number)
            NSYMM = NSYMM*ANH(I)          ! symmetry number
         ENDIF
         IF (IDOF(I) .EQ. TRN) THEN
            NTRN = I
            NTRNFLAG = NTRNFLAG + 1       ! Flag for presence of relative translation
            IF (NTRNFLAG .GT. 1) THEN
              STOP '**** ONLY ONE 3-D TRANSLATION ALLOWED ***'
            ENDIF
         ENDIF
      END DO
 
       IF ( (NR.GT.0) .AND. (NTRN .EQ. 0) ) THEN              ! CLASS. ROTORS BUT NO TRANSLATION
         R2 = DBLE(NR)/2.D+00             ! in Eq. 5.15
         FAC = FAC/GAM( R2+1.D+00 )       ! in Eq. 5.15
         DO K = 2 , JMAX              ! Initialize T(K)
            E = (K-1)*DELE
            T(K) = FAC*( SQRT(E**NR) - SQRT( (E-DELE)**NR ) )           !  number of classical free rotor states in grain
         END DO
         T(1) = 0.0d+00
      ELSEIF ( (NTRN .GT. 0) .AND. (NR .EQ. 0) ) THEN         ! TRANSLATION BUT NO CLASS. ROTORS
         FACTR = 2.43980d+20*( AMOM(NTRN) )**1.5d+00                      ! states per cc
         DO K = 1 , JMAX               
            E = (K-1)*DELE
            T(K) = FACTR*( E**1.5d+00 - (E-DELE)**1.5d+00 )               ! number of classical translational states (per cc) in grain
         END DO
      ELSEIF ( (NTRN .GT. 0) .AND. (NR .GT. 0) ) THEN         ! TRANSLATION WITH CLASS. ROTORS
         atr = (NR + 3.0D+00) / 2.D+00
         FACTR = 3.24332d+20*( AMOM(NTRN) )**1.5d+00                      ! states per cc
         FACTR =  FACTR*FAC / ( atr * GAM(atr) )
         DO K = 2 , JMAX               
            E = (K-1)*DELE
            T(K) = FACTR*( E**atr - (E-DELE)**atr )                       ! number of classical rot-translational states (per cc) in grain
         END DO
         T(1) = 0.0d+00
      ENDIF
C
C
C
      DO I = 1 , N
         WRITE (*,99001) I , WE(I)
 
         IF ( IDOF(I).EQ.KRO ) THEN                             ! ********  1-Dim. free K-rotor *********
            B2 = 16.85763D+00/AMOM(I)   ! 2D rotational constant (cm-1)
            B1 = ANH(I)                 ! 1D rotational constant (cm-1)
            JK    = NG(I)               ! quantum number for total angular momentum
            CALL KROTLEV( T, DELE, JMAX, B2, B1, JK )      ! EJmin is the lowest energy K-level for specified JK

         ELSEIF ( IDOF(I).EQ.TOP ) THEN                         ! ********  Symmetric Top (special: for testing purposes only)  *********
            NSYMM = NG(I)              ! Symmetry number
            CALL STOPLEV(T,DELE,JMAX,B2,B1,NSYMM)

         ELSEIF ( IDOF(I).EQ.VIB ) THEN                         ! ********  Vibrations **********************
            CALL MORLEV(T,DELE,JMAX,WE(I),ANH(I),NG(I),zap)
            zzpe(I) = zap                                       ! Zero Point Energy
            ZPPE = ZPPE + zzpe(I)

         ELSEIF ( IDOF(I).EQ.BOX ) THEN                         ! ********  Particle-in-a-Box ***************
            CALL BOXLEV(T,DELE,JMAX,WE(I),NG(I),zap)
            zzpe(I) = zap                                       ! Zero Point Energy
            ZPPE = ZPPE + zzpe(I)
 
         ELSEIF ( IDOF(I).EQ.QRO ) THEN                         ! ********  Quantized free rotors *********
             BB = 16.85763D+00/AMOM(I)   ! rotational constant (cm-1)
            NDIM = NG(I)                ! Rotor dimension
            NSYMM = ANH(I)              ! Symmetry number
            CALL ROTLEV(T,DELE,JMAX,BB,NDIM,NSYMM)

         ELSEIF ( (IDOF(I).EQ.HRA) .OR. (IDOF(I).EQ.HRB) .OR. 
     &            (IDOF(I).EQ.HRC) ) THEN                       ! ******  Quantized symmetrical hindered rotors ********
            BB = B(I)           ! rotational constant (cm-1)
            VVV = VV(I)         ! hindrance barrier (cm-1)
            NGG = NG(I)         ! Potential energy symmetry (foldedness)
	       CALL SHRLEV(T,DELE,JMAX,BB,VVV,NGG,1,IMAX,zap,
     &          'VHRD1','BHRD1',0.0d0,0.0d0)		! rigid, symmetrical hindered rotor
            zzpe(I) = zap                                   ! Zero Point Energy
            ZPPE = ZPPE + zzpe(I)
            HRZPE = HRZPE + zap                             ! total zpe of all hindered rotors

         ELSEIF ( IDOF(I).EQ.HRD ) THEN                         ! ******  Quantized GENERAL hindered rotors ********
	       DO L=1, NCV(I)
	         CVt(L)=CV(I,L)
	       ENDDO
            DO L=1, NCB(I)
              CBt(L)=CB(I,L)
            ENDDO
	       NGG = NG(I)         ! Hinder rotation symmetry No. (foldedness)
	       CALL UHRLEV(T,DELE,JMAX,NCB(I),NCV(I),CBt,CVt,
     &	    NSV(I),NSB(I),IMAX,zap,Vhr(I),Bhr(I),Phav(I),Phab(I),NGG)	! non-rigid, unsymmetrical hindered rotor
	       zzpe(I) = zap                                   ! Zero Point Energy
            ZPPE = ZPPE + zzpe(I)
            HRZPE = HRZPE + zap                             ! total zpe of all hindered rotors
         ENDIF
         
      END DO   ! I  loop quantized degrees of freedom
c
c    Reminder: T(I) is the number of states in the Ith grain.
c
      SUM(1) = T(1)                                    ! Energy at E = 0 (Energy at TOP of grain)
      DENS(1) = T(1) / DELE
      
      DO I = 2 , JMAX
         SUM(I)  = T(I) + SUM(I-1)                     ! Sum of states (Energy at TOP of energy grain)
         DENS(I) = T(I)/DELE
      END DO  ! I

      RETURN
99001 FORMAT (10X,' Mode no.',I5,F10.2)
      END SUBROUTINE STERAB


