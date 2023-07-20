      MODULE SEPSUBS
!
!    Module for sharing data between main program DENSUM, and subroutines 
!    used for sums and densities of separable degrees of freedom.
!
c
c A code for sums and densities of states calculations.
c Copyright (C) 2017 John R. Barker
c
c John R. Barker
c University of Michigan
c Ann Arbor, MI 48109-2143
c (734) 763 6239
c jrbarker@umich.edu
c
c This program is free software; you can redistribute it and/or
c modify it under the terms of the GNU General Public License (version 2)
c as published by the Free Software Foundation.
c
c This program is distributed in the hope that it will be useful,
c but WITHOUT ANY WARRANTY; without even the implied warranty of
c MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
c GNU General Public License for more details:
c
c     Free Software Foundation, Inc.
c     59 Temple Place - Suite 330
c     Boston, MA 02111-1307, USA.
c

      IMPLICIT NONE
      SAVE
      
      INTEGER, PARAMETER :: lo = 12               ! levels-file unit number
 
      INTEGER, PARAMETER :: NMAX = 50001          ! maximum no. of energy grains
      INTEGER, PARAMETER :: MaxDof = 500          ! maximum no. of degrees of freedom entries
      INTEGER, PARAMETER :: Nchr = 20             ! maximum no. of coefficients for hindered rotor Fourier series
      INTEGER, PARAMETER :: Levels = 20001        ! maximum no. of computed energy levels 
 
      CHARACTER(3), DIMENSION(MaxDof) :: IDOF(MaxDof)
      CHARACTER(3), PARAMETER :: VIB = 'VIB'      ! d.o.f. type vibration
      CHARACTER(3), PARAMETER :: ROT = 'ROT'      ! d.o.f. type classical rotation
      CHARACTER(3), PARAMETER :: QRO = 'QRO'      ! d.o.f. type quantum rotation
      CHARACTER(3), PARAMETER :: HRA = 'HRA'      ! d.o.f. type symmetric hindered rotor
      CHARACTER(3), PARAMETER :: HRB = 'HRB'      ! d.o.f. type symmetric hindered rotor
      CHARACTER(3), PARAMETER :: HRC = 'HRC'      ! d.o.f. type symmetric hindered rotor
      CHARACTER(3), PARAMETER :: HRD = 'HRD'      ! d.o.f. type general hindered rotor
      CHARACTER(3), PARAMETER :: TRN = 'TRN'      ! d.o.f. type translation
      CHARACTER(3), PARAMETER :: KRO = 'KRO'      ! d.o.f. type K-rotor (special uses, only)
      CHARACTER(3), PARAMETER :: TOP = 'TOP'      ! d.o.f. type symmetric top
      CHARACTER(3), PARAMETER :: BOX = 'BOX'      ! d.o.f. type particle-in-a-box
 
c      REAL(8) :: EJmin                            ! K-rotor zero of energy = VZPE - EJmin
      REAL(8) :: HRZPE                            ! total zpe of all hindered rotors

      INTEGER, DIMENSION(MaxDof) :: MODE
      CHARACTER(5), DIMENSION(MaxDof) :: Vhr
      CHARACTER(5), DIMENSION(MaxDof) :: Bhr 
 
      REAL(8), DIMENSION(MaxDof) :: MOD
      REAL(8), DIMENSION(MaxDof) :: WE
      REAL(8), DIMENSION(MaxDof) :: ANH
      INTEGER, DIMENSION(MaxDof) :: NG

      REAL(8), DIMENSION(MaxDof) :: B             ! rotational constant (cm-1)
      REAL(8), DIMENSION(MaxDof) :: VV
      REAL(8), DIMENSION(MaxDof) :: W
      REAL(8), DIMENSION(MaxDof) :: AMOM
      REAL(8), DIMENSION(MaxDof) :: zzpe
      REAL(8) :: ZPPE
      REAL(8) :: B1
      REAL(8) :: B2
      INTEGER, DIMENSION(MaxDof) :: NSIG
      REAL(8), DIMENSION(MaxDof) :: CV(MaxDof,Nchr)
      REAL(8), DIMENSION(MaxDof) :: CB(MaxDof,Nchr)
      REAL(8), DIMENSION(MaxDof) :: CVt(Nchr)
      REAL(8) :: CBt(Nchr)
      REAL(8), DIMENSION(MaxDof) :: Phav(MaxDof)
      REAL(8), DIMENSION(MaxDof) :: Phab(MaxDof)
      INTEGER, DIMENSION(MaxDof) :: NCV(MaxDof)
      INTEGER, DIMENSION(MaxDof) :: NCB(MaxDof)
      INTEGER, DIMENSION(MaxDof) :: NSV(MaxDof)
      INTEGER, DIMENSION(MaxDof) :: NSB(MaxDof)

!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
      CONTAINS
      
      SUBROUTINE MORLEV(T,DELE,JMAX,WE,XE,NG,zpe)
!
C	CALCULATES ENERGY LEVELS FOR MORSE OSCILLATOR TO BE USED
C          WITH STERAB
C
c	DELE	= energy grain
c	JMAX	= index of Emax (top of energy space)
c	WE	= harmonic vibration frequency (cm-1)
c	XE	= anharmonicity (cm-1)
c         NG          = degeneracy of vibration
c	IR	= vector of energy level indices
c	IMAX	= ceiling of array (e.g. dissociation energy for Morse Osc.)
c         zpe         = zero point energy (cm-1), including anharmonicity
c
C	Input OBSERVED WEo = WE + 2*XE for 0-1 transition
c	Note sign of XE:  negative for usual Morse Oscillator
c
c	E = WE*(v+1/2) + XE*(v+1/2)^2
C
C     ALL ENERGIES AND FREQUENCIES IN WAVENUMBERS, relative to ZPE = WE/2 + XE/4
C
c    2/02 Modification: now, hindered rotors are counted as vibs in Whitten-Rabinovitch
c         parameter calculations. Zero point energy is calculated based on anharmonic
c         vibrations and from the hindered rotor subroutine.
C
c    9/07 fixed small bug in calculation of RMAX 
c
      IMPLICIT NONE
      INTEGER, INTENT(IN) :: JMAX
      REAL(8), DIMENSION(JMAX), INTENT(INOUT) :: T
      REAL(8), INTENT(IN) :: DELE
      REAL(8), INTENT(IN) :: WE
      REAL(8), INTENT(IN) :: XE
      INTEGER, INTENT(IN) :: NG
      REAL(8), INTENT(OUT) :: zpe

      REAL(8), DIMENSION(JMAX) :: AT
      REAL(8) :: w0
      REAL(8) :: RMAX
      REAL(8) :: vi
      REAL(8) :: R
      INTEGER :: I
      INTEGER :: IMAX
      INTEGER, DIMENSION(Levels) :: IR
      INTEGER :: J
      INTEGER :: K
      INTEGER :: KARG
      INTEGER :: L
 
      DO I = 1, JMAX
        AT(I) = 0.0d+00
      END DO

      zpe = 0.5d+00*WE + 0.25d+00*XE                    ! Zero point energy at v=0
      write(lo,*) '  '
      write(lo,*) "Morse (or Harmonic) oscillator (vib), zpe = ", zpe
      write(lo,*) "          I     E-zpe                        IR(I)"

      IF ( XE .LT. 0.0 ) THEN
         w0 = WE + XE
         RMAX = -0.5d+00*w0/XE                          ! Highest bound state, E relative to ZPE
      ELSE
         RMAX = JMAX*DELE/WE                            ! Highest bound state
      ENDIF
      IMAX = INT(RMAX) 
      IF ( IMAX.GT.Levels ) IMAX = Levels

      DO I = 1 , IMAX                               ! Start at v=1
        vi = I + 0.5d+00
        R = (WE + XE*vi)*vi - zpe                   ! state energy relative to zpe
        IR(I) = 1 + NINT( R/DELE )                  ! Nearest integer number of grains
c        IR(I) = 1 + CEILING( R/DELE )		       ! CEILING integer number of grains
        write(lo,*) I, R, IR(I)
      END DO ! I

      DO L = 1 , NG                                         ! Loop for degenerate vibrations
         DO J = 1 , IMAX                                    ! IMAX is the number of energy states
            DO K = IR(J) , JMAX                             ! JMAX is the number of energy grains
               KARG = K - IR(J) + 1
               IF ( KARG.LE.JMAX ) AT(K) = AT(K) + T(KARG)  ! AT(K) is the number of states in the Kth grain
            END DO  ! K
         END DO  ! J
 
         DO J = 1 , JMAX
            T(J) = AT(J) + T(J)
            AT(J) = 0.0D+00
         END DO
      END DO   ! L

      write(lo,*) '---------------------------------------------------'
 
      RETURN
      END SUBROUTINE MORLEV
      
!-----------------------------------------------------------------------
      SUBROUTINE CROTLEV(T,DELE,JMAX,B,NDIM,NSYMM)
c
C      CALCULATES ENERGY LEVELS FOR **CLASSICAL** FREE ROTOR
c	  from W. Forst, "Unimolecular Reactions. A concise Introduction", 
c           Cambridge, 2003, Appendix I, p. 277ff
C
c      DELE      = energy grain size
c      JMAX      = index of Emax
c      B         = rotational constant (cm-1)
c      NDIM      = rotor dimension
c      NMAX      = ceiling of array
c
c      Eq. numbers refer to Forst 2003
c
      IMPLICIT NONE
      INTEGER, INTENT(IN) :: JMAX                      ! number of energy grains corresponding to Emax2
      REAL(8), DIMENSION(JMAX), INTENT(INOUT) :: T     ! array used by Beyer-Swinehart algorithm
      REAL(8), INTENT(IN) :: DELE                      ! energy grain size (units of cm-1)
      REAL(8), INTENT(IN) :: B                         ! rotational constant (cm-1)
      INTEGER, INTENT(IN) :: NDIM                      ! rotor dimension
      INTEGER, INTENT(IN) :: NSYMM                     ! rotor symmetry number

      REAL(8), DIMENSION(JMAX) :: AT                   ! array used by Beyer-Swinehart algorithm
      REAL(8), DIMENSION(JMAX) :: x
      REAL(8) :: qr
      REAL(8) :: R2
      REAL(8) :: E
      REAL(8) :: SE
      REAL(8) :: SELAST
      INTEGER :: I
      INTEGER :: J
      INTEGER :: K
      REAL(8) :: FAC
      REAL(8) :: PI

      PI = 4.0d+00*DATAN(1.0d+00)

       IF ( NDIM .EQ. 1 ) THEN
         qr =  SQRT( PI / B ) / NSYMM				! Eq. A1.12
       ELSEIF ( NDIM .EQ. 2 ) THEN
         qr = 1.d+00 / (NSYMM * B)					! Eq. A1.10
       ELSE
         WRITE(*,*) '***FATAL ERROR***: execution terminated'
         WRITE(*,*) 'Classical rotor dimension can be 1 or 2, only.'
         WRITE(*,*) 'Classical spherical tops must be computed using'
         WRITE(*,*) 'two separable rotors: 1-D and 2-D with same B.'
         STOP
       ENDIF
      
      FAC = qr / GAM( 1.d+00 + NDIM/2.d+00 )			! Eq. A1.13
      
      x(1) = 1.0D+00
      SELAST = 0.0d+00
      AT(1) = 0.0d+00
      DO  K = 1 , JMAX                                 ! Convolution using Beyer-Swinehart
         E = (K - 1)*DELE + 0.5*DELE
         SE = SQRT( E**NDIM)
         x(K) = FAC*( SE - SELAST )   ! number of classical free rotor states in grain
         write(*,*) SE, SELAST
         SELAST = SE
         AT(K) = 0.0d+00
         DO I = 1 , K
             AT(K) = AT(K) + x(K)*T(I)
         END DO         
      END DO
 
      DO K = 1 , JMAX
        T(K) = AT(K)
        AT(K) = 0.0D+00
      END DO

      RETURN
      END

!-----------------------------------------------------------------------
      SUBROUTINE ROTLEV(T,DELE,JMAX,B,NDIM,NSYMM)
!
!	Calculates energy levels for QUANTUM free rotations and prepares  
!    the T(i) array for computing sums and densities of states via the
!    Beyer-Swinehart algorithm.
!
!    Called by SUBROUTINE STERAB, which is called by PROGRAM DENSUM
!
!    Revised to modernized Fortran by John R. Barker (July 26, 2017)
!
      IMPLICIT NONE
      INTEGER, INTENT(IN) :: JMAX                      ! number of energy grains corresponding to Emax2
      REAL(8), DIMENSION(JMAX), INTENT(INOUT) :: T     ! array used by Beyer-Swinehart algorithm
      REAL(8), INTENT(IN) :: DELE                      ! energy grain size (units of cm-1)
      REAL(8), INTENT(IN) :: B                         ! rotational constant (cm-1)
      INTEGER, INTENT(IN) :: NDIM                      ! rotor dimension
      INTEGER, INTENT(IN) :: NSYMM                     ! rotor symmetry number

      REAL(8), DIMENSION(JMAX) :: AT                   ! array used by Beyer-Swinehart algorithm
      INTEGER :: IR                                    ! index of energy level
      REAL(8) :: R
      REAL(8) :: RMAX
      INTEGER :: IMAX
      INTEGER :: I
      INTEGER :: J
      INTEGER :: K
      INTEGER :: KARG
      REAL(8) :: F
 
      DO I = 1 , JMAX
        AT(I) = 0.0d+00
      END DO

      RMAX = SQRT(JMAX*DELE/B)   ! Maximum level
      IMAX = INT(RMAX) + 2
      IF ( IMAX.GT.20001 ) IMAX = 20001
c      write(*,*) 'JMAX=',JMAX,'  DELE=',DELE,'  B=',B,RMAX,IMAX
      write(lo,*) '  '
      write(lo,*) "Quantum Free-Rotor (qro)"
      write(lo,*) '          J   E                               IR   F'
      write(lo,*) '          0   0.0'
 
      DO  J = 1 , IMAX               ! J = rot quantum number; IMAX is the number of energy states
         IF ( NDIM.EQ.1 )     THEN   ! 1-D rotor (e.g. single-axis internal rotor)
            R = B*J*J
            F = 2
         ELSEIF ( NDIM.EQ.2 ) THEN   ! 2-D rotor
            R = B*J*( J + 1 )
            F = 2*J + 1
         ELSEIF ( NDIM.EQ.3 ) THEN   ! 3-D rotor (shperical top)
            R = B*J*( J + 1 )
            F = ( 2*J + 1 )**2
         ELSE
            WRITE(*,*) 'FATAL ERROR: execution terminated'
            WRITE(*,*) 'quantum rotor dimension must be 1, 2, or 3.'
            STOP
         ENDIF
         IR = NINT( R/DELE ) + 1              ! Nearest integer: index of energy grain
c         IR = CEILING( R/DELE ) + 1           ! CEILING integer: index of energy grain

         write(lo,*) J , R , IR , F

         DO K = IR , JMAX                     ! Jmax = number of energy grains (i.e. corresponding to Emax2)
            KARG = K - IR + 1
            IF ( KARG.LE.JMAX) AT(K) = AT(K) + F*T(KARG)
         END DO  ! K
      END DO  ! J
 
      DO J = 1 , JMAX     ! over energy grains
         T(J) = AT(J) + T(J)
         T(J) = T(J) / NSYMM
         AT(J) = 0.0D+00
      END DO  !  J
 
      write(lo,*) '---------------------------------------------------'

      RETURN
      END SUBROUTINE ROTLEV
      
!-----------------------------------------------------------------------
      SUBROUTINE STOPLEV(T,DELE,JMAX,B2,B1,NSYMM)
!
!	Calculates energy levels for a symmetric top and prepares the 
!    T(i) array for computing sums and densities of states via the
!    Beyer-Swinehart algorithm.
!
!    Called by SUBROUTINE STERAB, which is called by PROGRAM DENSUM
!
!    Revised to modernized Fortran by John R. Barker (July 26, 2017)
!
      IMPLICIT NONE
      INTEGER, INTENT(IN) :: JMAX                      ! number of energy grains corresponding to Emax2
      REAL(8), DIMENSION(JMAX), INTENT(INOUT) :: T     ! Beyer-Swinehart array
      REAL(8), INTENT(IN) :: DELE                      ! energy grain size
      REAL(8), INTENT(IN) :: B2                        ! 2D J-rotor rotational constant (cm-1)
      REAL(8), INTENT(IN) :: B1                        ! 1D K-rotor rotational constant (cm-1)
      INTEGER, INTENT(IN) :: NSYMM                     ! rotor symmetry

      REAL(8), DIMENSION(JMAX) :: AT                   ! Beyer-Swinehart array
      REAL(8) :: R
      REAL(8) :: EMAX
      REAL(8) :: Ej                                    ! J-rotor energy
      REAL(8) :: Ek                                    ! K-rotor energy
      REAL(8) :: F, Fj                                 ! multiplicity factor 
      REAL(8):: Emin                                   ! lowest J,K rotational energy for given J
      INTEGER :: I
      INTEGER :: IR                                    ! bin number for an energy level
      INTEGER :: j
      INTEGER :: jj
      INTEGER :: jjmax
      INTEGER :: k
      INTEGER :: kk
      INTEGER :: KARG
 
      DO I = 1, JMAX     ! Initialize AT array
        AT(I) = 0.0d+00
      END DO

      EMAX = DELE*(JMAX-1)
      IF ( B2 .GT. B1) THEN
        jjmax = CEILING( (-B2 + SQRT(B2**2 + 4.0*B1*EMAX) )/(2.0*B1) )   ! OBLATE
      ELSE
        jjmax = CEILING( (-B2 + SQRT(B2**2 + 4.0*B2*EMAX) )/(2.0*B2) )   ! PROLATE AND SPHERICAL
      ENDIF
      
      write(lo,*) '  '
      write(lo,*) 'Symmetric Top (top) '
      write(lo,*) '*** Only the levels ≤ 5000 cm-1 are printed ***'
      write(lo,*) '  '
      write(lo,*) 'Jmax consistent with Emax: ', jjmax
      write(lo,*) '  '
      write(lo,*) '          ',
     & 'J           K         cm-1                degen              IR'
 
      jj = 0
      R = 0.0d+00
      DO J = 0 , jjmax                                    ! start at qantum number J=0
          Ej = B2*J*(J + 1.d+00)
          DO K = 0 , J                                    ! K quantum number
             Ek = ( B1 - B2 )*K*K
             R = ( Ej + Ek )
             IF ( R.GE.0.0 .AND. R.LE.EMAX) THEN
                jj = jj + 1                               ! = level number
                IR = NINT( R/DELE ) + 1                   ! Nearest integer: number of grains
                IF ( k .EQ. 0 ) THEN
                   F = ( 2.d+00*J + 1.d+00 )                 ! (2*J+1) degeneracy when K=0
                 ELSE
                   F = 2.d+00 * ( 2.d+00*J + 1.d+00 )        ! 2*(2*J+1) when K >0
                ENDIF
                IF ( R .LE. 5000. ) write(lo,*) J , K , R , NINT(F) , IR      ! write out the energy level and degeneracy for R ≥ 5000 cm-1                ENDIF
                DO kk = 1 , JMAX                          ! Jmax = number of energy grains (i.e. corresponding to Emax2)
                   KARG = kk - IR + 1.d+0
                   IF ( KARG.GT.0 .AND. KARG.LE.JMAX) THEN
                      AT(kk) = AT(kk) + F*T(KARG)
                   ENDIF
                END DO  ! kk over energy range
             ENDIF   ! if R in energy range
          ENDDO ! K
       ENDDO   ! J
                 
      DO I = 1 , JMAX     ! over energy grains
         T(I) = AT(I) + T(I)
         T(I) = T(I) / NSYMM
         AT(I) = 0.0D+00
      END DO  !  I
       
      write(lo,*) '---------------------------------------------------'

      RETURN
      END SUBROUTINE STOPLEV

!-----------------------------------------------------------------------
      FUNCTION GAM(YY) 
c
c	Gamma FUNCTION of YY
c
      IMPLICIT NONE
      REAL(8), INTENT(IN) :: YY
      REAL(8) :: GAM
      REAL(8) :: Z
      REAL(8) :: X
      REAL(8) :: X2
      REAL(8) :: X3
      REAL(8) :: X5
      REAL(8) :: X7
      REAL(8) :: PI
      REAL(8) :: TWOPI

      PI = 4.0d+00*DATAN(1.0d+00)
      TWOPI = 2.0d+00*PI
      X = YY
      Z = 0.0
      IF ( X <= 0.d+00 ) THEN
        WRITE(*,*) '*** FATAL: argument of GAM(X) is less than 0  ****'
        STOP
      ENDIF
      
      DO
        IF ( X  >= 1.5 ) EXIT
        Z = Z + LOG(X)
        X = X + 1.0d+00
      END DO

      X2 = X*X
      X3 = X2*X
      X5 = X2*X3
      X7 = X2*X5
      GAM = EXP( -Z + (X-0.5d+00)*LOG(X) - X + 0.5d+00*LOG(TWOPI)
     &   + 1.0d+00/(12.d+00*X) - 1.0d+00/(360.d+00*X3) 
     &   + 1.0d+00/(1260.d+00*X5) - 1.0d+00/(1680.d+00*X7) )

! Stirling's formula
!      GAM = EXP( -Z + (X-0.5d+00)*LOG(X) - X + 0.5d+00*LOG(TWOPI)
!     &    + LOG( 1.d+00 + 1.d+00/(12.d+00*X) + 1.d+00/(288.d+00*X**2)
!     &    -139.d+00/(51840.d+00*X**3) - 571.d+00/(2488320.d+00*X**4) ) )

      RETURN
      END
      

!-----------------------------------------------------------------------
      SUBROUTINE BOXLEV( T , DELE , JMAX , WE , NG , zpe )
!
!	Calculates energy levels for a particle-in-a-box and prepares the 
!    T(i) array for computing sums and densities of states via the
!    Beyer-Swinehart algorithm.
!
!	E = WE*n^2  for n=1,2,3,...
!
!    OUTPUT ENERGIES IN WAVENUMBERS, relative to ZPE = WE
!
!    Called by SUBROUTINE STERAB, which is called by PROGRAM DENSUM
!
!    Revised to modernized Fortran by John R. Barker (July 26, 2017)
!
      IMPLICIT NONE

!  Calling parameters and types
      REAL(8), INTENT(INOUT), DIMENSION(JMAX) :: T      ! Beyer-Swinehart array
      REAL(8), INTENT(IN) :: DELE                       ! Energy grain size
      INTEGER, INTENT(IN) :: JMAX                       ! No. of energy grains
      REAL(8), INTENT(IN) :: WE                         ! 0-1 energy difference
      INTEGER, INTENT(IN) :: NG                         ! degeneracy
      REAL(8), INTENT(OUT) :: zpe                       ! energy of lowest state (n=1)
      
! Local variables
      INTEGER ::  I , J , K , L                   ! do-loop index
      INTEGER :: KARG                             ! argument
      INTEGER :: IMAX                             ! maximum number of states
      REAL(8) :: v2                               ! real value of n*n
      INTEGER :: IR                               
      REAL(8) :: R                                ! Energy relative to zpe
      REAL(8) :: RMAX                             ! energy of highest energy grain
      REAL(8), DIMENSION(JMAX) :: AT              ! Beyer-Swinehart temporary array
      
      RMAX = SQRT(JMAX*DELE/WE)                           ! Maximum bound level
      IMAX = INT(RMAX) + 1
      IF ( IMAX.GT.Levels ) IMAX = Levels
      zpe = WE                                            ! Zero point energy
      write(lo,*) '   '
      write(lo,*) "Particle-in-a-box (box), zpe = ", zpe
      write(lo,*) "          I     E-zpe                        IR(I)"

      DO I = 1 , JMAX
        AT(I) = 0.0d+00
      END DO

      DO L = 1 , NG                                          ! Loop for degenerate vibrations
         DO I = 2 , IMAX                                     ! IMAX is the number of energy states; start at 2 because "zero" level at I = 1
           v2 = DBLE(I)*DBLE(I)                              
           R = (WE*v2 - zpe)                                 ! Energy relative to zpe = WE
!           IR = INT( R/DELE ) + 2                            ! Truncate: number of grains
           IR = NINT( R/DELE ) + 1                           ! Nearest integer: number of grains
!           IR = CEILING( R/DELE ) + 1                        ! CEILING integer: index of energy grain
           write(lo,*)  I, R  , IR 
           DO K = IR , JMAX                                  ! JMAX is the number of energy grains
              KARG = K - IR + 1
              IF ( KARG.LE.JMAX ) AT(K) = AT(K) + T(KARG)    ! AT(K) is the number of states in the Kth grain
           END DO  ! K
         END DO  ! J

         DO J = 1 , JMAX
            T(J) = AT(J) + T(J)
            AT(J) = 0.0D+00
         END DO
      END DO   ! L

      write(lo,*) '---------------------------------------------------'

      RETURN
      END SUBROUTINE BOXLEV

!-----------------------------------------------------------------------

      Subroutine chkdens(Egrain1,Imax1,Isize,DENS,Emax2,bestImax)
!
! Tests fluctuations in densities of states to determine whether energy parameters are suitable
!
      IMPLICIT NONE
      
      REAL(8), INTENT(IN) :: Egrain1
      INTEGER, INTENT(IN) :: Imax1
      INTEGER, INTENT(IN) :: Isize
      REAL(8), INTENT(IN), DIMENSION(50001) :: DENS
      REAL(8), INTENT(IN) :: Emax2
      INTEGER, DIMENSION(6), INTENT(OUT) :: bestImax
      INTEGER :: x
      INTEGER :: i
      REAL(8) :: err
      
      x=1
      DO i=1,6
       bestImax(i)=0
      ENDDO

      DO i=INT(Emax2/Egrain1), 2, -1                                    ! Isize= INT(Emax2/Egrain1)+1
       if(Dens(i-1).lt.1E-2) goto 100
       err=ABS((Dens(i)-Dens(i-1))/Dens(i-1))*100
!      write(*,*) i,err
       IF(i.eq.INT(Emax2/Egrain1)) then
         if(err.gt.1.and.err.lt.2) then
           WRITE(*,*) " WARNING: Density fluctuation at Emax2 > 1%."
           bestImax(1)=INT(Emax2/Egrain1)
           x=2
         elseif(err.gt.2.and.err.lt.3) then
           WRITE(*,*) " WARNING: Density fluctuation at Emax2 > 2%."
           bestImax(1)=INT(Emax2/Egrain1)
           bestImax(2)=INT(Emax2/Egrain1)
           x=3
         elseif(err.gt.3.and.err.lt.4) then
           WRITE(*,*) " WARNING: Density fluctuation at Emax2 > 4%."
           bestImax(1)=INT(Emax2/Egrain1)
           bestImax(2)=INT(Emax2/Egrain1)
           bestImax(3)=INT(Emax2/Egrain1)
           x=4
         elseif(err.gt.4.and.err.lt.5) then
           WRITE(*,*) " WARNING: Density fluctuation at Emax2 > 4%."
           bestImax(1)=INT(Emax2/Egrain1)
           bestImax(2)=INT(Emax2/Egrain1)
           bestImax(3)=INT(Emax2/Egrain1)
           bestImax(4)=INT(Emax2/Egrain1)
           x=5
         elseif(err.gt.5) then
           WRITE(*,*) " WARNING: Density fluctuation at Emax2 > 5%."
           bestImax(1)=INT(Emax2/Egrain1)
           bestImax(2)=INT(Emax2/Egrain1)
           bestImax(3)=INT(Emax2/Egrain1)
           bestImax(4)=INT(Emax2/Egrain1)
           bestImax(5)=INT(Emax2/Egrain1)
           return
         endif

       ENDIF

       IF(err.gt.x) then              ! Error > 1%-5%
        bestImax(x)=i
!     bestImax(1): error< 1%, bestImax(2): error< 2%, bestImax(3): error< 3%
!     bestImax(4): error< 4%, bestImax(5): error< 5%
        x=x+1
        if(x.eq.6) goto 100
       ENDIF
      

      ENDDO
c      WRITE(*,*) "Wrong definition of energy parameters"
      WRITE(*,*)"Cannot evaluate fluctuations using these E parameters"
      STOP

100   IF(Imax1.lt.bestImax(5)) THEN
       WRITE(*,'(A45,I6)')" WARNING: Miminum suggested value for Imax1:"
     & , bestImax(5)
      ENDIF
      return
      END SUBROUTINE chkdens

!-----------------------------------------------------------------------
      SUBROUTINE KROTLEV( T, DELE, JMAX, B2, B1, J )
!
!    Calculates energy levels for a symmetric top K-rotor for a 
!    specified value of the total angular momentum quantum number J. 
!
!    THIS VERSION FOR PROLATE SYMMETRIC TOPS, ONLY!!
!
!   It prepares the T(i) array for computing sums and densities of 
!   states via the Beyer-Swinehart algorithm.
!
!    Called by SUBROUTINE STERAB0
!
!    Revised to modernized Fortran by John R. Barker (July 26, 2017)
!    Revised by John R. Barker (June 2, 2022)
!
      IMPLICIT NONE
      INTEGER, INTENT(IN) :: JMAX                      ! number of energy grains corresponding to Emax
      REAL(8), DIMENSION(JMAX), INTENT(INOUT) :: T     ! Beyer-Swinehart array
      REAL(8), INTENT(IN) :: DELE                      ! energy grain size
      REAL(8), INTENT(IN) :: B2                        ! 2D rotor rotational constant (cm-1)
      REAL(8), INTENT(IN) :: B1                        ! 1D rotor rotational constant (cm-1)
      INTEGER, INTENT(IN) :: J                         ! specified quantum number for total angular momentum

      REAL(8), DIMENSION(JMAX) :: AT                   ! Beyer-Swinehart array
      REAL(8) :: R                                     ! energy of a level
      REAL(8) :: Ejmax                                 ! highest rotational energy for given J
      REAL(8) :: EMAX
      REAL(8) :: F, Fj                                 ! multiplicity factors 
      INTEGER :: I
      INTEGER :: IR                                    ! bin number for an energy level
      INTEGER :: jj
      INTEGER :: jjmax
      INTEGER :: K
      INTEGER :: kk
      INTEGER :: KARG

      IF ( J.EQ. 0 ) THEN              ! RETURN IF J = 0
        RETURN
      END IF

      EMAX = DELE*(JMAX - 1)                              ! Highest energy

      DO I = 1, JMAX     ! Initialize AT array
        AT(I) = 0.d+00
      END DO
 
      write(lo,9901)
9901  FORMAT( 11X, 'J',11X, 'K',10X, 'E/cm-1',
     &        15x, 'degen', 10x, 'IR' )

      jj = 0
      R = 0.0d+00
      DO K = 1 , J                                    ! K quantum number
         R = ( B1 - B2 )*K*K
         IF ( R .LE. EMAX) THEN
            jj = jj + 1                               ! = level number
            IR = NINT( R/DELE ) + 1                   ! Nearest integer: number of grains
            IF ( k .EQ. 0 ) THEN
               F = 1.d+00
             ELSE
               F = 2.d+00 
            ENDIF
            IF ( R .LE. 5000. ) write(lo,*) J , K , R , NINT(F) , IR      ! write out the energy level and degeneracy for R ≥ 5000 cm-1
            DO kk = 1 , JMAX                                              ! Jmax = number of energy grains (i.e. corresponding to Emax2)
               KARG = kk - IR + 1.d+0
               IF ( (KARG .GT. 0) .AND. (KARG .LE. JMAX) ) THEN
                  AT(kk) = AT(kk) + F*T(KARG)
               ENDIF
            END DO  ! kk over energy range
         ENDIF   ! if R in energy range
      ENDDO ! K
          
      Fj = 2.d+00*J + 1.d+00                        ! 2J+1 degeneracy
      DO I = 1 , JMAX     ! over energy grains
         T(I) = AT(I) + T(I)
         T(I) = Fj*T(I)                             ! 2J+1 degeneracy
      END DO  !  I

      write(lo,*) '---------------------------------------------------'

      RETURN
      END SUBROUTINE KROTLEV

!-----------------------------------------------------------------------
c%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
c
c Thermo: a code for thermochemical calculations.
c Copyright (C) 2009 John R. Barker
c
c John R. Barker
c jrbarker@umich.edu
c University of Michigan
c Ann Arbor, MI 48109-2143
c (734) 763 6239
c
c This program is free software; you can redistribute it and/or
c modify it under the terms of the GNU General Public License (version 2)
c as published by the Free Software Foundation.
c
c This program is distributed in the hope that it will be useful,
c but WITHOUT ANY WARRANTY; without even the implied warranty of
c MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
c GNU General Public License for more details.
c
c See the 'ReadMe' file for a copy of the GNU General Public License,
c or contact:
c
c Free Software Foundation, Inc.
c 59 Temple Place - Suite 330
c Boston, MA 02111-1307, USA.
c
c%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

      SUBROUTINE rotunits( VROTIN , X , VROTOUT )
      
      IMPLICIT NONE
      CHARACTER(len=4), INTENT(IN) :: VROTIN      ! name of input units
      REAL(8), INTENT(INOUT) :: X                 ! numerical value
      CHARACTER(len=4), INTENT(IN) :: VROTOUT     ! name of output units
     
c      CONVERT UNITS FROM VROTIN TO VROTOUT
c
c       VROTIN or VROTOUT
c               = 'AMUA' for moment of inertia in units of amu*Ang^2
c       or      = 'GMCM' for moment of inertia units of gram*cm^2
c       or      = 'CM-1' for rotational constant in units of cm^-1
c       or      = 'MHZ' for rotational constant in units of MHz
c       or      = 'GHZ' for rotational constant in units of GHz
c         
      CALL ucase ( VROTIN )    ! convert to upper case
      CALL ucase ( VROTOUT )   ! convert to upper case
c
c         Convert VROTIN to *** AMUA ***
c
      IF ( VROTIN .EQ. 'AMUA' ) THEN
      ELSEIF ( VROTIN .EQ. 'GMCM' ) THEN
        X = X / 1.660538782d-040
      ELSEIF ( VROTIN .EQ. 'CM-1') THEN
        X = 16.85763D+00 / X
      ELSEIF ( VROTIN .EQ. 'MHZ' ) THEN
        X = 5.05379D+005 / X
      ELSEIF ( VROTIN .EQ. 'GHZ' ) THEN
        X = 5.05379D+002 / X
      ELSE
        write (*,*) 'FATAL: Rotation units (VROTIN) not recongized: ',
     &      VROTIN
        STOP
      ENDIF
c
c         Convert *** AMUA *** to VROTOUT and RETURN
c
      IF ( VROTOUT .EQ. 'AMUA' ) THEN
        RETURN
      ELSEIF ( VROTOUT .EQ. 'GMCM' ) THEN
        X = X * 1.660538782d-040
        RETURN
      ELSEIF ( VROTOUT .EQ. 'CM-1' ) THEN
        X = 16.85763D+00 / X
        RETURN
      ELSEIF ( VROTOUT .EQ. 'MHZ' ) THEN
        X = 5.05379D+005 / X
        RETURN
      ELSEIF ( VROTOUT .EQ. 'GHZ' ) THEN
        X = 5.05379D+002 / X
        RETURN
      ELSE
        write (*,*) 'FATAL: Rotation units (VROTOUT) not recongized: ',
     &      VROTOUT
        STOP
      ENDIF

      RETURN

      END SUBROUTINE rotunits

!-----------------------------------------------------------------------

      END MODULE SEPSUBS
      
      
