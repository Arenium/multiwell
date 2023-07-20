c%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
c
c DenSum: a code for calculating sums and densities of states.
c Copyright (C) 2001 John R. Barker
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
c See the "ReadMe" file for a copy of the GNU General Public License,
c or contact:
c
c Free Software Foundation, Inc.
c 59 Temple Place - Suite 330
c Boston, MA 02111-1307, USA.
c
c%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
C
C	CALCULATES ENERGY LEVELS FOR MORSE OSCILLATOR
C
c	DELE	= energy grain
c	JMAX	= index of Emax
c	WE	= harmonic vibration frequency (cm-1)
c	XE	= anharmonicity (cm-1)
c	IR	= vector of energy level indices
c	NMAX	= ceiling of array (e.g. dissociation energy for Morse Osc.)
c	zpe     = zero point energy (cm-1), including anharmonicity
c
C	Input OBSERVED WEo = WE + 2*XE for 0-1 transition
c	Note sign of XE:  negative for usual Morse Oscillator
c
c	E = WE*(v+1/2) + XE*(v+1/2)^2
C
C     ALL ENERGIES AND FREQUENCIES IN WAVENUMBERS, relative to ZPE = WE/2 + XE/4
C
c--------------------------------------------------------------------------------
      SUBROUTINE MORLEV(T,AT,DELE,JMAX,WE,XE,zpe)
 
      IMPLICIT REAL(8) (A-H,O-Z)
      EXTERNAL STERAB , GAM , WRAB , WRDEN , ALNGAM
      DIMENSION IR(20001)
      DIMENSION T(JMAX), AT(JMAX)
      SAVE 
 
      w0 = WE + XE
      IF ( XE.LT.0.0 ) THEN
         RMAX = -0.5*w0/XE                              ! Maximum bound level, E relative to ZPE
      ELSE
         RMAX = JMAX*DELE/WE                            ! Energy Ceiling level
      ENDIF
      NMAX = INT(RMAX) 

      zpe = 0.5d+00*WE + 0.25d+00*XE                    ! Zero point energy at v=0

      DO 100 I = 1 , NMAX                               ! Start at v=1
        vi = DBLE(I)
        R = (w0 + XE*vi)*vi/DELE
      IR(I) = IDNINT(R)                                 ! Nearest integer: number of grains
 100  CONTINUE
 
      DO J = 1 , NMAX                                   ! NMAX is the number of energy levels
        DO K = 1 , JMAX
          KARG = K - IR(J)
          IF ( KARG.GT.0 .AND. KARG.LE.JMAX ) AT(K) = AT(K) + T(KARG)
        END DO
      END DO

      DO J = 1 , JMAX
         T(J) = AT(J) + T(J)
         AT(J) = 0.0D+00
      END DO

      RETURN
      END
