      SUBROUTINE BOXLEV(T,AT,DELE,JMAX,WE,NSYMM,zpe)
 
c%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
c
c DenSum: a code for calculating sums and densities of states.
c Copyright (C) 2002 John R. Barker
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
C
C
C      CALCULATES ENERGY LEVELS FOR PARTICLE-IN-A-BOX TO BE USED
C          WITH STERAB
C
c      DELE      = energy grain
c      JMAX      = index of Emax
c      WE      = vibration frequency (cm-1)
c      IR      = vector of energy level indices
c      NMAX      = ceiling of array
c       zpe     = zero point energy (cm-1)
c
c      E = WE*n^2  for n=1,2,3,...
C
C     OUTPUT ENERGIES IN WAVENUMBERS, relative to ZPE = WE
C
c    2/02 counted as vibs in Whitten-Rabinovitch
c         parameter calculations.
C
      IMPLICIT REAL(8)(A-H,O-Z)
      DIMENSION IR(20001)
      DIMENSION T(JMAX), AT(JMAX)

      SAVE 
 
      Emax=(JMAX-1)*DELE
      RMAX = SQRT(Emax/WE)                           ! Maximum bound level
      NMAX = INT(RMAX) 
      IF ( NMAX.GT.20001 ) NMAX = 20001
      zpe = WE                                            ! Zero point energy
      DO 100 I = 1 , NMAX
        v2 = DBLE(I+1)*DBLE(I+1)                          ! "zero" level = 1, start loop at I+1=2
        R = (WE*v2 - zpe)/DELE                            ! Relative to zpe = WE
        IR(I) = NINT(R)                                 ! Nearest integer: number of grains
 100  CONTINUE

            DO 390 L = 1 , NSYMM                                ! Loop for degenerate vibrations
               DO 355 J = 1 , NMAX                              ! NMAX is the number of energy levels
                  DO 365 K = 1 , JMAX
                     KARG = K - IR(J)
                     IF ( KARG.GT.0 .AND. KARG.LE.JMAX )
     &                   AT(K) = AT(K) + T(KARG)
 365              CONTINUE
 355           CONTINUE

               DO 370 J = 1 , JMAX
                  T(J) = AT(J) + T(J)
                  AT(J) = 0.0D+00
 370           CONTINUE
 390        CONTINUE


      RETURN
      END
