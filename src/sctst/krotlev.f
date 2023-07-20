      SUBROUTINE KROTLEV(T,AT,DELE,JMAX,B,JK,NSYMM)
 
c%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
c
c DenSum: a code for calculating sums and densities of states.
c Copyright (C) 2007 John R. Barker
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
C      CALCULATES ENERGY LEVELS FOR 1-DIMENSIONAL K-ROTOR TO BE USED
C          WITH STERAB
C
c      DELE      = energy grain size
c      JMAX      = index of Emax
c      B      = rotational constant (cm-1)
c      JK      = quantum number for total angular momentum
c      IR      = vector of energy level indices
c      NMAX      = ceiling of array
c
c
      IMPLICIT REAL(8)(A-H,O-Z)
      DIMENSION IR(20001)
      DIMENSION T(JMAX), AT(JMAX)
      SAVE 


      Emax=(JMAX-1)*DELE 
      RMAX = SQRT(Emax/B)          ! Highest rotational quantum number consistent with energy ceiling
      NMAX = INT(RMAX)             ! INTEGER(4) version
      IF ( NMAX .GT. JK ) NMAX = JK
      IF ( NMAX .GT. 20001 ) NMAX = 20001
       
      IF ( NMAX .GT. 0 ) THEN
      DO 100 K = 1 , NMAX               ! K = quantum number of 1-dimensional free rotor
         R = B*K*K/DELE                 ! index number starting at 0
         IR(K) = NINT(R)              ! Nearest INTEGER(4): number of grains
 100  CONTINUE

CCCCCCCCCCCCCCCCCCCCCCC

            DO 281 J = 1 , NMAX ! J = quantum number for K-rotor;  NMAX is the number of energy levels
               IF ( J.EQ.0 ) THEN
                  F = 1.d+00         ! F is the multiplicity of the rotational level (see J.L. McHale, Molecular Spectroscopy (Prentice-Hall, 1999), 216f.
               ELSE
                  F = 2.d+00         ! 1-D rotor (e.g. single-axis internal rotor)
               ENDIF
               DO 261 K = 1 , JMAX
                  KARG = K - IR(J)
                  IF ( KARG.GT.0 .AND. KARG.LE.JMAX)
     &                   AT(K) = AT(K) + F*T(KARG)
 261           CONTINUE
 281        CONTINUE

            DO 301 J = 1 , JMAX
               T(J) = AT(J) + T(J)
               T(J) = T(J) / NSYMM
               AT(J) = 0.0D+00
 301        CONTINUE

CCCCCCCCCCCCCCCCCCCCCC

      ENDIF
 
      RETURN
      END
