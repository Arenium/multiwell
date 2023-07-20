c%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
c
c Thermo: a code for thermochemical calculations.
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
      FUNCTION Qvib(T,w,x,g,C,S,H,zpe)
c
c      Last modified: 
c       11/02 Fixed yet another bug related to vibrational degeneracies
c 
c
c      Vibrational partition function
c
c      T      = temperature
c      w      = harmonic vibrational frequency (cm-1)
c      x      = anharmonicity
c      g      = degeneracy
c      C      = heat capacity at constant pressure (units of R)
c      S      = entropy (units of R)
c      H      = [H¡(T)-H¡(0)]/RT
c       zpe     = zero point energy
c
 
      IMPLICIT NONE
      REAL(8) Q0 , Q1 , Q2 , T , w , C , S , H , Ered , E , 
     &                 hor , Fac , Qvib , x , zpe
      PARAMETER (hor = 1.4387752D+00)
      INTEGER(4) g , v , vmax
      SAVE 
 
      zpe = w*0.5 + x*0.25
      vmax = 10000
      IF ( x .LT. 0.0d+00) vmax = (-w/(2.0d+00*x) - 0.5d+00)
      Q0 = 1.0D+00
      Q1 = 0.0D+00
      Q2 = 0.0D+00
      Fac = 1.0d+00
      v = 1
 
      DO WHILE ( Fac.GE.1.0D-12 .AND. v .LE. vmax )
         E = w*(v+0.5d+00) + x*(v+0.5d+00)*(v+0.5d+00) - zpe
         Ered = hor*E/T
         Fac = EXP(-Ered)
         Q0 = Q0 + Fac
         Q1 = Q1 + Ered*Fac
         Q2 = Q2 + Ered*Ered*Fac
         v = v + 1
      ENDDO
 
c   Standard Thermo formulae taken from Pitzer and Gwinn, J. Chem. Phys. 10, 428 (1942)
 
      C = Q2/Q0 - (Q1/Q0)**2
      H = Q1/Q0
      S = (H+log(Q0))

      C = g*C
      H = g*H
      S = g*S
 
      Qvib = Q0**g
      zpe = g*zpe
 
      RETURN
 
      END
