**==qvib.spg  processed by SPAG 5.11R  at 20:12 on 26 Mar 2001
      SUBROUTINE qboxmod(T,w,g,Qbox,C,S,H,zpe)
 
c%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
c
c Thermo: a code for thermochemical calculations.
c Copyright (C) 20162 John R. Barker
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
c      Last modified: 
c       10/2016
c 
c
c      Particle-in-a-Box partition function
c
c         INPUT
c      T      = temperature
c      g       = degeneracy
c      w      = vibrational frequency (cm-1)
c
c         OUTPUT
c      Qbox   = partition function
c      C      = heat capacity at constant pressure (units of R)
c      S      = entropy (units of R)
c      H      = [H�(T)-H�(0)]/RT
c      zpe    = zero point energy
c
        IMPLICIT NONE
       REAL(8) Q0 , Q1 , Q2 , T , w , C , S , H , Ered , E , 
     &                 hor , Fac , Qbox , zpe
       PARAMETER (hor = 1.4387752D+00)
       INTEGER(4) g , v
       SAVE 
 
      zpe = w
      Qbox = 0.0D+00
      Q0 = 0.0D+00
      Q1 = 0.0D+00
      Q2 = 0.0D+00
      Fac = 1.0D+00
      v = 1
 
      DO WHILE ( Fac.GE.1.0D-8 )
         E = w*v*v - zpe             ! Energy zero at zero pint energy
         IF ( Ered .LT. 20. ) THEN
           Fac = EXP(-Ered)
          ELSE
           Fac = 0.0
         ENDIF
         Fac = EXP(-Ered)
         Q0 = Q0 + Fac
         Q1 = Q1 + Ered*Fac
         Q2 = Q2 + Ered*Ered*Fac
         v = v + 1
      ENDDO
 
c      Thermo formulae from Pitzer and Gwinn, J. Chem. Phys. 10, 428 (1942)
 
      C = Q2/Q0 - (Q1/Q0)**2
      H = Q1/Q0
      S = (H+log(Q0))

      C = g*C
      H = g*H
      S = g*S
 
      Qbox = Q0**g
      zpe = g*zpe
 
      RETURN
 
      END SUBROUTINE
