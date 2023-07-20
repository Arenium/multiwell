      SUBROUTINE quantrot(T,B,d,sig,Qroq,C,S,H)
 
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
 
c
c      Rotational partition function for quantum rotor
c
c         INPUT
c      T      = temperature
c      B      = rotational constant (cm-1)
c      d      = rotor dimension
c      sig      = rotor symmetry (foldedness)
C
C         OUTPUT
c      Qroq   = partition function
c      C      = constant p heat capacity (units of R)
c      S      = entropy (units of R)
c      H      = [H(T)-H(0)]/RT
c
 
      IMPLICIT NONE
      REAL(8) Q0 , Q1 , Q2 , T , B , C , S , H , Ered , E , 
     &                 hor , F , Fac , Qroq
      INTEGER(4) d , sig , j
      PARAMETER (hor = 1.4387752D+00)
      SAVE 
 
 
c     At j = 0:
c
      Q0 = 1.0D+00
      Q1 = 0.0D+00
      Q2 = 0.0D+00
      Fac = 1.0D+00
c
c     Now start loop

      j = 1
 
      DO WHILE ( Fac.GE.1.0D-12 )
         IF ( d.EQ.1 ) THEN
            F = 2.d+00
            E = B*j*j
         ELSEIF ( d.EQ.2 ) THEN
            F = 2.d+00*j + 1.d+00
            E = B*j*(j+1.d+00)
         ELSEIF ( d.EQ.3 ) THEN
            F = (2.d+00*j+1.d+00)*(2.d+00*j+1.d+00)
            E = B*j*(j+1.d+00)
         ENDIF
         Ered = hor*E/T
          IF ( Ered .LT. 20. ) THEN
           Fac = F*EXP(-Ered)
          ELSE
           Fac = 0.0
         ENDIF
         Q0 = Q0 + Fac
         Q1 = Q1 + Ered*Fac
         Q2 = Q2 + Ered*Ered*Fac
         j = j + 1
      ENDDO
 
c      Thermo formulae from Pitzer and Gwinn, J. Chem. Phys. 10, 428 (1942)
 
      Q0 = Q0/sig
      Q1 = Q1/sig
      Q2 = Q2/sig
      C = Q2/Q0 - (Q1/Q0)**2
      H = Q1/Q0
      S = (H+log(Q0))
 
      Qroq = Q0
 
      RETURN
 
      END SUBROUTINE
