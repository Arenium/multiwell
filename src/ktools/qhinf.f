c%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
c
c Qhin: a subroutine for hindered rotor partition fxn calculations.
c Copyright (C) 2009 Lam T. Nguyen and John R. Barker
c
c Lam T. Nguyen
c NGUYENLT@umich.edu
c University of Michigan
c Ann Arbor, MI 48109-2143
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
c
c      Quantum Hindered Rotor partition function 
c
c      T      = temperature
c      EV      = vector of eigenvalues (cm-1) with INMAX elements [from ghrlev.f]
c      INMAX= number of elements in vector EV 
c      nsig      = symmetry number of internal rotor
c      C      = constant p heat capacity (units of R)
c      SE      = entropy (units of R)
c      H      = [H(T)-H(0)]/RT
c    zpe  = zero point energy (cm-1)
c
c
c%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

      FUNCTION Qhin(T,EV,INMAX,nsig,C,SE,H,zpe)
 
      IMPLICIT NONE
      REAL(8) Qhin , T, C , SE , H, zpe
      REAL(8) Q0 , Q1 , Q2 , Ered , Fac, hor, Temp 
      INTEGER(4) INMAX 
      REAL(8) EV(INMAX)
      INTEGER(4)  L, nsig  
      PARAMETER (hor = 1.4387752D+00)
      SAVE 
 
      Temp = T

      Q0 = 0.0d+00
      Q1 = 0.0D+00
      Q2 = 0.0D+00
      Fac=1.0d0
      L=1

      zpe=EV(1)

      DO WHILE ((Fac.GE.1.0d-30).and.(L.LE.INMAX))
            Ered = hor*(EV(L)-zpe)/Temp
            Fac=exp(-Ered)
            Q0 = Q0 + Fac
            Q1 = Q1 + Ered*Fac
            Q2 = Q2 + Ered*Ered*Fac
            L=L+1
      ENDDO

      Q0 = Q0/nsig
      Q1 = Q1/nsig
      Q2 = Q2/nsig

c      Thermo formulae from Pitzer and Gwinn, J. Chem. Phys. 10, 428 (1942)
 
      Qhin = Q0
      C = Q2/Q0 - (Q1/Q0)**2
      H = Q1/Q0
      SE = (H+log(Q0))
 
      RETURN
 
      END



