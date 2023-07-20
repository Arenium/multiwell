**==qrot.spg  processed by SPAG 5.11R  at 20:12 on 26 Mar 2001
      FUNCTION Qele(T , ns , n , Elev , g , datpts , C , S , H)
c%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
 
 
c
c Thermo: a code for thermochemical calculations.
c Copyright (C) 2004 John R. Barker
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
c      ELECTRONIC PARTITION FUNCTION
c
c      T           = temperature
c      ns          = species index number
c      n           = number of electronic levels
c      Elev(ns,n)  = array of level energies (lowest level at Elev = 0)
c      g(ns,n)     = array of level degeneracies
c	   datpts	   = max a;;oed data points (see DECLARE_INC.f)
c      C           = constant p heat capacity (units of R)
c      S           = entropy (units of R)
c      H           = [H(T)-H(0)]/RT
c
 
      IMPLICIT NONE
      REAL(8) Qele , T , Elev , C , S , H , 
     &                 Ered , Fac , hor , Q0 , Q1 , Q2
      INTEGER(4) n , ns , i , g , datpts
      DIMENSION Elev(datpts,10) , g(datpts,10)		! Explicit shape
      PARAMETER (hor = 1.4387752D+00)
      SAVE
         Q0 = 0.0d+00
         Q1 = 0.0d+00
         Q2 = 0.00d+00
      
      DO 100 i = 1 , n
         Ered = hor*Elev(ns,i)/T
         Fac = g(ns,i)*EXP(-Ered)
         Q0 = Q0 + Fac
         Q1 = Q1 + Ered*Fac
         Q2 = Q2 + Ered*Ered*Fac
 100  CONTINUE

      Qele = Q0
      C = Q2/Q0 - (Q1/Q0)**2
      H = Q1/Q0
      S = ( Q1/Q0 + log(Q0) )
      RETURN
 
      END
