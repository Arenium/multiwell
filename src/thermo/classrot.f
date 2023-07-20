      SUBROUTINE classrot(T,B,d,sig,Qrot,C,S,H)
 
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
c      Rotational partition function (see Robinson & Holbrook, Appendix 5)
c
c         INPUT
c      T      = temperature
c      B      = rotational constant (cm-1)
c      d      = rotor dimension
c      sig      = rotor symmetry (foldedness)
c
c         OUTPUT
c      Qrot   = partition function
c      C      = constant p heat capacity (units of R)
c      S      = entropy (units of R)
c      H      = [H(T)-H(0)]/RT
c
 
      IMPLICIT NONE
      REAL(8) Qrot , T , B , ALNGAM , GAM , C , S , H , hor
      INTEGER(4) d , sig
      PARAMETER (hor = 1.4387752D+00)
      EXTERNAL ALNGAM
      SAVE 
 
      C = 0.5D+00*d
      GAM = EXP( ALNGAM( (REAL(d,8) / 2.0D+00) ) )
      Qrot = ( 1.0D+00/sig )*GAM*( SQRT( T/(hor*B) ) )**d
      S = log(Qrot) + 0.5D+00*d
      H = 0.5D+00*d
 
      RETURN
 
      END SUBROUTINE
