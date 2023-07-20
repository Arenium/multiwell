**==gam.spg  processed by SPAG 5.11R  at 20:12 on 26 Mar 2001
      DOUBLE PRECISION FUNCTION GAM(X)
 
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
 
c
c    gamma function of X, from FUNCTION ALNGAM
c
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
      EXTERNAL ALNGAM
      SAVE 
      XX = X
      GAM = EXP(ALNGAM(XX))
 100  RETURN
      END
