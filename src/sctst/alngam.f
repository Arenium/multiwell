**==alngam.spg  processed by SPAG 5.11R  at 20:12 on 26 Mar 2001
 
      DOUBLE PRECISION FUNCTION ALNGAM(YY)
 
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
c	Natural logarithm of the gamma FUNCTION of X
c
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
      SAVE 
      PI = 4.0D+00*DATAN(1.0D+00)
      TWOPI = 2.0D+00*PI
      X = YY
      Z = 0.0
 100  IF ( X.GE.1.5 ) THEN
         X2 = X*X
         X3 = X2*X
         X5 = X2*X3
         X7 = X2*X5
         ALNGAM = -Z + (X-0.5)*LOG(X) - X + 0.5*LOG(TWOPI)
     &            + 1.0/(12.0*X) - 1.0/(360.0*X3) + 1.0/(1260.0*X5)
     &            - 1.0/(1680.0*X7)
      ELSE
         ALNGAM = 0.0
         IF ( X.GT.0.0 ) THEN
            Z = Z + LOG(X)
            X = X + 1.0
            GOTO 100
         ENDIF
      ENDIF
 200  RETURN
      END
