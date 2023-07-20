**==rotate.spg  processed by SPAG 5.11R  at 20:12 on 26 Mar 2001
 
 
      SUBROUTINE ROTATE(XCOORD,YCOORD,ZCOORD,NATOMS)
 
c%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
 
 
c
c MomInert: a code for calculating moments of inertia
c Copyright (C) 2001 Nicolas F. Ortiz and John R. Barker
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
 
C     Modified by Amity Andersen and William Schneider (July 2000)
C
C     Rotates the cartesian coordinates with the three Euler angle
C     transformation matrices; the Euler angles are chosen randomly from 0 to
C     2*pi.  This procedure allows the cartesian coordinates of each atom to
C     change while leaving the geometry of the molecule intact, thereby
C     allowing the definition of planes in the main code.
 
      INTEGER I , NATOMS
      DOUBLE PRECISION XCOORD(150) , YCOORD(150) , ZCOORD(150)
      DOUBLE PRECISION X(150) , Y(150) , Z(150)
      DOUBLE PRECISION PHI , THETA , PSI , TWOPI , FAC
      DOUBLE PRECISION LAMBDA(3,3)
 
 
 
 
c%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
 
C..rotation algorithm using pseudorandom euler angles
C..define pi...
      TWOPI = 8.0D+00*ATAN(1.0D+00)
      FAC = 0.05D+00
      PHI = FAC*TWOPI
      THETA = FAC*TWOPI
      PSI = FAC*TWOPI

C..the phi, theta, and psi combined tranformation matrix, lambda
      LAMBDA(1,1) = COS(PSI)*COS(PHI) - SIN(PSI)*COS(THETA)*SIN(PHI)
      LAMBDA(1,2) = COS(PSI)*SIN(PHI) + SIN(PSI)*COS(THETA)*COS(PHI)
      LAMBDA(1,3) = SIN(PSI)*SIN(THETA)
      LAMBDA(2,1) = -SIN(PSI)*COS(PHI) - COS(PSI)*COS(THETA)*SIN(PHI)
      LAMBDA(2,2) = -SIN(PSI)*SIN(PHI) + COS(PSI)*COS(THETA)*COS(PHI)
      LAMBDA(2,3) = COS(PSI)*SIN(THETA)
      LAMBDA(3,1) = SIN(THETA)*SIN(PHI)
      LAMBDA(3,2) = -SIN(THETA)*COS(PHI)
      LAMBDA(3,3) = COS(THETA)
 
C..use matrix to rotate axes
      DO 100 I = 1 , NATOMS
         X(I) = LAMBDA(1,1)*XCOORD(I) + LAMBDA(1,2)*YCOORD(I)
     &          + LAMBDA(1,3)*ZCOORD(I)
         Y(I) = LAMBDA(2,1)*XCOORD(I) + LAMBDA(2,2)*YCOORD(I)
     &          + LAMBDA(2,3)*ZCOORD(I)
         Z(I) = LAMBDA(3,1)*XCOORD(I) + LAMBDA(3,2)*YCOORD(I)
     &          + LAMBDA(3,3)*ZCOORD(I)
         XCOORD(I) = X(I)
         YCOORD(I) = Y(I)
         ZCOORD(I) = Z(I)
 100  CONTINUE
 
      RETURN
      END
 
 
