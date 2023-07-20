 
      SUBROUTINE CMdist(XCOORD,YCOORD,ZCOORD,MASS,FRGM1,FRGM2,NATOMS,
     &   re,XCM1,XCM,YCM1,YCM,ZCM1,ZCM)
 
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
 
c     Calculate distance between centers of mass
 
      IMPLICIT NONE
      INTEGER I , J, NUMROT, NATOMS, FRGM1, FRGM2
      DOUBLE PRECISION XCOORD(100) , YCOORD(100) , ZCOORD(100)
      DOUBLE PRECISION MASS(100)
      DOUBLE PRECISION XCM1, XCM, YCM1, YCM, ZCM1, ZCM
      DOUBLE PRECISION MASSFRG1, MASSFRG2, re , xx

        XCM1=0  
        YCM1=0  
        ZCM1=0  
        XCM=0
        YCM=0
        ZCM=0
        MASSFRG1=0
        MASSFRG2=0

        do I=1,FRGM1     !   CM fragment1
         XCM=XCM+XCOORD(I)*MASS(I)
         YCM=YCM+YCOORD(I)*MASS(I)
         ZCM=ZCM+ZCOORD(I)*MASS(I)
         MASSFRG1=MASSFRG1+MASS(I)
        enddo
        XCM=XCM/MASSFRG1
        YCM=YCM/MASSFRG1
        ZCM=ZCM/MASSFRG1


        do I=FRGM1+1,NATOMS    !   CM fragment2
         XCM1=XCM1+XCOORD(I)*MASS(I)
         YCM1=YCM1+YCOORD(I)*MASS(I)
         ZCM1=ZCM1+ZCOORD(I)*MASS(I)
         MASSFRG2=MASSFRG2+MASS(I)
        enddo
        XCM1=XCM1/MASSFRG2
        YCM1=YCM1/MASSFRG2
        ZCM1=ZCM1/MASSFRG2


       re=((XCM1-XCM)**2+(YCM1-YCM)**2+(ZCM1-ZCM)**2)**0.5  !distance

      RETURN
      END
 
 
