c%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
c
c Copyright (C) 2009 Thanh Lam Nguyen
c
c Authors: Thanh Lam Nguyen   
c          nguyenlt@umich.edu       
c          December 2009
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
c See the 'ReadMe' file for a copy of the GNU General Public License,
c or contact:
c
c Free Software Foundation, Inc.
c 59 Temple Place - Suite 330
c Boston, MA 02111-1307, USA.
c
c%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
c
c    PROGRAM Large-Amplitude Motion (LAMM) 
c
c    by
c
c    Thanh Lam Nguyen    
c    nguyenlt@umich.edu        
c
c
c
c@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
c    Citations:
c
c    1. M. A. Harthcock and J. Laane, J. Phys. Chem. 89, 4231-4240 (2001).
c
c    2. T. L. Nguyen and J. R. Barker, to be published (2010).
c
c@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

C
C DEFINITION for VARIABLES used in the Red-I program
C

CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C
C N       = number of atoms in the molecule
C M       = number of optimized points on the potential energy surface
C AM(I)   = mass of ith atom, in amu
C ATM     = total mass of molecule, in g/mol or amu
C X, Y, Z = coordinates
C ANG     = a vector of torsional/inversion angle
C Amin    = minimum of ANGLE
C Amax    = maximum of ANGLE
C DELA    = stepsize of ANGLE
C GI      = GI matrix
C G44     = 1/I
C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC


      Program LAMM 

        Integer NMAX, MMAX
        Parameter (NMAX=100, MMAX=500)

        Integer N, M, I, J
        Double precision AMT, DELA, ANG(MMAX), AM(NMAX)
        Double precision X(MMAX,NMAX), Y(MMAX,NMAX), Z(MMAX,NMAX)
        Double precision GI(MMAX,4,4)
        Double precision G44(MMAX), E(MMAX)

      Character(len=100) C1 , C2 , C3 

CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C
C To compute reduced moment of inertia as a function of torsional/inversion/ring-puckering angle
C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC


C
C To read data in 
C

      CALL READ_IN(N,M,ANG,X,Y,Z,AMT,AM,DELA,C1,C2,C3,E)

C
C To identify center of mass and to transform the origin to the center of mass 
C

      CALL CMXYZ(N,M,X,Y,Z,AMT,AM)

C
C To construct ro-vibrational matrix GI
C

      CALL CONGI(N,M,GI,X,Y,Z,AM,DELA)

C
C To invert GI matrix to gain G44 = 1/I
C

      CALL INVG(G44,GI,M)

C
C To write out results for the effective reduced I
C

      CALL WRITE_OUT(G44,ANG,M,C1,C2,C3,DELA,E)

      END


