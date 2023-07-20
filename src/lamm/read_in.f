c%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
c
c Copyright (C) 2009 Thanh Lam Nguyen
c
c Authors: Thanh Lam Nguyen 
c          nguyenlt@umich.edu
c          Dec. 23, 2009
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

C
C To read data for the LAMM program
C


      SUBROUTINE READ_IN(N,M,ANG,X,Y,Z,AMT,AM,DELA,C1,C2,C3,E)
        Integer NMAX, MMAX
      Parameter (NMAX=100, MMAX=500)

      Integer N, M, I, J
        Double precision AMT, Amin, Amax, DELA, ANG(MMAX), AM(NMAX)
        Double precision X(MMAX,NMAX), Y(MMAX,NMAX), Z(MMAX,NMAX)
        Double precision E(MMAX), TANG

      Character(len=100) C1 , C2 , C3 

CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
!
!	N       = number of atoms in the molecule 
!	M       = number of optimized points on the potential energy surface
! 	AM(I)   = mass of ith atom, in amu 
! 	ATM     = total mass of molecule, in g/mol or amu
! 	X, Y, Z = coordinates
! 	ANG     = a vector of torsional/inversion angle
! 	Amin    = minimum of ANG
! 	Amax    = maximum of ANG
! 	DELA    = stepsize of ANG 
!
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC

        open (unit=1,file='lamm.dat',status='old',
     &  access='sequential')

      read(1,99) C1
      read(1,99) C2

      read(1,*) N
      read(1,*) Amin, M, DELA
        AMT=0.0d0
        DO I=1, N
                read(1,*) AM(I)
                AMT=AMT + AM(I)
        ENDDO

!!!      M=INT((Amax-Amin)/DELA) + 1

      Amax = (M-1)*DELA + Amin

      DO I=1, M
                ANG(I)=Amin + (I-1)*DELA
      ENDDO
      
      DO I=1, M
            DO J=1, N
                  read(1,*) X(I,J), Y(I,J), Z(I,J)
            ENDDO
      ENDDO

      read(1,*) C3

      DO I=1, M
            read(1,*) TANG, E(I)
      ENDDO

        close (unit=1, status='keep')

99      FORMAT(A100)

      RETURN
      END

CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC

