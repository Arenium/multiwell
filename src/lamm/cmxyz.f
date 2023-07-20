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


      SUBROUTINE CMXYZ(N,M,X,Y,Z,AMT,AM)

        Integer NMAX, MMAX
        Parameter (NMAX=100, MMAX=500)

        Integer N, M, I, J
        Double precision AM(NMAX), AMT
        Double precision X(MMAX,NMAX), Y(MMAX,NMAX), Z(MMAX,NMAX)
        Double precision CMX, CMY, CMZ

CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C
C To transform X,Y,Z to the central mass of the molecule
C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC

      DO I=1, M
            CMX=0.0d0
            CMY=0.0d0
            CMZ=0.0d0
            DO J=1, N
                  CMX=CMX+AM(J)*X(I,J)
                  CMY=CMY+AM(J)*Y(I,J)
                  CMZ=CMZ+AM(J)*Z(I,J)
            ENDDO

            CMX=CMX/AMT
            CMY=CMY/AMT
            CMZ=CMZ/AMT
            DO J=1, N
                  X(I,J)=X(I,J)-CMX
                  Y(I,J)=Y(I,J)-CMY
                  Z(I,J)=Z(I,J)-CMZ
            ENDDO
      ENDDO
      RETURN
      END


