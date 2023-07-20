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


      SUBROUTINE CONGI(N,M,GI,X,Y,Z,AM,DELA)

        Integer NMAX, MMAX
        Parameter (NMAX=100, MMAX=500)

        Integer N, M, I, J
        Double precision DELA, AM(NMAX), PI, DEL
        Double precision X(MMAX,NMAX), Y(MMAX,NMAX), Z(MMAX,NMAX)
        Double precision A(3,3), GI(MMAX,4,4)
        Double precision ABX(NMAX), ABY(NMAX), ABZ(NMAX)
        Double precision DX(NMAX), DY(NMAX), DZ(NMAX)

CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C
C To construct G(-1) matrix
C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC

      PI=acos(-1.0d0)
      DEL=DELA*PI/180.0d0

      DO J=2, M-1

      DO I=1, N
            DX(I)=(X(J+1,I)-X(J-1,I))/2/DEL   ! The first derivative of X using the central different approach
            DY(I)=(Y(J+1,I)-Y(J-1,I))/2/DEL
            DZ(I)=(Z(J+1,I)-Z(J-1,I))/2/DEL
      ENDDO

      DO I=1, N
            ABX(I)=Y(J,I)*DZ(I)-Z(J,I)*DY(I) ! The product of two vectors: X and DX
            ABY(I)=Z(J,I)*DX(I)-X(J,I)*DZ(I)
            ABZ(I)=X(J,I)*DY(I)-Y(J,I)*DX(I)
      ENDDO

      DO I=1, 3
            DO K=1, 3
                  A(I,K)=0.0d0
            ENDDO
      ENDDO

      DO I=1, N
            A(1,1)=A(1,1)+AM(I)*( Y(J,I)**2 + Z(J,I)**2 )
            A(2,2)=A(2,2)+AM(I)*( X(J,I)**2 + Z(J,I)**2 )
            A(3,3)=A(3,3)+AM(I)*( X(J,I)**2 + Y(J,I)**2 )
            A(1,2)=A(1,2)-AM(I)*X(J,I)*Y(J,I)
            A(2,3)=A(2,3)-AM(I)*Y(J,I)*Z(J,I)
            A(1,3)=A(1,3)-AM(I)*X(J,I)*Z(J,I)
      ENDDO
      A(2,1)=A(1,2)
      A(3,1)=A(1,3)
      A(3,2)=A(2,3)

      DO I=1, 3
            DO K=1, 3
                  GI(J,I,K)=A(I,K)
            ENDDO
      ENDDO

      DO I=1, 4
            GI(J,4,I)=0.0d0
      ENDDO

      DO I=1, N
            GI(J,4,1)=GI(J,4,1)+AM(I)*ABX(I)
            GI(J,4,2)=GI(J,4,2)+AM(I)*ABY(I)
            GI(J,4,3)=GI(J,4,3)+AM(I)*ABZ(I)
            GI(J,4,4)=GI(J,4,4)+AM(I)*( DX(I)**2 + DY(I)**2 + DZ(I)**2 )
      ENDDO
      GI(J,1,4)=GI(J,4,1)
      GI(J,2,4)=GI(J,4,2)
      GI(J,3,4)=GI(J,4,3)

      ENDDO

      RETURN
      END


