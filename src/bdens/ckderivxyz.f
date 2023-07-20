c%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
c
c Copyright (C) 2009 John R. Barker
c
c Authors: Thanh Lam Nguyen and John R. Barker
c          nguyenlt@umich.edu
c          Aug. 9, 2009
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


CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
c
c To compute the first derivative of Ev energy function 
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

      SUBROUTINE CKDERIVXYZ(ns, nv, wa, xa, ya, za, ntest)

      IMPLICIT NONE
      INTEGER(4) ns , i , j , k , l, ntest , nv(100)
c      REAL(8) wa(100) , xa(100,100), deriv , sum 
      REAL(8) wa(100) , xa(100,100), deriv 
      REAL(8) ya(100,100,100) , za(100,100,100,100) 

        REAL(8) XX, A1, A2, A3, A4, tp

      SAVE


      ntest = 0
      k = 1
      DO WHILE ( (k .LE. ns) .AND. (ntest .EQ. 0) )    ! start checking

        A1=4.0d0*za(k,k,k,k)

        A2=ya(k,k,k)
        DO j=1, ns
                IF(j.ne.k) A2=A2+za(k,k,k,j)*(nv(j)+0.5d0)
        ENDDO
        A2=3.0d0*A2

        A3=xa(k,k)
        DO j=1, ns
                IF(j.ne.k) A3=A3+ya(k,k,j)*(nv(j)+0.5d0)
        ENDDO
        DO i=1, ns
                IF(i.ne.k) THEN
                        DO j=1, i
        IF(j.ne.k) A3=A3+za(k,k,i,j)*(nv(i)+0.5d0)*(nv(j)+0.5d0)
                        ENDDO
                ENDIF
        ENDDO
        A3=2.0d0*A3

        A4=wa(k)
        DO i=1, ns
                IF(i.ne.k) A4=A4+xa(k,i)*(nv(i)+0.5d0)
        ENDDO
        DO i=1, ns
                IF(i.ne.k) THEN
                        DO j=1, i
        IF(j.ne.k) A4=A4+ya(k,i,j)*(nv(i)+0.5d0)*(nv(j)+0.5d0)
                        ENDDO
                ENDIF
        ENDDO
        DO i=1, ns
                IF(i.ne.k) THEN
                DO j=1, i
                        IF(j.ne.k) THEN
                                DO l=1, j
                                IF(l.ne.k) THEN
        tp=(nv(i)+0.5d0)*(nv(j)+0.5d0)*(nv(l)+0.5d0)
        A4=A4+za(k,i,j,l)*tp
                                ENDIF
                                ENDDO
                        ENDIF
                ENDDO
                ENDIF
        ENDDO

        XX = 0.5d0 + nv(k) 
        deriv = A1*(XX**3) + A2*(XX**2) + A3*XX + A4

         IF ( deriv .LE. 0.0d+00 ) ntest = 1           ! failed test
         k = k + 1
      END DO

      RETURN
      END

CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC


