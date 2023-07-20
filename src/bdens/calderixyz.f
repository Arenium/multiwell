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

      SUBROUTINE CALDERIXYZ(F,NVK,ns,k,nv,w,x,y,z)
        INTEGER(4) ns, k, NVK, No
      PARAMETER (No=100)
        REAL(8) A1, A2, A3, A4
        REAL(8) x(No,No), y(No,No,No), z(No,No,No,No)
        REAL(8)  w(No)
        REAL(8) F
      INTEGER(4) nv(No)

C
C Definition: 
C      F   : return a value of partitial derivative of Ev at the mode k
C      NVK : NVK = nv(k) 
C

      CALL CALA1234XYZ(A1,A2,A3,A4,ns,k,nv,w,x,y,z)
      CALL CALFXYZ(F,NVK,A1,A2,A3,A4)

      RETURN
      END

CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC

      SUBROUTINE CALFXYZ(F,NVK,A1,A2,A3,A4)
      INTEGER(4) NVK
      REAL(8) X, F, A1, A2, A3, A4
      
      X = 0.5d0 + NVK
      F = A1*(X**3) + A2*(X**2) + A3*X + A4
      RETURN
      END

CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC

      SUBROUTINE CALA1234XYZ(A1,A2,A3,A4,ns,k,nv,w,x,y,z)
      INTEGER(4) ns, i, j, k, l, No
      PARAMETER (No=100)
        REAL(8) A1, A2, A3, A4
        REAL(8) x(No,No), y(No,No,No), z(No,No,No,No)
        REAL(8)  w(No)
      INTEGER(4) nv(No)

        A1=4.0d0*z(k,k,k,k)

        A2=y(k,k,k)
        DO j=1, ns
                IF(j.ne.k) A2=A2+z(k,k,k,j)*(nv(j)+0.5d0)
        ENDDO
        A2=3.0d0*A2

        A3=x(k,k)
        DO j=1, ns
                IF(j.ne.k) A3=A3+y(k,k,j)*(nv(j)+0.5d0)
        ENDDO
        DO i=1, ns
                IF(i.ne.k) THEN 
                        DO j=1, i
        IF(j.ne.k) A3=A3+z(k,k,i,j)*(nv(i)+0.5d0)*(nv(j)+0.5d0)
                        ENDDO
                ENDIF
        ENDDO
        A3=2.0d0*A3
        
        A4=w(k) 
        DO i=1, ns
                IF(i.ne.k) A4=A4+x(k,i)*(nv(i)+0.5d0)
        ENDDO
        DO i=1, ns
                IF(i.ne.k) THEN 
                        DO j=1, i
        IF(j.ne.k) A4=A4+y(k,i,j)*(nv(i)+0.5d0)*(nv(j)+0.5d0)
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
      A4=A4+z(k,i,j,l)*tp
                        ENDIF
                              ENDDO
                      ENDIF
            ENDDO
                ENDIF
        ENDDO

      RETURN
      END

CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC


