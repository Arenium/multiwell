c%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
c
c Copyright (C) 2009 John R. Barker
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

      SUBROUTINE shuffle ( ns, w0, x0, w, x, nx)
c
c	W0 = Inital frequency vector
c	x0 = Initial anharmonicity coefficient matrix
c	w  = shuffled frequency vector
c	x  = shuffled anharmonicity coefficient matrix
c	nx = shuffled vector of initial index numbers

      IMPLICIT NONE
      INTEGER(4) ns, index(100), idum, i, j, k, m, nx(100), nt(100)
      REAL(8) r(100), RAN1, rtemp, 
     &                 w(100), x(100,100), w0(100), x0(100,100)
      COMMON/I/ idum
      EXTERNAL RAN1
      SAVE

      DO i = 1, ns
        r(i) = RAN1(idum)
      END DO

      DO j = 1 , ns
        rtemp = r(1)
        k = 1
        DO i = 2, ns
          IF ( r(i) .gt. rtemp ) THEN
            k = i
            rtemp = r(i)
          ENDIF
        END DO
        index(k) = j
        r(k) = 0.0
      END DO

      DO i = 1, ns
        nt(i) = nx( index(i) )
        w(i)  = w0( index(i) )
        DO j = 1, ns
          x(i,j) = x0( index(i),index(j) )
        END DO
      END DO
      
      DO i = 1, ns
        nx(i) = nt(i)
      END DO

      RETURN
      END

       

      
