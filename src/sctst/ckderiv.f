c%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
c
c Copyright (C) 2009 John R. Barker and Thanh Lam Nguyen
c
c Authors: John R. Barker and Thanh Lam Nguyen
c          nguyenlt@umich.edu
c          September 2009
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
c-----------------------------------------------------------------------------------------
c    Test partial derivatives of Evib with respect to each quantum number
c
c    for bound states, all partial derivatives should be >0
c    ntest = 0 for PASS
c          = 1 for FAIL
c-----------------------------------------------------------------------------------------

      SUBROUTINE ckderiv( ns , nv , wa, xa, ntest ) 
      
      IMPLICIT NONE
c      INTEGER(4) ns , i , j , k , ntest , nv(100)
      INTEGER(4) ns , i , k , ntest , nv(100)
      REAL(8) wa(100) , xa(100,100)
      REAL(8) deriv, sum

      SAVE
      
      ntest = 0
      k = 1
      DO WHILE ( (k .LE. ns) .AND. (ntest .EQ. 0) )    ! start checking
         sum = 0.0d+00
         DO i = 1 , ns
           sum = sum + xa(i,k)*( nv(i)+0.5d0 ) 
         END DO
         deriv = wa(k) + ( nv(k)+0.5d0 )*xa(k,k) + sum
         IF ( deriv .LE. 0.0d+00 ) ntest = 1           ! failed test
         k = k + 1
      END DO
        
      RETURN
      END


