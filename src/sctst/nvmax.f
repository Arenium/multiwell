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

       INTEGER(4) FUNCTION nvmax( ns, w, x, nv, k, Eu )

c     Select highest quantum number to be sampled

      IMPLICIT NONE
c      INTEGER(4) ns, index(100), nv(100), i, k, nvd
c      INTEGER(4) ns, nv(100), i, k, nvd
      INTEGER(4) ns, nv(100), i, k
c      REAL(8) r(100), w(100), x(100,100), Eu, D , vd
      REAL(8) w(100), x(100,100), Eu, D , vd
c      REAL(8) vmax, sum, energy, enrg, zpe
      REAL(8) vmax, sum, zpe
      SAVE

      IF ( Eu .LT. 0.0d+00 ) THEN
         nvmax = -1
         RETURN
      ELSE
         nvmax = 0
      ENDIF

      IF ( x(k,k) .NE. 0.0 ) THEN
         sum = 0.0d+00
         DO i = 1, ns
           IF ( i .NE. k ) sum = sum + ( nv(i) + 0.5d0 )*x(i,k)
         END DO
        zpe = x(k,k)/4.0d0 + w(k)/2.0d0 + sum/2.0d0

         vd = -(w(k) + sum)/( 2.0d+00*x(k,k) ) - 0.5d0             ! Eq. 15: quantum no. of highest bound energy level
         D = -( ( w(k) + sum )**2)/( 4.0d+00*x(k,k) ) - zpe       ! Do (not De)
         vmax = vd*( 1.0d+00 - SQRT( 1.0d+00 - Eu/D) )      ! Eq. 16

         IF ( (w(k) + sum) .GT. 0.0 .AND. x(k,k) .LT. 0.0) THEN
           IF ( Eu .GT. D ) THEN
             nvmax = INT( vd )
           ELSE
             nvmax = INT( vmax )
           ENDIF
         ELSEIF ( (w(k) + sum) .GT. 0.0 .AND. x(k,k) .GT. 0.0) THEN
           nvmax = INT( vmax )
         ELSE                                ! (w(k) + sum) < 0
           nvmax = -1
         ENDIF
      ELSEIF ( x(k,k) .EQ. 0.0 ) THEN
           nvmax = INT( Eu/w(k) )                              ! harmonic when x(k,k) = 0
      ENDIF   ! Eu < 0

      RETURN
      END


