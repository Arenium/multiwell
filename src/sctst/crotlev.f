 
c%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
c
c DenSum: a code for calculating sums and densities of states.
c Copyright (C) 2001 John R. Barker
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

      SUBROUTINE CROTLEV(T,DELE,JMAX,B,NDIM,NSYMM)
c
C      CALCULATES ENERGY LEVELS FOR CLASSICAL FREE ROTOR
C
c      DELE      = energy grain size
c      JMAX      = index of Emax
c      B         = rotational constant (cm-1)
c      NDIM      = rotor dimension
c      NMAX      = ceiling of array
c
c      Eq. numbers refer to Robinson & Holbrook
c
      IMPLICIT REAL(8)(A-H,O-Z)
      DIMENSION T(JMAX) , x(JMAX) , AT(JMAX)
      SAVE 

      R2 = DBLE(NDIM)/2.D+00                          ! in Eq. 5.15
      FAC = GAM( R2 )/( NSYMM*GAM( R2 + 1.D+00 ) )    ! in Eq. 5.15
      x(1) = 1.0D+00
      DO  K = 2 , JMAX                                ! Convolution
         E = (K-1)*DELE
         x(K) = FAC*( SQRT((E**NDIM)/B) - SQRT( ((E-DELE)**NDIM)/B ) )      !  number of classical free rotor states in grain
         AT(K) = 0.0d+00
         DO I = 1 , K
             AT(K) = AT(K) + x(K)*T(I)
         END DO         
      END DO
 
      DO K = 1 , JMAX
        T(K) = AT(K)
        AT(K) = 0.0D+00
      END DO

      RETURN
      END

