c%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
c
c DenSum: a code for calculating sums and densities of states.
c Copyright (C) 20161 John R. Barker
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
 
c    Last modified:
c
c    5/2016 Programmed by J. R. Barker
CC
C	CALCULATES ENERGY LEVELS FOR SYMMETRIC TOP TO BE USED
C          WITH STERAB
C
c	DELE	= energy grain size
c	JMAX	= number of energy grains corresponding to Emax2
c	B1	= 1D rotor rotational constant (cm-1)
c	B2	= 2D rotor rotational constant (cm-1)
c         NSYMM     = rotor symmetry
c	IR	= vector of energy level indices
c
c
      SUBROUTINE STOPLEV(T,AT,DELE,JMAX,B2,B1,NSYMM)
 
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
      Double precision AT(JMAX) , T(JMAX)  
      SAVE 
 
      EMAX = DELE*(JMAX-1)
c      write(12,*) '  '
c      write(12,*) 'Symmetric Top (top)'
c      write(12,*) 
c     & '   Index           J           K         cm-1            degen'
 
      jjmax = SQRT( EMAX / B2 ) + 1
      jj = 0
      R = 0.0d+00
      DO j = 1 , jjmax                                    ! start at j=1
          Ej = B2*j*(j + 1.d+00)
          DO k = 0 , j
             Ek = ( B1 - B2 )*k*k
             R = ( Ek + Ej)
             IF ( R.GE.0.0 .AND. R.LE.EMAX) THEN
                jj = jj + 1                               ! = level number
                IR = NINT( R/DELE ) + 1                   ! Nearest integer: number of grains
                IF ( k .EQ. 0 ) THEN
                   F = ( 2.d+0*j + 1.d+0 )                ! 2J+1 degeneracy when K=0
                 ELSE
                   F = 2.d+0 * ( 2.d+0*j + 1.d+0 )        ! 2(2J+1) when K >0
                ENDIF
c                write(12,*) jj , j , k , R , F   ! write out the energy level and degeneracy
                DO kk = 1 , JMAX                     ! Jmax = number of energy grains (i.e. corresponding to Emax2)
                   KARG = kk - IR + 1.d+0
                   IF ( KARG.GT.0 .AND. KARG.LE.JMAX) THEN
                      AT(kk) = AT(kk) + F*T(KARG)
                   ENDIF
                END DO  ! k over energy range
             ENDIF   ! if R in energy range
          ENDDO ! k
       ENDDO   ! j
                 
      DO J = 1 , JMAX     ! over energy grains
         T(J) = AT(J) + T(J)
         T(J) = T(J) / NSYMM
         AT(J) = 0.0D+00
      END DO  !  J
 
      RETURN
      END
