**==rotlev.spg  processed by SPAG 5.11R  at 20:12 on 26 Mar 2001
 
 
      SUBROUTINE ROTLEV(T,AT,DELE,JMAX,B,NDIM,NSYMM)
 
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
C
C      CALCULATES ENERGY LEVELS FOR FREE ROTOR TO BE USED
C          WITH STERAB
C
c      DELE      = energy grain size
c      JMAX      = index of Emax
c      B         = rotational constant (cm-1)
c      NDIM      = rotor dimension
c      IR        = vector of energy level indices
c      NMAX      = ceiling of array
c
c
      IMPLICIT REAL(8)(A-H,O-Z)
      DIMENSION IR(20001)
      DIMENSION T(JMAX), AT(JMAX)
      SAVE 

      Emax=(JMAX-1)*DELE 
      RMAX = SQRT(Emax/B)  ! Maximum bound level
      NMAX = INT(RMAX)  
      IF ( NMAX.GT.20001 ) NMAX = 20001
 
      DO 100 J = 1 , NMAX
         IF ( NDIM.EQ.1 ) THEN
             R = B*J*J/DELE
           ELSE
             R = B*J*(J+1)/DELE
         ENDIF
         IR(J) = NINT(R)              ! Nearest INTEGER(4): number of grains
 100  CONTINUE

CCCCCCCCCCCCCCCCCCCCCCCCC

      DO J = 1 , NMAX ! J = rot quantum number;  NMAX is the number of energy levels
        IF ( J.EQ.0 ) THEN
           F = 1d+00             ! F is the multiplicity of the rotational level (see J.L. McHale, Molecular Spectroscopy (Prentice-Hall, 1999), 216f.
        ELSEIF ( NDIM.EQ.1 ) THEN   ! 1-D rotor (e.g. single-axis internal rotor)
           F = 2.0d+00
        ELSEIF ( NDIM.EQ.2 ) THEN   ! 2-D rotor
           F = 2.0d+00*J + 1.0d+00
        ELSEIF ( NDIM.EQ.3 ) THEN   ! 3-D rotor (shperical top)
           F = (2.0d+00*J+1.0d+00)*(2.0d+00*J+1.0d+00)
        ENDIF
        DO K = 1 , JMAX
           KARG = K - IR(J)
           IF ( KARG.GT.0 .AND. KARG.LE.JMAX) THEN
               AT(K) = AT(K) + F*T(KARG)
           ENDIF
        END DO
      END DO

      DO K = 1 , JMAX
          T(K) = AT(K) + T(K)
          T(K) = T(K) / NSYMM
          AT(K) = 0.0D+00
      END DO

CCCCCCCCCCCCCCCCCCCCCCCC
 
      RETURN
      END
