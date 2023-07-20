      FUNCTION tQstop(T,B2,B1,sig,C,S,H)
 
c%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
c
c Thermo: a code for thermochemical calculations.
c Copyright (C) 2016 John R. Barker
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
 
c
c      Rotational partition function for symmetric top rotor
c
c      NOTE: given rotational constants A, B=C, then B1=A and B2=B=C
c
c      B1     = 1D rotational constant (cm-1) 
c      B2     = 2D rotational constant (cm-1)
c      sig    = symmetry number
c      T      = temperature
c      C      = heat capacity at constant p (units of R)
c      S      = entropy (units of R)
c      H      = [H(T)-H(0)]/RT
c
 
      IMPLICIT NONE
      REAL(8) T, B2, B1, C , S , H , R , Q0 , Q1 , Q2 , Emax , Ered ,
     &        Ej , Ek , hor , H2, S2 , term , tQstop
      INTEGER(4)  sig, j , jj , k , kmax , jmax
      PARAMETER (hor = 1.4387752D+00)
      SAVE 

      Emax = 20 * T / hor
      jmax = SQRT( Emax/B2 ) + 1
       
       Q0 = 0.d+0
       Q1 = 0.d+0
       Q2 = 0.d+0
           
      R = 0.0d+0
      jj = 0
      DO j = 0 , jmax
          Ej = B2*j*(j + 1.d+00)
          DO k = 0 , j
             Ek = ( B1 - B2 )*k*k
             R = ( Ek + Ej)
             Ered = R*hor/T
             IF ( Ered .GE. 0.0) THEN
                jj = jj + 1                               ! = level number
c                write(12,*) jj , j , k , R   ! write out the energy levels
                IF ( k .EQ. 0 ) THEN
                   term = ( 2.d+0*j + 1.d+0 ) * exp( -Ered ) 
                 ELSE
                   term = 2.d+0 * ( 2.d+0*j + 1.d+0 ) * exp( -Ered ) 
                ENDIF
                Q0 = Q0 + term
                Q1 = Q1 + Ered*term
                Q2 = Q2 + Ered*Ered*term
             ENDIF   ! if Ered≥0
          ENDDO ! k
       ENDDO   ! j

c      Thermo formulae from Pitzer and Gwinn, J. Chem. Phys. 10, 428 (1942)
 
      Q0 = Q0/sig
      Q1 = Q1/sig
      Q2 = Q2/sig
      C = Q2/Q0 - (Q1/Q0)**2
      H = Q1/Q0
      S = (H+log(Q0))
 
      tQstop = Q0

      RETURN
 
      END
