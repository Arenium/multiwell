!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!    
!    LICENSE NOTICE
!
!    Copyright (C) 2019 Michele Ceotto, Chiara Aieta, Fabio Gabas, 
!                       Thanh Lam Nguyen, and John R. Barker
!
!    Contact:
!
!    Michele Ceotto  (email: michele.ceotto@unimi.it)
!    Dipartimento di Chimica
!    Università degli Studi di Milano
!    via Golgi 19, Milano 20133, ITALY
!
!    or:
!
!    John R. Barker   (email: jrbarker@umich.edu)
!    Department of Atmospheric, Oceanic, and Space Sciences
!    College of Engineering
!    University of Michigan
!    Ann Arbor, MI 48109
!
!    http://clasp-research.engin.umich.edu/multiwell
!    This program is free software; you can redistribute it and/or
!    modify it under the terms of the GNU General Public License (version 2)
!    as published by the Free Software Foundation.
!   
!    This program is distributed in the hope that it will be useful,
!    but WITHOUT ANY WARRANTY; without even the implied warranty of
!    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
!    GNU General Public License for more details.
!   
!    See the 'ReadMe' file for a copy of the GNU General Public License,
!    or contact:
!   
!    Free Software Foundation, Inc.
!    59 Temple Place - Suite 330
!    Boston, MA 02111-1307, USA.
!
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!                                                                 
!                         PROGRAM paradensum    
!
!                               by      
!
!           by Michele Ceotto, Chiara Aieta, Fabio Gabas,
!               Thanh Lam Nguyen, and John R. Barker    
!                                                                  
!           ***PARallel Anharmonic DENsities and SUM of states***
!
!                             based on
!
!          The parallel implementation* of the Wang-Landau algorithms 
!          for densities of states**
!
!    Literature Citations:                                                    
!    *parallel implementation of anaharmonic density of states
!    C. Aieta, F. Gabas and M. Ceotto, J. Phys. Chem. A., 120(27), 
!    4853-4862 (2016).
!                                                                  
!    **Density of states algorithms
!    F. Wang and D. P. Landau, Phys. Rev. Letters 86, 2050-2053 (2001).           
!    M. Basire, P. Parneix, and F. Calvo, J. Chem.Phys. 129, 081101 (2008).    
!    T. L. Nguyen and J.R. Barker, J. Phys. Chem. A., 114, 3718–3730 (2010).
!    
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

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
c      SUBROUTINE STOPLEV(T,AT,DELE,JMAX,B2,B1,NSYMM)
      SUBROUTINE STOPLEV(DELE,JMAX,B2,B1,NSYMM)

      USE decl_alloc
 
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
!      Double precision AT(JMAX) , T(JMAX)  
!      SAVE 
 
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
                      AT(kk) = AT(kk) + F*Ttot(KARG)
                   ENDIF
                END DO  ! k over energy range
             ENDIF   ! if R in energy range
          ENDDO ! k
       ENDDO   ! j
                 
      DO J = 1 , JMAX     ! over energy grains
         Ttot(J) = AT(J) + Ttot(J)
         Ttot(J) = Ttot(J) / NSYMM
         AT(J) = 0.0D+00
      END DO  !  J
 
      RETURN
      END
