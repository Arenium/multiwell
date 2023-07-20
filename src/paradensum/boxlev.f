      !SUBROUTINE BOXLEV(Ttot,AT,DELE,JMAX,WE,NSYMM)
      SUBROUTINE BOXLEV(DELE,JMAX,WE,NSYMM)

      USE decl_alloc
C
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
C 
C
C
C
C      CALCULATES ENERGY LEVELS FOR PARTICLE-IN-A-BOX TO BE USED
C          WITH STERAB
C
c      DELE      = energy grain
c      JMAX      = index of Emax
c      WE      = vibration frequency (cm-1)
c      IR      = vector of energy level indices
c      NMAX      = ceiling of array
c       zpe     = zero point energy (cm-1)
c
c      E = WE*n^2  for n=1,2,3,...
C
C     OUTPUT ENERGIES IN WAVENUMBERS, relative to ZPE = WE
C
c    2/02 counted as vibs in Whitten-Rabinovitch
c         parameter calculations.
C
      IMPLICIT REAL(8)(A-H,O-Z)
      DIMENSION IR(20001)
!      DIMENSION T(JMAX), AT(JMAX)

!      SAVE 
 
      Emax=(JMAX-1)*DELE
      RMAX = SQRT(Emax/WE)                           ! Maximum bound level
      NMAX = INT(RMAX) 
      IF ( NMAX.GT.20001 ) NMAX = 20001
      zpe = WE                                            ! Zero point energy
      DO 100 I = 1 , NMAX
        v2 = DBLE(I+1)*DBLE(I+1)                          ! "zero" level = 1, start loop at I+1=2
        R = (WE*v2 - zpe)/DELE                            ! Relative to zpe = WE
        IR(I) = NINT(R)                                 ! Nearest integer: number of grains
 100  CONTINUE

            DO 390 L = 1 , NSYMM                                ! Loop for degenerate vibrations
               DO 355 J = 1 , NMAX                              ! NMAX is the number of energy levels
                  DO 365 K = 1 , JMAX
                     KARG = K - IR(J)
                     IF ( KARG.GT.0 .AND. KARG.LE.JMAX )
     &                   AT(K) = AT(K) + Ttot(KARG)
 365              CONTINUE
 355           CONTINUE

               DO 370 J = 1 , JMAX
                  Ttot(J) = AT(J) + Ttot(J)
                  AT(J) = 0.0D+00
 370           CONTINUE
 390        CONTINUE


      RETURN
      END
