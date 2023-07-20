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
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@c
c%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
C
C	CALCULATES ENERGY LEVELS FOR MORSE OSCILLATOR
C
c	DELE	= energy grain
c	JMAX	= index of Emax
c	WE	= harmonic vibration frequency (cm-1)
c	XE	= anharmonicity (cm-1)
c	IR	= vector of energy level indices
c	NMAX	= ceiling of array (e.g. dissociation energy for Morse Osc.)
c	zpe     = zero point energy (cm-1), including anharmonicity
c
C	Input OBSERVED WEo = WE + 2*XE for 0-1 transition
c	Note sign of XE:  negative for usual Morse Oscillator
c
c	E = WE*(v+1/2) + XE*(v+1/2)^2
C
C     ALL ENERGIES AND FREQUENCIES IN WAVENUMBERS, relative to ZPE = WE/2 + XE/4
C
c--------------------------------------------------------------------------------
      !SUBROUTINE MORLEV(Ttot,AT,DELE,JMAX,WE,XE,zpe)
      SUBROUTINE MORLEV(DELE,JMAX,WE,XE,zpe_1)

      USE decl_alloc

      IMPLICIT REAL(8) (A-H,O-Z)
      EXTERNAL STERAB , GAM , WRAB , WRDEN , ALNGAM
      DIMENSION IR(20001)
!      DIMENSION T(JMAX), AT(JMAX)
!      SAVE 
 
      w_0 = WE + XE
      
      IF ( XE.LT.0.0 ) THEN
         RMAX = -0.5*w_0/XE                              ! Maximum bound level, E relative to ZPE
      ELSE
         RMAX = JMAX*DELE/WE                            ! Energy Ceiling level
      ENDIF
      NMAX = INT(RMAX) 
      
      WRITE(100,*) 'NMAX=', NMAX 

      zpe_1 = 0.5d+00*WE + 0.25d+00*XE                    ! Zero point energy at v=0

      DO 100 I = 1 , NMAX                               ! Start at v=1
        vi = DBLE(I)
        R = (w_0 + XE*vi)*vi/DELE
      IR(I) = IDNINT(R)                                 ! Nearest integer: number of grains
 100  CONTINUE
 
      DO J = 1 , NMAX                                   ! NMAX is the number of energy levels
        DO K = 1 , JMAX
          KARG = K - IR(J)
          IF ( KARG.GT.0 .AND. KARG.LE.JMAX ) AT(K) = AT(K) + Ttot(KARG)
        END DO
      END DO

      DO J = 1 , JMAX
         Ttot(J) = AT(J) + Ttot(J)
         AT(J) = 0.0D+00
      END DO

      RETURN
      END
