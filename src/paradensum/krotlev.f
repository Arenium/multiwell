      !SUBROUTINE KROTLEV(Ttot,AT,DELE,JMAX,B,JK,NSYMM)
      SUBROUTINE KROTLEV(DELE,JMAX,B,JK,NSYMM)
 
      USE decl_alloc
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
c%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
C
C      CALCULATES ENERGY LEVELS FOR 1-DIMENSIONAL K-ROTOR TO BE USED
C          WITH STERAB
C
c      DELE      = energy grain size
c      JMAX      = index of Emax
c      B      = rotational constant (cm-1)
c      JK      = quantum number for total angular momentum
c      IR      = vector of energy level indices
c      NMAX      = ceiling of array
c
c
      IMPLICIT REAL(8)(A-H,O-Z)
      DIMENSION IR(20001)
!      DIMENSION T(JMAX), AT(JMAX)
!      SAVE 


      Emax=(JMAX-1)*DELE 
      RMAX = SQRT(Emax/B)          ! Highest rotational quantum number consistent with energy ceiling
      NMAX = INT(RMAX)             ! INTEGER(4) version
      IF ( NMAX .GT. JK ) NMAX = JK
      IF ( NMAX .GT. 20001 ) NMAX = 20001
       
      IF ( NMAX .GT. 0 ) THEN
      DO 100 K = 1 , NMAX               ! K = quantum number of 1-dimensional free rotor
         R = B*K*K/DELE                 ! index number starting at 0
         IR(K) = NINT(R)              ! Nearest INTEGER(4): number of grains
 100  CONTINUE

CCCCCCCCCCCCCCCCCCCCCCC

            DO 281 J = 1 , NMAX ! J = quantum number for K-rotor;  NMAX is the number of energy levels
               IF ( J.EQ.0 ) THEN
                  F = 1.d+00         ! F is the multiplicity of the rotational level (see J.L. McHale, Molecular Spectroscopy (Prentice-Hall, 1999), 216f.
               ELSE
                  F = 2.d+00         ! 1-D rotor (e.g. single-axis internal rotor)
               ENDIF
               DO 261 K = 1 , JMAX
                  KARG = K - IR(J)
                  IF ( KARG.GT.0 .AND. KARG.LE.JMAX)
     &                   AT(K) = AT(K) + F*Ttot(KARG)
 261           CONTINUE
 281        CONTINUE

            DO 301 J = 1 , JMAX
               Ttot(J) = AT(J) + Ttot(J)
               Ttot(J) = Ttot(J) / NSYMM
               AT(J) = 0.0D+00
 301        CONTINUE

CCCCCCCCCCCCCCCCCCCCCC

      ENDIF
 
      RETURN
      END
