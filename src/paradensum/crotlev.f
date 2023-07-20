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

      !SUBROUTINE CROTLEV(Ttot,DELE,JMAX,B,NDIM,NSYMM)
      SUBROUTINE CROTLEV(DELE,JMAX,B,NDIM,NSYMM)
 
      USE decl_alloc
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
!      DIMENSION T(JMAX) , x(JMAX) , AT(JMAX)
      DIMENSION x(JMAX) 
!      SAVE 

      R2 = DBLE(NDIM)/2.D+00                          ! in Eq. 5.15
      FAC = GAM( R2 )/( NSYMM*GAM( R2 + 1.D+00 ) )    ! in Eq. 5.15
      x(1) = 1.0D+00
      DO  K = 2 , JMAX                                ! Convolution
         E = (K-1)*DELE
         x(K) = FAC*( SQRT((E**NDIM)/B) - SQRT( ((E-DELE)**NDIM)/B ) )      !  number of classical free rotor states in grain
         AT(K) = 0.0d+00
         DO I = 1 , K
             AT(K) = AT(K) + x(K)*Ttot(I)
         END DO         
      END DO
 
      DO K = 1 , JMAX
        Ttot(K) = AT(K)
        AT(K) = 0.0D+00
      END DO

      RETURN
      END

